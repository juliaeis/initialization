import os
import sys
import pandas as pd
import geopandas as gpd
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from oggm import cfg, utils, workflow, tasks
from oggm.utils._downloads import get_demo_file
from oggm.core.flowline import FluxBasedModel, FileModel
import oggm
from leclercq.leclercq_plots import *
sys.path.append('../')
sys.path.append('../../')
from  initialization.core import *
from paper.plots_paper import *


def fitness_function(ye, model1, model2):
    """
    calculates the objective value (difference in geometry)
    :param model1: oggm.flowline.FluxBasedModel
    :param model2: oggm.flowline.FluxBasedModel
    :return:       float
    """

    model1 = model1
    model2 = model2
    model2.run_until(0)
    model1.run_until(ye)

    fls1 = model1.fls
    fls2 = model2.fls

    fitness = 0
    m = 0
    for i in range(len(model1.fls)):
        fitness = fitness + np.sum(
            abs(fls1[i].surface_h - fls2[i].surface_h)**2) + \
                    np.sum(abs(fls1[i].widths - fls2[i].widths)**2)
        m = m + fls1[i].nx
    fitness = fitness / m

    return fitness


def _run_experiment(gdir, temp_bias, bias, ys,ye):
    """
    Creates the synthetic experiment for one glacier. model_run_experiment.nc
    will be saved in working directory.
    """

    # check, if this experiment already exists
    try:
        rp = gdir.get_filepath('model_run', filesuffix='_advanced_experiment_' +str(temp_bias)+'_'+ str(bias))
        model = FileModel(rp)

    # otherwise create experiment
    except:

        fls = gdir.read_pickle('model_flowlines')
        try:
            model = tasks.run_random_climate(gdir, nyears=400, y0=ys, bias=bias, seed=1,
                                             temperature_bias=temp_bias,
                                             init_model_fls=fls)

            # construct observed glacier, previous glacier will be run forward from
            # 1917 - 2000 with past climate file

            fls = deepcopy(model.fls)
            model = tasks.run_from_climate_data(gdir, ys=ys, ye=ye, init_model_fls=fls,bias=bias,
                                        output_filesuffix='_advanced_experiment_'+str(temp_bias)+'_'+str(bias))
        except:
            pass
    return model


def find_residual(gdir, temp_bias_list, ys,a=-2000,b=2000):

    best_df = pd.DataFrame()

    fls = gdir.read_pickle('model_flowlines')
    mod = FluxBasedModel(flowlines=fls)

    for temp_bias in temp_bias_list:
        print(temp_bias)
        try:
            ye = gdir.rgi_date
            max_it = 15
            i = 0
            bounds = [a,b]

            df = pd.DataFrame()

            while i < max_it:
                bias = round((bounds[0] + bounds[1]) / 2,1)
                ex_mod2 = _run_experiment(gdir, temp_bias, bias, ys, ye)
                fit = fitness_function(ye,ex_mod2,mod)
                df = df.append(pd.Series({'bias':bias,'fitness':fit}),ignore_index=True)
                if  (abs(mod.area_km2-ex_mod2.area_km2)<1e-4 and fit<125) or bounds[1]-bounds[0]<=1:
                    break

                elif ex_mod2.area_km2 > mod.area_km2:
                    bounds[0] = bias
                else:
                    bounds[1] = bias
                i +=1

            # best bias found
            bias = df.iloc[df.fitness.idxmin()].bias
            rp = gdir.get_filepath('model_run', filesuffix='_advanced_experiment_'+str(temp_bias)+'_'+str(bias))
            model = FileModel(rp)
            model.run_until(ye)

            rp = gdir.get_filepath('model_run', filesuffix='_advanced_experiment_'+str(temp_bias)+'_'+str(0.0))
            ex_mod = FileModel(rp)
            ex_mod.run_until(ye)
            '''
            plt.figure(figsize=(15,10))
            plt.plot(model.fls[-1].surface_h,'r',label='best')
            plt.plot(mod.fls[-1].surface_h, 'orange', label='original')
            plt.plot(ex_mod.fls[-1].surface_h, 'r:', label='old experiment')
            plt.plot(model.fls[-1].bed_h,'k', label='bed')
            plt.legend()
            utils.mkdir(os.path.join(cfg.PATHS['plot_dir'],'bias_test'))
            plt.savefig(os.path.join(cfg.PATHS['plot_dir'],'bias_test',gdir.rgi_id+'.png'),dpi=200)
            '''
            diff = mod.area_km2 - model.area_km2_ts()[gdir.rgi_date]
            model.reset_y0(ys)

            series = pd.Series({'rgi_id':gdir.rgi_id,'bias':bias,'iterations':i, 'fitness':df.fitness.min(), 'area_diff':diff, 'model':model, 'temp_bias':temp_bias})
        except:
            series =  pd.Series({'rgi_id':gdir.rgi_id, 'temp_bias':temp_bias})
        best_df = best_df.append(series, ignore_index=True)

    plt.figure()

    x = np.arange(mod.fls[-1].nx) * mod.fls[-1].dx * mod.fls[-1].map_dx

    for temp_bias, model in zip(best_df.temp_bias, best_df.model):
        model.run_until(model.length_m_ts().index[-1])
        plt.plot(x,model.fls[-1].surface_h, label=str(temp_bias))
        #model.volume_km3_ts().plot()
    plt.plot(x,mod.fls[-1].surface_h,'r:')
    plt.plot(x,mod.fls[-1].bed_h, 'k')

    plt.legend(title='temp_bias')
    plt.xlabel('Distance along the main flowline (m)')
    plt.ylabel('Altitude (m)')
    plt.title(gdir.rgi_id)
    utils.mkdir(os.path.join(cfg.PATHS['plot_dir'], 'bias_test'))
    plt.savefig(os.path.join(cfg.PATHS['plot_dir'], 'bias_test', gdir.rgi_id + '.png'), dpi=200)
    plt.show()

    return best_df


def advanced_experiments(gdirs, temp_bias_list ,ys , region):

    exp_df = pd.DataFrame()

    pool = Pool()
    list = pool.map(partial(find_residual,temp_bias_list=temp_bias_list,ys=ys),gdirs)
    pool.close()
    pool.join()

    exp_df = exp_df.append(list, ignore_index=True)
    p = os.path.join(cfg.PATHS['working_dir'], str(region)+'_advanced_experiments.pkl')
    exp_df.to_pickle(p, compression='gzip')
    return p

if __name__ == '__main__':

    cfg.initialize()

    ON_CLUSTER = False
    REGION = '11'

    # Local paths
    if ON_CLUSTER:
        WORKING_DIR = os.environ.get("S_WORKDIR")
        cfg.PATHS['working_dir'] = WORKING_DIR
        REGION = str(os.environ.get('I')).zfill(2)
    else:
        WORKING_DIR = '/home/juliaeis/Dokumente/OGGM/work_dir/reconstruction/leclercq_experiments'
        cfg.PATHS['working_dir'] = WORKING_DIR
        utils.mkdir(WORKING_DIR, reset=False)

    cfg.PATHS['plot_dir'] = os.path.join(cfg.PATHS['working_dir'], 'plots')
    utils.mkdir(cfg.PATHS['plot_dir'], reset=False)

    # Use multiprocessing?
    cfg.PARAMS['use_multiprocessing'] = True

    # How many grid points around the glacier?
    cfg.PARAMS['border'] = 200

    # Set to True for operational runs
    cfg.PARAMS['continue_on_error'] = True

    cfg.PARAMS['run_mb_calibration'] = False
    cfg.PARAMS['optimize_inversion_params'] = False
    cfg.PARAMS['dl_verify'] = False

    '''
    # Use HISTALP climate file
    cfg.PARAMS['baseline_climate'] = 'HISTALP'
    cfg.PARAMS['prcp_scaling_factor'] = 1.75
    cfg.PARAMS['temp_all_liq'] = 2.0
    cfg.PARAMS['temp_default_gradient'] = -0.0065
    cfg.PARAMS['temp_melt'] = -1.75
    cfg.PARAMS['temp_all_solid'] = 0.0
    '''

    # add to BASENAMES
    _doc = 'contains observed and searched glacier from synthetic experiment to find intial state'
    cfg.BASENAMES['synthetic_experiment'] = ('synthetic_experiment.pkl', _doc)

    diff = pd.DataFrame()

    # read leclercq links
    lec = pd.read_csv('rgi_leclercq_links_2014_RGIV6.csv')
    lec['REGION'] = lec.RGI_ID.apply(lambda x: x.split('-')[-1].split('.')[0])

    # We use intersects
    db = utils.get_rgi_intersects_region_file(version='61', region=REGION)
    cfg.set_intersects_db(db)

    # RGI file
    path = utils.get_rgi_region_file(REGION, version='61')
    rgidf = gpd.read_file(path)

    # only the ones with leclercq observation
    rgidf = rgidf[rgidf.RGIId.isin(lec[lec.REGION==REGION].RGI_ID.values)].head(5)
    #rgidf = rgidf[rgidf.RGIId.isin(['RGI60-05.02112'])]

    # exclude non-landterminating glaciers
    rgidf = rgidf[rgidf.TermType==0]
    rgidf = rgidf[rgidf.Connect !=2]
    # sort for efficient using
    rgidf = rgidf.sort_values('Area', ascending=True)

    gdirs = workflow.init_glacier_regions(rgidf)

    #preprocessing(gdirs)
    temp_bias_list = [0, -0.25, -0.5, -0.75, -1, -1.25, -1.5, -1.75, -2]
    temp_bias_list = [0,-1]
    #p = advanced_experiments(gdirs, temp_bias_list, 1917,REGION)
    '''
    best_df = pd.DataFrame()
    for REGION in [str(x).zfill(2) for x in range(0,19)]:
        try:
            p = os.path.join(cfg.PATHS['working_dir'],
                             REGION + '_advanced_experiments.pkl')
            df = pd.read_pickle(p, compression='gzip')

            for rgi in df.rgi_id.unique():
                best = df.iloc[df[df.rgi_id==rgi].fitness.idxmin()]
                best_df = best_df.append(pd.Series({'rgi_id':best.rgi_id, 'temp_bias':best.temp_bias, 'bias':best.bias, 'fitness':best.fitness/125}),ignore_index=True)
        except:
            pass

    best_df = best_df.set_index('rgi_id')
    best_df.to_csv(os.path.join(cfg.PATHS['working_dir'],'best_experiment.csv'))
    '''

    df = pd.read_csv(os.path.join(cfg.PATHS['working_dir'],'best_experiment.csv'))
    df = df.set_index('rgi_id')
    for gdir in gdirs:
        bias = df.loc[gdir.rgi_id].bias
        temp_bias = df.loc[gdir.rgi_id].temp_bias
        print(bias, temp_bias)
        ex_mod = _run_experiment(gdir, temp_bias, bias, 1917,gdir.rgi_date)
        plot_experiment(gdir, ex_mod, 1917, cfg.PATHS['plot_dir'])
        print(ex_mod)
