
import os
import sys
import pandas as pd
import geopandas as gpd
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from oggm import cfg, utils, workflow, tasks
from oggm.utils._downloads import get_demo_file
from oggm.core.flowline import FluxBasedModel
import oggm
from leclercq.leclercq_plots import *
sys.path.append('../')
sys.path.append('../../')
from  initialization.core import *




def fitness_function(model1, model2):
    """
    calculates the objective value (difference in geometry)
    :param model1: oggm.flowline.FluxBasedModel
    :param model2: oggm.flowline.FluxBasedModel
    :return:       float
    """

    model1 = model1
    model2 = model2
    model2.run_until(0)
    model1.run_until(2016)

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


def _run_experiment(gdir, bias):
    """
    Creates the synthetic experiment for one glacier. model_run_experiment.nc
    will be saved in working directory.
    """

    # check, if this experiment already exists
    try:
        rp = gdir.get_filepath('model_run', filesuffix='_advanced_experiment_' + str(bias))
        model = FileModel(rp)

    # otherwise create experiment
    except:

        fls = gdir.read_pickle('model_flowlines')
        try:
            model = tasks.run_random_climate(gdir, nyears=400, y0=1917, bias=bias, seed=1,
                                             temperature_bias=-0.5,
                                             init_model_fls=fls)

            # construct observed glacier, previous glacier will be run forward from
            # 1917 - 2000 with past climate file

            fls = copy.deepcopy(model.fls)
            model = tasks.run_from_climate_data(gdir, ys=1917, ye=2016, init_model_fls=fls,bias=bias,
                                        output_filesuffix='_advanced_experiment_'+str(bias))
        except:
            pass
    return model


def find_residual(gdir, a=-2000,b=2000):

    try:

        max_it = 15
        i = 0
        bounds = [a,b]

        df = pd.DataFrame()

        fls = gdir.read_pickle('model_flowlines')
        mod = FluxBasedModel(flowlines=fls)

        while i < max_it:
            bias = round((bounds[0] + bounds[1]) / 2,1)
            ex_mod2 = _run_experiment(gdir, bias)
            fit = fitness_function(ex_mod2,mod)
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
        rp = gdir.get_filepath('model_run', filesuffix='_advanced_experiment_' + str(bias))
        model = FileModel(rp)
        model.run_until(2016)

        rp = gdir.get_filepath('model_run', filesuffix='_advanced_experiment_' + str(0.0))
        ex_mod = FileModel(rp)
        ex_mod.run_until(2016)


        plt.figure(figsize=(15,10))
        plt.plot(model.fls[-1].surface_h,'r',label='best')
        plt.plot(mod.fls[-1].surface_h, 'orange', label='original')
        plt.plot(ex_mod.fls[-1].surface_h, 'r:', label='old experiment')
        plt.plot(model.fls[-1].bed_h,'k', label='bed')
        plt.legend()
        utils.mkdir(os.path.join(cfg.PATHS['plot_dir'],'bias_test'))
        plt.savefig(os.path.join(cfg.PATHS['plot_dir'],'bias_test',gdir.rgi_id+'.png'),dpi=200)


        diff = mod.area_km2 - model.area_km2_ts()[2016]
        model.reset_y0(1917)
        print(gdir.rgi_id, i, bias, df.fitness.min())

        series = pd.Series({'rgi_id':gdir.rgi_id,'bias':bias,'iterations':i, 'fitness':df.fitness.min(), 'area_diff':diff, 'model':model})
    except:
        series =  pd.Series({'rgi_id':gdir.rgi_id})

    return series


def get_ref_length_data(gdir):
    """Get the glacier lenght data from P. Leclercq's data base.

     https://folk.uio.no/paulwl/data.php

     For some glaciers only!
     """

    df = pd.read_csv('rgi_leclercq_links_2014_RGIV6.csv')
    df = df.loc[df.RGI_ID == gdir.rgi_id]
    if len(df) == 0:
        raise RuntimeError('No length data found for this glacier!')
    ide = df.LID.values[0]

    f = get_demo_file('Glacier_Lengths_Leclercq.nc')
    with xr.open_dataset(f) as dsg:
        # The database is not sorted by ID. Don't ask me...
        grp_id = np.argwhere(dsg['index'].values == ide)[0][0] + 1
    with xr.open_dataset(f, group=str(grp_id)) as ds:
        df = ds.to_dataframe()
        df.name = ds.glacier_name
    return df


def advanced_experiments(gdirs):
    exp_df = pd.DataFrame()

    pool = Pool()
    list = pool.map(find_residual, gdirs)
    pool.close()
    pool.join()

    exp_df = exp_df.append(list, ignore_index=True)
    exp_df.to_pickle(
        os.path.join(cfg.PATHS['working_dir'], '11_advanced_experiments.pkl'))


if __name__ == '__main__':

    cfg.initialize()

    ON_CLUSTER = False

    # Local paths
    if ON_CLUSTER:
        WORKING_DIR = os.environ.get("S_WORKDIR")
        cfg.PATHS['working_dir'] = WORKING_DIR
    else:
        WORKING_DIR = '/home/juliaeis/Dokumente/OGGM/work_dir/reconstruction/leclercq'
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

    # We use intersects
    db = utils.get_rgi_intersects_region_file(version='61', region='11')
    cfg.set_intersects_db(db)

    cfg.PARAMS['run_mb_calibration'] = False
    cfg.PARAMS['optimize_inversion_params'] = False

    # RGI file
    path = utils.get_rgi_region_file('11', version='61')
    rgidf = gpd.read_file(path)

    # read leclercq links
    lec = pd.read_csv('rgi_leclercq_links_2014_RGIV6.csv')
    lec['REGION'] = lec.RGI_ID.apply(lambda x: x.split('-')[-1].split('.')[0])

    # only the ones with leclercq observation
    rgidf = rgidf[rgidf.RGIId.isin(lec[lec.REGION=='11'].RGI_ID.values)]
    #rgidf = rgidf[rgidf.RGIId.isin(['RGI60-11.03229'])]
    # sort for efficient using
    rgidf = rgidf.sort_values('Area', ascending=True)
    gdirs = workflow.init_glacier_regions(rgidf)

    preprocessing(gdirs)
    advanced_experiments(gdirs)

    if ON_CLUSTER:
        df = pd.read_pickle('11_advanced_experiments.pkl')
    else:
        df = pd.read_pickle('11_advanced_experiments.pkl')
    df = df.set_index('rgi_id')
    df.fitness = df.fitness/125


    for gdir in gdirs:
        if df.loc[gdir.rgi_id].fitness<1:
            bias = df.loc[gdir.rgi_id].bias
            ex_mod = df.loc[gdir.rgi_id].model

            ini_df = find_possible_glaciers(gdir, 1917, 2016, 200, ex_mod, bias)
            ini_df.fitness = ini_df.fitness/125

            lec = get_ref_length_data(gdir).dL
            lec = (lec -lec.iloc[-1]) + ex_mod.length_m_ts()[lec.index[-1]]

            plot_fitness_values(gdir, df, ex_mod, 1917, 2016, lec, cfg.PATHS['plot_dir'])
            plot_median(gdir, df, 125, ex_mod, 1917, 2016,lec, cfg.PATHS['plot_dir'])


