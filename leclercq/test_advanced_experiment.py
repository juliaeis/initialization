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
import advanced_experiments
from leclercq.leclercq_plots import *
import copy

mpl.rcParams['axes.linewidth'] = 3
mpl.rcParams['xtick.major.width'] = 3
mpl.rcParams['ytick.major.width'] = 3
mpl.rcParams['font.size'] =25
mpl.rcParams['font.weight'] = 'medium'
mpl.rcParams['axes.labelweight'] = 'medium'
mpl.rcParams['legend.fontsize'] = 20 #30
mpl.rcParams['lines.linewidth'] = 3


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
            model = tasks.run_random_climate(gdir, nyears=600, y0=ys, bias=bias, seed=1,
                                             temperature_bias=temp_bias,
                                             init_model_fls=fls,output_filesuffix='_random_experiment_'+str(temp_bias)+'_'+str(bias) )

            # construct observed glacier, previous glacier will be run forward from
            # 1917 - rgi_date with past climate file

            fls = deepcopy(model.fls)
            tasks.run_from_climate_data(gdir, ys=ys, ye=ye, init_model_fls=fls,bias=bias,
                                        output_filesuffix='_advanced_experiment_'+str(temp_bias)+'_'+str(bias))
            # to return FileModel
            rp = gdir.get_filepath('model_run',filesuffix='_advanced_experiment_' + str(
                                       temp_bias) + '_' + str(bias))
            model = FileModel(rp)


        except:
            pass

    return model


def find_residual(gdir, temp_bias_list, ys,a=-2000,b=2000):

    best_df = pd.DataFrame()

    fls = gdir.read_pickle('model_flowlines')
    mod = FluxBasedModel(flowlines=fls)

    for temp_bias in temp_bias_list:

        try:
            ye = gdir.rgi_date
            max_it = 15
            i = 0
            bounds = [a,b]

            df = pd.DataFrame()

            while i < max_it:
                bias = round((bounds[0] + bounds[1]) / 2,1)

                ex_mod2 = _run_experiment(gdir, temp_bias, bias, ys, ye)

                diff = mod.area_km2 - ex_mod2.area_km2_ts()[ye]

                df = df.append(pd.Series({'bias':bias,'area_diff':diff}),ignore_index=True)

                if  (abs(diff)<1e-4) or bounds[1]-bounds[0]<=1:
                    break

                elif ex_mod2.area_km2_ts()[ye] > mod.area_km2:
                    bounds[0] = bias
                else:
                    bounds[1] = bias
                i +=1

            # best bias found
            bias = df.iloc[df.area_diff.abs().idxmin()].bias

            rp = gdir.get_filepath('model_run', filesuffix='_advanced_experiment_'+str(temp_bias)+'_'+str(bias))
            model = FileModel(rp)

            diff = gdir.rgi_area_km2 - model.area_km2_ts()[gdir.rgi_date]

            series = pd.Series({'rgi_id':gdir.rgi_id,'bias':bias,'iterations':i,  'area_diff':diff, 'model':model, 'temp_bias':temp_bias})

        except:
            series =  pd.Series({'rgi_id':gdir.rgi_id, 'temp_bias':temp_bias})
        best_df = best_df.append(series, ignore_index=True)


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
        WORKING_DIR =  WORKING_DIR = '/home/juliaeis/Dokumente/OGGM/work_dir/reconstruction/leclercq_diff/temp_0'
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

    # add to BASENAMES
    _doc = 'contains observed and searched glacier from synthetic experiment to find intial state'
    cfg.BASENAMES['synthetic_experiment'] = ('synthetic_experiment.pkl', _doc)


    # We use intersects
    db = utils.get_rgi_intersects_region_file(version='61', region=REGION)
    cfg.set_intersects_db(db)

    # RGI file
    path = utils.get_rgi_region_file(REGION, version='61')
    rgidf = gpd.read_file(path)

    # only the ones with leclercq observation
    rgidf = rgidf[rgidf.RGIId.isin(['RGI60-11.00002'])]


    gdir = workflow.init_glacier_regions(rgidf)

    p =advanced_experiments(gdir, [0] ,1917 , REGION)
    print(p)
    plot_advanced_experiment(gdir[0])
    plt.show()