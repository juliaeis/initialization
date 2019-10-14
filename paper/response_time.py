import os
import sys
import pandas as pd
import geopandas as gpd
import numpy as np
import matplotlib.pyplot as plt
from oggm import cfg, workflow, utils


sys.path.append('../')
from initialization.core import *


def _response(gdir, model_df):

    fig,ax = plt.subplots(1,1)
    ex_mod2 = model_df.loc[gdir.rgi_id].experiment
    ax.plot(0,copy.deepcopy(ex_mod2.volume_km3),'o', label=r'experiment$_{1850}$')
    tasks.run_constant_climate(gdir, nyears=600, y0=1850, halfsize=15,
                               temperature_bias=-1,
                               store_monthly_step=False,
                               output_filesuffix='experiment_equilibrium',
                               init_model_fls=copy.deepcopy(ex_mod2.fls))
    rp = gdir.get_filepath('model_run', filesuffix='experiment_equilibrium')
    ex_mod = FileModel(rp)
    ex_mod.run_until(600)

    try:
        # scenario a
        tasks.run_constant_climate(gdir, nyears=600, y0=1865, halfsize=15,
                                   temperature_bias=-1.05,
                                   store_monthly_step=False,
                                   output_filesuffix='_-1.05',
                                   init_model_fls=copy.deepcopy(ex_mod.fls))
    except:
        pass

    try:

        # senario b
        tasks.run_constant_climate(gdir, nyears=600, y0=1865, halfsize=15,
                                   store_monthly_step=False, temperature_bias=-1.1,
                                   output_filesuffix='_-1.1',
                                   init_model_fls=copy.deepcopy(ex_mod.fls))

    except:
        pass

    try:

        # senario c
        tasks.run_constant_climate(gdir, nyears=600, y0=1865, halfsize=15,
                                   temperature_bias=-1,
                                   store_monthly_step=False,
                                   output_filesuffix='_-1',
                                   init_model_fls=copy.deepcopy(ex_mod.fls))
    except:
        pass

    try:
        # senario d
        tasks.run_constant_climate(gdir, nyears=600, y0=1865, halfsize=15,
                                   store_monthly_step=False, temperature_bias=-0.9,
                                   output_filesuffix='_-0.9',
                                   init_model_fls=copy.deepcopy(ex_mod.fls))
    except:
        pass

    try:
        # senario_e
        tasks.run_constant_climate(gdir, nyears=600, y0=1865, halfsize=15,temperature_bias=-0.95,
                                   store_monthly_step=False,
                                   output_filesuffix='_-0.95',
                                   init_model_fls=copy.deepcopy(ex_mod.fls))
    except:
        pass

    response = pd.DataFrame()

    colors = ['C0','C1','C2','C3','C4','C5',]
    for i, s in enumerate(['_-1.1', '_-1.05', '_-1', '_-0.95', '_-0.9']):

        #try:
        rp = gdir.get_filepath('model_run', filesuffix=s)
        mod2 = FileModel(rp)
        if mod2.volume_km3_ts()[600] > 0:

            diff = mod2.volume_km3_ts()[600] - ((
                                                mod2.volume_km3_ts()[600] -
                                                ex_mod.volume_km3_ts()[
                                                    600]) / np.exp(1))

            t = abs(mod2.volume_km3_ts() - diff).idxmin()
            response.at[gdir.rgi_id, s.split('_')[-1]] = t
            mod2.volume_km3_ts().plot(ax=ax, color=colors[i], label=s.split('_')[-1])
            ax.axvline(t,color=colors[i],linestyle=':')
        #except:
        #    pass
    plt.legend(loc='best')
    p = os.path.join(cfg.PATHS['plot_dir'],'response_time')
    utils.mkdir(p)
    plt.savefig(os.path.join(p,str(gdir.rgi_id)+'.png'))
    return response


def response_time(gdirs, model_df, job_nr):

    response_df = pd.DataFrame()

    pool = Pool()
    list = pool.map(partial(_response, model_df=model_df), gdirs)
    pool.close()
    pool.join()

    response_df = response_df.append(list, ignore_index=False)
    p = os.path.join(cfg.PATHS['working_dir'], 'response_time_'+str(job_nr)+'.pkl')
    response_df.to_pickle(p)
    return p

if __name__ == '__main__':
    cfg.initialize()

    ON_CLUSTER = False

    # Local paths
    if ON_CLUSTER:
        WORKING_DIR = os.environ.get("S_WORKDIR")
        cfg.PATHS['working_dir'] = WORKING_DIR
        job_nr = int(os.environ.get('I'))
    else:
        WORKING_DIR = '/home/juliaeis/Dokumente/OGGM/work_dir/reconstruction/600_paper_correction/'
        cfg.PATHS['working_dir'] = WORKING_DIR
        utils.mkdir(WORKING_DIR, reset=False)
        job_nr=0

    cfg.PATHS['plot_dir'] = os.path.join(cfg.PATHS['working_dir'], 'plots')
    utils.mkdir(cfg.PATHS['plot_dir'], reset=False)

    # Use multiprocessing?
    cfg.PARAMS['use_multiprocessing'] = True

    # How many grid points around the glacier?
    cfg.PARAMS['border'] = 200

    # Set to True for operational runs
    cfg.PARAMS['continue_on_error'] = True

    # Use HISTALP climate file
    cfg.PARAMS['baseline_climate'] = 'HISTALP'
    cfg.PARAMS['prcp_scaling_factor'] = 1.75
    cfg.PARAMS['temp_all_liq'] = 2.0
    cfg.PARAMS['temp_default_gradient'] = -0.0065
    cfg.PARAMS['temp_melt'] = -1.75
    cfg.PARAMS['temp_all_solid'] = 0.0

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
    if ON_CLUSTER:
         p = os.path.join('/home/users/julia/initialization/out/paper_correction/paper_600', 'models_' + str(job_nr) + '.pkl')
         model_df = pd.read_pickle(p,compression='gzip')
    else:
        p = os.path.join(cfg.PATHS['working_dir'], 'models_0.pkl')
        model_df = pd.read_pickle(p, compression='gzip')

    rgidf = rgidf[rgidf.RGIId.isin(model_df.index)].sort_values(by='Area', ascending=False)

    gdirs = workflow.init_glacier_regions(rgidf)

    preprocessing(gdirs)
    p = response_time(gdirs, model_df, job_nr)
    print(pd.read_pickle(p))


