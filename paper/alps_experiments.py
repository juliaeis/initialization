
import os
import sys
from copy import deepcopy
sys.path.append('../')
from initialization.core import *
from plots_paper import *

import matplotlib.pyplot as plt
import salem
import geopandas as gpd
from oggm import cfg, workflow, utils
from oggm.core.flowline import FluxBasedModel
pd.options.mode.chained_assignment = None
import time

def find_median(df):


    try:
        accept_df = df[df.fitness <= 1]

        quant_df = accept_df[accept_df.fitness <= accept_df.fitness.quantile(0.05)]

        # median state
        quant_df.loc[:, 'length'] = quant_df.model.apply(lambda x: x.length_m)

        quant_df = quant_df.sort_values('length', ascending=False)

        l = len(quant_df)

        if l % 2:
            index = int((l - 1) / 2)
        else:
            index = int(l / 2)

        return deepcopy(quant_df.iloc[index].model), quant_df.at[quant_df.length.idxmin(),'model'], quant_df.at[quant_df.length.idxmax(),'model']

    except:

        return deepcopy(df.iloc[df.fitness.idxmin()].model), None, None


if __name__ == '__main__':
    cfg.initialize()

    ON_CLUSTER = False

    # Local paths
    if ON_CLUSTER:
        WORKING_DIR = os.environ.get("S_WORKDIR")
        WORKING_DIR = '/home/users/julia/initialization/out/paper_correction/paper_600'
        cfg.PATHS['working_dir'] = WORKING_DIR
        job_nr = int(os.environ.get('I'))
    else:
        WORKING_DIR = '/home/juliaeis/Dokumente/OGGM/work_dir/reconstruction/600_paper_correction/'
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
    rgidf = rgidf[rgidf.RGIId == 'RGI60-11.00897']

    # sort for efficient using
    rgidf = rgidf.sort_values('Area', ascending=False)

    if ON_CLUSTER:
        rgidf = rgidf[job_nr:len(rgidf):80]

    gdirs = workflow.init_glacier_regions(rgidf)

    #preprocessing(gdirs)

    # experiments.py
    #synthetic_experiments_parallel(gdirs[:1])

    t_0 = 1850
    t_e = 2000
    epsilon = 125

    model_df = pd.DataFrame()
    time_df = pd.DataFrame()


    for gdir in gdirs:

        if os.path.isfile(os.path.join(gdir.dir, 'model_run_experiment.nc')):
            start = time.time()


            rp = gdir.get_filepath('model_run', filesuffix='_experiment')
            ex_mod = FileModel(rp)

            if ex_mod.area_km2_ts()[2000] > 0.01:

                    df = find_possible_glaciers(gdir, t_0, t_e, 200)
                    df['fitness'] = df.fitness/125
                    med_mod, perc_min, perc_max = find_median(df)
                    min_mod = deepcopy(df.loc[df.fitness.idxmin(), 'model'])

                    print(med_mod.length_m)
                    print(min_mod.length_m)
            print(time.time()-start)

    '''

                        med_mod, perc_min, perc_max = find_median(df)
                        min_mod = deepcopy(df.loc[df.fitness.idxmin(), 'model'])

                        df['fitness2'] = df.model.apply(lambda x: abs(x.area_km2_ts()[2000] - ex_mod.area_km2_ts()[2000]) ** 2)
                        df['fitness3'] = df.model.apply(lambda x: abs(x.length_m_ts()[2000] - ex_mod.length_m_ts()[2000]) ** 2)

                        # saves median state, minimum state and experiment model
                        model_df.loc[gdir.rgi_id, 'median'] = deepcopy(med_mod)
                        model_df.loc[gdir.rgi_id, 'minimum'] = deepcopy(min_mod)
                        model_df.loc[gdir.rgi_id, 'experiment'] = deepcopy(ex_mod)
                        model_df.loc[gdir.rgi_id, 'flowline'] = FluxBasedModel(flowlines=gdir.read_pickle('model_flowlines'))
                        model_df.loc[gdir.rgi_id, 'perc_min'] = deepcopy(perc_min)
                        model_df.loc[gdir.rgi_id, 'perc_max'] = deepcopy(perc_max)
                        model_df.loc[gdir.rgi_id, 'fit2'] = deepcopy(df.loc[df.fitness2.idxmin(), 'model'])
                        model_df.loc[gdir.rgi_id, 'fit3'] = deepcopy(df.loc[df.fitness3.idxmin(), 'model'])


                        # time_df
                        time_df.loc[gdir.rgi_id, 'time'] = time.time()-start

                        try:
                            # plots
                            plot_experiment(gdir, ex_mod, t_0, cfg.PATHS['plot_dir'])
                            plot_candidates(gdir, df, t_0, 'step3', cfg.PATHS['plot_dir'])
                            plot_fitness_values(gdir, df, ex_mod, t_0, cfg.PATHS['plot_dir'])
                            plot_median(gdir, df, epsilon, ex_mod, t_0, t_e, cfg.PATHS['plot_dir'])

                        except:
                            pass

            except:
                print(gdir.rgi_id+' failed')

    if ON_CLUSTER:
        model_df.to_pickle(os.path.join(cfg.PATHS['working_dir'], 'models_'+str(job_nr)+'.pkl'), compression='gzip')
        time_df.to_pickle(os.path.join(cfg.PATHS['working_dir'], 'time_'+str(job_nr)+'.pkl'), compression='gzip')

    model_df = pd.read_pickle(os.path.join(cfg.PATHS['working_dir'], 'models.pkl'), compression='gzip')
    time_df = pd.read_pickle(os.path.join(cfg.PATHS['working_dir'], 'time.pkl'), compression='gzip')

    median = model_df['median'].apply(lambda x: x.volume_km3_ts())
    minimum = model_df['minimum'].apply(lambda x: x.volume_km3_ts())
    experiment = model_df['experiment'].apply(lambda x: x.volume_km3_ts())

    error1 = median-experiment
    error2 = minimum-experiment

    error3 = (median-experiment)/experiment
    error4 = (minimum-experiment)/experiment

    error5 = np.log(median/experiment)
    error6 = np.log(minimum/experiment)

    plot_relative_error(error1, error2, 'abs', cfg.PATHS['plot_dir'], all=True)
    plot_relative_error(error3, error4, 'rel', cfg.PATHS['plot_dir'], all=True)
    plot_relative_error(error5, error6, 'log', cfg.PATHS['plot_dir'], all=True)

    '''