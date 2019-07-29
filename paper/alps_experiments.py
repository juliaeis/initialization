
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
import matplotlib as mpl

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

def read_results(gdirs):

    model_df = pd.DataFrame()
    pool = mp.Pool()
    list = pool.map(read_result_parallel, gdirs)
    pool.close()
    pool.join()
    model_df = model_df.append(list, ignore_index=True)
    return model_df


def plot_ratio_volume(df,ex_mod,gdir,ratio1, ratio2, ratio3):

    fig, ax = plt.subplots(1,1, figsize=(10,7))
    norm = mpl.colors.LogNorm(vmin=0.01 / 125, vmax=10)
    cmap = mpl.cm.get_cmap('viridis')

    # df = df.sort_values('objective', ascending=False)
    df = df.sort_values('fitness', ascending=False)

    for i, model in df['model'].iteritems():
        try:

            model = deepcopy(model)
            color = cmap(norm(df.loc[i, 'fitness']))
            model.volume_km3_ts().plot(ax=ax, color=[color], label='')
        except:
            pass
    ex_mod.volume_km3_ts().plot(ax=ax, color='red', linestyle=':',
                                linewidth=3,
                                label='')
    plt.title(gdir.rgi_id)
    plt.text(1950,0.95*df.volume.max(), 'ratio1: '+str(ratio1.round(4)))
    plt.text(1950, 0.85 * df.volume.max(), 'ratio2: ' + str(ratio2.round(4)))
    plt.text(1950, 0.75 * df.volume.max(), 'ratio3: ' + str(ratio3.round(4)))
    plt.ylabel('Volume (km$^3$)')
    plt.ylabel('Volume (km$^3$)')
    plt.xlabel('Time (years)')
    utils.mkdir(os.path.join(cfg.PATHS['plot_dir'],'ratio'))
    plt.savefig(os.path.join(cfg.PATHS['plot_dir'],'ratio',gdir.rgi_id+'.png'))
    #plt.show()

def read_result_parallel(gdir):


    if os.path.isfile(os.path.join(gdir.dir, 'model_run_experiment.nc')):
        print(gdir.rgi_id)
        start = time.time()

        rp = gdir.get_filepath('model_run', filesuffix='_experiment')
        ex_mod = FileModel(rp)

        if ex_mod.area_km2_ts()[2000] > 0.01:

            df = pd.read_pickle(os.path.join(gdir.dir, 'result1850.pkl'),
                                compression='gzip')
            df.fitness = df.fitness / 125
            # replace all values smaller than 1e-4
            df.fitness[df.fitness < 1e-4] = 1e-4

            m = df.fitness.min()
            median = df.fitness.median()
            q5 = df.fitness.quantile(0.1)

            ratio1 = np.exp((1-(m / q5)) + (1-(q5 / median)))/np.exp(2)
            ratio2 = 0.5*(1-(m / q5)) + 0.5*(1-(q5 / median))
            ratio3 = 1-(m/q5)

            plot_ratio_volume(df,ex_mod,gdir,ratio1, ratio2, ratio3)



    return pd.Series({'rgi':gdir.rgi_id, 'min':m,'median':median,'q5':q5,
                      'min/q5':m/q5,'q5/median':q5/median,
                      'area':gdir.rgi_area_km2,'ratio1':ratio1,'ratio2':ratio2,
                      'ratio3': ratio3})
    #return pd.Series({'rgi_id':gdir.rgi_id, 'v_min':v_min,'v_mac':v_max, 'range':v_max-v_min})
    #except:
    #   pass
if __name__ == '__main__':
    cfg.initialize()

    ON_CLUSTER = False

    # Local paths
    if ON_CLUSTER:
        WORKING_DIR = os.environ.get("S_WORKDIR")
        WORKING_DIR = '/home/users/julia/initialization/out/paper_correction/paper_600'
        cfg.PATHS['working_dir'] = WORKING_DIR
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
    #rgidf = rgidf[rgidf.RGIId.isin(['RGI60-11.00897','RGI60-11.00779', 'RGI60-11.00029', 'RGI60-11.00036', 'RGI60-11.00001','RGI60-11.00026'])]



    # sort for efficient using
    rgidf = rgidf.sort_values('Area', ascending=False)


    gdirs = workflow.init_glacier_regions(rgidf)
    df = read_results(gdirs)
    print(df)
    '''

        print(gdir.rgi_id)
        print(len(df[df.fitness<1e-4]))
        print(len(df[(df.fitness > 1e-4) & (df.fitness< 1e-3)]))
        print(len(df[(df.fitness > 1e-3) & (df.fitness < 1e-2)]))
        print(len(df[(df.fitness > 1e-2) & (df.fitness < 1e-1)]))
        print(len(df[(df.fitness > 1e-1) & (df.fitness < 1e-0)]))
        print(len(df[(df.fitness > 1e0) & (df.fitness < 1e1)]))
        print(len(df[df.fitness > 1e1]))
    plt.show()

    df = pd.read_pickle(os.path.join(cfg.PATHS['working_dir'], 'ratio.pkl'),
                 compression='gzip')
    ratio = df.dropna().reconstructability
    ratio.plot.hist(bins=30)
    plt.yscale('log')
    plt.show()

    df = read_results(gdirs)
    if ON_CLUSTER:
        df.to_pickle(os.path.join(os.environ.get("S_WORKDIR"),'ratio.pkl'),compression='gzip')
    else:
        df.to_pickle(os.path.join(cfg.PATHS['working_dir'],'ratio.pkl'),compression='gzip')
    print(df)




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
    '''