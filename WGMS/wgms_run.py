
import sys
sys.path.append('../')
from initialization.core import *
from paper.plots_paper import *


import geopandas as gpd
from oggm import cfg, utils
pd.options.mode.chained_assignment = None


if __name__ == '__main__':
    cfg.initialize()

    ON_CLUSTER = True

    # Local paths
    if ON_CLUSTER:
        WORKING_DIR = os.environ.get("S_WORKDIR")
        cfg.PATHS['working_dir'] = WORKING_DIR
        OUT_DIR = os.environ.get("OUTDIR")
        REGION = str(os.environ.get('REGION')).zfill(2)
        JOB_NR =  int(os.environ.get("I"))
    else:
        WORKING_DIR = '/home/juliaeis/Dokumente/OGGM/work_dir/reconstruction/global/'
        OUT_DIR = WORKING_DIR
        cfg.PATHS['working_dir'] = WORKING_DIR
        utils.mkdir(WORKING_DIR, reset=False)
        REGION='05'

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
    cfg.PARAMS['dl_verify'] = True

    # add to BASENAMES
    _doc = 'contains observed and searched glacier from synthetic experiment to find intial state'
    cfg.BASENAMES['synthetic_experiment'] = ('synthetic_experiment.pkl', _doc)

    # We use intersects
    db = utils.get_rgi_intersects_region_file(version='61', region=REGION)
    cfg.set_intersects_db(db)

    # RGI file
    path = utils.get_rgi_region_file(REGION, version='61')
    rgidf = gpd.read_file(path)
    rgidf = rgidf.sort_values('Area', ascending=False)

    # exclude non-landterminating glaciers
    rgidf = rgidf[rgidf.TermType == 0]
    rgidf = rgidf[rgidf.Connect != 2]

    # gdirs with all glaciers
    gdirs = workflow.init_glacier_regions(rgidf)

    # here you will need to select the gdirs which are wgms glaciers



    # runs only a quarter of glaciers per job array
    gdirs = gdirs[JOB_NR:len(gdirs):4]

    t_0 = 1917
    epsilon = 125
    diff = pd.DataFrame()

    for gdir in gdirs:
        try:
            # copy previous files to gdir.dir
            dir = os.path.join(OUT_DIR,'per_glacier',gdir.dir.split('per_glacier/')[-1])
            os.system('cp -rf '+dir+'/* '+ gdir.dir)

            t_e = gdir.rgi_date
            ex = [f for f in os.listdir(gdir.dir) if f.startswith('model_run_ad')]
            if len(ex)==1 :
                dst = os.path.join(gdir.dir,ex[0])
                ex_mod = FileModel(dst)

                bias = float(ex[0].split('_')[-1].split('.nc')[0])

                # run initialization
                df = find_possible_glaciers(gdir, t_0, t_e, 200, ex_mod, bias, delete=False)

                df.fitness = pd.to_numeric(df.fitness / 125)
                df = df.dropna(subset=['fitness'])

                med_mod, perc_min, perc_max = find_median(df)


                # from her onwards: how I have calcuated the rmse, ... to leclercq
                lec.loc[1917] = np.nan
                lec = lec.sort_index().interpolate()[lec.index >= 1917]

                rmse = np.sqrt(((lec - med_mod.length_m_ts(rollmin=5)[lec.index]) ** 2).mean())
                rmspe = np.sqrt((((lec - med_mod.length_m_ts(rollmin=5)[lec.index]) / lec) ** 2).mean()) * 100
                error = (lec - med_mod.length_m_ts(rollmin=5)[lec.index]).mean()
                perc_error = ((lec - med_mod.length_m_ts(rollmin=5)[lec.index]) / lec).mean()

                max = (lec - med_mod.length_m_ts(rollmin=5)[lec.index]).max()
                min = (lec - med_mod.length_m_ts(rollmin=5)[lec.index]).min()
                if abs(max) > abs(min):
                    max_diff = max
                else:
                    max_diff = min

                temp_bias = cfg.PATHS['working_dir'].split('_')[-1]

                diff.loc[gdir.rgi_id, 'region'] = REGION
                diff.loc[gdir.rgi_id, 'rmse'] = rmse
                diff.loc[gdir.rgi_id, 'rmspe'] = rmspe
                diff.loc[gdir.rgi_id, 'error'] = error
                diff.loc[gdir.rgi_id, 'perc_error'] = perc_error
                diff.loc[gdir.rgi_id, 'max_diff'] = max_diff


        except Exception as e:
            print(e)

    diff.to_pickle(os.path.join(cfg.PATHS['working_dir'], REGION + '_leclercq_difference.pkl'))
