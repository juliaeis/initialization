
import sys
sys.path.append('../')
from initialization.core import *
from shutil import copyfile



import geopandas as gpd
from oggm import cfg, utils
pd.options.mode.chained_assignment = None


if __name__ == '__main__':
    cfg.initialize()

    ON_CLUSTER = False

    # Local paths
    if ON_CLUSTER:
        WORKING_DIR = os.environ.get("S_WORKDIR")
        cfg.PATHS['working_dir'] = WORKING_DIR
        OUT_DIR = os.path.join(os.environ.get("OUTDIR"),'global')
        REGION = str(os.environ.get('REGION')).zfill(2)
        JOB_NR = int(os.environ.get('I'))
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
    rgidf = rgidf.sort_values('Area', ascending=False)
    #rgidf = rgidf[rgidf.RGIId == 'RGI60-05.00789']

    #select 24th-48th largest glaciers
    gdirs = workflow.init_glacier_regions(rgidf[-1:])

    t_0 = 1917

    epsilon = 125
    #preprocessing(gdirs)
    #advanced_experiments(gdirs, [0], 1917, REGION)


    for gdir in gdirs:
        # split command works different on Cluster and localy
        dir = os.path.join(OUT_DIR, gdir.dir.split('/global/')[-1])
        ex = [f for f in os.listdir(dir) if f.startswith('model_run_ad')]
        t_e = gdir.rgi_date
        if len(ex)==1 :
            if dir!=gdir.dir:
                # copy advanced experiment to gdir.dir
                src = os.path.join(dir,ex[0])
                dst = os.path.join(gdir.dir,ex[0])
                copyfile(src, dst)
            else:
                dst = os.path.join(gdir.dir,ex[0])

            ex_mod = FileModel(dst)
            bias = float(ex[0].split('_')[-1].split('.nc')[0])
            print(ex_mod.area_km2_ts()[t_e])
    '''
            df = find_possible_glaciers(gdir, t_0, t_e, 200, ex_mod, bias)

    '''

