
import sys
sys.path.append('../')
from initialization.core import *
import shutil
import geopandas as gpd
from oggm import cfg, utils
pd.options.mode.chained_assignment = None
from oggm.workflow import execute_entity_task



if __name__ == '__main__':
    cfg.initialize()

    ON_CLUSTER = False

    # Local paths
    if ON_CLUSTER:
        WORKING_DIR = os.environ.get("S_WORKDIR")
        cfg.PATHS['working_dir'] = WORKING_DIR
        OUT_DIR = os.path.join(os.environ.get("OUTDIR"), WORKING_DIR.split('/')[-1])

        # get this information from bash script
        REGION = str(os.environ.get('REGION')).zfill(2)
        # should be 0-3
        JOB_NR = int(os.environ.get('I'))
        TEMP_BIAS = float(os.environ.get('TEMP_BIAS'))
    else:
        WORKING_DIR = '/home/juliaeis/Dokumente/OGGM/work_dir/reconstruction/wgms/'
        OUT_DIR = WORKING_DIR
        cfg.PATHS['working_dir'] = WORKING_DIR
        utils.mkdir(WORKING_DIR, reset=False)
        REGION='11'
        JOB_NR = 0

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

    # exclude non-landterminating glaciers
    rgidf = rgidf[rgidf.TermType == 0]
    rgidf = rgidf[rgidf.Connect != 2]

    wgms = utils.get_ref_mb_glaciers_candidates()

    # Keep only the wgms reference glaciers
    rgidf = rgidf.loc[rgidf.RGIId.isin(wgms)]

    # initialize glaciers
    gdirs = workflow.init_glacier_regions(rgidf)

    # runs only a quarter of glaciers per job array
    #gdirs= gdirs[JOB_NR:len(gdirs):4]

    t_0 = 1917
    epsilon = 125

    preprocessing(gdirs)

    advanced_experiments(gdirs, [TEMP_BIAS], t_0)


