import sys
sys.path.append('../')
from initialization.core import *
from paper.plots_paper import *
from oggm.core.massbalance import MultipleFlowlineMassBalance, PastMassBalance

import geopandas as gpd
from oggm import cfg, utils
pd.options.mode.chained_assignment = None
from copy import deepcopy

if __name__ == '__main__':
    cfg.initialize()

    ON_CLUSTER = True

    # Local paths
    if ON_CLUSTER:
        WORKING_DIR = os.environ.get("S_WORKDIR")
        cfg.PATHS['working_dir'] = WORKING_DIR
        OUT_DIR = os.environ.get("OUTDIR")
        REGION = str(os.environ.get('REGION')).zfill(2)

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

    wgms = utils.get_ref_mb_glaciers_candidates()

    # Keep only the wgms reference glaciers
    rgidf = rgidf.loc[rgidf.RGIId.isin(wgms)]

    # initialize glaciers
    gdirs = workflow.init_glacier_regions(rgidf)

    t_0 = 1917
    epsilon = 125
    exp_df = pd.DataFrame()


    for gdir in gdirs:
        df = pd.DataFrame()
        try:
            # copy previous files to gdir.dir
            dir = os.path.join(OUT_DIR,'per_glacier',gdir.dir.split('per_glacier/')[-1])
            os.system('cp -rf '+dir+'/* '+ gdir.dir)

            t_e = gdir.rgi_date
            ex = [f for f in os.listdir(gdir.dir) if f.startswith('model_run_ad')]
            if len(ex)==1 :
                # read experiment
                dst = os.path.join(gdir.dir,ex[0])
                ex_mod = FileModel(dst)

                # get mb bias and temp_bias
                bias = float(ex[0].split('_')[-1].split('.nc')[0])
                temp_bias = cfg.PATHS['working_dir'].split('_')[-1]

                exp_df.loc[gdir.rgi_id, 'region'] = REGION
                exp_df.loc[gdir.rgi_id, 'rgi_date'] = gdir.rgi_date
                exp_df.loc[gdir.rgi_id, 'mb_bias'] = bias
                exp_df.loc[gdir.rgi_id, 'temp_bias'] = temp_bias


        except:
            exp_df.loc[gdir.rgi_id, 'region'] = REGION
    exp_df.to_csv(os.path.join(cfg.PATHS['working_dir'], REGION + '_experiment.pkl'))