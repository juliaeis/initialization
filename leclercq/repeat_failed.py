import os
import pandas as pd
import geopandas as gpd
from oggm import cfg, utils

if __name__ == '__main__':

    cfg.initialize()

    ON_CLUSTER = False
    REGION = '13'

    # Local paths
    if ON_CLUSTER:
        WORKING_DIR = os.environ.get("S_WORKDIR")
        cfg.PATHS['working_dir'] = WORKING_DIR
        TEMP_BIAS = float(os.environ.get('TEMP_BIAS'))
        REGION = str(os.environ.get('I')).zfill(2)

    else:
        WORKING_DIR = '/home/juliaeis/Dokumente/OGGM/work_dir/reconstruction/leclercq_diff2'
        cfg.PATHS['working_dir'] = WORKING_DIR
        utils.mkdir(WORKING_DIR, reset=False)
        TEMP_BIAS = 0

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

    # RGI file
    path = utils.get_rgi_region_file(REGION, version='61')
    rgidf = gpd.read_file(path)

    # exclude non-landterminating glaciers
    rgidf = rgidf[rgidf.TermType == 0]
    rgidf = rgidf[rgidf.Connect != 2]

    # We use intersects
    db = utils.get_rgi_intersects_region_file(version='61', region=REGION)
    cfg.set_intersects_db(db)

    # read leclercq links
    lec = pd.read_csv('rgi_leclercq_links_2014_RGIV6.csv')
    lec.loc[:,'REGION'] = lec.RGI_ID.apply(lambda x: x.split('-')[-1].split('.')[0])
    lec = lec[lec.REGION==str(REGION)]

    # find missing glaciers
    exp_df = pd.read_pickle(os.path.join(cfg.PATHS['working_dir'], 'quick_experiment_df.pkl'))
    exp_df = exp_df[exp_df.temp_bias==TEMP_BIAS]
    exp_df.loc[:,'region'] = exp_df.rgi_id.apply(lambda x: x.split('RGI60-')[-1].split('.')[0])
    exp_df = exp_df[exp_df.region==REGION]
    id = pd.concat([lec.RGI_ID, exp_df.rgi_id]).drop_duplicates(keep=False)
    rgidf = rgidf[rgidf.RGIId.isin(id)]

    print(rgidf)




