import sys
sys.path.append('../')
from initialization.core import *
import shutil
import geopandas as gpd
from oggm import cfg, utils
from oggm.core.massbalance import MultipleFlowlineMassBalance, PastMassBalance
pd.options.mode.chained_assignment = None
from oggm.workflow import execute_entity_task
import matplotlib.pyplot as plt



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
        WORKING_DIR = '/home/juliaeis/Dokumente/OGGM/work_dir/reconstruction/wgms/temp_0'
        OUT_DIR = WORKING_DIR
        cfg.PATHS['working_dir'] = WORKING_DIR
        utils.mkdir(WORKING_DIR, reset=False)
        REGION='03'
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
    rgidf = rgidf.loc[rgidf.RGIId.isin(['RGI60-03.00840'])]

    # initialize glaciers
    gdirs = workflow.init_glacier_regions(rgidf)
    t_0 = 1917
    diff = pd.DataFrame()
    delta_diff = pd.DataFrame()

    for gdir in gdirs:
        df = pd.DataFrame()
        temp_bias = cfg.PATHS['working_dir'].split('_')[-1]
        refmb = gdir.get_ref_mb_data().copy()

        # get experiment
        ex = [f for f in os.listdir(gdir.dir) if f.startswith('model_run_ad')]
        if len(ex) == 1:
            dst = os.path.join(gdir.dir, ex[0])
            mod = FileModel(dst)
            # get optimal mb bias
            bias = float(ex[0].split('_')[-1].split('.nc')[0])
            years = mod.volume_km3_ts().index



            # if observation record longer than rgi_date: create new model which can be run until observation record
            if refmb.index[-1]>years[-1]:
                mod.run_until(t_0)
                tasks.run_from_climate_data(gdir, ys=t_0, ye=refmb.index[-1],
                                            init_model_fls=copy.deepcopy(mod.fls),
                                            output_filesuffix='_until_refmb', bias=bias)
                mod = FileModel(gdir.get_filepath('model_run', filesuffix='_until_refmb'))

            # get mass balance from volume difference
            df.loc[:-1, 'OGGM_dv'] = mod.volume_m3_ts().diff() * cfg.PARAMS['ice_density'] / mod.area_m2_ts()
            df = df.shift(-1)

            for yr in  mod.volume_km3_ts().index:
                mod.run_until(yr)
                mb = MultipleFlowlineMassBalance(gdir,fls=copy.deepcopy( mod.fls),
                                                 mb_model_class=PastMassBalance, bias=bias)
                df.loc[yr, 'OGGM_mb']=mb.get_specific_mb(year=[mod.yr])

        df.loc[:, 'WGMS'] = refmb.ANNUAL_BALANCE
        df.index = df.index.astype(int)

        # difference between Mass Balance and volume delta
        rmse_d = np.sqrt(((df.OGGM_mb-df.OGGM_dv)**2).mean())
        max_d = (df.OGGM_mb-df.OGGM_dv).abs().max()
        delta_diff.loc[gdir.rgi_id, 'region'] = REGION
        delta_diff.loc[gdir.rgi_id, 'rmse'] = rmse_d
        delta_diff.loc[gdir.rgi_id, 'max_diff'] = max_d
        delta_diff.loc[gdir.rgi_id, 'temp_bias'] = temp_bias

        # difference between modelled and observed mass balance
        df = df.dropna(subset=['WGMS'])
        rmse = np.sqrt(((df.WGMS - df.OGGM_mb) ** 2).mean())
        error = (df.WGMS - df.OGGM_mb).mean()
        max = (df.WGMS - df.OGGM_mb).max()
        min = (df.WGMS - df.OGGM_mb).min()
        if abs(max) > abs(min):
            max_diff = max
        else:
            max_diff = min

        diff.loc[gdir.rgi_id, 'region'] = REGION
        diff.loc[gdir.rgi_id, 'rmse'] = rmse
        diff.loc[gdir.rgi_id, 'error'] = error
        diff.loc[gdir.rgi_id, 'max_diff'] = max_diff
        diff.loc[gdir.rgi_id, 'temp_bias'] = temp_bias

    delta_diff.to_csv(os.path.join(cfg.PATHS['working_dir'], REGION+'_'+ str(JOB_NR) + 'OGGM_instablity.csv'))
    diff.to_csv(os.path.join(cfg.PATHS['working_dir'], REGION+'_'+ str(JOB_NR) + '_leclercq_difference.csv'))


