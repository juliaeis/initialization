
import sys
sys.path.append('../')
from initialization.core import *
from paper.plots_paper import *
from oggm.core.massbalance import MultipleFlowlineMassBalance, PastMassBalance

import geopandas as gpd
from oggm import cfg, utils
pd.options.mode.chained_assignment = None
from copy import deepcopy

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
        cfg.PATHS['working_dir'] = WORKING_DIR
        OUT_DIR = os.environ.get("OUTDIR")
        REGION = str(os.environ.get('REGION')).zfill(2)
        JOB_NR =  int(os.environ.get("I"))
    else:
        WORKING_DIR = '/home/juliaeis/Dokumente/OGGM/work_dir/reconstruction/wgms/temp_-0.25'
        OUT_DIR = WORKING_DIR
        cfg.PATHS['working_dir'] = WORKING_DIR
        utils.mkdir(WORKING_DIR, reset=False)
        REGION='01'
        JOB_NR=0

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
    rgidf = rgidf.loc[rgidf.RGIId.isin(['RGI60-01.04591'])]

    # initialize glaciers
    gdirs = workflow.init_glacier_regions(rgidf)

    # runs only a quarter of glaciers per job array
    gdirs = gdirs[JOB_NR:len(gdirs):4]

    t_0 = 1917
    epsilon = 125
    diff = pd.DataFrame()
    delta_diff = pd.DataFrame()

    for gdir in gdirs:
        df = pd.DataFrame()
        try:
            # copy previous files to gdir.dir
            dir = os.path.join(OUT_DIR,'per_glacier',gdir.dir.split('per_glacier/')[-1])
            #os.system('cp -rf '+dir+'/* '+ gdir.dir)
            refmb = gdir.get_ref_mb_data().copy()
            t_e = gdir.rgi_date
            ex = [f for f in os.listdir(gdir.dir) if f.startswith('model_run_ad')]
            if len(ex)==1 :
                # read experiment
                dst = os.path.join(gdir.dir,ex[0])
                ex_mod = FileModel(dst)

                # get mb bias and temp_bias
                bias = float(ex[0].split('_')[-1].split('.nc')[0])
                temp_bias = cfg.PATHS['working_dir'].split('_')[-1]

                # run initialization
                res_df = find_possible_glaciers(gdir, t_0, t_e, 200, ex_mod, bias, delete=False)

                res_df.fitness = pd.to_numeric(res_df.fitness / 125)
                res_df = res_df.dropna(subset=['fitness'])

                # get median and percentile states
                mod, perc_min, perc_max = find_median(res_df)

                # if observation record longer than rgi_date: create new model which can be run until last observation record
                if refmb.index[-1] > gdir.rgi_date:
                    mod.run_until(t_0)
                    tasks.run_from_climate_data(gdir, ys=t_0,
                                                ye=refmb.index[-1],
                                                init_model_fls=deepcopy(mod.fls),
                                                output_filesuffix='_until_refmb',
                                                bias=bias)
                    mod = FileModel(gdir.get_filepath('model_run',
                                                      filesuffix='_until_refmb'))

                # get modelled mass balance from volume difference
                df.loc[:-1, 'OGGM_dv'] = mod.volume_m3_ts().diff() * cfg.PARAMS['ice_density'] / mod.area_m2_ts()
                df = df.shift(-1)

                # get mass balance from MassBalanceModel
                for yr in mod.volume_km3_ts().index:
                    mod.run_until(yr)
                    mb = MultipleFlowlineMassBalance(gdir, fls=deepcopy(mod.fls),
                                                     mb_model_class=PastMassBalance,
                                                     bias=bias)
                    df.loc[yr, 'OGGM_mb'] = mb.get_specific_mb(year=[mod.yr])


                # set WGMS data
                df.loc[:, 'WGMS'] = refmb.ANNUAL_BALANCE
                df.index = df.index.astype(int)

                # difference between Mass Balance and volume delta
                rmse_d = np.sqrt(((df.OGGM_mb - df.OGGM_dv) ** 2).mean())
                max_d = (df.OGGM_mb - df.OGGM_dv).abs().max()
                delta_diff.loc[gdir.rgi_id, 'region'] = REGION
                delta_diff.loc[gdir.rgi_id, 'rmse'] = rmse_d
                delta_diff.loc[gdir.rgi_id, 'max_diff'] = max_d
                delta_diff.loc[gdir.rgi_id, 'temp_bias'] = temp_bias

                # difference between modelled and observed mass balance at WGMS years
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


        except Exception as e:
            print(e)
    diff.to_csv(os.path.join(cfg.PATHS['working_dir'], REGION + '_' + str(JOB_NR)+ '_leclercq_difference.csv'))
    delta_diff.to_csv(os.path.join(cfg.PATHS['working_dir'], REGION + '_' + str(JOB_NR) + '_OGGM_instablity.csv'))
