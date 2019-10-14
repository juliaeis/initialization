
import os
import sys
from copy import deepcopy
sys.path.append('../')
from initialization.core import *
from plots_paper import *

import matplotlib.pyplot as plt
import salem
import geopandas as gpd
from oggm import cfg, workflow, utils, tasks
from oggm.core.flowline import FluxBasedModel
pd.options.mode.chained_assignment = None
import time
import matplotlib as mpl
import xarray as xr

mpl.rcParams['axes.linewidth'] = 3
mpl.rcParams['xtick.major.width'] = 3
mpl.rcParams['ytick.major.width'] = 3
mpl.rcParams['font.size'] =20
mpl.rcParams['font.weight'] = 'medium'
mpl.rcParams['axes.labelweight'] = 'medium'
mpl.rcParams['legend.fontsize'] = 20 #30
mpl.rcParams['lines.linewidth'] = 3

from scipy.stats import pearsonr



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


def read_result_parallel(gdir):

    rp = gdir.get_filepath('model_run', filesuffix='_experiment')
    ex_mod = FileModel(rp)


    df = pd.read_pickle(os.path.join(gdir.dir, 'result1850.pkl'),
                        compression='gzip')

    df.fitness = df.fitness / 125

    plot_supplement(gdir, df, ex_mod, 1850,2000,125, cfg.PATHS['plot_dir'])
    try:
        pass
        '''
        acc_df = df[df.fitness<1]
        acc_df.loc[:, 'length'] = acc_df.model.apply(lambda x: x.length_m)
        med_mod, perc_min, perc_max = find_median(df)
        min_mod = deepcopy(df.loc[df.fitness.idxmin(), 'model'])
        df['fitness2'] = df.model.apply(
            lambda x: abs(x.area_km2_ts()[2000] - ex_mod.area_km2_ts()[2000]) ** 2)
        df['fitness3'] = df.model.apply(
            lambda x: abs(x.length_m_ts()[2000] - ex_mod.length_m_ts()[2000]) ** 2)

        # saves median state, minimum state and experiment model
        return pd.Series({'rgi':gdir.rgi_id,'median':med_mod,
                          'minimum':min_mod, 'experiment':ex_mod,
                          'flowline':FluxBasedModel(flowlines=gdir.read_pickle('model_flowlines')),
                          'fit2':df.loc[df.fitness2.idxmin(), 'model'],
                          'fit3':df.loc[df.fitness3.idxmin(), 'model']})
        '''
        return pd.Series({'rgi': gdir.rgi_id})
    except:
        return pd.Series({'rgi':gdir.rgi_id})

def mean_thick(mod):
    try:
        mod.run_until(2000)
    except:
        pass
    thick = []
    for t in [mod.fls[i].thick for i in range(len(mod.fls))]:
        thick.extend(t[t>0])
    return np.mean(thick)

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
    rgidf = rgidf[rgidf.RGIId.isin(['RGI60-11.01450', 'RGI60-11.00003',
                                    'RGI60-11.00026', 'RGI60-11.00049',
                                    'RGI60-11.00068', 'RGI60-11.00080',
                                    'RGI60-11.00106', 'RGI60-11.00251',
                                    'RGI60-11.00289', 'RGI60-11.00300',
                                    'RGI60-11.00402', 'RGI60-11.00719',
                                    'RGI60-11.00787', 'RGI60-11.00887',
                                    'RGI60-11.00918', 'RGI60-11.02791',
                                    'RGI60-11.02715', 'RGI60-11.01246'])]
    #rgidf = rgidf[rgidf.RGIId.isin(['RGI60-11.00074'])]
    # sort for efficient using
    rgidf = rgidf.sort_values('Area', ascending=False)
    gdirs = workflow.init_glacier_regions(rgidf)

    df = read_results(gdirs).dropna().set_index('rgi')
    print(df)
    '''
    df = df[['flowline','experiment','median']]
    df['area']=df['experiment'].apply(lambda x: x.area_km2_ts()[2000])
    df['volume'] = df['experiment'].apply(lambda x: x.volume_km3_ts()[2000])
    df['length'] = df['experiment'].apply(lambda x: x.length_m_ts()[2000]/1000)

    print(df[['area','volume','length']])

    print(df['flowline'].apply(lambda x: x.length_m/1000))
    print(df['flowline'].apply(lambda x: x.area_km2))
    print(df['flowline'].apply(lambda x: x.volume_km3))

    print(df['median'].apply(mean_thick))
    print(df['flowline'].apply(mean_thick))
    '''