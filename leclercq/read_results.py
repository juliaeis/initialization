import os
import sys
import copy
import numpy as np
import xarray as xr
import pandas as pd
import multiprocessing as mp
from oggm import utils, cfg
import geopandas as gpd
from oggm.utils._downloads import get_demo_file
from oggm.core.flowline import FluxBasedModel, FileModel

sys.path.append('../')
sys.path.append('../../')
from  initialization.core import *

def find_median(df):

    try:
        accept_df = df[df.fitness <= 1]
        quant_df = accept_df[accept_df.fitness <= accept_df.fitness.quantile(0.05)]

        # median state
        quant_df.loc[:,'length']= quant_df.model.apply(lambda x: x.length_m)
        quant_df = quant_df.sort_values('length', ascending=False)
        l = len(quant_df)
        if l % 2:
            index = int((l - 1) / 2)
        else:
            index = int(l / 2)
        return copy.deepcopy(quant_df.iloc[index].model), quant_df.at[quant_df.length.idxmin(),'model'], quant_df.at[quant_df.length.idxmax(),'model']

    except:
        return copy.deepcopy(df.loc[df.fitness.idxmin()].model), None, None


def get_ref_length_data(gdir):
    """Get the glacier lenght data from P. Leclercq's data base.

     https://folk.uio.no/paulwl/data.php

     For some glaciers only!
     """

    df = pd.read_csv('rgi_leclercq_links_2014_RGIV6.csv')
    df = df.loc[df.RGI_ID == gdir.rgi_id]
    if len(df) == 0:
        raise RuntimeError('No length data found for this glacier!')
    ide = df.LID.values[0]

    f = get_demo_file('Glacier_Lengths_Leclercq.nc')
    with xr.open_dataset(f) as dsg:
        # The database is not sorted by ID. Don't ask me...
        grp_id = np.argwhere(dsg['index'].values == ide)[0][0] + 1
    with xr.open_dataset(f, group=str(grp_id)) as ds:
        df = ds.to_dataframe()
        df.name = ds.glacier_name
    return df

def read_results(gdirs):

    model_df = pd.DataFrame()
    pool = mp.Pool()
    list = pool.map(read_result_parallel, gdirs)
    pool.close()
    pool.join()
    model_df = model_df.append(list, ignore_index=True)
    return model_df

def read_result_parallel(gdir):


    ye = gdir.rgi_date
    fls = gdir.read_pickle('model_flowlines')
    mod = FluxBasedModel(flowlines=fls)

    lec = get_ref_length_data(gdir).dL
    lec = lec[lec.index <= ye]

    lec = (lec - lec.iloc[-1]) + mod.length_m

    ini_df = pd.read_pickle(os.path.join(gdir.dir,'result1917.pkl'), compression='gzip')
    ini_df.loc[:,'fitness'] = pd.to_numeric(ini_df.fitness / 125)
    ini_df = ini_df.dropna(subset=['fitness'])

    med_mod, perc_min, perc_max = find_median(ini_df)
    temp_bias = cfg.PATHS['working_dir'].split('_')[-1]
    return pd.Series({'rgi':gdir.rgi_id, 'temp_bias':temp_bias, 'median':med_mod})

    '''
    lec.loc[1917] = np.nan
    lec = lec.sort_index().interpolate()[lec.index >= 1917]
    med_mod.length_m_ts(rollmin=5)[lec.index]


    rmse = np.sqrt(((lec - med_mod.length_m_ts(rollmin=5)[lec.index]) ** 2).mean())
    error = (lec - med_mod.length_m_ts(rollmin=5)[lec.index]).mean()
    max = (lec - med_mod.length_m_ts(rollmin=5)[lec.index]).max()
    min = (lec - med_mod.length_m_ts(rollmin=5)[lec.index]).min()
    if abs(max)>abs(min):
        max_diff = max
    else:
        max_diff = min


    # saves median state, minimum state and experiment model
    return pd.Series({'rgi':gdir.rgi_id,'region':gdir.rgi_region,
                      'length':mod.length_m, 'temp_bias':temp_bias,
                      'rmse':rmse, 'max_diff':max_diff, 'error':error})

    '''

    #except:
    #    return pd.Series({'rgi':gdir.rgi_id})


if __name__ == '__main__':

    cfg.initialize()

    ON_CLUSTER = False
    REGION = '01'

    # Local paths
    if ON_CLUSTER:
        WORKING_DIR = '/home/users/julia/initialization/out/temp_0'
        cfg.PATHS['working_dir'] = WORKING_DIR
        #REGION = str(os.environ.get('I')).zfill(2)
    else:
        WORKING_DIR = '/home/juliaeis/Dokumente/OGGM/work_dir/reconstruction/leclercq_diff/temp_0'
        cfg.PATHS['working_dir'] = WORKING_DIR
        utils.mkdir(WORKING_DIR, reset=False)

    cfg.PATHS['plot_dir'] = os.path.join(cfg.PATHS['working_dir'], 'plots')
    utils.mkdir(cfg.PATHS['plot_dir'], reset=False)

    for REGION in range(11,12):
        REGION = str(REGION).zfill(2)

        # read leclercq links
        lec = pd.read_csv('rgi_leclercq_links_2014_RGIV6.csv')
        lec['REGION'] = lec.RGI_ID.apply(lambda x: x.split('-')[-1].split('.')[0])

        # We use intersects
        db = utils.get_rgi_intersects_region_file(version='61', region=REGION)
        cfg.set_intersects_db(db)

        # RGI file
        path = utils.get_rgi_region_file(REGION, version='61')
        rgidf = gpd.read_file(path)

        # only the ones with leclercq observation
        rgidf = rgidf[rgidf.RGIId.isin(lec[lec.REGION==REGION].RGI_ID.values)]
        rgidf = rgidf[rgidf.RGIId.isin(['RGI60-11.00897'])]

        # exclude non-landterminating glaciers
        rgidf = rgidf[rgidf.TermType==0]
        rgidf = rgidf[rgidf.Connect !=2]
        # sort for efficient using
        rgidf = rgidf.sort_values('Area', ascending=True)

        gdirs = workflow.init_glacier_regions(rgidf)
        df = read_results(gdirs)
        df.to_pickle(os.path.join(cfg.PATHS['working_dir'],'lec_median_0.pkl'))