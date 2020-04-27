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
import matplotlib.pyplot as plt
from oggm.utils._workflow import GlacierDirectory
from oggm.graphics import plot_googlemap, plot_centerlines, plot_modeloutput_section, plot_modeloutput_map


sys.path.append('../paper')
from plots_paper import *

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
    if len(list)>0:
        model_df = model_df.append(list, ignore_index=True)
        return model_df

def read_result_parallel(gdir):

    try:
        ye = gdir.rgi_date
        fls = gdir.read_pickle('model_flowlines')
        mod = FluxBasedModel(flowlines=fls)

        lec = get_ref_length_data(gdir).dL

        if lec.index[0] < 1917:
            lec.loc[1917] = np.nan
            lec = lec.sort_index().interpolate()
            lec = lec - lec.loc[1917]
            lec = lec[lec.index >=1917]

        '''
        lec = lec[lec.index <= ye]

        lec = (lec - lec.iloc[-1]) + mod.length_m

        ini_df = pd.read_pickle(os.path.join(gdir.dir,'result1917.pkl'), compression='gzip')
        ini_df.loc[:,'fitness'] = pd.to_numeric(ini_df.fitness / 125)
        ini_df = ini_df.dropna(subset=['fitness'])

        med_mod, perc_min, perc_max = find_median(ini_df)

        temp_bias = cfg.PATHS['working_dir'].split('_')[-1]
        #return pd.Series({'rgi':gdir.rgi_id, 'temp_bias':temp_bias, 'median':med_mod})


        lec.loc[1917] = np.nan
        lec = lec.sort_index().interpolate()[lec.index >= 1917]
        lec.loc['rgi_id']= gdir.rgi_id
        print(med_mod.length_m_ts()[lec.index])
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
        lec.loc['rgi'] = gdir.rgi_id
        return lec
    except:
        return pd.Series({'rgi':gdir.rgi_id})



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
        WORKING_DIR = '/home/juliaeis/Dokumente/OGGM/work_dir/reconstruction/new_leclercq/temp_0'
        cfg.PATHS['working_dir'] = WORKING_DIR
        utils.mkdir(WORKING_DIR, reset=False)

    cfg.PATHS['plot_dir'] = os.path.join(cfg.PATHS['working_dir'], 'plots')
    utils.mkdir(cfg.PATHS['plot_dir'], reset=False)
    lec_df = pd.DataFrame()
    n1 = 0
    n2 = 0
    n3 = 0
    for REGION in range(1,19):
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
        #rgidf = rgidf[rgidf.RGIId.isin(['RGI60-09.00144'])]
        #rgidf = rgidf[rgidf.RGIId.isin(['RGI60-11.00897'])]

        # exclude non-landterminating glaciers
        rgidf = rgidf[rgidf.TermType==0]
        rgidf = rgidf[rgidf.Connect !=2]
        # sort for efficient using
        rgidf = rgidf.sort_values('Area', ascending=True)

        gdirs = workflow.init_glacier_regions(rgidf)
        #workflow.gis_prepro_tasks(gdirs)
        #workflow.climate_tasks(gdirs)
        #workflow.inversion_tasks(gdirs)
        #workflow.execute_entity_task(tasks.init_present_time_glacier, gdirs)

        '''
        try:
            df = read_results(gdirs).set_index('rgi')
            print(df)
            lec_df = lec_df.append(df, ignore_index=False)

        except:
            pass
        lec_df.to_csv(os.path.join(cfg.PATHS['working_dir'],'leclercq_dl.csv'))
        #df.to_pickle(os.path.join(cfg.PATHS['working_dir'],'lec_median_0.pkl'))

        print(lec_df)
        '''
        lec_df = pd.DataFrame(columns=np.arange(1917,2018))
        n500 = 0
        n1000 = 0
        n2000 = 0

        for gdir in gdirs:

            if gdir.was_tidewater(l_min=500):
                n500=n500+1
                ex = [f for f in os.listdir(gdir.dir) if
                      f.startswith('model_run_ad')]
                dst = os.path.join(gdir.dir, ex[0])
                ex_mod = FileModel(dst)


                #print(get_ref_length_data(gdir))
            if gdir.was_tidewater(l_min=1000):
                n1000 = n1000 + 1
            if gdir.was_tidewater(l_min=2000):
                n2000 = n2000 + 1


            '''
                fig, (ax1, ax2) = plt.subplots(1, 2,figsize=(20,10))
                plot_googlemap(gdir, ax=ax1)
                plot_centerlines(gdir,use_flowlines=True,add_downstream=True, ax=ax2)
                plt.savefig(os.path.join(cfg.PATHS['plot_dir'], '1000m', gdir.rgi_id+'.png'))
                #plt.show()


            ex_df = pd.DataFrame()

            t_e = gdir.rgi_date
            ex = [f for f in os.listdir(gdir.dir) if
                  f.startswith('model_run_ad')]

            lec = get_ref_length_data(gdir).dL

            if (lec.index[0]<1917) and (not 1917 in lec.index):
                lec.loc[1917] = np.nan
                lec = lec.sort_index().interpolate(method='slinear')[lec.index >= 1917]
            lec = lec-lec.iloc[0]
            lec.name=gdir.rgi_id
            lec_df = lec_df.append(lec, ignore_index=False)

        print(lec_df)
        '''
        if n500!= 0:
            print(REGION, len(gdirs),  n500, n1000, n2000)
        n1 = n1 + n500
        n2 = n2 + n1000
        n3 = n3 + n2000
    print(n1, n2, n3)
