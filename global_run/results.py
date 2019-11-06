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


sys.path.append('../')
sys.path.append('../../')
from  initialization.core import *

def read_results(gdirs):

    model_df = pd.DataFrame()
    pool = mp.Pool()
    list = pool.map(read_result_parallel, gdirs)
    pool.close()
    pool.join()
    model_df = model_df.append(list, ignore_index=True)
    return model_df

def read_result_parallel(gdir):

    t_e = gdir.rgi_date
    ex = [f for f in os.listdir(gdir.dir) if f.startswith('model_run_ad')]

    if len(ex) == 1:
        path = os.path.join(gdir.dir, ex[0])
        ex_mod = FileModel(path)
        ye = ex_mod.volume_km3_ts().index[-1]
        bias = float(ex[0].split('_')[-1].split('.nc')[0])
        ts = ex_mod.volume_km3_ts()
        if ye <2016:
            try:
                ex_mod.run_until(ye)
                tasks.run_from_climate_data(gdir, ys=ye, ye=2016, bias=bias,
                                            output_filesuffix='_to_2016',
                                            init_model_fls=copy.deepcopy(ex_mod.fls))

                res_mod = FileModel(gdir.get_filepath('model_run', filesuffix='_to_2016'))
                ts2 = res_mod.volume_km3_ts()
                ts2 = ts2[ts2.index[1:]]
                ts = pd.concat([ts,ts2])
            except:
                pass
        ts['rgi_id'] = gdir.rgi_id
        return ts
    else:
        return pd.Series({'rgi_id':gdir.rgi_id})


if __name__ == '__main__':
    cfg.initialize()

    ON_CLUSTER = False
    REG = '01'

    # Local paths
    if ON_CLUSTER:
        OUT_DIR = os.environ.get("OUTDIR")
        cfg.PATHS['working_dir'] = OUT_DIR
        REG = os.environ.get("I")
    else:
        WORKING_DIR = '/home/juliaeis/Dokumente/OGGM/work_dir/reconstruction/global'
        OUT_DIR = WORKING_DIR
        cfg.PATHS['working_dir'] = WORKING_DIR
        utils.mkdir(WORKING_DIR, reset=False)

    cfg.PATHS['plot_dir'] = os.path.join(cfg.PATHS['working_dir'], 'plots')
    utils.mkdir(cfg.PATHS['plot_dir'], reset=False)

    for dir in os.listdir(OUT_DIR):
        if dir.startswith('reg1'):
            cfg.PATHS['working_dir'] = os.path.join(OUT_DIR,dir)
            REGION = dir.split('reg')[-1].split('-')[0].zfill(2)
            if REGION ==  REG:

                # RGI file
                path = utils.get_rgi_region_file(REGION, version='61')
                rgidf = gpd.read_file(path)
                #rgidf = rgidf.sort_values('Area', ascending=False)

                # exclude non-landterminating glaciers
                rgidf = rgidf[rgidf.TermType == 0]
                rgidf = rgidf[rgidf.Connect != 2]

                gdirs = workflow.init_glacier_regions(rgidf.head(9))
                df = read_results(gdirs)
                print(df)