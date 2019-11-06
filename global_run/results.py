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
        bias = float(ex[0].split('_')[-1].split('.nc')[0])
        print(ex_mod.length_m_ts())

    #except:
    return pd.Series({'rgi':gdir.rgi_id})


if __name__ == '__main__':
    cfg.initialize()

    ON_CLUSTER = False
    REGION = '02'

    # Local paths
    if ON_CLUSTER:
        OUT_DIR = os.environ.get("OUTDIR")
        cfg.PATHS['working_dir'] = OUT_DIR
    else:
        WORKING_DIR = '/home/juliaeis/Dokumente/OGGM/work_dir/reconstruction/leclercq_diff/temp_0'
        cfg.PATHS['working_dir'] = WORKING_DIR
        utils.mkdir(WORKING_DIR, reset=False)

    cfg.PATHS['plot_dir'] = os.path.join(cfg.PATHS['working_dir'], 'plots')
    utils.mkdir(cfg.PATHS['plot_dir'], reset=False)

    for dir in os.listdir(OUT_DIR):
        if dir.startswith('reg1'):
            REGION = dir.split('reg')[-1].split('-')[0].zfill(2)

            # RGI file
            path = utils.get_rgi_region_file(REGION, version='61')
            rgidf = gpd.read_file(path)
            rgidf = rgidf.sort_values('Area', ascending=False)

            # exclude non-landterminating glaciers
            rgidf = rgidf[rgidf.TermType == 0]
            rgidf = rgidf[rgidf.Connect != 2]

            gdirs = workflow.init_glacier_regions(rgidf)
            df = read_results(gdirs)