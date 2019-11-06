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



if __name__ == '__main__':
    cfg.initialize()

    ON_CLUSTER = False
    REGION = '01'

    # Local paths
    if ON_CLUSTER:
        WORKING_DIR = os.environ.get("S_WORKDIR")
        cfg.PATHS['working_dir'] = WORKING_DIR
        OUT_DIR = os.environ.get("OUTDIR")

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
            for gdir in gdirs:
                print(gdir.dir)