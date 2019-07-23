
import os
import sys
from copy import deepcopy
sys.path.append('../')
from initialization.core import *
from plots_paper import *

import matplotlib.pyplot as plt
import salem
import geopandas as gpd
from oggm import cfg, workflow, utils
from oggm.core.flowline import FluxBasedModel
pd.options.mode.chained_assignment = None
import time


def experiments(gdirs):
    """
    creates searched and observed glacier to test the method, need only to
    be run once

    :param gdirs: list of oggm.GlacierDirectories
    :return:
    """
    df = pd.DataFrame()
    pool = mp.Pool()
    list = pool.map(_run_parallel_experiment, gdirs)
    pool.close()
    pool.join()
    df = df.append(list, ignore_index=True)
    return df

def _run_parallel_experiment(gdir):
    """
    Creates the synthetic experiment for one glacier. model_run_experiment.nc
    will be saved in working directory.
    """
    series = pd.Series({'rgi_id':gdir.rgi_id})

    for temp_bias in np.arange(-2.5,1.25,0.25):
        try:
            rp = gdir.get_filepath('model_run', filesuffix='_experiment_' + str(temp_bias))
            model = FileModel(rp)
        except:

            try:
                fls = gdir.read_pickle('model_flowlines')
                # try to run random climate with temperature bias -1
                model = tasks.run_random_climate(gdir, nyears=600, y0=1850, bias=0, seed=1,
                                                 temperature_bias=temp_bias,
                                                 init_model_fls=fls)
                # construct observed glacier, previous glacier will be run forward from
                # 1850 - 2000 with past climate file
                b = fls[-1].bed_h

                fls = deepcopy(model.fls)
                model = tasks.run_from_climate_data(gdir, ys=1850, ye=2000, init_model_fls=fls,
                                            output_filesuffix='_experiment_'+str(temp_bias))
            except:
                model=None

        series.at[str(temp_bias)]=model
    return series


if __name__ == '__main__':
    cfg.initialize()

    ON_CLUSTER = False

    # Local paths
    if ON_CLUSTER:
        WORKING_DIR = os.environ.get("S_WORKDIR")
        cfg.PATHS['working_dir'] = WORKING_DIR
        job_nr = int(os.environ.get('I'))
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
    #rgidf = rgidf.loc[rgidf.RGIId.isin(['RGI60-11.00013','RGI60-11.00062'])]

    # sort for efficient using
    rgidf = rgidf.sort_values('Area', ascending=False)

    gdirs = workflow.init_glacier_regions(rgidf)

    preprocessing(gdirs)

    # experiments
    df =experiments(gdirs).set_index('rgi_id')
    df.to_pickle(os.path.join(cfg.PATHS['working_dir'],'experiments.py.pkl'))
