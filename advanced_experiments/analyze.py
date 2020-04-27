import os
import sys
from copy import deepcopy
sys.path.append('../../')
from  initialization.core import *
from paper.plots_paper import *
from mpl_toolkits.axes_grid.anchored_artists import AnchoredText
import geoviews as gv
import geoviews.tile_sources as gts

matplotlib.use('TkAgg')
import holoviews as hv
hv.extension('matplotlib')
import matplotlib
import matplotlib.pyplot as plt
import time
import geopandas as gpd
from oggm import cfg, workflow, utils, graphics
from oggm.core.flowline import FluxBasedModel
pd.options.mode.chained_assignment = None

mpl.rcParams['axes.linewidth'] = 3
mpl.rcParams['xtick.major.width'] = 3
mpl.rcParams['ytick.major.width'] = 3
mpl.rcParams['font.size'] =25 #35
mpl.rcParams['font.weight'] = 'medium'
mpl.rcParams['axes.labelweight'] = 'medium'
mpl.rcParams['legend.fontsize'] = 40 #30
mpl.rcParams['lines.linewidth'] = 3 #5

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
        return deepcopy(quant_df.iloc[index].model), quant_df.iloc[quant_df.length.idxmin()].model, quant_df.iloc[quant_df.length.idxmax()].model

    except:

        return deepcopy(df.iloc[df.fitness.idxmin()].model), None, None


if __name__ == '__main__':

    cfg.initialize()

    ON_CLUSTER = False

    # Local paths
    if ON_CLUSTER:
        WORKING_DIR = os.environ.get("S_WORKDIR")
        cfg.PATHS['working_dir'] = WORKING_DIR
        job_nr = int(os.environ.get('I'))
    else:
        WORKING_DIR = '/home/juliaeis/Dokumente/OGGM/work_dir/reconstruction/paper_correction'
        cfg.PATHS['working_dir'] = WORKING_DIR
        utils.mkdir(WORKING_DIR, reset=False)

    print(WORKING_DIR)
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
    # cfg.PARAMS['baseline_y0'] = 1850
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
    rgidf = rgidf[rgidf.RGIId == 'RGI60-11.02400']

    gdirs = workflow.init_glacier_regions(rgidf)
    #synthetic_experiments_parallel(gdirs)
    model_df = pd.DataFrame()
    for gdir in gdirs:

        rp = gdir.get_filepath('model_run', filesuffix='_experiment')
        ex_mod = FileModel(rp)

        df = find_possible_glaciers(gdir, 1850, 2000, 200)
        df['fitness'] = df.fitness / 125


        plot_candidates(gdir, df, 1850, 'step1', cfg.PATHS['plot_dir'])
        plot_candidates(gdir, df, 1850, 'step2', cfg.PATHS['plot_dir'])
        plot_candidates(gdir, df, 1850, 'step3', cfg.PATHS['plot_dir'])


    '''
    df = gpd.GeoDataFrame()
    l = [os.path.join(WORKING_DIR,p) for p in os.listdir(WORKING_DIR)
         if p.endswith('.pkl')]

    for ele in l:
        df = pd.concat([df,pd.read_pickle(ele)],ignore_index=True)

    plt.figure()
    df.area_diff.plot.hist(bins=20)
    plt.xlabel('Area')
    plt.figure()
    df.fitness = df.fitness/125

    df.fitness.plot.hist(bins=50)
    plt.xlabel('Fitness')
    plt.figure()
    df.bias.plot.hist(bins=20)
    plt.xlabel('bias')

    df = df.dropna()
    df = df.set_index('rgi_id')

    print(df.sort_values('fitness').fitness)

    rgidf = rgidf.set_index('RGIId')

    for id in df.index:
        df.loc[id,'geometry'] = rgidf.loc[id,'geometry']
    df = df.set_geometry('geometry')
    df.fitness = df.fitness/125

    layout = gv.Layout([ts.relabel(name) for name, ts in gts.tile_sources.items()])
    layout.opts('WMTS', xaxis=None, yaxis=None)
    (gv.Polygons(df))
    plt.show()
    '''


