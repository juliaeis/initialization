
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
from mpl_toolkits.axes_grid.anchored_artists import AnchoredText
import time


mpl.rcParams['axes.linewidth'] = 3
mpl.rcParams['xtick.major.width'] = 3
mpl.rcParams['ytick.major.width'] = 3
mpl.rcParams['font.size'] =20
mpl.rcParams['font.weight'] = 'medium'
mpl.rcParams['axes.labelweight'] = 'medium'
mpl.rcParams['legend.fontsize'] = 17 #30
mpl.rcParams['lines.linewidth'] = 3

def add_at(ax, t, loc=2):
    fp = dict(size=20)
    _at = AnchoredText(t, loc=loc, prop=fp ,borderpad=0)
    _at.patch.set_linewidth(3)
    ax.add_artist(_at)
    return _at

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

    #for temp_bias in np.arange(-2.5,1.25,0.25):
    for temp_bias in [0,-1]:
        try:
            rp = gdir.get_filepath('model_run', filesuffix='_experiment_' + str(temp_bias))
            model = FileModel(rp)
        except:

            try:
                fls = gdir.read_pickle('model_flowlines')
                # try to run random climate with temperature bias -1
                model = tasks.run_random_climate(gdir, nyears=600, y0=1850, bias=0, seed=1,
                                                 temperature_bias=temp_bias,
                                                 init_model_fls=fls, output_filesuffix='random_experiment_'+str(temp_bias))
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


class MidpointNormalize(mpl.colors.Normalize):
    def __init__(self, vmin=None, vmax=None, vcenter=None, clip=False):
        self.vcenter = vcenter
        mpl.colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.vcenter, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))



if __name__ == '__main__':
    cfg.initialize()

    ON_CLUSTER = False

    # Local paths
    if ON_CLUSTER:
        WORKING_DIR = os.environ.get("S_WORKDIR")
        cfg.PATHS['working_dir'] = WORKING_DIR
        job_nr = int(os.environ.get('I'))
    else:
        WORKING_DIR = '/home/juliaeis/Dokumente/OGGM/work_dir/reconstruction/600_paper_correction/experiment'
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
    rgidf = rgidf.loc[rgidf.RGIId.isin(['RGI60-11.00779'])]

    #sort for efficient using
    #rgidf = rgidf.sort_values('Area', ascending=False)

    gdirs = workflow.init_glacier_regions(rgidf)

    #preprocessing(gdirs)
    #df = experiments(gdirs)

    p = gdirs[0].get_filepath('model_run', filesuffix='random_experiment_0')
    ex_mod1 = FileModel(p)
    p = gdirs[0].get_filepath('model_run', filesuffix='random_experiment_-1')
    ex_mod2 = FileModel(p)
    p = gdirs[0].get_filepath('model_run', filesuffix='_experiment_0')
    ex_mod3 = FileModel(p)
    p = gdirs[0].get_filepath('model_run', filesuffix='_experiment_-1')
    ex_mod4 = FileModel(p)

    x = np.arange(ex_mod4.fls[-1].nx) * ex_mod4.fls[-1].dx * ex_mod4.fls[-1].map_dx
    x = x/1000
    fig = plt.figure(figsize=(15,10))
    grid1 = plt.GridSpec(1, 2, bottom=0.57, top=0.91)
    ax3 = plt.subplot(grid1[0, 0],)
    ax4 = plt.subplot(grid1[0, 1], sharey=ax3)
    grid2 = plt.GridSpec(1, 2, bottom=0.09, top=0.43)
    grid2.update( wspace=0, hspace=0)
    ax1 = plt.subplot(grid2[0, 0])
    ax2 = plt.subplot(grid2[0, 1], sharey=ax1)

    ax3.plot(x, deepcopy(ex_mod4.fls[-1].surface_h), color='C3',label=r'$z_{1850}^{exp}$')
    ex_mod4.run_until(2000)
    ax4.plot(x, ex_mod4.fls[-1].surface_h, color='C3',label=r'$z_{2000}^{exp}$')
    ax3.plot(x, ex_mod4.fls[-1].bed_h,'k', label='b')
    ax4.plot(x, ex_mod4.fls[-1].bed_h, 'k', label='b')
    ax3.set_xlabel('Distance along the main flowline (km)')
    ax4.set_xlabel('Distance along the main flowline (km)')

    ax3.set_ylabel('Altitude (m)')
    ax4.set_ylabel('Altitude (m)')

    ex_mod2.volume_km3_ts().plot(ax=ax1,color='C3')
    #ax1.axhline(y=ex_mod1.volume_km3_ts()[0],color='grey', linestyle=':', label='volume in 2000')
    ax1.set_xticks([0,100,200,300,400,500])
    ax1.set_xlabel('Time (years)')
    ax1.set_ylabel(r'Volume ($km^3$)')
    ax1.set_title('Random climate forcing', fontsize=20)

    ex_mod4.volume_km3_ts().plot(ax=ax2,color='C3', label=r'$s_{1850-2000}^{exp}$')
    #ax2.axhline(y=ex_mod1.volume_km3_ts()[0], color='grey', linestyle=':',
    #            label='volume in 2000')
    ax2.set_xlabel('Time (years)')
    ax2.set_title('HISTALP 1850-2000', fontsize=20)
    ax3.set_yticks([2000,2500,3000,3500])

    plt.suptitle('RGI60-11.00779: Guslarferner',fontsize=20)
    add_at(ax3, r"a", loc=3)
    add_at(ax4, r"b", loc=3)
    add_at(ax1, r"c", loc=2)
    add_at(ax2, r"d", loc=2)

    ax2.legend(loc='best')
    ax3.legend(loc='best')
    ax4.legend(loc='best')

    p = '/home/juliaeis/Dropbox/Apps/Overleaf/reconstruction_paper/plots/experiment2.pdf'
    plt.savefig(p,dpi=300)
    plt.show()




    '''
    a = pd.read_pickle(os.path.join(cfg.PATHS['working_dir'], 'area_experiments.pkl'))
    norm = MidpointNormalize(vmin=-2.5, vmax=1.25, vcenter=0.25)

    fig = plt.figure(figsize=(12, 8))
    grid = plt.GridSpec(2,4, hspace=0.2, wspace=0.2)
    ax = plt.subplot(grid[:,:3])

    a.plot(ax=ax, kind='scatter', x='1850', y='diff', c='bias', cmap='RdBu_r', s=a['size'].values, edgecolors='k', zorder=5,
           norm=norm, colorbar=False)
    ax.axhline(color='k')
    ax.axvline(x=4475, color='grey', linestyle=':', label='Zemp et al.(2006)')
    for i in a[['s', 'class']].sort_values(by='s',
                                           ascending=False).drop_duplicates().index:
        plt.scatter([], [], c='k', alpha=0.3, s=a.loc[i, 's'],
                    label=str(a.loc[i, 'class']) + ' - ' + str(
                        (a.loc[i, 'class']) + 100))
    handles, labels = ax.get_legend_handles_labels()
    handles = np.append(handles[1::], handles[0])
    labels = np.append(labels[1::], labels[0])

    lg = plt.legend(handles, labels, scatterpoints=1, frameon=True,
                    labelspacing=1, title=r'sample size', bbox_to_anchor=(1.34, 1))
    lg.get_title().set_fontsize(20)

    # add colorbar

    sm = plt.cm.ScalarMappable(cmap='RdBu_r', norm=norm)
    sm.set_array([])
    labels = np.append(a.bias.values[::2], [1.5])

    cbar = fig.colorbar(sm, boundaries=np.append(a.bias.values, [1.25]),
                        label=' Temperature bias')
    loc = labels - 0.375
    cbar.set_ticks(loc)
    cbar.set_ticklabels(labels)

    ax.set_ylabel(r'Difference in 2000 to RGI (km$^2$)')
    ax.set_xlabel('Glacierized area in 1850 (km$^2$)')
    plt.savefig(os.path.join(cfg.PATHS['working_dir'], 'bias.pdf'), dpi=300)
    plt.savefig(os.path.join(cfg.PATHS['working_dir'], 'bias.png'), dpi=300)
    plt.show()
    '''

