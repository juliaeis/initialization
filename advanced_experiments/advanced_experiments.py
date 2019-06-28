import os
import sys
from copy import deepcopy
sys.path.append('../')
sys.path.append('../../')
from  initialization.core import *
from paper.plots_paper import *
from mpl_toolkits.axes_grid.anchored_artists import AnchoredText

import matplotlib
import matplotlib.pyplot as plt
import time
import geopandas as gpd
from oggm import cfg, workflow, utils
from oggm.core.flowline import FluxBasedModel
pd.options.mode.chained_assignment = None

mpl.rcParams['axes.linewidth'] = 3
mpl.rcParams['xtick.major.width'] = 3
mpl.rcParams['ytick.major.width'] = 3
mpl.rcParams['font.size'] =35
mpl.rcParams['font.weight'] = 'medium'
mpl.rcParams['axes.labelweight'] = 'medium'
mpl.rcParams['legend.fontsize'] = 40 #30
mpl.rcParams['lines.linewidth'] = 5

def add_at(ax, t, loc=2):
    fp = dict(size=30)
    _at = AnchoredText(t, loc=loc, prop=fp ,borderpad=0)
    _at.patch.set_linewidth(3)
    ax.add_artist(_at)
    return _at

def fitness_function(model1, model2):
    """
    calculates the objective value (difference in geometry)
    :param model1: oggm.flowline.FluxBasedModel
    :param model2: oggm.flowline.FluxBasedModel
    :return:       float
    """

    model1 = model1
    model2 = model2
    model2.run_until(0)
    model1.run_until(2000)

    fls1 = model1.fls
    fls2 = model2.fls

    fitness = 0
    m = 0
    for i in range(len(model1.fls)):
        fitness = fitness + np.sum(
            abs(fls1[i].surface_h - fls2[i].surface_h)**2) + \
                    np.sum(abs(fls1[i].widths - fls2[i].widths)**2)
        m = m + fls1[i].nx
    fitness = fitness / m

    return fitness

def plot_compare_median(gdir, df, ex_mod, mod,plot_dir):

    v = pd.DataFrame()


    #fig, (ax1,ax2) = plt.subplots(2,1, sharex=True, figsize=(10,10))

    fig = plt.figure(figsize=(27, 18))
    grid = plt.GridSpec(2, 3, hspace=0.3, wspace=0.35)
    ax2 = plt.subplot(grid[1,:2])
    ax1 = plt.subplot(grid[0,:2],sharex=ax2)
    ax3 = plt.subplot(grid[1,2])

    plt.suptitle(gdir.rgi_id + ': Untersulzbach Kees')
    ax1.plot(2000, mod.length_m / 1000, 'o', color='C1',
             label=r'OGGM$_{init}$')
    ax1.plot(ex_mod.length_m_ts() / 1000, 'k:', label='experiment')

    for i in df.index:
        model = deepcopy(df.loc[i, 'model'])
        #model.run_until(2000)
        v = v.append(model.length_m_ts()/1000, ignore_index=True)
        v.loc[v.index[-1],'fitness'] = df.loc[i,'fitness']
        v.loc[v.index[-1], 'fls_fitness'] = df.loc[i, 'fls_fitness']


    # acceptable states
    acc1 = v[v.fitness<1]
    acc2 = v[v.fls_fitness < 1]

    ax1.fill_between(model.length_m_ts().index, acc1.iloc[:,:-2].min().values,acc1.iloc[:,:-2].max().values, color='grey', alpha=0.3, label='accepted')
    ax2.fill_between(model.length_m_ts().index, acc2.iloc[:,:-2].min().values,acc2.iloc[:,:-2].max().values, color='grey', alpha=0.3, label='accepted')

    # 5th percentile
    acc3 = acc1[acc1.fitness < acc1.fitness.quantile(0.05)]
    acc4 = acc2[acc2.fls_fitness < acc2.fls_fitness.quantile(0.05)]

    ax1.fill_between(model.length_m_ts().index, acc3.iloc[:, :-2].min().values,
                     acc3.iloc[:, :-2].max().values, color='C0', alpha=0.5,
                     label='best 5%')
    ax2.fill_between(model.length_m_ts().index, acc4.iloc[:, :-2].min().values,
                     acc4.iloc[:, :-2].max().values, color='C0', alpha=0.5, label='best 5%')

    # median state
    acc3 = acc3.sort_values(2000)
    acc4 = acc4.sort_values(2000)

    idx1 = acc3.index[int(len(acc3)/2)]
    idx2 = acc4.index[int(len(acc4)/2)]

    ax1.plot(model.length_m_ts().index, acc3.loc[idx1][:-2].values,color='C0',label='median')
    ax2.plot(model.length_m_ts().index, acc4.loc[idx2][:-2].values,color='C0', label='median')
    ax1.plot(ex_mod.length_m_ts() / 1000, 'k:',label='')

    x = (np.arange(ex_mod.fls[-1].nx) * ex_mod.fls[-1].dx * ex_mod.fls[-1].map_dx) / 1000
    n = len(x[x<(ex_mod.length_m_ts()[2000]/1000)+0.5])
    ex_mod.run_until(2000)
    ax3.plot(x[:n], mod.fls[-1].surface_h[:n] / 1000, color='C1')
    ax3.plot(x[:n], ex_mod.fls[-1].surface_h[:n]/1000,'k:')
    ax3.plot(x[:n], df.loc[idx1].model.fls[-1].bed_h[:n]/1000, 'k')

    # print observation
    gdir.rgi_id = 'RGI50-' + gdir.rgi_id.split('-')[-1]
    df_l = gdir.get_ref_length_data()

    df_l = df_l.loc[1848:2000]['dL']
    df_l = df_l - df_l.iloc[-1]
    df_l = (df_l + mod.length_m_ts()[0]) / 1000
    df_l.plot(ax=ax1, color='C3', label='Leclercq')
    df_l.plot(ax=ax2, color='C3', label='Leclercq')

    ax1.plot(2000, mod.length_m / 1000, 'o', markersize=12, color='C1')
    ax2.plot(2000, mod.length_m / 1000, 'o', markersize=12, color='C1')

    add_at(ax1, r"a", loc=1)
    add_at(ax2, r"b", loc=1)
    add_at(ax3, r"c", loc=1)

    ax1.set_ylabel('Length (km)')
    ax2.set_ylabel('Length (km)')
    ax2.set_xlabel('Time (years)')
    ax3.set_ylabel('Altitude (km)')
    ax3.set_xlabel('Distance along flowline (km)')
    ax3.set_title('t=2000')

    ax1.legend(bbox_to_anchor=(1.04,1), loc="upper left",prop={'size': 25})
    ax2.set_ylim(ax1.get_ylim())
    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir,'example_1b.png'))
    plt.show()


def plot_fitness_values(gdir, df, ex_mod, mod,plot_dir):

    x = (np.arange(ex_mod.fls[-1].nx) * ex_mod.fls[-1].dx * ex_mod.fls[-1].map_dx)/1000
    fig = plt.figure(figsize=(30,18))
    grid = plt.GridSpec(2, 3, hspace=0.3, wspace=0.5)
    ax2 = plt.subplot(grid[1,:2])
    ax1 = plt.subplot(grid[0,:2],sharex=ax2)
    ax3 = plt.subplot(grid[1,2])

    plt.suptitle(gdir.rgi_id + ': '+gdir.name)
    plt.suptitle(gdir.rgi_id + ': ' + 'Unteraargletscher')
    ax2.plot(2000, mod.length_m / 1000, 'o', color='C1',
             label=r'OGGM$_{init}$')
    ax1.plot(2000, mod.length_m / 1000, 'o', color='C1',
             label=r'OGGM$_{init}$')
    ax1.plot(ex_mod.length_m_ts() / 1000, 'r:', label='experiment',zorder=3)
    ax2.plot(ex_mod.length_m_ts() / 1000, 'r:', label='experiment', zorder=3)
    max = df.loc[df.volume.idxmax(),'model']
    n = len(x[x < (max.length_m_ts()[2000] / 1000) + 0.5])

    norm = mpl.colors.LogNorm(vmin=0.001, vmax=10)
    cmap = matplotlib.cm.get_cmap('viridis')

    df = df.sort_values(by='fitness',ascending=False)

    for i in df.index:
        model = deepcopy(df.loc[i, 'model'])
        model.run_until(2000)
        c = cmap(norm(df.loc[i,'fitness']))
        l = model.length_m_ts()/1000
        ax1.plot(l, color=c,label='')


    df = df.sort_values(by='fls_fitness', ascending=False)
    for i in df.index:
        model = deepcopy(df.loc[i, 'model'])
        model.run_until(2000)
        c = cmap(norm(df.loc[i, 'fls_fitness']))
        l = model.length_m_ts()/1000
        ax2.plot(model.length_m_ts()/1000, color=c,label='')


    try:
        # print observation
        gdir.rgi_id = 'RGI50-' + gdir.rgi_id.split('-')[-1]
        df_l = gdir.get_ref_length_data()
        df_l = df_l.loc[1848:2000]['dL']
        df_l = df_l - df_l.iloc[-1]
        df_l = (df_l + mod.length_m_ts()[0])/1000
        df_l.plot(ax=ax1, color='C3', label='Leclercq')
        df_l.plot(ax=ax2, color='C3', label='Leclercq')

    except:
        pass

    ax1.plot(2000, mod.length_m / 1000, 'o', markersize=12, color='C1')
    ax2.plot(2000, mod.length_m / 1000, 'o', markersize=12, color='C1')
    e=int(np.max([ex_mod.length_m_ts()[2000],mod.length_m_ts()[0]])/1000 +0.5)
    n = len(x[x<e])
    ex_mod.run_until(2000)
    ax3.plot(x[:n], mod.fls[-1].surface_h[:n] / 1000, color='C1')
    ax3.plot(x[:n], ex_mod.fls[-1].surface_h[:n] / 1000, 'r:')
    ax3.plot(x[:n], ex_mod.fls[-1].bed_h[:n] / 1000, 'k')

    # add colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cax, kw = mpl.colorbar.make_axes([ax1, ax2])
    cbar = fig.colorbar(sm,cax=cax, extend='both')
    cbar.ax.tick_params(labelsize=30)
    cbar.set_label('Fitness value', fontsize=30)

    add_at(ax1, r"a", loc=4)
    add_at(ax2, r"b", loc=4)
    add_at(ax3, r"c", loc=3)

    ax1.set_ylabel('Length (km)')
    ax2.set_ylabel('Length (km)')
    ax2.set_xlabel('Time (years)')
    ax3.set_ylabel('Altitude (km)')
    ax3.set_xlabel('Distance along flowline (km)')
    ax3.set_title('t=2000')

    ax1.legend(bbox_to_anchor=(1.5, 1), loc="upper left", prop={'size': 30})
    ax2.set_ylim(ax1.get_ylim())
    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir,'example_1.png'))
    plt.show()


def plot_compare_2_exp(diff,exp_df, cfg):
    fig, (ax3, ax2, ax1) = plt.subplots(1, 3, sharey=True, figsize=(23, 15))
    diff.volume.plot.hist(bins=30, ax=ax1)
    textstr = '\n'.join((
        r'$\min=%.2f  km^3$' % (diff.volume.min()),
        r'$\mathrm{mean}=%.2f  km^3$' % (diff.volume.mean(),),
        r'$max=%.2f  km^3$' % (diff.volume.max()),))

    # these are matplotlib.patch.Patch properties
    props = dict(boxstyle='round', facecolor='grey', alpha=0.3)

    # place a text box in upper left in axes coords
    ax1.text(0.05, 0.95, textstr, transform=ax1.transAxes, fontsize=25,
             verticalalignment='top', bbox=props)
    # ax.set_yscale('log')
    ax1.set_xlabel(r'Difference in volume (km$^3$)')
    ax1.set_ylabel('Frequency')

    diff.area.plot.hist(bins=30, ax=ax2)
    # ax.set_yscale('log')
    textstr = '\n'.join((
        r'$\min=%.2f  km^2$' % (diff.area.min()),
        r'$\mathrm{mean}=%.2f  km^2$' % (diff.area.mean(),),
        r'$max=%.2f  km^2$' % (diff.area.max()),))

    # these are matplotlib.patch.Patch properties
    props = dict(boxstyle='round', facecolor='grey', alpha=0.3)

    # place a text box in upper left in axes coords
    ax2.text(0.05, 0.95, textstr, transform=ax2.transAxes, fontsize=25,
             verticalalignment='top', bbox=props)
    ax2.set_xlabel(r'Difference in area (km$^2$)')
    ax2.set_ylabel('Frequency')

    diff.length.plot.hist(bins=30, ax=ax3)
    textstr = '\n'.join((
        r'$\min=%.2f  km$' % (diff.length.min()),
        r'$\mathrm{mean}=%.2f  km$' % (diff.length.mean(),),
        r'$max=%.2f  km$' % (diff.length.max()),))

    # these are matplotlib.patch.Patch properties
    props = dict(boxstyle='round', facecolor='grey', alpha=0.3)

    # place a text box in upper left in axes coords
    ax3.text(0.05, 0.95, textstr, transform=ax3.transAxes, fontsize=25,
             verticalalignment='top', bbox=props)
    ax3.set_xlabel(r'Difference in length (km)')
    ax3.set_ylabel('Frequency')
    # ax.set_yscale('log')
    plt.suptitle(' OGGM$_{init}$ - experiment')
    plt.savefig(os.path.join(cfg.PATHS['plot_dir'], 'fls_vs_exp.png'), dpi=200)
    plt.show()

    plt.figure(figsize=(20, 15))
    grid = plt.GridSpec(1, 3, wspace=0.4, hspace=0.3)
    ax2 = plt.subplot(grid[1:])
    ax1 = plt.subplot(grid[0], sharey=ax2)

    diff['fls_area'] = exp_df.OGGM.apply(lambda x: x.area_km2)
    diff.plot.scatter(y='fitness', x='fls_area', ax=ax2, alpha=0.9,s=30)
    diff[diff.fitness < 1].plot.scatter(y='fitness', x='fls_area', ax=ax2,
                                        color='r', label='accepted',s=30)
    ax2.legend(loc='best')
    N, bins, patches = ax1.hist(x=diff.fitness.values,
                                bins=range(int(diff.fitness.max()) + 1),
                                orientation='horizontal')
    ax1.invert_xaxis()
    textstr = '\n'.join((
        r'$\min=%.2f$' % (diff.fitness.min()),
        r'$\mathrm{median}=%.2f$' % (diff.fitness.median(),),
        r'$max=%.2f$' % (diff.fitness.max()),))

    # these are matplotlib.patch.Patch properties
    props = dict(boxstyle='round', facecolor='grey', alpha=0.3)

    # place a text box in upper left in axes coords
    ax1.text(0.05, 0.97, textstr, transform=ax1.transAxes, fontsize=25,
             verticalalignment='top', bbox=props)
    t = str(format(len(diff[diff.fitness < 1]) / len(diff) * 100,
                   '.2f')) + ' % of all tested glaciers'
    ax2.text(0.35, 0.85, t, transform=ax2.transAxes, fontsize=30,
             verticalalignment='top', color='r')
    ax1.set_ylabel(r'Fitness value ')
    ax1.set_xlabel('Frequency')

    ax2.set_xlabel(r'OGGM$_{init}$ area ($km^3$)')
    ax1.yaxis.tick_right()
    ax1.yaxis.set_label_position("right")
    for i in range(0, 1):
        patches[i].set_facecolor('r')
    # ax.set_yscale('log')
    plt.suptitle(' OGGM$_{init}$ - experiment')
    plt.savefig(os.path.join(cfg.PATHS['plot_dir'], 'fls_fitness.png'),
                dpi=200)

    plt.show()


def _run_experiment(gdir, bias):
    """
    Creates the synthetic experiment for one glacier. model_run_experiment.nc
    will be saved in working directory.
    """

    # check, if this experiment already exists
    try:
        rp = gdir.get_filepath('model_run', filesuffix='_advanced_experiment_' + str(bias))
        model = FileModel(rp)

    # otherwise create experiment
    except:

        fls = gdir.read_pickle('model_flowlines')
        # try to run random climate with temperature bias -1
        try:
            model = tasks.run_random_climate(gdir, nyears=400, y0=1866, bias=bias, seed=1,
                                             temperature_bias=-1,
                                             init_model_fls=fls)

            # construct observed glacier, previous glacier will be run forward from
            # 1850 - 2000 with past climate file

            fls = deepcopy(model.fls)
            model = tasks.run_from_climate_data(gdir, ys=1866, ye=2000, init_model_fls=fls,bias=bias,
                                        output_filesuffix='_advanced_experiment_'+str(bias))
        except:
            pass

    return model

def find_residual(gdir, a=-2000,b=2000):

    try:

        max_it = 15
        i = 0
        bounds = [a,b]

        df = pd.DataFrame()

        fls = gdir.read_pickle('model_flowlines')
        mod = FluxBasedModel(flowlines=fls)

        while i < max_it:
            bias = round((bounds[0] + bounds[1]) / 2,1)
            ex_mod2 = _run_experiment(gdir, bias)
            fit = fitness_function(ex_mod2,mod)
            df = df.append(pd.Series({'bias':bias,'fitness':fit}),ignore_index=True)

            if bounds[1]-bounds[0]<=1:
                break

            elif ex_mod2.area_km2 > mod.area_km2:
                bounds[0] = bias
            else:
                bounds[1] = bias
            i +=1

        # best bias found
        bias = df.iloc[df.fitness.idxmin()].bias
        rp = gdir.get_filepath('model_run', filesuffix='_advanced_experiment_' + str(bias))
        model = FileModel(rp)
        model.run_until(2000)

        rp = gdir.get_filepath('model_run', filesuffix='_advanced_experiment_' + str(0.0))
        ex_mod = FileModel(rp)
        ex_mod.run_until(2000)

        plt.figure(figsize=(15,10))
        plt.plot(model.fls[-1].surface_h,'r',label='best')
        plt.plot(mod.fls[-1].surface_h, 'orange', label='original')
        plt.plot(ex_mod.fls[-1].surface_h, 'r:', label='old experiment')
        plt.plot(model.fls[-1].bed_h,'k', label='bed')
        plt.legend()
        utils.mkdir(os.path.join(cfg.PATHS['plot_dir'],'bias_test'))
        plt.savefig(os.path.join(cfg.PATHS['plot_dir'],'bias_test',gdir.rgi_id+'.png'),dpi=200)
        diff = mod.area_km2 - model.area_km2_ts()[2000]
        model.reset_y0(1865)

        series = pd.Series({'rgi_id':gdir.rgi_id,'bias':bias,'iterations':i, 'fitness':fit, 'area_diff':diff, 'model':model})
    except:
        series =  pd.Series({'rgi_id':gdir.rgi_id})

    return series


if __name__ == '__main__':
    cfg.initialize()

    ON_CLUSTER = True

    # Local paths
    if ON_CLUSTER:
        WORKING_DIR = os.environ.get("S_WORKDIR")
        cfg.PATHS['working_dir'] = WORKING_DIR
    else:
        WORKING_DIR = '/home/juliaeis/Dokumente/OGGM/work_dir/reconstruction/advanced_experiments2/'
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
    cfg.PARAMS['baseline_y0'] = 1850
    cfg.PARAMS['prcp_scaling_factor'] = 1.75
    cfg.PARAMS['temp_all_liq'] = 2.0
    cfg.PARAMS['temp_default_gradient'] = -0.0065
    cfg.PARAMS['temp_melt']= -1.75
    cfg.PARAMS['temp_all_solid'] = 0.0


    # We use intersects
    db = utils.get_rgi_intersects_region_file(version='61', region='11')
    cfg.set_intersects_db(db)

    cfg.PARAMS['run_mb_calibration'] = False
    cfg.PARAMS['optimize_inversion_params'] = False

    # RGI file
    path = utils.get_rgi_region_file('11', version='61')
    rgidf = gpd.read_file(path)

    # sort for efficient using
    rgidf = rgidf.sort_values('Area', ascending=False)
    gdirs = workflow.init_glacier_regions(rgidf)

    preprocessing(gdirs)

    # experiments
    #synthetic_experiments_parallel(gdirs)

    t_0 = 1865
    t_e = 2000
    epsilon = 125

    exp_df = pd.DataFrame()

    pool = Pool()
    list = pool.map(find_residual, gdirs)
    pool.close()
    pool.join()

    exp_df = exp_df.append(list, ignore_index=True)
    exp_df.to_pickle(os.path.join(cfg.PATHS['working_dir'],'advanced_experiments.pkl'))

    #exp_df = pd.read_pickle(os.path.join(cfg.PATHS['working_dir'],'advanced_experiments.pkl'))

    '''
    for gdir in gdirs:
        print(gdir.rgi_id)
        try:
            rp = gdir.get_filepath('model_run', filesuffix='_experiment' )
            ex_mod = FileModel(rp)
            exp_df.loc[gdir.rgi_id, 'experiment'] = ex_mod

            rp = gdir.get_filepath('model_run', )
            mod = FileModel(rp)
            exp_df.loc[gdir.rgi_id,'OGGM'] = mod
        except:
            pass

    print(exp_df)
    exp_df.to_pickle(os.path.join(cfg.PATHS['working_dir'],'models.pkl'),compression='gzip')


    exp_df = pd.read_pickle(os.path.join(cfg.PATHS['working_dir'],'models.pkl'),compression='gzip')

    diff = pd.read_pickle(os.path.join(cfg.PATHS['working_dir'], 'diff.pkl'))


    diff['fitness'] = (exp_df[['experiment','OGGM']].apply(lambda x: fitness_function(*x), axis=1))/125

    diff['length'] = (exp_df.OGGM.apply(lambda x: x.length_m)-exp_df.experiment.apply(lambda x: x.length_m))/1000
    diff['area'] = exp_df.OGGM.apply(lambda x: x.area_km2)-exp_df.experiment.apply(lambda x: x.area_km2)
    diff['volume'] = exp_df.OGGM.apply(lambda x: x.volume_km3)-exp_df.experiment.apply(lambda x: x.volume_km3)

    plot_compare_2_exp(diff,exp_df,cfg)


    for gdir in gdirs:
        if gdir.rgi_id in diff.index:
            gdir.rgi_id = 'RGI50-' + gdir.rgi_id.split('-')[-1]
            try:
                df_l = gdir.get_ref_length_data()
                gdir.rgi_id = 'RGI60-' + gdir.rgi_id.split('-')[-1]
                diff.loc[gdir.rgi_id, 'dL'] = True

            except:
                gdir.rgi_id = 'RGI60-' + gdir.rgi_id.split('-')[-1]
                diff.loc[gdir.rgi_id, 'dL'] = False

    #diff.to_pickle(os.path.join(cfg.PATHS['working_dir'], 'diff.pkl'))

    diff = pd.read_pickle(os.path.join(cfg.PATHS['working_dir'], 'diff.pkl'))

    diff = diff.sort_values('fitness')
    print(diff[(diff.dL == True) & (diff.fitness>1)])
    print(diff.loc['RGI60-11.00002'])

    for gdir in gdirs:
        fls = gdir.read_pickle('model_flowlines')
        print(fls)

        rp = gdir.get_filepath('model_run', filesuffix='_experiment')

        ex_mod = FileModel(rp)
        ex_mod.reset_y0(1850)

        rp = gdir.get_filepath('model_run', )
        mod = FileModel(rp)
        mod.reset_y0(0)

        p = os.path.join(gdir.dir,'result1850.pkl')
        df = pd.read_pickle(p, compression='gzip')
        df = df.sort_values('fitness', ascending=False)

        df['fls_fitness'] = df.model.apply(fitness_function,
                                           model2=deepcopy(mod))

        df.fls_fitness = df.fls_fitness / 125
        df.fitness = df.fitness / 125

        #plot_compare_median(gdir, df, ex_mod, mod, cfg.PATHS['plot_dir'])
        plot_fitness_values(gdir, df, ex_mod, mod, cfg.PATHS['plot_dir'])

    '''
