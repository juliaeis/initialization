import os
from functools import partial
from pylab import *
from oggm.core.flowline import FluxBasedModel, FileModel
import matplotlib.patches as mpatches
from oggm import graphics, tasks,utils
from matplotlib import cm
import xarray as xr
import pandas as pd
from multiprocessing import Pool
from copy import deepcopy
FlowlineModel = partial(FluxBasedModel, inplace=False)
from mpl_toolkits.axes_grid.anchored_artists import AnchoredText
import matplotlib as mpl
from matplotlib.collections import LineCollection
from matplotlib.legend_handler import HandlerLineCollection
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MaxNLocator
import random
pd.options.mode.chained_assignment = None
warnings.filterwarnings("ignore")

mpl.rcParams['axes.linewidth'] = 3
mpl.rcParams['xtick.major.width'] = 3
mpl.rcParams['ytick.major.width'] = 3
mpl.rcParams['ytick.minor.width'] = 2
mpl.rcParams['font.size'] =25
mpl.rcParams['font.weight'] = 'medium'
mpl.rcParams['axes.labelweight'] = 'medium'
mpl.rcParams['legend.fontsize'] = 30 #30
mpl.rcParams['lines.linewidth'] = 3

from matplotlib.ticker import Locator

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



class MinorSymLogLocator(Locator):
    """
    Dynamically find minor tick positions based on the positions of
    major ticks for a symlog scaling.
    """
    def __init__(self, linthresh):
        """
        Ticks will be placed between the major ticks.
        The placement is linear for x between -linthresh and linthresh,
        otherwise its logarithmically
        """
        self.linthresh = linthresh

    def __call__(self):
        'Return the locations of the ticks'
        majorlocs = self.axis.get_majorticklocs()

        # iterate through minor locs
        minorlocs = []

        # handle the lowest part
        for i in range(1, len(majorlocs)):
            majorstep = majorlocs[i] - majorlocs[i-1]
            if abs(majorlocs[i-1] + majorstep/2) < self.linthresh:
                ndivs = 10
            else:
                ndivs = 9
            minorstep = majorstep/ndivs
            locs = np.arange(majorlocs[i-1], majorlocs[i], minorstep)[1:]
            minorlocs.extend(locs)

        return self.raise_if_exceeds(np.array(minorlocs))

    def tick_values(self, vmin, vmax):
        raise NotImplementedError('Cannot get tick locations for a '
                                  '%s type.' % type(self))


def add_at(ax, t, loc=2):
    fp = dict(size=20)
    _at = AnchoredText(t, loc=loc, prop=fp ,borderpad=0)
    _at.patch.set_linewidth(3)
    ax.add_artist(_at)
    return _at


class HandlerColorLineCollection(HandlerLineCollection):
    def create_artists(self, legend, artist, xdescent, ydescent,
                       width, height, fontsize, trans):
        x = np.linspace(0, width, self.get_numpoints(legend) + 1)
        y = np.zeros(
            self.get_numpoints(legend) + 1) + height / 2. - ydescent
        points = np.array([x, y]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)
        lc = LineCollection(segments, cmap=artist.cmap,
                            transform=trans)
        lc.set_array(x)
        lc.set_linewidth(artist.get_linewidth())
        return [lc]


def plot_experiment(gdir, ex_mod, ys, plot_dir):

    x = np.arange(ex_mod.fls[-1].nx) * ex_mod.fls[-1].dx * \
        ex_mod.fls[-1].map_dx

    fig = plt.figure(figsize=(15, 14))
    grid = plt.GridSpec(2, 1, hspace=0.2, wspace=0.2)
    ax1 = plt.subplot(grid[0, 0])
    ax2 = plt.subplot(grid[1, 0], sharex=ax1)

    if gdir.name != '':
        ax1.set_title(gdir.rgi_id+':'+gdir.name)
    else:
        ax1.set_title(gdir.rgi_id)

    # plot experiments.py, run until ys
    ex_mod = deepcopy(ex_mod)
    ex_mod.reset_y0(ys)
    ex_mod.run_until(ys)
    i = np.where(ex_mod.fls[-1].thick > 0)[0][-1] + 10
    i = ex_mod.fls[-1].nx
    ax1.plot(x[:i], ex_mod.fls[-1].surface_h[:i], 'k:',
             label=r'$z_{'+str(ys)+'}^{exp}$', linewidth=3)
    ax1.plot(x[:i], ex_mod.fls[-1].bed_h[:i], 'k', label=r'$b$', linewidth=3)

    ex_mod.run_until(2000)

    ax2.plot(x[:i], ex_mod.fls[-1].surface_h[:i], 'k:',
             label=r'$z_{2000}^{exp = obs} $', linewidth=3)
    ax2.plot(x[:i], ex_mod.fls[-1].bed_h[:i], 'k', label=r'$b$', linewidth=3)

    # add figure names and legends
    add_at(ax1, r"a", loc=3)
    add_at(ax2, r"b", loc=3)

    ax1.legend(loc=1)
    ax2.legend(loc=1)

    ax1.set_ylabel('Altitude (m)')
    ax1.set_xlabel('Distance along the main flowline (m)')
    ax2.set_ylabel('Altitude (m)')
    ax2.set_xlabel('Distance along the main flowline (m)')

    ax1.tick_params(axis='both', which='major')
    ax2.tick_params(axis='both', which='major')

    plot_dir = os.path.join(plot_dir, '00_experiment')
    utils.mkdir(plot_dir)
    fig_name = 'experiment_'+str(ys)+'_'+gdir.rgi_id
    plt.savefig(os.path.join(plot_dir, fig_name+'.pdf'), dpi=300)
    plt.savefig(os.path.join(plot_dir, fig_name+'.png'), dpi=300)
    #plt.show()
    plt.close()


def plot_candidates(gdir, df, yr, step, plot_dir):
    plot_dir = os.path.join(plot_dir, '06_candidates')
    utils.mkdir(plot_dir, reset=False)
    fig, ax = plt.subplots(figsize=(10, 10))

    for file in os.listdir(os.path.join(gdir.dir, str(yr))):

        if file.startswith('model_run'+str(yr)+'_random'):
            suffix = file.split('model_run')[1].split('.nc')[0]
            rp = os.path.join(gdir.dir, str(yr), 'model_run'+suffix+'.nc')
            try:
                fmod = FileModel(rp)
                fmod.volume_km3_ts().plot(ax=ax, color='grey', label='',
                                         zorder=1)

            except:
                pass

    # last one again for labeling
    label = r'temperature bias $\in [$' + str(
        df['temp_bias'].min()) + ',' + str(df['temp_bias'].max()) + '$]$'
    df.time = df.time.apply(lambda x: int(x))
    t_eq = df['time'].sort_values().iloc[0]

    df['Fitness value'] = df.fitness

    plt.title(gdir.rgi_id)

    if step == 'step1':
        fmod.volume_km3_ts().plot(ax=ax, color='grey', label=label, zorder=1)
        plt.legend(loc=0, fontsize=28)
        plt.xlabel('Time (years)')
        plt.ylabel(r'Volume $(km^3)$')
        plt.savefig(os.path.join(plot_dir, 'candidates1_' + str(yr) + '_' +
                                 str(gdir.rgi_id) + '.png'), dpi=300)
    elif step == 'step2':
        ax.axvline(x=int(t_eq), color='k', zorder=1, label=r'$t_{stag}$')
        fmod.volume_km3_ts().plot(ax=ax, color='grey', label='', zorder=1)
        # black points
        df.plot.scatter(x='time', y='volume', ax=ax, color='k',
                        label='candidates', s=250, zorder=2)
        plt.legend(loc=0, fontsize=27.5)
        plt.xlabel('Time (years)')
        plt.ylabel(r'Volume $(km^3)$')
        plt.xlim((int(t_eq)-10,605))
        plt.savefig(os.path.join(plot_dir, 'candidates2_' + str(yr) + '_' +
                                 str(gdir.rgi_id) + '.png'), dpi=300)
    elif step == 'step3':
        fmod.volume_km3_ts().plot(ax=ax, color='grey', label=None, zorder=1)
        ax.axvline(x=int(t_eq), color='k', zorder=1)

        cmap = matplotlib.cm.get_cmap('viridis')
        norm = mpl.colors.LogNorm(vmin=0.01 / 125, vmax=10)

        im=df.plot.scatter(x='time', y='volume', ax=ax, c='Fitness value',
                        colormap='viridis',
                        norm=mpl.colors.LogNorm(vmin=0.01/125, vmax=10, clip=True),
                        s=250, edgecolors='k', zorder=2,colorbar=False)
        # plot again points with objective == 0, without norm
        if len(df[df.fitness == 0]) > 0:
            df[df.fitness == 0].plot.scatter(x='time', y='volume', ax=ax,
                                               c=cmap(0), s=250,
                                               edgecolors='k', zorder=2, colorbar=False)

        plt.xlim(int(t_eq)-10,605)
        plt.xlabel('Time (years)')
        plt.ylabel(r'Volume $(km^3)$')

        # add colorbar
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        #cax, kw = mpl.colorbar.make_axes([ax1, ax2, ax3])
        cbar = fig.colorbar(sm, extend='both')
        cbar.ax.tick_params(labelsize=30)
        cbar.set_label('Fitness value', fontsize=30)

        plt.savefig(os.path.join(plot_dir, 'candidates3_' + str(yr) + '_' +
                                 str(gdir.rgi_id) + '.png'), dpi=300)
        plt.show()

    plt.close()

    plt.figure(figsize=(15, 14))
    plt.hist(df.volume.values, bins=20)
    plt.xlabel(r'Volume $(km^3)$')
    plt.ylabel(r'Frequency')
    plt.title(gdir.rgi_id)
    plt.savefig(os.path.join(plot_dir, 'hist_candidates' + str(yr) + '_' +
                             str(gdir.rgi_id) + '.png'), dpi=300)
    plt.close()


def plot_fitness_values(gdir, df, ex_mod, ys, plot_dir):

    plot_dir = os.path.join(plot_dir, '03_surface_by_fitness')
    utils.mkdir(plot_dir)
    x = (np.arange(ex_mod.fls[-1].nx) * ex_mod.fls[-1].dx * \
        ex_mod.fls[-1].map_dx)/1000
    fig = plt.figure(figsize=(25, 18))
    grid = plt.GridSpec(2, 2, hspace=0.2, wspace=0.2)
    ax1 = plt.subplot(grid[0, 0])
    ax2 = plt.subplot(grid[0, 1], sharey=ax1)
    ax3 = plt.subplot(grid[1, :])

    if gdir.name != '':
        plt.suptitle(gdir.rgi_id + ': ' + gdir.name, fontsize=30)
    elif gdir.rgi_id.endswith('779'):
        plt.suptitle(gdir.rgi_id + ': Guslarferner', fontsize=30)
    else:
        plt.suptitle(gdir.rgi_id, fontsize=30)

    norm = mpl.colors.LogNorm(vmin=0.01/125, vmax=10)
    cmap = matplotlib.cm.get_cmap('viridis')

    # df = df.sort_values('objective', ascending=False)
    df = df.sort_values('fitness', ascending=False)

    for i, model in df['model'].iteritems():
        model = deepcopy(model)
        model.reset_y0(ys)
        # color = cmap(norm(df.loc[i, 'objective']))
        color = cmap(norm(df.loc[i, 'fitness']))
        ax1.plot(x, deepcopy(model.fls[-1].surface_h), color=color,
                 label='')
        model.volume_km3_ts().plot(ax=ax3, color=[color], label='')
        model.run_until(2000)

        ax2.plot(x, model.fls[-1].surface_h, color=color, label='')

    # plot experiments.py
    ex_mod = deepcopy(ex_mod)
    ex_mod.volume_km3_ts().plot(ax=ax3, color='red', linestyle=':',
                               linewidth=3,
                               label='')
    ex_mod.reset_y0(1850)
    ex_mod.run_until(ys)

    ax1.plot(x, ex_mod.fls[-1].surface_h, ':', color='red', label='',
             linewidth=3)
    ax1.plot(x, ex_mod.fls[-1].bed_h, 'k', label='', linewidth=3)

    ex_mod.run_until(2000)

    ax2.plot(x, ex_mod.fls[-1].surface_h, ':', color='red', label='',
             linewidth=3)
    ax2.plot(x, ex_mod.fls[-1].bed_h, 'k', label='',
             linewidth=3)

    # add colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cax, kw = mpl.colorbar.make_axes([ax1, ax2, ax3])
    cbar = fig.colorbar(sm, cax=cax, extend='both', **kw)
    cbar.ax.tick_params(labelsize=30)
    cbar.set_label('Fitness value', fontsize=30)

    # add figure names and x-/ylabels
    add_at(ax1, r"a", loc=3)
    add_at(ax2, r"b", loc=3)
    add_at(ax3, r"c", loc=3)

    ax1.set_ylabel('Altitude (m)', fontsize=30)
    ax1.set_xlabel('Distance along the main flowline (km)', fontsize=30)
    ax2.set_ylabel('Altitude (m)', fontsize=30)
    ax2.set_xlabel('Distance along the main flowline (km)', fontsize=30)
    ax3.set_ylabel(r'Volume ($km^3$)', fontsize=30)
    ax3.set_xlabel('Time (years)', fontsize=30)

    ax1.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax2.xaxis.set_major_locator(MaxNLocator(integer=True))

    ax1.tick_params(axis='both', which='major', labelsize=30)
    ax2.tick_params(axis='both', which='major', labelsize=30)
    ax3.tick_params(axis='both', which='major', labelsize=30)
    ax3.yaxis.offsetText.set_fontsize(30)

    # add legend
    # legend

    t = np.linspace(0, 10, 200)
    x = np.cos(np.pi * t)
    y = np.sin(t)
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    lc = LineCollection(segments, cmap=cmap,
                        norm=plt.Normalize(0, 10), linewidth=3)
    lc2 = LineCollection(segments, color='k',
                         norm=plt.Normalize(0, 10), linewidth=3)
    lc3 = LineCollection(segments, color='r', linestyle=':',
                         norm=plt.Normalize(0, 10), linewidth=3)

    l1 = ax1.legend(handles=[lc3, lc, lc2], handler_map={
        lc: HandlerColorLineCollection(numpoints=100)},
                    labels=[r'$z_{' + str(ys) + '}^{exp}$',
                            r'$z_{' + str(ys) + '}$',
                            r'$b$'], loc=1)

    l2 = ax2.legend(handles=[lc3, lc, lc2],
                    handler_map={
                        lc: HandlerColorLineCollection(numpoints=100)},
                    labels=[r'$z_{2000}^{obs}$', r'$z_{2000}$',
                            r'$b$'], loc=1)

    l3 = ax3.legend(handles=[lc3, lc],
                    handler_map={
                        lc: HandlerColorLineCollection(numpoints=100)},
                    labels=[r'$s_{' + str(ys) + '-2000}^{exp}$',
                            r'$s_{' + str(ys) +
                            '-2000}$'], loc=1)

    l1.set_zorder(1)
    l2.set_zorder(1)
    l3.set_zorder(1)


    ax3.set_xlim(xmin=1847, xmax=2003)
    fig_name = 'surface_' + str(ys) + '_' + gdir.rgi_id
    #plt.savefig(os.path.join(plot_dir, fig_name + '.pdf'), dpi=300)
    plt.savefig(os.path.join(plot_dir, fig_name + '.png'), dpi=300)
    plt.show()
    plt.close()


def plot_median(gdir, df, eps, ex_mod, ys, ye, plot_dir):
    plot_dir = os.path.join(plot_dir, '04_median')
    utils.mkdir(plot_dir)
    x = (np.arange(ex_mod.fls[-1].nx) * ex_mod.fls[-1].dx * ex_mod.fls[
        -1].map_dx)/1000
    fig = plt.figure(figsize=(25, 18))
    grid = plt.GridSpec(2, 2, hspace=0.2, wspace=0.2)
    ax1 = plt.subplot(grid[0, 0])
    ax2 = plt.subplot(grid[0, 1], sharey=ax1)
    ax3 = plt.subplot(grid[1, :])

    if gdir.name != '':
        plt.suptitle(gdir.rgi_id + ': ' + gdir.name, fontsize=30)
    elif gdir.rgi_id.endswith('779'):
        plt.suptitle(gdir.rgi_id + ': Guslarferner', fontsize=30)
    else:
        plt.suptitle(gdir.rgi_id, fontsize=30)

    df = df.sort_values('fitness', ascending=False)

    # acceptable glacier states
    df = df[df.fitness <=1]
    s_t0 = pd.DataFrame()
    s_te = pd.DataFrame()
    v = pd.DataFrame()
    for i in df.index:
        model = deepcopy(df.loc[i, 'model'])
        s_t0 = s_t0.append(pd.Series(model.fls[-1].surface_h),
                                 ignore_index=True)
        model.run_until(2000)
        s_te = s_te.append(pd.Series(model.fls[-1].surface_h),
                                 ignore_index=True)
        v = v.append(model.volume_km3_ts(), ignore_index=True)

    ax1.fill_between(x, s_t0.max().values, s_t0.min().values, alpha=0.3, color='grey',
                     label=r'$\mathcal{S}_{' + str(ys) + '}^{' + str(
                         eps) + '}$')
    ax2.fill_between(x, s_te.max().values, s_te.min().values, alpha=0.3, color='grey',
                     label=r'$\mathcal{S}_{' + str(ye) + '}^{' + str(
                         eps) + '}$')
    ax3.fill_between(model.volume_m3_ts().index,
                     v.max().values,
                     v.min().values, alpha=0.3, color='grey',
                     label=r'$\mathcal{S}_{' +str(ys)+'-'+str(ye) +'}^{' + str(eps) + '}$')

    # 5% quantile and median of 5% quantile

    df = df[df.fitness <= df.fitness.quantile(0.05)]
    s_t0 = pd.DataFrame()
    s_te = pd.DataFrame()
    v = pd.DataFrame()
    for i in df.index:
        model = deepcopy(df.loc[i, 'model'])
        s_t0 = s_t0.append(pd.Series(model.fls[-1].surface_h),
                           ignore_index=True)
        model.run_until(2000)
        s_te = s_te.append(pd.Series(model.fls[-1].surface_h),
                           ignore_index=True)
        v = v.append(model.volume_km3_ts(), ignore_index=True)

    ax1.fill_between(x, s_t0.max().values, s_t0.min().values, alpha=0.5,
                     label=r'$Q_{0.05}(\mathcal{S}_{'+str(ys)+'}^{'+str(eps)+'})$')

    ax2.fill_between(x, s_te.max().values, s_te.min().values, alpha=0.5,
                     label=r'$Q_{0.05}(\mathcal{S}_{2000}^{'+str(eps)+'})$')

    ax3.fill_between(v.columns, v.max().values, v.min().values, alpha=0.5, linewidth=3,
                     label=r'$Q_{0.05}(\mathcal{S}_{'+str(ys)+'-2000}^{'+str(eps)+'})$')

    # median of 5% quantile
    df.loc[:, 'length'] = df.model.apply(lambda x: x.length_m)
    df = df.sort_values('length', ascending=False)
    l = len(df)
    if l % 2:
        index = int((l - 1) / 2)
    else:
        index = int(l / 2)

    median_model = deepcopy(df.iloc[index].model)
    median_model.volume_km3_ts().plot(ax=ax3, linewidth=3, label=r'$s_{1850-2000}^{med}$')
    median_model.reset_y0(1850)
    median_model.run_until(ys)

    ax1.plot(x, median_model.fls[-1].surface_h, label=r'$z_{1850}^{med}$',
             linewidth=3)
    median_model.run_until(2000)
    ax2.plot(x, median_model.fls[-1].surface_h, label=r'$z_{2000}^{med}$',
             linewidth=3)

    # min model
    min_mod = deepcopy(df.loc[df.fitness.idxmin(), 'model'])
    min_mod.volume_km3_ts().plot(ax=ax3, color='C1',
                                linewidth=3, label=r'$s_{1850-2000}^{min}$')
    min_mod.reset_y0(ys)
    min_mod.run_until(ys)

    ax1.plot(x, min_mod.fls[-1].surface_h, 'C1', label=r'$z_{1850}^{min}$',
             linewidth=3)

    min_mod.run_until(ye)

    ax2.plot(x, min_mod.fls[-1].surface_h, 'C1', label=r'$z_{2000}^{min}$',
             linewidth=3)

    # experiment
    ex_mod = deepcopy(ex_mod)
    ex_mod.volume_km3_ts().plot(ax=ax3, color='k', linestyle=':',
                               linewidth=3, label=r'$s_{1850-2000}^{exp}$')
    ex_mod.reset_y0(1850)
    ex_mod.run_until(ys)

    ax1.plot(x, ex_mod.fls[-1].surface_h, 'k:', label=r'$z_{1850}^{exp}$',
             linewidth=3)
    ax1.plot(x, ex_mod.fls[-1].bed_h, 'k', label=r'$b$', linewidth=3)

    ex_mod.run_until(2000)

    ax2.plot(x, ex_mod.fls[-1].surface_h, 'k:', label=r'$z_{1850}^{exp}$',
             linewidth=3)
    ax2.plot(x, ex_mod.fls[-1].bed_h, 'k', label=r'$b$', linewidth=3)

    # add figure names and x-/ylabels
    add_at(ax1, r"a", loc=3)
    add_at(ax2, r"b", loc=3)
    add_at(ax3, r"c", loc=3)

    ax1.set_ylabel('Altitude (m)', fontsize=30)
    ax1.set_xlabel('Distance along the main flowline (km)', fontsize=30)
    ax2.set_ylabel('Altitude (m)', fontsize=30)
    ax2.set_xlabel('Distance along the main flowline (km)', fontsize=30)
    ax3.set_ylabel(r'Volume ($km^3$)', fontsize=30)
    ax3.set_xlabel('Time (years)', fontsize=30)

    ax1.tick_params(axis='both', which='major', labelsize=30)
    ax2.tick_params(axis='both', which='major', labelsize=30)
    ax3.tick_params(axis='both', which='major', labelsize=30)
    ax3.yaxis.offsetText.set_fontsize(30)
    ax3.set_xlim(xmin=1847, xmax=2003)

    ax1.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax2.xaxis.set_major_locator(MaxNLocator(integer=True))

    l1 = ax1.legend(loc=1, fontsize=25)
    l1.set_zorder(1)

    l2 = ax2.legend(loc=1, fontsize=25)
    l2.set_zorder(1)

    l3 = ax3.legend(loc=1, fontsize=25)
    l3.set_zorder(1)

    fig_name = 'median_'+str(ys)+'_'+gdir.rgi_id
    plt.savefig(os.path.join(plot_dir, fig_name+'.pdf'), dpi=300)
    #plt.savefig(os.path.join(plot_dir, fig_name+'.png'), dpi=300)
    #plt.show()
    plt.close()

    return median_model


def plot_relative_error_min_vs_med(df1, df2, plot_dir):

    plot_dir = os.path.join(plot_dir, 'errors')
    utils.mkdir(plot_dir)

    df1 = df1.replace([np.inf, -np.inf], np.nan).dropna()
    df2 = df2.replace([np.inf, -np.inf], np.nan).dropna()

    f, (ax1, ax2) = plt.subplots(1, 2, sharey=True, figsize=(15, 10), gridspec_kw = {'width_ratios':[1, 3],'wspace':0.05})
    plt.subplots_adjust(hspace=.001)
    plt.suptitle('Relative error (Alps)' + ', n=' + str(len(df1.index)), y=0.99)
    df1.median().plot(ax=ax2, zorder=4, label=r'median($e_{'+'rel'+'}^{med}$)')
    ax2.fill_between(df1.columns.values, df1.quantile(0.95),
                     df1.quantile(0.05), label=r'$Q_{0.05-0.95}(e_{'+'rel'+'}^{med})$',
                     facecolor='C0', alpha=0.5, interpolate=True)

    df2.median().plot(ax=ax2, zorder=4, label=r'median($e_{'+'rel'+'}^{min})$')
    ax2.fill_between(df2.columns.values,
                     df2.quantile(0.95),
                     df2.quantile(0.05), label=r'$Q_{0.05-0.95}(e_{'+'rel'+'}^{min})$',
                     facecolor='C1', alpha=0.5, interpolate=True)

    df1 = df1.loc[:, 1850]
    df1 = df1[(df1 > df1.quantile(0.05)) & (df1 < df1.quantile(0.95))]

    df2 = df2.loc[:, 1850]
    df2 = df2[(df2 > df2.quantile(0.05)) & (df2 < df2.quantile(0.95))]

    bins = np.arange(np.min([df1.min(), df2.min()]), np.max([df1.max(), df2.max()]), 5)

    df1.plot.hist(ax=ax1, bins=bins, color='C0', alpha=0.5, orientation='horizontal')
    df2.plot.hist(ax=ax1, bins=bins, color='C1', alpha=0.5,
                  orientation='horizontal')

    ax1.invert_xaxis()
    ax1.set_xticks([1000, 500, 0])

    ax1.set_title('1850', fontsize=20)
    ax2.set_title('1850-2000', fontsize=20)

    add_at(ax1, r"a", loc=1)
    add_at(ax2, r"b", loc=1)

    ax1.grid()
    ax2.grid()

    ax1.set_ylabel(r'Volume error ($\%$)')
    ax2.set_xlabel('Time (years)')
    plt.legend(loc='best', bbox_to_anchor=(0.94, 1))
    plt.savefig(os.path.join(plot_dir, 'relative_median.png'), dpi=150)
    plt.savefig(os.path.join(plot_dir, 'relative_median.pdf'), dpi=150)
    plt.close()


def plot_absolute_error(df1, df2, plot_dir, all=False):

    plot_dir = os.path.join(plot_dir, 'errors')
    utils.mkdir(plot_dir)

    df1 = df1.replace([np.inf, -np.inf], np.nan).dropna()
    df2 = df2.replace([np.inf, -np.inf], np.nan).dropna()

    title = 'Absolute error (Alps)'
    name = 'absolute'
    unit = '$(km^3)$'

    #plt.figure(figsize=(18, 10))
    f, (ax1, ax2) = plt.subplots(2, 1, figsize=(18, 12), sharex=True)

    ax1.set_title(title+', n=' + str(len(df1.index)))


    if all:
        #plt.figure(figsize=(18, 10))
        #plt.title(title + ', n=' + str(len(df1.index)))

        ax1.fill_between(df1.columns.values,
                         df1.max(),
                         df1.min(),
                         facecolor='C0', alpha=0.2, interpolate=True, zorder=1)

        ax1.fill_between(df2.columns.values,
                         df2.max(),
                         df2.min(),
                         facecolor='C1', alpha=0.2, interpolate=True, zorder=1)

        for yr in df1.columns.values:
            if yr % 10 == 0:
                y1 = df1.loc[:, yr].values
                y2 = df2.loc[:, yr].values

                x1 = np.zeros(len(y1)) + yr - [random.uniform(0.5, 1.5) for x
                                               in range(len(y1))]
                x2 = np.zeros(len(y2)) + yr + [random.uniform(0.5, 1.5) for x
                                               in range(len(y2))]
                ax1.plot(x1, y1, 'o', color='C0', markersize=3)
                ax1.plot(x2, y2, 'o', color='C1', markersize=3)
        ax1.plot(x1, y1, 'o', color='C0', markersize=3,
                 label=r'$e_{abs}^{med}$')
        ax1.plot(x2, y2, 'o', color='C1', markersize=3,
                 label=r'$e_{abs}^{min}$')


    df1.median().plot(zorder=4, label=r'median($e_{abs}^{med}$)', ax=ax2)
    ax2.fill_between(df1.columns.values, df1.quantile(0.95),
                     df1.quantile(0.05), label=r'$Q_{0.05-0.95}(e_{abs}^{med})$',
                     facecolor='C0', alpha=0.5, interpolate=True)

    df2.median().plot(ax=ax2,zorder=4, label=r'median($e_{abs}^{min})$')
    ax2.fill_between(df2.columns.values,
                     df2.quantile(0.95),
                     df2.quantile(0.05), label=r'$Q_{0.05-0.95}(e_{abs}^{min})$',
                     facecolor='C1', alpha=0.5, interpolate=True)
    add_at(ax1, r"a", loc=3)
    add_at(ax2, r"b", loc=3)

    ax1.grid()
    ax2.grid()
    plt.xlabel(r' Time (years)')
    ax1.set_ylabel(r'Volume error ' + unit,labelpad=55)
    ax2.set_ylabel(r'Volume error ' + unit)
    plt.legend(loc='best')
    #plt.savefig(os.path.join(plot_dir, name+'_median.png'), dpi=300)
    #plt.savefig(os.path.join(plot_dir, name + '_median.pdf'), dpi=300)

    #plt.yscale('symlog', linthreshy=1e-2)

    #yaxis = plt.gca().yaxis
    #yaxis.set_minor_locator(MinorSymLogLocator(1e-2))
    ax1.tick_params(which='minor', length=4)
    ax2.tick_params(which='minor', length=4)
    ax1.tick_params(which='major', length=6)
    ax2.tick_params(which='major', length=6)
    from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                                   AutoMinorLocator)
    ax2.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax1.yaxis.set_minor_locator(AutoMinorLocator(2))
    #ax2.set_yticks([-0.015,-0.010, -0.005, 0, 0.005, 0.01])

    ax2.set_xlim(1845, 2005)
    #plt.ylim(-10,10)

    plt.xlabel(r'Time (years)')
    plt.ylabel(r'Volume error ' + unit)
    ax2.legend(loc=4, fontsize=20)
    ax1.legend( fontsize=20)

    plt.savefig(os.path.join(plot_dir, name + '_median_all.png'), dpi=150)
    plt.savefig(os.path.join(plot_dir, name + '_median_all.pdf'), dpi=150)
    plt.show()
    plt.close()


def plot_compare_fitness_functions(df1,df2,df3, plot_dir):

    plot_dir = os.path.join(plot_dir, 'errors')
    utils.mkdir(plot_dir)

    f, (ax1, ax2) = plt.subplots(1, 2, sharey=True, figsize=(15, 10),
                                 gridspec_kw={'width_ratios': [1, 3],
                                              'wspace': 0.05})
    plt.subplots_adjust(hspace=.001)
    plt.suptitle(r'Relative error (Alps), n=' + str(len(df1.index)), y=0.99)

    p1, = ax2.plot(df1.columns.values, df1.median().values, zorder=5)
    ax2.fill_between(df1.columns.values, df1.quantile(0.95),
                     df1.quantile(0.05), alpha=0.5, zorder=3)
    p3, = ax2.plot(df2.columns.values, df2.median().values, zorder=4)
    ax2.fill_between(df2.columns.values, df2.quantile(0.95),
                     df2.quantile(0.05), alpha=0.4, zorder=2)
    p5, = ax2.plot(df3.columns.values, df3.median().values, zorder=4)
    ax2.fill_between(df3.columns.values, df3.quantile(0.95),
                     df3.quantile(0.05), alpha=0.4, zorder=1)

    # plt.grid()
    ax2.set_xlabel('Time (years)')
    ax1.set_ylabel(r'Volume error ($\%$) ')
    p2 = mpatches.Patch(color='C0', alpha=0.5, linewidth=0)
    p4 = mpatches.Patch(color='C1', alpha=0.4, linewidth=0)
    p6 = mpatches.Patch(color='C2', alpha=0.4, linewidth=0)

    df1 = df1.loc[:, 1850]
    df1 = df1[(df1 > df1.quantile(0.05)) & (df1 < df1.quantile(0.95))]

    df2 = df2.loc[:, 1850]
    df2 = df2[(df2 > df2.quantile(0.05)) & (df2 < df2.quantile(0.95))]

    df3 = df3.loc[:, 1850]
    df3 = df3[(df3 > df3.quantile(0.05)) & (df3 < df3.quantile(0.95))]

    bins = np.arange(np.min([df1.min(), df2.min(), df3.min()]),
                     np.max([df1.max(), df2.max(), df3.max()]), 30)

    df1.plot.hist(ax=ax1, bins=bins, color='C0', alpha=0.5,
                  orientation='horizontal')
    df2.plot.hist(ax=ax1, bins=bins, color='C1', alpha=0.5,
                  orientation='horizontal')
    df3.plot.hist(ax=ax1, bins=bins, color='C2', alpha=0.5,
                  orientation='horizontal')
    ax1.invert_xaxis()
    ax1.set_xticks([1000, 500, 0])

    ax1.set_title('1850', fontsize=20)
    ax2.set_title('1850-2000', fontsize=20)

    add_at(ax1, r"a", loc=1)
    add_at(ax2, r"b", loc=1)

    ax2.legend(((p1, p2), (p3, p4), (p5, p6)), ('Geometry', 'Area', 'Length'),
               loc='best', bbox_to_anchor=(0.95, 1))

    ax1.grid()
    ax2.grid()

    plt.savefig(os.path.join(plot_dir, 'compare_median.png'), dpi=150)
    plt.savefig(os.path.join(plot_dir, 'compare_median.pdf'), dpi=150)
    #plt.show()
    plt.close()

def plot_supplement(gdir, df, ex_mod, ys,ye,eps, plot_dir):

    plot_dir = os.path.join(plot_dir, 'supplement')
    utils.mkdir(plot_dir)
    x = (np.arange(ex_mod.fls[-1].nx) * ex_mod.fls[-1].dx * \
         ex_mod.fls[-1].map_dx) / 1000
    f = plt.figure(figsize=(35, 20))

    gs0 = gridspec.GridSpec(1, 2, figure=f, hspace=0.2, wspace=0.2, right=1.05,left=0.1)

    gs00 = gridspec.GridSpecFromSubplotSpec(2, 2, subplot_spec=gs0[0])

    ax1 = f.add_subplot(gs00[0, 0])
    ax2 = f.add_subplot(gs00[0, 1],sharey=ax1)
    ax3 = f.add_subplot(gs00[1, :])

    # the following syntax does the same as the GridSpecFromSubplotSpec call above:
    gs01 = gs0[1].subgridspec(2, 2)

    ax4 = f.add_subplot(gs01[0, 0])
    ax5 = f.add_subplot(gs01[0, 1],sharey=ax4)
    ax6 = f.add_subplot(gs01[1, :])

    if gdir.name != '':
        plt.suptitle(gdir.rgi_id + ': ' + gdir.name, fontsize=30)
    elif gdir.rgi_id.endswith('779'):
        plt.suptitle(gdir.rgi_id + ': Guslarferner', fontsize=30)
    else:
        plt.suptitle(gdir.rgi_id, fontsize=30)

    norm = mpl.colors.LogNorm(vmin=0.01 / 125, vmax=10)
    cmap = matplotlib.cm.get_cmap('viridis')

    # df = df.sort_values('objective', ascending=False)
    df = df.sort_values('fitness', ascending=False)


    for i, model in df['model'].iteritems():
        model = deepcopy(model)
        model.reset_y0(ys)

        # color = cmap(norm(df.loc[i, 'objective']))
        color = cmap(norm(df.loc[i, 'fitness']))
        ax4.plot(x, deepcopy(model.fls[-1].surface_h), color=color,
                 label='')
        model.volume_km3_ts().plot(ax=ax6, color=[color], label='')
        model.run_until(2000)
        ax5.plot(x, model.fls[-1].surface_h, color=color, label='')

    # acceptable glacier states
    df = df[df.fitness <= 1]
    s_t0 = pd.DataFrame()
    s_te = pd.DataFrame()
    v = pd.DataFrame()

    for i in df.index:
        model = deepcopy(df.loc[i, 'model'])
        s_t0 = s_t0.append(pd.Series(model.fls[-1].surface_h),
                                 ignore_index=True)
        model.run_until(2000)
        s_te = s_te.append(pd.Series(model.fls[-1].surface_h),
                                 ignore_index=True)
        v = v.append(model.volume_km3_ts(), ignore_index=True)

    ax1.fill_between(x, s_t0.max().values, s_t0.min().values, alpha=0.3,
                     color='grey',
                     label=r'$\mathcal{S}_{' + str(ys) + '}^{' + str(
                         eps) + '}$')
    ax2.fill_between(x, s_te.max().values, s_te.min().values, alpha=0.3,
                     color='grey',
                     label=r'$\mathcal{S}_{' + str(ye) + '}^{' + str(
                         eps) + '}$')
    ax3.fill_between(model.volume_m3_ts().index,
                     v.max().values,
                     v.min().values, alpha=0.3, color='grey',
                     label=r'$\mathcal{S}_{' + str(ys) + '-' + str(
                         ye) + '}^{' + str(eps) + '}$')

    # median of 5% quantile

    median_model, min_mod, qu = find_median(df)
    print(median_model)
    median_model.volume_km3_ts().plot(ax=ax3, linewidth=4,
                                      label=r'$s_{1850-2000}^{med}$')
    median_model.reset_y0(1850)
    median_model.run_until(ys)

    ax1.plot(x, median_model.fls[-1].surface_h, label=r'$z_{1850}^{med}$',
             linewidth=4)
    median_model.run_until(2000)
    ax2.plot(x, median_model.fls[-1].surface_h, label=r'$z_{2000}^{med}$',
             linewidth=4)

    # min model
    min_mod = deepcopy(df.loc[df.fitness.idxmin(), 'model'])
    min_mod.volume_km3_ts().plot(ax=ax3, color='C1',
                                 linewidth=3, label=r'$s_{1850-2000}^{min}$')
    min_mod.reset_y0(ys)
    min_mod.run_until(ys)

    ax1.plot(x, min_mod.fls[-1].surface_h, 'C1', label=r'$z_{1850}^{min}$',
             linewidth=3)

    min_mod.run_until(ye)

    ax2.plot(x, min_mod.fls[-1].surface_h, 'C1', label=r'$z_{2000}^{min}$',
             linewidth=3)



    # 5% quantile and median of 5% quantile
    df.fitness[df.fitness < 1e-4] = 1e-4
    df = df[df.fitness <= round(df.fitness.quantile(0.05),4)]

    s_t0 = pd.DataFrame()
    s_te = pd.DataFrame()
    v = pd.DataFrame()
    for i in df.index:
        model = deepcopy(df.loc[i, 'model'])
        s_t0 = s_t0.append(pd.Series(model.fls[-1].surface_h),
                           ignore_index=True)
        model.run_until(2000)
        s_te = s_te.append(pd.Series(model.fls[-1].surface_h),
                           ignore_index=True)
        v = v.append(model.volume_km3_ts(), ignore_index=True)


    ax1.fill_between(x, s_t0.max().values, s_t0.min().values, alpha=0.5,
                     label=r'$Q_{0.05}$')

    ax2.fill_between(x, s_te.max().values, s_te.min().values, alpha=0.5,
                     label=r'$Q_{0.05}$')

    ax3.fill_between(v.columns, v.max().values, v.min().values, alpha=0.5,
                     linewidth=3,
                     label=r'$Q_{0.05}(\mathcal{S}_{' + str(
                         ys) + '-2000}^{' + str(eps) + '})$')



    # plot experiment
    ex_mod = deepcopy(ex_mod)
    ex_mod.volume_km3_ts().plot(ax=ax6, color='red', linestyle=':',
                                linewidth=3,
                                label='')
    ex_mod.volume_km3_ts().plot(ax=ax3, color='red', linestyle=':',
                                linewidth=3,
                                label='')
    ex_mod.reset_y0(1850)
    ex_mod.run_until(ys)

    ax4.plot(x, ex_mod.fls[-1].surface_h, ':', color='red', label='',
             linewidth=3)
    ax4.plot(x, ex_mod.fls[-1].bed_h, 'k', label='', linewidth=3)
    ax1.plot(x, ex_mod.fls[-1].surface_h, ':', color='red', label='',
             linewidth=3)
    ax1.plot(x, ex_mod.fls[-1].bed_h, 'k', label='', linewidth=3)

    ex_mod.run_until(2000)

    ax5.plot(x, ex_mod.fls[-1].surface_h, ':', color='red', label='',
             linewidth=3)
    ax5.plot(x, ex_mod.fls[-1].bed_h, 'k', label='',
             linewidth=3)
    ax2.plot(x, ex_mod.fls[-1].surface_h, ':', color='red', label='',
             linewidth=3)
    ax2.plot(x, ex_mod.fls[-1].bed_h, 'k', label='',
             linewidth=3)

    # add colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cax, kw = mpl.colorbar.make_axes([ax1,ax2,ax3,ax4,ax5,ax6])
    cbar = f.colorbar(sm, cax=cax, extend='both', **kw)
    cbar.ax.tick_params(labelsize=30)
    cbar.set_label('Fitness value', fontsize=30)

    # add figure names and x-/ylabels
    add_at(ax1, r"a", loc=3)
    add_at(ax2, r"b", loc=3)
    add_at(ax3, r"c", loc=3)

    ax1.set_ylabel('Altitude (m)', fontsize=30)
    ax3.set_title('Distance along the main flowline (km)', fontsize=30, pad=40)
    ax2.set_ylabel('Altitude (m)', fontsize=30)
    ax3.set_ylabel(r'Volume ($km^3$)', fontsize=30)
    ax3.set_xlabel('Time (years)', fontsize=30)

    ax4.set_ylabel('Altitude (m)', fontsize=30)
    ax6.set_title('Distance along the main flowline (km)', fontsize=30, pad=40)
    ax5.set_ylabel('Altitude (m)', fontsize=30)
    ax6.set_ylabel(r'Volume ($km^3$)', fontsize=30)
    ax6.set_xlabel('Time (years)', fontsize=30)


    ax1.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax2.xaxis.set_major_locator(MaxNLocator(integer=True))

    ax1.tick_params(axis='both', which='major', labelsize=30)
    ax2.tick_params(axis='both', which='major', labelsize=30)
    ax3.tick_params(axis='both', which='major', labelsize=30)
    ax3.yaxis.offsetText.set_fontsize(30)

    ax4.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax5.xaxis.set_major_locator(MaxNLocator(integer=True))

    ax4.tick_params(axis='both', which='major', labelsize=30)
    ax5.tick_params(axis='both', which='major', labelsize=30)
    ax6.tick_params(axis='both', which='major', labelsize=30)
    ax6.yaxis.offsetText.set_fontsize(30)

    # add legend
    # legend

    t = np.linspace(0, 10, 200)
    x = np.cos(np.pi * t)
    y = np.sin(t)
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    lc = LineCollection(segments, cmap=cmap,
                        norm=plt.Normalize(0, 10), linewidth=3)
    lc2 = LineCollection(segments, color='k',
                         norm=plt.Normalize(0, 10), linewidth=3)
    lc3 = LineCollection(segments, color='r', linestyle=':',
                         norm=plt.Normalize(0, 10), linewidth=3)

    l1 = ax1.legend(loc='best', fontsize=30)
    l1.set_zorder(1)

    l2 = ax2.legend(loc='best', fontsize=30)
    l2.set_zorder(1)

    l3 = ax3.legend(loc='best', fontsize=30)
    l3.set_zorder(1)

    l4 = ax4.legend(handles=[lc3, lc, lc2], handler_map={
        lc: HandlerColorLineCollection(numpoints=100)},
                    labels=[r'$z_{' + str(ys) + '}^{exp}$',
                            r'$z_{' + str(ys) + '}$',
                            r'$b$'], loc='best', fontsize=30)

    l5 = ax5.legend(handles=[lc3, lc, lc2],
                    handler_map={
                        lc: HandlerColorLineCollection(numpoints=100)},
                    labels=[r'$z_{2000}^{obs}$', r'$z_{2000}$',
                            r'$b$'], loc='best', fontsize=30)

    l6 = ax6.legend(handles=[lc3, lc],
                    handler_map={
                        lc: HandlerColorLineCollection(numpoints=100)},
                    labels=[r'$s_{' + str(ys) + '-2000}^{exp}$',
                            r'$s_{' + str(ys) +
                            '-2000}$'], loc='best', fontsize=30)

    l4.set_zorder(1)
    l5.set_zorder(1)
    l6.set_zorder(1)

    ax6.set_xlim(xmin=1847, xmax=2003)
    ax3.set_xlim(xmin=1847, xmax=2003)
    ax3.set_ylim(ax6.get_ylim())

    fig_name = 'supplement_' + str(ys) + '_' + gdir.rgi_id.split('.')[0]+'_'+gdir.rgi_id.split('.')[1]
    plt.savefig(os.path.join(plot_dir, fig_name + '.pdf'), dpi=150)
    # plt.savefig(os.path.join(plot_dir, fig_name + '.png'), dpi=300)
    #plt.show()
    plt.close()