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
#FlowlineModel = partial(FluxBasedModel, inplace=False)
from mpl_toolkits.axes_grid.anchored_artists import AnchoredText
import matplotlib as mpl
from matplotlib.collections import LineCollection
from matplotlib.legend_handler import HandlerLineCollection
import random
pd.options.mode.chained_assignment = None
warnings.filterwarnings("ignore")

mpl.rcParams['axes.linewidth'] = 3
mpl.rcParams['xtick.major.width'] = 3
mpl.rcParams['ytick.major.width'] = 3
mpl.rcParams['font.size'] =25
mpl.rcParams['font.weight'] = 'medium'
mpl.rcParams['axes.labelweight'] = 'medium'
mpl.rcParams['legend.fontsize'] = 20 #30
mpl.rcParams['lines.linewidth'] = 3


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

def plot_advanced_experiment(gdir):

    fig = plt.figure(figsize=(15, 14))
    grid = plt.GridSpec(1, 2, hspace=0.2, wspace=0.2)
    ax1 = plt.subplot(grid[0])
    ax2 = plt.subplot(grid[1], sharey=ax1)
    ax2.plot(gdir.rgi_date,gdir.rgi_area_km2,'o', label='RGI area '+ str(gdir.rgi_date))
    mod = FluxBasedModel(flowlines=gdir.read_pickle('model_flowlines'))
    ax2.plot(gdir.rgi_date, mod.area_km2, 'o',
             label='RGI area ' + str(gdir.rgi_date))

    for f in os.listdir(gdir.dir):
        if f.startswith('model_run_'):
            rp = os.path.join(gdir.dir,f)
            model = FileModel(rp)
            model.area_km2_ts().plot(ax=ax2, label=f)

    #ax2.set_xlim((1915,2005))
    #ax2.legend()
    plt.show()


def plot_fitness_values(gdir, df, ex_mod, ys, ye, lec, plot_dir):

    plot_dir = os.path.join(plot_dir, '03_surface_by_fitness')
    utils.mkdir(plot_dir)
    x = np.arange(ex_mod.fls[-1].nx) * ex_mod.fls[-1].dx * \
        ex_mod.fls[-1].map_dx
    fig, ax = plt.subplots(1,1,figsize=(25, 18))

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

        model.length_m_ts(rollmin=5).plot(ax=ax, color=[color], label='')

    # plot experiments.py
    lec.plot(ax=ax, color='red',linewidth=3, label='')

    min_mod = df.loc[df.fitness.idxmin(),'model']

    # add colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cax, kw = mpl.colorbar.make_axes([ax])
    cbar = fig.colorbar(sm, cax=cax, extend='both', **kw)
    cbar.ax.tick_params(labelsize=30)
    cbar.set_label('Fitness value', fontsize=30)


    ax.set_ylabel(r'Volume ($km^3$)', fontsize=30)
    ax.set_xlabel('Time (years)', fontsize=30)

    ax.tick_params(axis='both', which='major', labelsize=30)
    ax.yaxis.offsetText.set_fontsize(30)

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

    l3 = ax.legend(handles=[lc3, lc],
                    handler_map={
                        lc: HandlerColorLineCollection(numpoints=100)},
                    labels=['Leclercq',
                            r'$s_{' + str(ys) +
                            '-'+str(ye)+'}$'], loc=1)

    l3.set_zorder(0)

    ax.set_xlim(xmin=ys-5, xmax=ye-5)
    fig_name = 'surface_' + str(ys) + '_' + gdir.rgi_id
    #plt.savefig(os.path.join(plot_dir, fig_name + '.pdf'), dpi=300)
    #plt.savefig(os.path.join(plot_dir, fig_name + '.png'), dpi=300)
    plt.show()
    plt.close()

def plot_fitness_old(gdir, df, ex_mod, ys, ye, lec,plot_dir):


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
    elif gdir.rgi_id.endswith('11.00779'):
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
        (model.length_m_ts()/1000).plot(ax=ax3, color=[color], label='')
        model.run_until(ye)

        ax2.plot(x, model.fls[-1].surface_h, color=color, label='')

    # plot experiments.py
    ex_mod = deepcopy(ex_mod)
    (lec/1000).plot(ax=ax3, color='red', linewidth=3, label='')
    #(ex_mod.length_m_ts()/1000).plot(ax=ax3, color='red', linestyle=':',
    #                           linewidth=3,
    #                           label='')
    ex_mod.reset_y0(ys)
    ex_mod.run_until(ys)

    ax1.plot(x, ex_mod.fls[-1].surface_h, ':', color='red', label='',
             linewidth=3)
    ax1.plot(x, ex_mod.fls[-1].bed_h, 'k', label='', linewidth=3)

    ex_mod.run_until(ye)

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
    ax1.set_xlabel('Distance along the main flowline (m)', fontsize=30)
    ax2.set_ylabel('Altitude (m)', fontsize=30)
    ax2.set_xlabel('Distance along the main flowline (m)', fontsize=30)
    ax3.set_ylabel(r'Length ($km$)', fontsize=30)
    ax3.set_xlabel('Time (years)', fontsize=30)

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

    l1.set_zorder(0)
    l2.set_zorder(0)
    l3.set_zorder(5)

    ax3.set_xlim(xmin=ys-5, xmax=ye+5)
    fig_name = 'surface_' + str(ys) + '_' + gdir.rgi_id
    plt.savefig(os.path.join(plot_dir, fig_name + '.pdf'), dpi=300)
    plt.savefig(os.path.join(plot_dir, fig_name + '.png'), dpi=300)
    #plt.show()
    plt.close()



def plot_median(gdir, df, eps, ex_mod, ys, ye, lec, plot_dir):
    plot_dir = os.path.join(plot_dir, '04_median')
    utils.mkdir(plot_dir)
    x = np.arange(ex_mod.fls[-1].nx) * ex_mod.fls[-1].dx * ex_mod.fls[
        -1].map_dx
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

    # min model

    min_mod = deepcopy(df.loc[df.fitness.idxmin(), 'model'])
    min_mod.length_m_ts(rollmin=5).plot(ax=ax3, color='C1',
                                        linewidth=3,
                                        label=r'$s_{' + str(ys) + '-' + str(
                                            ye) + '}^{min}$')
    min_mod.reset_y0(ys)
    min_mod.run_until(ys)
    # real flowlines
    fls = gdir.read_pickle('model_flowlines')
    mod = FluxBasedModel(flowlines=fls)
    ax2.plot(x, mod.fls[-1].surface_h, 'C2', label=r'OGGM$_{init}$')

    ax1.plot(x, min_mod.fls[-1].surface_h, 'C1',
             label=r'$z_{' + str(ys) + '}^{min}$',
             linewidth=3)

    min_mod.run_until(ye)

    ax2.plot(x, min_mod.fls[-1].surface_h, 'C1', label=r'$z_{2000}^{min}$',
             linewidth=3)

    # acceptable glacier states
    df = df[df.fitness <=1]
    s_t0 = pd.DataFrame()
    s_te = pd.DataFrame()
    v = pd.DataFrame()
    for i in df.index:
        model = deepcopy(df.loc[i, 'model'])
        s_t0 = s_t0.append(pd.Series(model.fls[-1].surface_h),
                                 ignore_index=True)
        model.run_until(ye)
        s_te = s_te.append(pd.Series(model.fls[-1].surface_h),
                                 ignore_index=True)
        v = v.append(model.length_m_ts(rollmin=5), ignore_index=True)

    ax1.fill_between(x, s_t0.max().values, s_t0.min().values, alpha=0.3, color='grey',
                     label=r'$\mathcal{S}_{' + str(ys) + '}^{' + str(
                         eps) + '}$')
    ax2.fill_between(x, s_te.max().values, s_te.min().values, alpha=0.3, color='grey',
                     label=r'$\mathcal{S}_{' + str(ye) + '}^{' + str(
                         eps) + '}$')
    ax3.fill_between(model.length_m_ts(rollmin=5).index,
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
        model.run_until(ye)
        s_te = s_te.append(pd.Series(model.fls[-1].surface_h),
                           ignore_index=True)
        v = v.append(model.length_m_ts(rollmin=5), ignore_index=True)

    ax1.fill_between(x, s_t0.max().values, s_t0.min().values, alpha=0.5,
                     label=r'$Q_{0.05}(\mathcal{S}_{'+str(ys)+'}^{'+str(eps)+'})$')

    ax2.fill_between(x, s_te.max().values, s_te.min().values, alpha=0.5,
                     label=r'$Q_{0.05}(\mathcal{S}_{'+str(ye)+'}^{'+str(eps)+'})$')

    ax3.fill_between(v.columns, v.max().values, v.min().values, alpha=0.5, linewidth=3,
                     label=r'$Q_{0.05}(\mathcal{S}_{'+str(ys)+'-'+str(ye)+'}^{'+str(eps)+'})$')

    # median of 5% quantile
    df.loc[:, 'length'] = df.model.apply(lambda x: x.length_m)
    df = df.sort_values('length', ascending=False)
    l = len(df)
    if l % 2:
        index = int((l - 1) / 2)
    else:
        index = int(l / 2)

    median_model = deepcopy(df.iloc[index].model)
    median_model.length_m_ts(rollmin=5).plot(ax=ax3, linewidth=3, label=r'$s_{'+str(ys)+'-'+str(ye)+'}^{med}$')
    median_model.reset_y0(ys)
    median_model.run_until(ys)

    ax1.plot(x, median_model.fls[-1].surface_h, label=r'$z_{'+str(ys)+'}^{med}$',
             linewidth=3)
    median_model.run_until(ye)
    ax2.plot(x, median_model.fls[-1].surface_h, label=r'$z_{'+str(ye)+'}^{med}$',
             linewidth=3)


    # experiment
    lec.plot(ax=ax3, color='r', linewidth=3, label='Leclercq')
    ax1.plot(x, ex_mod.fls[-1].bed_h, 'k', label=r'$b$', linewidth=3)
    ax2.plot(x, ex_mod.fls[-1].bed_h, 'k', label=r'$b$', linewidth=3)

    #ex_mod.length_m_ts(rollmin=5).plot(ax=ax3, color='k')

    # add figure names and x-/ylabels
    add_at(ax1, r"a", loc=3)
    add_at(ax2, r"b", loc=3)
    add_at(ax3, r"c", loc=3)

    ax1.set_ylabel('Altitude (m)', fontsize=30)
    ax1.set_xlabel('Distance along the main flowline (m)', fontsize=30)
    ax2.set_ylabel('Altitude (m)', fontsize=30)
    ax2.set_xlabel('Distance along the main flowline (m)', fontsize=30)
    ax3.set_ylabel(r'Length ($m$)', fontsize=30)
    ax3.set_xlabel('Time (years)', fontsize=30)

    ax1.tick_params(axis='both', which='major', labelsize=30)
    ax2.tick_params(axis='both', which='major', labelsize=30)
    ax3.tick_params(axis='both', which='major', labelsize=30)
    ax3.yaxis.offsetText.set_fontsize(30)
    ax3.set_xlim(xmin=1915, xmax=2019)

    l1 = ax1.legend(loc=1, fontsize=25)
    l1.set_zorder(1)

    l2 = ax2.legend(loc=1, fontsize=25)
    l2.set_zorder(1)

    l3 = ax3.legend(loc=1, fontsize=25)
    l3.set_zorder(1)

    fig_name = 'median_'+str(ys)+'_'+gdir.rgi_id
    plt.savefig(os.path.join(plot_dir, fig_name+'.pdf'), dpi=300)
    plt.savefig(os.path.join(plot_dir, fig_name+'.png'), dpi=300)
    plt.show()
    plt.close()