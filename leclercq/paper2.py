
import os
import sys
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid.anchored_artists import AnchoredText
from oggm.utils._downloads import get_demo_file
import matplotlib as mpl
from oggm import cfg, utils
from oggm.core.flowline import FluxBasedModel, FileModel
from matplotlib.legend_handler import HandlerLineCollection
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from matplotlib.collections import LineCollection
import oggm
import copy
import seaborn as sns
import numpy as np
from initialization.core import fitness_value_fls
pd.set_option('display.max_columns', None)

mpl.rcParams['axes.linewidth'] = 3
mpl.rcParams['xtick.major.width'] = 3
mpl.rcParams['ytick.major.width'] = 3
mpl.rcParams['font.size'] =25
mpl.rcParams['font.weight'] = 'medium'
mpl.rcParams['axes.labelweight'] = 'medium'
mpl.rcParams['legend.fontsize'] = 15 #30
mpl.rcParams['legend.title_fontsize'] = 20
mpl.rcParams['lines.linewidth'] = 3

def add_at(ax, t, loc=2):
    fp = dict(size=20)
    _at = AnchoredText(t, loc=loc, prop=fp, borderpad=0)
    _at.patch.set_linewidth(3)
    ax.add_artist(_at)
    return _at

def find_median(df):
    try:
        df.loc[:,'length'] = df.model.apply(lambda x: x.length_m)
        df = df.sort_values(by='length')
        accept_df = df[df.fitness <= 1]
        quant_df = accept_df[accept_df.fitness <= accept_df.fitness.quantile(0.05)]
        # median state
        #quant_df.loc[:, 'length'] = quant_df.model.apply(lambda x: x.length_m)
        #quant_df = quant_df.sort_values('length', ascending=False)
        l = len(quant_df)
        if l % 2:
            index = int((l - 1) / 2)
        else:
            index = int(l / 2)
        print(quant_df.iloc[index].fitness)

        return quant_df.iloc[index].model, quant_df.at[quant_df.length.idxmin(),'model'], quant_df.at[quant_df.length.idxmax(),'model']
    except:
        return df.iloc[df.fitness.idxmin()].model, None, None


def plot_bias_density(exp_df, plot_dir):
    exp_df = exp_df.replace([np.inf, -np.inf], np.nan)
    exp_df = exp_df.dropna()

    plt.figure(figsize=(15,15))
    grid = plt.GridSpec(2, 1, hspace=0.3, wspace=0.2)
    ax2 = plt.subplot(grid[0])
    ax1 = plt.subplot(grid[1])

    cmap = plt.get_cmap("tab10")

    for i, temp_bias in enumerate(np.sort(exp_df.temp_bias.unique())):
        j = i
        if i >= 3:
            j = j + 1
        # Subset to the airline
        subset = exp_df[exp_df['temp_bias'] == temp_bias]

        # Draw the density plot
        sns.distplot(subset['bias'], hist=False, kde=True,
                     kde_kws={'shade': True, 'linewidth': 3, 'alpha':0.15},
                     bins=30, ax=ax2)
        ax2.plot(ax2.lines[-1].get_xydata()[:,0],ax2.lines[-1].get_xydata()[:,1], color=cmap(j), label=temp_bias)

        ax2.axvline(x=subset.bias.mean(), color=cmap(j), linestyle=':')

    exp_df.area_diff.plot.hist(ax=ax1, bins=50)

    ax1.set_yscale('log')
    # ax1.set_xscale('symlog',linthreshx=1e-2)
    ax1.set_ylabel('Frequency', labelpad=30)
    ax1.set_xlabel(r'Area difference (km$^2$)')
    ax1.set_xlim(-1.05, 1.05)
    ax2.set_ylim(0,None)
    ax2.set_yticks([0,1e-3,2e-3])
    ax2.set_xticks(np.arange(-1500,2000,500))
    from matplotlib.ticker import LogLocator
    ax1.yaxis.set_minor_locator(LogLocator(base=10, subs='auto'))
    ax1.tick_params(which='major', length=6, width=3)
    ax1.tick_params(which='minor', length=4, color='k', width=2)




    add_at(ax2, 'a')
    add_at(ax1, 'b')

    ax1.grid()
    ax2.grid()
    ax2.set_xlabel(r'Optimal mass balance balance bias $\beta^*_{mb}$ (mm w.e.)')
    ax2.set_ylabel('Density')
    ax2.legend(title=r'temp. bias $\widetilde{\beta}$')
    plt.savefig(os.path.join(plot_dir,'mb_bias.png'),dpi=300)
    plt.show()


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

def read_multi_pickle(path):
    df = pd.DataFrame()
    for file in [file for file in os.listdir(path) if file.startswith('result')]:
        p = os.path.join(path,file)
        result = pd.read_pickle(p,compression='gzip')
        temp_bias = float(file.split('_')[-1].split('.pkl')[0])
        result.loc[:,'temp_bias']= temp_bias
        df = df.append(result, ignore_index=True)
    df.fitness = df.fitness/125
    return df


def plot_multi_fitness(df, exp_df,leclercq, dir):


    fig = plt.figure(figsize=(35, 15))

    grid = fig.add_gridspec(ncols=5, nrows=2, wspace=0.1, left=0.05, right=1.05)

    ye = exp_df.iloc[0].model.length_m_ts().index[-1]

    leclercq = leclercq[leclercq.index <= ye]

    norm = mpl.colors.LogNorm(vmin=0.01 / 125, vmax=10)
    cmap = mpl.cm.get_cmap('viridis')
    axes = []
    num1 = ['a','b','c','d','e']
    num2 = ['f','g','h','i','j']
    for col,tb in enumerate([0,-0.25,-0.5,-0.75,-1]):
        if col != 0:
            ax = fig.add_subplot(grid[0, col], sharey=ax)
            ax2 = fig.add_subplot(grid[1, col], sharey=ax2)
        else:
            ax = fig.add_subplot(grid[0, col])
            ax.set_ylabel('Length (km)')
            ax2 = fig.add_subplot(grid[1, col])
            ax2.set_ylabel('Length (km)')
        axes.append(ax)
        axes.append(ax2)

        add_at(ax,num1[col], loc=1)
        add_at(ax2, num2[col], loc=1)

        exp = exp_df[exp_df.temp_bias==tb].iloc[0].model

        # get fls model
        fls =  pd.read_pickle(os.path.join(dir,'model_flowlines.pkl'),compression='gzip')
        fls_mod = FluxBasedModel(flowlines=fls)

        lec = (leclercq - leclercq.iloc[-1]) + fls_mod.length_m
        res = df[df.temp_bias==tb].sort_values(by='fitness', ascending=False)
        for i, model in res['model'].iteritems():
            color = cmap(norm(res.loc[i, 'fitness']))
            (model.length_m_ts(rollmin=5)/1000).plot(ax=ax, color=[color], label='', linewidth=4)
            fitness_fls = fitness_value_fls(model, fls_mod, ye)
            color2 = cmap(norm(fitness_fls))
            (model.length_m_ts(rollmin=5) / 1000).plot(ax=ax2, color=[color2], label='', linewidth=4)

        (exp.length_m_ts(rollmin=5)/1000).plot(ax=ax, color='red',linestyle=':', linewidth=4, label='')

        (lec/1000).plot(ax=ax, color='red', linewidth=4, label='Leclercq (2014)')
        (lec / 1000).plot(ax=ax2, color='red', linewidth=4, label='Leclercq (2014)')

        ax.plot(ye, fls_mod.length_m / 1000,'ro', markersize=15)
        ax2.plot(ye, fls_mod.length_m / 1000, 'ro', markersize=15)

        mb_bias = exp_df[exp_df.temp_bias==tb].bias.values[0]
        ax.set_title(r'$\widetilde{\beta}$= '+str(tb)+r', $\beta^*_{mb}$=' + str(mb_bias), pad=10)
        ax.set_xlim(1917,None)
        ax.set_xlabel('Time (years)')

        ax2.set_xlim(1917, None)
        ax2.set_xlabel('Time (years)')

    l1 = Line2D([0],[0],color='r',lw=4)
    l2 = Line2D([0], [0], marker='o', color='w', label='Scatter',
           markerfacecolor='r', markersize=15)
    l3 = Line2D([0],[0],color='r',linestyle=':',lw=4)
    t = np.linspace(0, 10, 200)
    x = np.cos(np.pi * t)
    y = np.sin(t)
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    lc  = LineCollection(segments, cmap=cmap, norm=plt.Normalize(0, 10), linewidth=4)
    fig.legend(handles=[l2,l1,l3,lc], handler_map={lc: HandlerColorLineCollection(numpoints=100)},
               labels=[r'$s^{OGGM}_{2003}$','Leclercq (2014)',r'$s_{1917-2003}^{exp}$',r'$s_{1917-2003}$'],
               loc = 'upper center', bbox_to_anchor=(0.5, 0.0625), ncol=4, fontsize=25)


    # colorbar

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cax, kw = mpl.colorbar.make_axes(axes)
    cbar = fig.colorbar(sm, cax=cax, extend='both',**kw)
    cbar.ax.tick_params(labelsize=25)
    cbar.set_label('Fitness value', fontsize=25)
    plt.suptitle('RGI60-11.00897, Hintereisferner',y=0.975)
    plt.savefig(os.path.join(dir,'multi_fitness.png'))
    plt.show()


def plot_mulit_median(df, exp_df, leclercq, dir):
    ye = exp_df.iloc[0].model.length_m_ts().index[-1]
    ex_mod = exp_df.iloc[0].model
    x = (np.arange(ex_mod.fls[-1].nx) * ex_mod.fls[-1].dx * ex_mod.fls[
        -1].map_dx) / 1000

    fig = plt.figure(figsize=(25, 18))
    grid = plt.GridSpec(2, 2, hspace=0.2, wspace=0.2,left=0.1, right=0.8)
    ax1 = plt.subplot(grid[0, 0])
    ax2 = plt.subplot(grid[0, 1], sharey=ax1)
    ax3 = plt.subplot(grid[1, :])

    cmap = plt.get_cmap("tab10")
    for i, temp_bias in enumerate(np.sort(exp_df.temp_bias.unique())):
        if i>=3:
            i=i+1

        res = df[df.temp_bias==temp_bias]
        median,q_min, q_max  = find_median(res)

        ax1.plot(x,copy.deepcopy(median.fls[-1].surface_h), color=cmap(i))
        ax1.fill_between(x,copy.deepcopy(q_min.fls[-1].surface_h),copy.deepcopy(q_max.fls[-1].surface_h), color=cmap(i),alpha=0.5)
        print(temp_bias, median.length_m / 1000, median.area_km2,
              median.volume_km3)
        median.run_until(2003)
        #print(temp_bias, median.length_m/1000, median.area_km2, median.volume_km3)
        q_min.run_until(2003)
        q_max.run_until(2003)
        ax2.plot(x, copy.deepcopy(median.fls[-1].surface_h), color=cmap(i))
        ax2.fill_between(x, copy.deepcopy(q_min.fls[-1].surface_h),
                         copy.deepcopy(q_max.fls[-1].surface_h), color=cmap(i),
                         alpha=0.5)
        ax3.plot(median.length_m_ts().index,(median.length_m_ts()/1000).values, color=cmap(i))
        ax3.fill_between(median.length_m_ts().index, (q_max.length_m_ts()/1000).values,(q_min.length_m_ts()/1000).values, color=cmap(i),
                         alpha=0.5)
    # get fls model
    fls = pd.read_pickle(os.path.join(dir, 'model_flowlines.pkl'), compression='gzip')
    fls_mod = FluxBasedModel(flowlines=fls)
    ax2.plot(x, fls_mod.fls[-1].surface_h,color='gold')

    print(fls_mod.length_m/1000, fls_mod.area_km2, fls_mod.volume_km3)
    lec = (leclercq - leclercq.iloc[-1]) + fls_mod.length_m
    (lec / 1000).plot(ax=ax3, color='red', linewidth=3, label='Leclercq (2014)')
    #ax3.plot(ye,(lec/1000)[ye],'o', color='gold', markersize=10)

    #ax2.plot(fls_mod.length_m/1000,ex_mod.fls[-1].bed_h[np.where(x == fls_mod.length_m/1000)],'ro', markersize=10, zorder=5)


    ax1.plot(x, ex_mod.fls[-1].bed_h, 'k')
    ax2.plot(x, ex_mod.fls[-1].bed_h, 'k')
    ax3.set_xlim((1917,2003))

    l1 = Line2D([0], [0], color='r', lw=4)
    l2 = Line2D([0], [0], color='gold', lw=4)
    l3 = (Line2D([0], [0], color=cmap(0), lw=4), mpatches.Patch(color=cmap(0), alpha=0.5, linewidth=0))
    l4 = (Line2D([0], [0], color=cmap(1), lw=4), mpatches.Patch(color=cmap(1), alpha=0.5, linewidth=0))
    l5 = (Line2D([0], [0], color=cmap(2), lw=4), mpatches.Patch(color=cmap(2), alpha=0.5, linewidth=0))
    l6 = (Line2D([0], [0], color=cmap(4), lw=4), mpatches.Patch(color=cmap(4), alpha=0.5, linewidth=0))
    l7 = (Line2D([0], [0], color=cmap(5), lw=4), mpatches.Patch(color=cmap(5), alpha=0.5, linewidth=0))


    ax3.legend(handles=[l2, l1],
               labels=[r'$s^{OGGM}_{2003}$', 'Leclercq (2014)'],
               loc='upper left', bbox_to_anchor=(1.02,1.6),
               fontsize=25)
    ax2.legend(handles=[l3,l4,l5,l6,l7], labels=[362.9,285.1,211.9,132.8,47.8],
               bbox_to_anchor=(1.04,1), loc='upper left',
               title=r'mb. bias $\beta^*_{mb}$', fontsize=25)

    add_at(ax1,'a', loc=1)
    add_at(ax2,'b', loc=1)
    add_at(ax3,'c', loc=1)

    ax1.set_ylabel('Altitude (m)')
    ax1.set_xlabel('Distance along the main flowline (km)')

    ax2.set_ylabel('Altitude (m)')
    ax2.set_xlabel('Distance along the main flowline (km)')

    ax3.set_ylabel('Length (km)')
    ax3.set_xlabel('Time (years)')

    ax1.set_title('t=1917')
    ax2.set_title('t=2003')
    plt.suptitle('RGI60-11.00897: Hintereisferner')

    plt.savefig(os.path.join(dir,'multi_median.png'))
    plt.show()


def get_ref_length_data(rgi_id):
    """Get the glacier lenght data from P. Leclercq's data base.

     https://folk.uio.no/paulwl/data.php

     For some glaciers only!
     """

    df = pd.read_csv('rgi_leclercq_links_2014_RGIV6.csv')
    df = df.loc[df.RGI_ID == rgi_id]
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
    df = df.dL

    return df



def merge_results(path):
    exp_df = pd.DataFrame()
    df = pd.DataFrame()
    for temp_bias in os.listdir(path):

        if temp_bias.startswith('temp'):

            dir = os.path.join(path, temp_bias)
            temp_bias = float(temp_bias.split('_')[-1])

            for f in os.listdir(dir):
                p = os.path.join(dir, f)

                if p.endswith('difference.pkl'):
                    res = pd.read_pickle(p)
                    res['temp_bias'] = temp_bias
                    df = df.append(res, sort=False)
                elif p.endswith('experiments.pkl'):
                    exp = pd.read_pickle(p, compression='gzip')
                    exp_df = exp_df.append(exp, ignore_index=True,sort=True)

    exp_df = exp_df.sort_values(by='bias')
    exp_df = exp_df[~exp_df.duplicated(subset=['rgi_id','temp_bias'])]
    exp_df.loc[:, 'region'] = exp_df.rgi_id.apply(lambda x: int(x.split('RGI60-')[-1].split('.')[0]))

    for rgi in exp_df.dropna().rgi_id.unique():
        i = exp_df[exp_df.rgi_id == rgi].area_diff.abs().idxmin()
        temp_bias = exp_df.loc[i, 'temp_bias']
        best = df[(df.index == rgi) & (df.temp_bias == temp_bias)]
        best.loc[:,'temp_bias'] = 'best'
        df = df.append(best)

    # save exp_df including "model" column
    exp_df.to_pickle(os.path.join(path,'experiment_df.pkl'))

    # remove "model" for quicker reading
    exp_df = exp_df.drop('model', axis=1)
    exp_df.to_pickle(os.path.join(path, 'quick_experiment_df.pkl'))

    # save df
    df.to_pickle(os.path.join(path,'df.pkl'))


if __name__ == '__main__':
    cfg.initialize()
    WORK_DIR = '/home/juliaeis/Dokumente/OGGM/work_dir/reconstruction/leclercq_diff2'
    PLOT_DIR= os.path.join(WORK_DIR, 'plots')

    #merge_results(WORK_DIR)

    exp_df = pd.read_pickle(os.path.join(WORK_DIR, 'experiment_df.pkl'))


    exp_df = exp_df.replace([np.inf, -np.inf], np.nan)
    exp_df = exp_df.dropna()
    exp_hef = exp_df[exp_df.rgi_id=='RGI60-11.00897']
    #plot_bias_density(exp_df, PLOT_DIR)

    # plot fitness with leclercq
    #read pickle for specific glaicer
    dir = os.path.join(PLOT_DIR,'HEF')
    df = read_multi_pickle(dir)
    lec = get_ref_length_data('RGI60-11.00897')

    #plot_multi_fitness(df, exp_hef, lec, dir)
    #plot_mulit_median(df, exp_hef, lec, dir)

    df = pd.read_pickle(os.path.join(WORK_DIR, 'df.pkl'))
    df = df.replace([np.inf, -np.inf], np.nan)
    df = df.dropna()
    print(df.columns)
    error_df = pd.DataFrame()

    df = df[df.temp_bias==0.0]

    error_df['n'] = (df.groupby('region')['rmse'].count())
    error_df['rmse_mean'] = df.groupby('region')['rmse'].mean().round(2)
    error_df['rmse_std'] = df.groupby('region')['rmse'].std().round(2)

    error_df['rmspe_mean'] = df.groupby('region')['rmspe'].mean().round(2)
    error_df['rmspe_std'] = df.groupby('region')['rmspe'].std().round(2)

    error_df['error_mean'] = df.groupby('region')['error'].mean().round(2)
    error_df['error_std'] = df.groupby('region')['error'].std().round(2)

    error_df['perc_mean'] = df.groupby('region')['perc_error'].mean().round(2)
    error_df['perc_std'] = df.groupby('region')['perc_error'].std().round(2)

    error_df['max_diff_mean'] = df.groupby('region')['max_diff'].mean().round(2)
    error_df['max_diff_std'] = df.groupby('region')['max_diff'].std().round(2)
    error_df = error_df.astype({'n': 'int32'})
    error_df.to_latex('error.tex')

    print(error_df.n.sum())
    print(df.rmse.mean(), df.rmse.std())
    print(df.rmspe.mean(), df.rmspe.std())
    print(df.error.mean(), df.error.std())
    print(df.perc_error.mean(), df.perc_error.std())
    print(df.max_diff.mean(), df.max_diff.std())