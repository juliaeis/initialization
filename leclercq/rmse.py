import os
import pandas as pd
import numpy as np
import seaborn as sns
from copy import deepcopy
import matplotlib.pyplot as plt
import matplotlib as mpl
from oggm import cfg, tasks, workflow, utils
from oggm.core.flowline import FileModel
import geopandas as gpd
#from  leclercq_plots import *
from matplotlib.pyplot import cm

mpl.rcParams['axes.linewidth'] = 3
mpl.rcParams['xtick.major.width'] = 3
mpl.rcParams['ytick.major.width'] = 3
mpl.rcParams['font.size'] =25
mpl.rcParams['font.weight'] = 'medium'
mpl.rcParams['axes.labelweight'] = 'medium'
mpl.rcParams['legend.fontsize'] = 15 #30
mpl.rcParams['legend.title_fontsize'] = 15
mpl.rcParams['lines.linewidth'] = 3


def add_at(ax, t, loc=2):
    fp = dict(size=20)
    _at = AnchoredText(t, loc=loc, prop=fp, borderpad=0)
    _at.patch.set_linewidth(3)
    ax.add_artist(_at)
    return _at

def read_temp_pkl(dir):

    try:
        df = pd.read_pickle(os.path.join(dir, 'error_df.pkl'))
    except:
        df = pd.DataFrame()
        ex = pd.DataFrame()
        for file in [f for f in os.listdir(dir) if
                     f.endswith('.pkl') and f.startswith('temp')]:
            df = df.append(pd.read_pickle(os.path.join(dir, file)),
                           ignore_index=True)


        for d in [f for f in os.listdir(dir) if not f.endswith('.pkl')]:
            temp_bias = float(d.split('_')[-1])
            for f in os.listdir(os.path.join(dir, d)):
                if f.endswith('experiments.pkl'):
                    ex_df = pd.read_pickle(os.path.join(dir, d, f),
                                           compression='gzip')

                    ex_df.loc[:, 'temp_bias'] = temp_bias
                    ex = ex.append(ex_df, ignore_index=True, sort=False)
                    pass

        df = df.dropna(subset=['region'])
        df.region = df.region.apply(lambda x: int(x))
        df.temp_bias = df.temp_bias.apply(lambda x: float(x))
        df.rmse = df.rmse / 1000
        df.error = df.error / 1000
        df.perc_error = df.perc_error * 100

        ex = ex.dropna()
        for gdir in ex.rgi_id.unique():
            temp_bias = ex.loc[ex[ex.rgi_id == gdir].fitness.idxmin(), 'temp_bias']
            df2 = df[(df.rgi == gdir) & (df.temp_bias == temp_bias)]
            df2.temp_bias = 'best'
            df = df.append(df2, ignore_index=True, sort=False)
        df.to_pickle(os.path.join(dir, 'error_df.pkl'))
    return df

def plot_rmse(df, dir):
    fig, [ax1, ax2] = plt.subplots(2, 1, figsize=(25,18))
    flierprops = dict(marker='o', markersize=4)
    meanprops = dict(marker='x', markeredgecolor="silver")
    sns.boxplot(x="region", y="rmse", hue="temp_bias",
                data=df[df.region <= 10], palette="Set1", ax=ax1,
                showmeans=True, flierprops=flierprops, meanprops=meanprops,
                hue_order=[-1, -0.75, -0.5, -0.25, 0, 'best'])
    sns.boxplot(x="region", y="rmse", hue="temp_bias",
                data=df[df.region > 10], palette="Set1", ax=ax2,
                showmeans=True, flierprops=flierprops, meanprops=meanprops,
                hue_order=[-1, -0.75, -0.5, -0.25, 0, 'best'])

    nobs = df[df.temp_bias == 'best'][
        'region'].value_counts().sort_index().values
    nobs = [str(x) for x in nobs.tolist()]
    nobs = ["n: " + i for i in nobs]

    # Add it to the plot
    pos = range(len(nobs))
    for tick, label in zip(pos, ax1.get_xticklabels()):
        ax1.text(pos[tick], -0.2, nobs[tick],
                 horizontalalignment='center', size='x-small', color='k',
                 weight='semibold')
    for tick, label in zip(pos, ax2.get_xticklabels()):
        ax2.text(pos[tick], -0.2, nobs[tick + 8],
                 horizontalalignment='center', size='x-small', color='k',
                 weight='semibold')

    #ax1.set_yticks(range(1, 5))
    #ax2.set_yticks(range(1, 5))

    ax1.set_xlabel('Region')
    ax2.set_xlabel('Region')
    ax1.set_ylabel('Root Mean Square Error (km)')
    ax2.set_ylabel('Root Mean Square Error (km)')

    ax1.grid()
    ax2.grid()


    ax1.legend(loc='center left', bbox_to_anchor=(1, 1), title='temperature bias')
    ax2.legend_.remove()


    ax1.set_title('Comparisons with Leclercq')
    plt.savefig(os.path.join(dir,'rmse.png'), dpi=300)
    plt.show()
    #plt.close()


def plot_rmspe(df, dir):
    fig, [ax1, ax2] = plt.subplots(2, 1, figsize=(25,18))
    flierprops = dict(marker='o', markersize=4)
    meanprops = dict(marker='x', markeredgecolor="silver")
    sns.boxplot(x="region", y="rmspe", hue="temp_bias",
                data=df[df.region <= 10], palette="Set1", ax=ax1,
                showmeans=True, flierprops=flierprops, meanprops=meanprops,
                hue_order=[-1, -0.75, -0.5, -0.25, 0, 'best'])
    sns.boxplot(x="region", y="rmspe", hue="temp_bias",
                data=df[df.region > 10], palette="Set1", ax=ax2,
                showmeans=True, flierprops=flierprops, meanprops=meanprops,
                hue_order=[-1, -0.75, -0.5, -0.25, 0, 'best'])

    nobs = df[df.temp_bias == 'best'][
        'region'].value_counts().sort_index().values
    nobs = [str(x) for x in nobs.tolist()]
    nobs = ["n: " + i for i in nobs]

    # Add it to the plot
    pos = range(len(nobs))
    for tick, label in zip(pos, ax1.get_xticklabels()):
        ax1.text(pos[tick], -7, nobs[tick],
                 horizontalalignment='center', size='x-small', color='k',
                 weight='semibold')
    for tick, label in zip(pos, ax2.get_xticklabels()):
        ax2.text(pos[tick], -7, nobs[tick + 8],
                 horizontalalignment='center', size='x-small', color='k',
                 weight='semibold')

    ax1.set_yticks(np.arange(0, 175, 25))
    ax2.set_yticks(np.arange(0, 175, 25))

    ax1.set_xlabel('Region')
    ax2.set_xlabel('Region')
    ax1.set_ylabel('Root Mean Square Percentage Error (%)')
    ax2.set_ylabel('Root Mean Square Percentage Error (%)')

    ax1.grid()
    ax2.grid()


    ax1.legend(loc='center left', bbox_to_anchor=(1, 1), title='temperature bias')
    ax2.legend_.remove()

    ax1.set_title('Comparisons with Leclercq')
    plt.savefig(os.path.join(dir,'rmspe.png'), dpi=300)
    #plt.show()
    plt.close()

def plot_mean_error(df, dir):

    fig, [ax1, ax2] = plt.subplots(2, 1, figsize=(25, 18))
    flierprops = dict(marker='o', markersize=4)
    meanprops = dict(marker='x', markeredgecolor="silver")
    sns.boxplot(x="region", y="error", hue="temp_bias",
                data=df[df.region <= 10], palette="Set1", ax=ax1,
                showmeans=True, flierprops=flierprops, meanprops=meanprops,
                hue_order=[-1, -0.75, -0.5, -0.25, 0, 'best'])
    sns.boxplot(x="region", y="error", hue="temp_bias",
                data=df[df.region > 10], palette="Set1", ax=ax2,
                showmeans=True, flierprops=flierprops, meanprops=meanprops,
                hue_order=[-1, -0.75, -0.5, -0.25, 0, 'best'])

    nobs = df[df.temp_bias == 'best'][
        'region'].value_counts().sort_index().values
    nobs = [str(x) for x in nobs.tolist()]
    nobs = ["n: " + i for i in nobs]

    # Add it to the plot
    pos = range(len(nobs))
    for tick, label in zip(pos, ax1.get_xticklabels()):
        ax1.text(pos[tick], -9, nobs[tick],
                 horizontalalignment='center', size='x-small', color='k',
                 weight='semibold')
    for tick, label in zip(pos, ax2.get_xticklabels()):
        ax2.text(pos[tick], -5.5, nobs[tick + 8],
                 horizontalalignment='center', size='x-small', color='k',
                 weight='semibold')

    #ax1.set_yticks(np.arange(-7, 8, 1))
    #ax2.set_yticks(np.arange(-7, 8, 1))

    ax1.set_xlabel('Region')
    ax2.set_xlabel('Region')
    ax1.set_ylabel('Mean Error (km)')
    ax2.set_ylabel('Mean Error (km)')

    ax1.legend(loc='center left', bbox_to_anchor=(1, 0.81), title='temp. bias')
    ax2.legend_.remove()

    ax1.grid()
    ax2.grid()

    ax1.set_title('Comparisons with Leclercq')
    plt.savefig(os.path.join(dir, 'mean_error.png'), dpi=300)
    plt.show()

def f(x):
    a = x[0]
    b = x[1]

    import math
    if math.isnan(b):
        b = 'Nan'
        return str(int(round(a)))+' '+b
    else:
        return str(int(round(a)))+ ' ' + str(int(round(b)))


if __name__ == '__main__':

    cfg.initialize()

    ON_CLUSTER = False
    REGION = '11'

    # Local paths
    if ON_CLUSTER:
        WORKING_DIR = os.environ.get("S_WORKDIR")
        cfg.PATHS['working_dir'] = WORKING_DIR
        REGION = str(os.environ.get('I')).zfill(2)
    else:
        WORKING_DIR = '/home/juliaeis/Dokumente/OGGM/work_dir/reconstruction/leclercq_diff2'
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

    cfg.PARAMS['run_mb_calibration'] = False
    cfg.PARAMS['optimize_inversion_params'] = False
    cfg.PARAMS['dl_verify'] = False


    # add to BASENAMES
    _doc = 'contains observed and searched glacier from synthetic experiment to find intial state'
    cfg.BASENAMES['synthetic_experiment'] = ('synthetic_experiment.pkl', _doc)

    # We use intersects
    db = utils.get_rgi_intersects_region_file(version='61', region=REGION)
    cfg.set_intersects_db(db)

    # RGI file
    path = utils.get_rgi_region_file(REGION, version='61')
    rgidf = gpd.read_file(path)

    # only the ones with leclercq observation

    #rgidf = rgidf[rgidf.RGIId.isin(['RGI60-11.00897'])]
    #gdir = workflow.init_glacier_regions(rgidf)[0]




    #print(os.listdir(hef))
    '''
    df = pd.DataFrame()
    exp_df = pd.DataFrame()
    for temp_bias in os.listdir(cfg.PATHS['working_dir']):

        if temp_bias.startswith('temp'):
            dir = os.path.join(cfg.PATHS['working_dir'],temp_bias)
            temp_bias = float(temp_bias.split('_')[-1])
            for f in os.listdir(dir):
                p = os.path.join(dir,f)
                if p.endswith('difference.pkl'):
                    res = pd.read_pickle(p)
                    res['temp_bias'] = temp_bias
                    df = df.append(res)
                elif p.endswith('experiments.pkl'):
                    exp = pd.read_pickle(p, compression='gzip')
                    exp_df = exp_df.append(exp, ignore_index=True)
    exp_df = exp_df.dropna()
    for rgi in exp_df.rgi_id.unique():
        i = exp_df[exp_df.rgi_id==rgi].area_diff.abs().idxmin()
        temp_bias = exp_df.loc[i,'temp_bias']
        best = df[(df.index==rgi) &(df.temp_bias==temp_bias)]
        best.temp_bias = 'best'
        df = df.append(best)

    df.to_pickle(os.path.join(cfg.PATHS['working_dir'], 'error_df.pkl'))
    exp_df.to_pickle(os.path.join(cfg.PATHS['working_dir'], 'experiment_df.pkl'))

    '''
    exp_df =  pd.read_pickle(os.path.join(cfg.PATHS['working_dir'], 'experiment_df.pkl'))
    exp_df.loc[:,'region'] = exp_df.rgi_id.apply(lambda x: int(x.split('RGI60-')[-1].split('.')[0]))
    '''
    grid = plt.GridSpec(2, 1, hspace=0.4, wspace=0.2)
    ax1 = plt.subplot(grid[0])
    ax2 = plt.subplot(grid[1])

    cmap = plt.get_cmap("tab10")

    for i, temp_bias in enumerate(np.sort(exp_df.temp_bias.unique())):

        # Subset to the airline
        subset = exp_df[exp_df['temp_bias'] == temp_bias]

        # Draw the density plot
        sns.distplot(subset['bias'], hist=False, kde=True,
                     kde_kws={'shade': True,'linewidth': 3},
                     label=temp_bias,bins=30, ax=ax2 )

        ax2.axvline(x= subset.bias.mean(), color=cmap(i), linestyle=':')


    exp_df.area_diff.plot.hist(ax=ax1, bins=40)
    ax1.set_yscale('log')
    #ax1.set_xscale('symlog',linthreshx=1e-2)
    ax1.set_ylabel('Frequency', labelpad=30)
    ax1.set_xlabel(r'Area difference (km$^2$)')
    ax1.set_xlim(-1.05,1.05)

    ax1.grid()


    ax2.grid()
    ax2.set_xlabel(r'Optimal mass balance balance bias $\beta^*_{mb}$')
    ax2.set_ylabel('Density')
    plt.legend(title='temp. bias')

    plt.show()


    df = pd.DataFrame()
    for temp_bias in exp_df.temp_bias.unique():

        df.loc[:,'mean'] = exp_df[exp_df.temp_bias==temp_bias].groupby(['region']).bias.mean()
        df.loc[:,'std'] = exp_df[exp_df.temp_bias==temp_bias].groupby('region').bias.std()

        df.loc['all', 'mean'] = exp_df[exp_df.temp_bias==temp_bias].bias.mean()
        df.loc['all', 'std'] = exp_df[exp_df.temp_bias == temp_bias].bias.std()
        df.loc[:, temp_bias] = df[['mean', 'std']].apply(f, axis=1)
        df = df.drop(['mean','std'], axis=1)

    df.loc[:, 'n'] = exp_df[exp_df.temp_bias == 0].groupby('region').bias.count()
    df.loc['all', 'n'] = len(exp_df[exp_df.temp_bias == 0])

    df = df[['n',-1.0,-0.75,-0.5,-0.25,0]]
    print(df)
    df.to_latex('table.tex')
   
    df = pd.read_pickle(os.path.join(cfg.PATHS['working_dir'], 'error_df.pkl'))
    df.region = pd.to_numeric(df.region)
    df.rmse = df.rmse/1000



    #plot_rmse(df, cfg.PATHS['plot_dir'])
    #plot_rmspe(df, cfg.PATHS['plot_dir'])
    #plot_mean_error(df, cfg.PATHS['plot_dir'])

    # read leclercq links
    lec = pd.read_csv('rgi_leclercq_links_2014_RGIV6.csv')
    lec['REGION'] = lec.RGI_ID.apply(lambda x: x.split('-')[-1].split('.')[0])
    print(lec)

    exp_df = pd.read_pickle(os.path.join(cfg.PATHS['working_dir'], 'experiment_df.pkl'))
    print(len(exp_df.rgi_id.unique()))
    '''
    exp_df = exp_df[exp_df.rgi_id == 'RGI60-11.00897'].sort_values(by='temp_bias')
    grid = plt.GridSpec(1, 2, hspace=0.2, wspace=0)
    plt.figure(figsize=(20,10))
    ax1 = plt.subplot(grid[0])
    ax2 = plt.subplot(grid[1], sharey=ax1)
    cmap = plt.get_cmap("tab10")
    for  j,i in enumerate(exp_df.index):
        if j >= 3:
            j = j + 1
        color = cmap(j)
        # plot past climate
        temp_bias = str(exp_df.loc[i,'temp_bias'])
        if exp_df.loc[i,'temp_bias']== -1.0:
            temp_bias='-1.00'
        if exp_df.loc[i,'temp_bias']== -0.5:
            temp_bias='-0.50'
        if exp_df.loc[i, 'temp_bias'] == 0.0:
            temp_bias = ' 0.00'
        exp_df.loc[i,'model'].area_km2_ts().plot(ax=ax2,color=color,label= temp_bias+r'$\qquad$'+str(exp_df.loc[i,'bias']))

        # plot_random climate
        if exp_df.loc[i, 'temp_bias'] == -1.0:
            temp_bias = '-1'
        if exp_df.loc[i, 'temp_bias'] == 0.0:
            temp_bias = '0'
        if exp_df.loc[i,'temp_bias']== -0.5:
            temp_bias='-0.5'
        name = 'model_run_random_experiment_'+str(exp_df.loc[i,'temp_bias'])+'_'+str(exp_df.loc[i,'bias'])+'.nc'
        file = os.path.join(cfg.PATHS['working_dir'],'temp_'+temp_bias, 'per_glacier','RGI60-11','RGI60-11.00','RGI60-11.00897',name)
        model = FileModel(file)
        model.area_km2_ts().plot(ax=ax1, color=color)
    add_at(ax1,'a')
    add_at(ax2, 'b')

    ax1.grid()
    ax2.grid()
    ax1.set_xticks([0,100,200,300,400,500])

    ax1.set_xlabel('Time(years)')
    ax2.set_xlabel('Time(years)')
    ax1.set_ylabel(r'Area (km$^2$)')
    ax1.set_title('Random climate forcing')
    ax2.set_title(r'CRU 1917-2003 (no $\widetilde{\beta}$)')

    plt.suptitle('RGI60-11.00897: Hintereisferner')
    plt.legend(title='\t'+r'temp.bias $\widetilde{\beta}$'+r'$\qquad$'+r'mb bias $\beta^*_{mb}$')
    plt.savefig(os.path.join(cfg.PATHS['plot_dir'],'experiment.png'),dpi=200)
    plt.show()

