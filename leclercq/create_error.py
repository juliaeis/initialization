import pandas as pd
import os
import numpy as np
import oggm
from oggm import utils
import seaborn as sns
import geopandas as gpd
import matplotlib.pyplot as plt

import sys
sys.path.append('../')
from initialization.core import *
from paper.plots_paper import *

mpl.rcParams['axes.linewidth'] = 3
mpl.rcParams['xtick.major.width'] = 3
mpl.rcParams['ytick.major.width'] = 3
mpl.rcParams['ytick.minor.width'] = 2
mpl.rcParams['font.size'] =25
mpl.rcParams['font.weight'] = 'medium'
mpl.rcParams['axes.labelweight'] = 'medium'
mpl.rcParams['legend.fontsize'] = 25 #30
mpl.rcParams['lines.linewidth'] = 3

def plot_mb_bias(df, oggm_df, plot_dir):

    df.loc[:, 'temp_bias'] = df.index.get_level_values(1)

    plt.figure(figsize=(15, 15))
    grid = plt.GridSpec(2, 1, hspace=0.3, wspace=0.2)
    ax2 = plt.subplot(grid[0])
    ax1 = plt.subplot(grid[1])

    cmap = plt.get_cmap("tab10")
    diff = []
    p1_list = []
    for j,col in enumerate(np.sort(df.temp_bias.unique())):
        if j>=3:
            j=j+1
        sub = df[df.temp_bias==col]
        label = col
        sns.distplot(sub.mb_bias, hist=False, kde=True,
                     kde_kws={'shade': True, 'linewidth': 3, 'alpha': 0.15},
                     bins=30, ax=ax2)
        ax2.plot(ax2.lines[-1].get_xydata()[:, 0],
                 ax2.lines[-1].get_xydata()[:, 1], color=cmap(j),
                 label=label)
        ax2.axvline(x=sub.mb_bias.mean(), color=cmap(j), linestyle=':')

    sns.distplot(oggm_df.bias, hist=False, kde=True, ax=ax2, color='C3',label='')
    ax2.axvline(x=oggm_df.bias.mean(), color=cmap(3), linestyle=':')
    ax1.text(-58, 1.1e3, 'n:' + str(len(diff)))
    ax1.hist(df.relative_diff)

    ax1.set_yscale('log')
    ax1.set_xscale('symlog',linthreshx=1e-3)

    ax1.set_ylabel('Frequency', labelpad=30)
    ax1.set_xlabel(r'Relative area difference ($\%$)')
    #ax1.set_xlim(-0.07, 0.07)
    ax2.set_ylim(0, None)
    ax2.set_yticks([0, 1e-3, 2e-3])
    ax2.set_xticks(np.arange(-1500, 2000, 500))
    from matplotlib.ticker import LogLocator

    ax1.yaxis.set_minor_locator(LogLocator(base=10, subs='auto'))
    ax1.tick_params(which='major', length=6, width=3)
    ax1.tick_params(which='minor', length=4, color='k', width=2)

    add_at(ax2, 'a')
    add_at(ax1, 'b')

    ax1.text(0.05,900,'n:'+str(len(df.relative_diff.dropna())))

    ax1.grid()
    ax2.grid()
    ax2.set_xlabel(  r'Optimal mass balance balance bias $\beta^*_{mb}$ (mm w.e. yr$^{-1}$)')
    ax2.set_ylabel('Density')
    leg1 = ax2.legend(title=r'temp. bias $\widetilde{\beta}$')
    legend_elements = [Line2D([0], [0], color='C3', label=' OGGM')]
    ax2.legend(handles=legend_elements, loc=7, bbox_to_anchor=(1,0.23))
    ax2.add_artist(leg1)
    plt.savefig(os.path.join(plot_dir, 'mb_bias.pdf'), dpi=300)
    plt.show()



def create_model_lec_df(home, repeat=False):

    model_df = pd.DataFrame()
    p = os.path.join(home,'lec_model_df.pkl')
    if repeat or not os.path.isfile(p):

        for d in [d for d in os.listdir(home) if d.startswith('temp')]:
            temp_bias = float(d.split('_')[-1])
            dir = os.path.join(home, d)
            if temp_bias <= 0:
                for file in [os.path.join(dir, d1) for d1 in os.listdir(dir) if
                             d1.endswith('.pkl')]:
                    df = pd.read_pickle(file, compression='gzip')
                    df = df.assign(temp_bias=temp_bias)
                    model_df = model_df.append(df, ignore_index=False, sort=True)

        model_df.index.name = 'rgi_id'
        model_df.loc[:, 'rgi_date'] = model_df.ex_mod.apply(
            lambda x: int(x.area_km2_ts().index[-1]))
        model_df.loc[:, 'region'] = model_df.index.map(
            lambda x: int(x.split('RGI60-')[-1].split('.')[0]))

        for region in model_df.region.unique():
            ids = model_df[model_df.region == region].index.get_level_values(0)

            # RGI file
            path = utils.get_rgi_region_file(str(region).zfill(2), version='61')
            rgidf = gpd.read_file(path).set_index('RGIId')
            rgidf = rgidf[rgidf.index.isin(ids)]

            model_df.loc[rgidf.index, 'rgi_area'] = rgidf.Area

        model_df = model_df.reset_index().set_index(
            ['rgi_id', 'temp_bias']).sort_index()
        model_df.loc[:, 'ex_area'] = model_df.ex_mod.apply(
            lambda x: x.area_km2_ts().values[-1])
        model_df.to_pickle(p,compression='gzip')
        return model_df
    else:
        return pd.read_pickle(p, compression='gzip')

def create_lec_df(home, repeat=False):

    lec_df = pd.DataFrame()
    p = os.path.join(home,'lec_df.csv')
    if repeat or not os.path.isfile(p):

        dir = os.path.join(home, 'temp_0')

        for file in [os.path.join(dir, d1) for d1 in os.listdir(dir) if
                     'leclercq' in d1]:
            df = pd.read_csv(file)

            lec_df = lec_df.append(df, ignore_index=False, sort=True)
        lec_df = lec_df.rename(columns={"Unnamed: 0": "rgi_id"}).set_index('rgi_id').rename(columns=int)
        lec_df = lec_df.sort_index()

        lec_df.to_csv(p)
        return lec_df

    else:
        return pd.read_csv(p)


def create_lec_median_df(home, repeat=False):
    lec_median_df = pd.DataFrame()
    p = os.path.join(home, 'lec_median_df.csv')
    if repeat or not os.path.isfile(p):

        for d in [d for d in os.listdir(home) if d.startswith('temp')]:
            temp_bias = float(d.split('_')[-1])
            print(temp_bias)
            dir = os.path.join(home, d)

            for file in [os.path.join(dir, d1) for d1 in os.listdir(dir) if
                         'median' in d1]:
                df = pd.read_csv(file)
                df = df.assign(temp_bias=temp_bias)

                lec_median_df = lec_median_df.append(df, ignore_index=False, sort=True)

        lec_median_df = lec_median_df.rename(columns={"Unnamed: 0": "rgi_id"}).set_index(
            ['rgi_id','temp_bias']).rename(columns=int)
        lec_median_df = lec_median_df.sort_index()

        lec_median_df.to_csv(p)
        return lec_median_df

    else:
        return pd.read_csv(p)


if __name__ == '__main__':

    all_ex = pd.DataFrame()

    home = '/home/juliaeis/Dokumente/OGGM/work_dir/reconstruction/new_leclercq'
    lec_model_df = create_model_lec_df(home,repeat=False).drop(['mod','ex_mod'], axis=1)
    lec_model_df.loc[:,'area_diff'] = lec_model_df.rgi_area-lec_model_df.ex_area
    lec_model_df.loc[:,'relative_diff'] = lec_model_df.area_diff/lec_model_df.rgi_area
    lec_model_df = lec_model_df.rename(columns={'bias':'mb_bias'})

    lec_df = create_lec_df(home,repeat=False)
    lec_median_df = create_lec_median_df(home, repeat=True)
    print(lec_median_df.head())


    all_ex = all_ex.append(lec_model_df[['mb_bias','area_diff','relative_diff']])

    wgms = pd.read_csv('/home/juliaeis/Dokumente/OGGM/work_dir/reconstruction/wgms/wgms_experiement_df.csv').set_index(['rgi_id','temp_bias'])

    all_ex = all_ex.append(wgms[['mb_bias', 'area_diff','relative_diff']])
    # filter all double glaciers
    all_ex = all_ex[~all_ex.index.duplicated(keep='last')]

    #all_ex = all_ex.drop_duplicates(keep='last')

    print(all_ex.mean(level=1).round({'mb_bias': 1, 'area_diff': 4, 'relative_diff':4}))
    print(all_ex.count(level=1))

    all_ex.to_csv(os.path.join(home,'all_experiments.csv'))

    oggm_df = pd.read_csv('oggm_ref_tstars_rgi6_cru4.csv')
    #plot_mb_bias(all_ex,oggm_df,os.path.join(home,'plots'))

