import sys
sys.path.append('../')
from initialization.core import *
from paper.plots_paper import *
from oggm.core.massbalance import MultipleFlowlineMassBalance, PastMassBalance
import seaborn as sns
import geopandas as gpd
from oggm import cfg, utils
pd.options.mode.chained_assignment = None
from leclercq import rmse
from copy import deepcopy

def create_wgms_df(dir):

    model_df = pd.DataFrame()
    for dir in os.listdir(cfg.PATHS['working_dir']):
        if dir.startswith('temp'):
            p = os.path.join(cfg.PATHS['working_dir'], dir)
            df = pd.DataFrame()
            for file in os.listdir(p):
                if file.endswith('.pkl'):
                    rp = os.path.join(p, file)
                    df = df.append(pd.read_pickle(rp), ignore_index=False)
            model_df = model_df.append(df)
    model_df = model_df.sort_index()

    model_df.loc[:,'relative_diff'] = model_df.area_diff/model_df.rgi_area
    model_df.index.name='rgi_id'
    model_df.to_pickle(os.path.join(cfg.PATHS['working_dir'],'wgms_experiement_df.pkl'))
    model_df = model_df.drop('ex_mod',axis=1)
    model_df.to_csv(os.path.join(cfg.PATHS['working_dir'],'wgms_experiement_df.csv'))

def wgms_get_mb_errors(dir):
    model_df = pd.DataFrame()
    for dir in os.listdir(cfg.PATHS['working_dir']):
        if dir.startswith('temp'):
            p = os.path.join(cfg.PATHS['working_dir'], dir)
            df = pd.DataFrame()
            for file in os.listdir(p):
                if file.endswith('difference.csv'):
                    rp = os.path.join(p, file)
                    df = df.append(pd.read_csv(rp), ignore_index=False)
            model_df = model_df.append(df)

    model_df = model_df.sort_index()

    return model_df

def merge_wgms_lec(wgms, lec, repeat=False):

    if repeat:

        lec.loc[:, 'region'] = lec.rgi_id.apply(lambda x: x.split('RGI60-')[-1].split('.')[0])

        # sort values by RGI Id and set index
        lec = lec.sort_values(by='rgi_id')
        lec = lec.set_index('rgi_id').drop('Unnamed: 0', axis=1)
        lec = lec[~lec.index.isin(wgms_df.index)]

        lec.loc[:, 'rgi_area'] = np.nan
        lec.loc[:, 'rgi_date'] = np.nan
        for region in lec.region.unique():
            ids = lec[lec.region == region].drop('region', axis=1)

            # RGI file
            path = utils.get_rgi_region_file(str(region).zfill(2), version='61')
            rgidf = gpd.read_file(path).set_index('RGIId')
            rgidf.index.name = 'rgi_id'
            df2 = pd.DataFrame()
            df2.loc[:, 'rgi_area'] = rgidf[rgidf.index.isin(ids.index)].Area

            date = rgidf[rgidf.index.isin(ids.index)].BgnDate.apply(lambda x:int(str(x)[:4]))
            df2.loc[:, 'rgi_date'] = date
            print(df2.rgi_date)
            lec.update(df2)

        lec.loc[:, 'relative_diff'] = lec.area_diff / lec.rgi_area
        lec = lec.rename(columns={"bias": "mb_bias"})

        all = wgms.append(lec)
        all = all[['region','rgi_date', 'rgi_area','temp_bias',
                   'iterations','mb_bias', 'area_diff',  'relative_diff', 'error']]
        all.to_csv(os.path.join(cfg.PATHS['working_dir'],'all_experiment_df.csv'))
    else:
        all = pd.read_csv(os.path.join(cfg.PATHS['working_dir'],'all_experiment_df.csv'))
    return all


def plot_mb_bias(df, oggm_df, plot_dir):

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
    ax1.hist(df.relative_diff, bins=30)

    ax1.set_yscale('log')
    # ax1.set_xscale('symlog',linthreshx=1e-2)
    ax1.set_ylabel('Frequency', labelpad=30)
    ax1.set_xlabel(r'Relative area difference ($\%$)')
    ax1.set_xlim(-0.07, 0.07)
    ax2.set_ylim(0, None)
    ax2.set_yticks([0, 1e-3, 2e-3])
    ax2.set_xticks(np.arange(-1500, 2000, 500))
    from matplotlib.ticker import LogLocator

    ax1.yaxis.set_minor_locator(LogLocator(base=10, subs='auto'))
    ax1.tick_params(which='major', length=6, width=3)
    ax1.tick_params(which='minor', length=4, color='k', width=2)

    add_at(ax2, 'a')
    add_at(ax1, 'b')

    ax1.text(0.054,1400,'n:'+str(len(df.relative_diff.dropna())))

    ax1.grid()
    ax2.grid()
    ax2.set_xlabel(  r'Optimal mass balance balance bias $\beta^*_{mb}$ (mm w.e. yr$^{-1}$)')
    ax2.set_ylabel('Density')
    leg1 = ax2.legend(title=r'temp. bias $\widetilde{\beta}$')
    legend_elements = [Line2D([0], [0], color='C3', label=' OGGM')]
    ax2.legend(handles=legend_elements, loc=7, bbox_to_anchor=(1,0.54))
    ax2.add_artist(leg1)
    plt.savefig(os.path.join(plot_dir, 'mb_bias.pdf'), dpi=300)
    plt.show()

def plot_error_density(df,plot_dir):
    cmap = plt.get_cmap("tab10")
    diff = []
    for j, col in enumerate(np.sort(df.temp_bias.unique())):
        if j >= 3:
            j = j + 1
        sub = df[df.temp_bias == col]
        label = col
        sns.distplot(sub.error, hist=False, kde=True,
                     kde_kws={'shade': True, 'linewidth': 3, 'alpha': 0.15},
                     bins=30, ax=ax2)
        ax2.plot(ax2.lines[-1].get_xydata()[:, 0],
                 ax2.lines[-1].get_xydata()[:, 1], color=cmap(j),
                 label=label)
        ax2.axvline(x=sub.error.mean(), color=cmap(j), linestyle=':')
    plt.xlabel('bias (m)')
    plt.ylabel('density')
    plt.grid()
    plt.xlim((-3, 3))
    plt.legend(title='temp. bias')
    plt.savefig(os.path.join(plot_dir, 'density.png'))
    plt.show()

def plot_error_hist(df,plot_dir):

    bins = np.arange(14) * 0.400 - 2.600
    plt.figure(figsize=(25, 8))
    grid = plt.GridSpec(1, 5, hspace=0.3, wspace=0.2)
    cmap = plt.get_cmap("tab10")
    for i, temp_bias in enumerate(np.sort(df.temp_bias.unique())):
        if i != 0:
            ax = plt.subplot(grid[i], sharey=ax)
        else:
            ax = plt.subplot(grid[i])
        if i >= 3:
            i = i + 1
        df[df.temp_bias == temp_bias].error.hist(bins=bins, ax=ax,
                                                 color=cmap(i))
        plt.axvline(x=df[df.temp_bias == temp_bias].error.mean(), color='k')
        plt.axvline(x=df[df.temp_bias == temp_bias].error.median(), color='k',
                    linestyle=':')
        plt.axvline(x=df[df.temp_bias == temp_bias].error.quantile(0.05),
                    color='grey')
        plt.axvline(x=df[df.temp_bias == temp_bias].error.quantile(0.95),
                    color='grey')
        ax.set_title(r'$\widetilde{\beta}$:' + str(temp_bias))
        ax.set_xlabel('model bias (m)')
        ax.set_xticks([-2,-1,0,1,2])
        ax.text(-2.7,85,'N='+str(len(df[df.temp_bias == temp_bias])))
    plt.savefig(os.path.join(plot_dir, 'histgrams.png'))
    plt.show()

def plot_error_box(df,col,plot_dir):
    plt.figure(figsize=(20, 15))
    grid = plt.GridSpec(2, 1, hspace=0.3, wspace=0.2)
    ax1 = plt.subplot(grid[0])
    ax2 = plt.subplot(grid[1], sharey=ax1)
    sns.boxplot(x="region", y=col, hue='temp_bias',
                data=df[df.region <= 10], ax=ax1, palette=cmap)

    sns.boxplot(x="region", y=col, hue='temp_bias',
                data=df[df.region > 10], ax=ax2, palette=cmap)

    # show number of glaciers in Region
    nobs1 = df[(df.region <= 10) & (
    df.temp_bias == 0.0)].region.value_counts().sort_index().values
    nobs1 = [str(x) for x in nobs1.tolist()]
    nobs1 = ["n: " + i for i in nobs1]
    pos = range(len(nobs1))
    for tick, label in zip(pos, ax1.get_xticklabels()):
        ax1.text(pos[tick], -0.05, nobs1[tick],
                 horizontalalignment='center', size='x-small', color='k',
                 weight='semibold')
    nobs1 = df[(df.region > 10) & (
    df.temp_bias == 0.0)].region.value_counts().sort_index().values
    nobs1 = [str(x) for x in nobs1.tolist()]
    nobs1 = ["n: " + i for i in nobs1]
    pos = range(len(nobs1))
    for tick, label in zip(pos, ax1.get_xticklabels()):
        ax2.text(pos[tick], -0.05, nobs1[tick],
                 horizontalalignment='center', size='x-small', color='k',
                 weight='semibold')
    df.loc[:,'region_lab'] = df.region.apply(lambda x: str(int(x)).zfill(2))

    ax1.set_xticklabels(df.region_lab.unique()[:11])
    ax2.set_xticklabels(df.region_lab.unique()[9:])
    print(df[df.region==18].rmse.values)

    ax1.set_ylabel('RMSE (m)')
    ax2.set_ylabel('RMSE (m)')
    ax1.set_xlabel('Region')
    ax2.set_xlabel('Region')

    ax1.legend(bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.,
               title='temp. bias')
    ax2.legend(bbox_to_anchor=(5.04, 1), loc=2, borderaxespad=0.)
    ax1.set_ylim((-0.1,None))
    ax2.set_ylim((-0.1,None))
    ax1.grid()
    ax2.grid()
    plt.savefig(os.path.join(plot_dir, col+'_box.png'))
    plt.show()

if __name__ == '__main__':
    cfg.initialize()

    ON_CLUSTER = False

    # Local paths
    if ON_CLUSTER:
        WORKING_DIR = os.environ.get("S_WORKDIR")
        cfg.PATHS['working_dir'] = WORKING_DIR
        OUT_DIR = os.environ.get("OUTDIR")
        REGION = str(os.environ.get('REGION')).zfill(2)

    else:
        WORKING_DIR = '/home/juliaeis/Dokumente/OGGM/work_dir/reconstruction/wgms'
        OUT_DIR = WORKING_DIR
        cfg.PATHS['working_dir'] = WORKING_DIR
        utils.mkdir(WORKING_DIR, reset=False)
        REGION='05'

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
    cfg.PARAMS['dl_verify'] = True

    # add to BASENAMES
    _doc = 'contains observed and searched glacier from synthetic experiment to find intial state'
    cfg.BASENAMES['synthetic_experiment'] = ('synthetic_experiment.pkl', _doc)

    # We use intersects
    db = utils.get_rgi_intersects_region_file(version='61', region=REGION)
    cfg.set_intersects_db(db)

    # RGI file
    path = utils.get_rgi_region_file(REGION, version='61')
    rgidf = gpd.read_file(path)
    rgidf = rgidf.sort_values('Area', ascending=False)

    # exclude non-landterminating glaciers
    rgidf = rgidf[rgidf.TermType == 0]
    rgidf = rgidf[rgidf.Connect != 2]

    wgms_id = utils.get_ref_mb_glaciers_candidates()

    # Keep only the wgms reference glaciers
    rgidf = rgidf.loc[rgidf.RGIId.isin(wgms_id)]

    # initialize glaciers
    #gdirs = workflow.init_glacier_regions(rgidf)

    t_0 = 1917
    epsilon = 125

    # merge all experiment.pkl from each temp_bias together
    #create_wgms_df(cfg.PATHS['working_dir'])


    # reads files from WGMS and Leclercq runs
    wgms_df = pd.read_csv(os.path.join(cfg.PATHS['working_dir'], 'wgms_experiement_df.csv')).set_index('rgi_id')
    lec_df = pd.read_csv('/home/juliaeis/Dokumente/OGGM/work_dir/reconstruction/leclercq_diff2/quick_experiment_df.csv')

    all_df = merge_wgms_lec(wgms_df,lec_df, repeat=False).sort_values(by='rgi_id')
    all_df = all_df.replace([np.inf, -np.inf], np.nan)
    all_df = all_df.dropna(subset=['mb_bias'])

    oggm_df = pd.read_csv('oggm_ref_tstars_rgi6_cru4.csv')

    plot_mb_bias(all_df, oggm_df, cfg.PATHS['plot_dir'])

    '''
    df = wgms_get_mb_errors(cfg.PATHS['working_dir'])
    df = df.rename(columns={'Unnamed: 0':'rgi_id'})
    df = df.sort_values(by=['rgi_id', 'temp_bias']).set_index('rgi_id')
    df.error= df.error/1000
    df.rmse = df.rmse/1000

    df.to_csv(os.path.join(cfg.PATHS['working_dir'],'wgms_error.csv'))

    cmap = plt.get_cmap("tab10")
    cmap = [cmap(0), cmap(1), cmap(2),cmap(4), cmap(5)]

    #sns.boxplot(df.error, groupby=df.temp_bias)
    plt.figure(figsize=(15,15))
    ax = sns.violinplot(x="temp_bias", y="rmse", data=df, palette=cmap)
    plt.ylabel('bias (m)')
    plt.xlabel(r'temperature bias $\widetilde{\beta}$ (K)')
    plt.grid()
    plt.savefig(os.path.join(cfg.PATHS['plot_dir'],'rmse_box_by_temp.png'))
    plt.close()

    #plt.show()
    plot_error_box(df,'rmse',cfg.PATHS['plot_dir'])

    #plot_error_density(df,cfg.PATHS['plot_dir'])
    #plot_error_hist(df, cfg.PATHS['plot_dir'])

    #rmse.plot_rmse(df, cfg.PATHS['plot_dir'])
    rmse.plot_mean_error(df, cfg.PATHS['plot_dir'])
    df = df.reset_index()
    best = pd.DataFrame()
    for id in df.rgi_id.unique():
        #print(df[df.rgi_id == id].error.idxmin())
        best.loc[id,'rmse'] = df[df.rgi_id==id].rmse.min()
        best.loc[id, 'error'] = df.iloc[df[df.rgi_id == id].rmse.idxmin()].error
        best.loc[id,'temp_bias'] = df.iloc[df[df.rgi_id==id].rmse.idxmin()].temp_bias

        best.loc[id,'rgi_area'] = round(wgms_df.loc[id].iloc[0].rgi_area,3)

    best.to_csv(os.path.join(cfg.PATHS['working_dir'], 'best_wgms_error_by_rmse.csv'))

    plt.figure(figsize=(10,10))
    bins = np.arange(14) * 0.400 - 2.600
    best.error.hist(bins=bins)
    plt.axvline(x=best.error.median(),color='k', label='median')
    plt.axvline(x=best.error.mean(), color='k', linestyle=':', label='mean')
    plt.axvline(x=best.error.quantile(0.05), color='grey', label='5% and 95% percentiles')
    plt.axvline(x=best.error.quantile(0.95), color='grey')
    plt.xlabel('Model bias (m)')
    plt.ylabel('Frequency')
    plt.text(-2.8,135,'N='+str(len(best)))
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(cfg.PATHS['plot_dir'],'best_hist.png'))
    #plt.show()

    '''


