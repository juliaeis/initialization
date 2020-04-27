import sys
sys.path.append('../')
from initialization.core import *
from paper.plots_paper import *
from oggm.core.massbalance import MultipleFlowlineMassBalance, PastMassBalance
import seaborn as sns
import geopandas as gpd
from oggm import cfg, utils
pd.options.mode.chained_assignment = None
from copy import deepcopy


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

    wgms = utils.get_ref_mb_glaciers_candidates()

    # Keep only the wgms reference glaciers
    rgidf = rgidf.loc[rgidf.RGIId.isin(wgms)]

    # initialize glaciers
    #gdirs = workflow.init_glacier_regions(rgidf)

    t_0 = 1917
    epsilon = 125
    exp_df = pd.DataFrame()

    '''

    for gdir in gdirs:
        df = pd.DataFrame()
        try:
            # copy previous files to gdir.dir
            dir = os.path.join(OUT_DIR,'per_glacier',gdir.dir.split('per_glacier/')[-1])
            os.system('cp -rf '+dir+'/* '+ gdir.dir)

            t_e = gdir.rgi_date
            ex = [f for f in os.listdir(gdir.dir) if f.startswith('model_run_ad')]
            if len(ex)==1 :
                # read experiment
                dst = os.path.join(gdir.dir,ex[0])
                ex_mod = FileModel(dst)

                # get mb bias and temp_bias
                bias = float(ex[0].split('_')[-1].split('.nc')[0])
                temp_bias = cfg.PATHS['working_dir'].split('_')[-1]

                exp_df.loc[gdir.rgi_id, 'region'] = REGION
                exp_df.loc[gdir.rgi_id, 'rgi_date'] = gdir.rgi_date
                exp_df.loc[gdir.rgi_id, 'mb_bias'] = bias
                exp_df.loc[gdir.rgi_id, 'temp_bias'] = temp_bias


        except:
            exp_df.loc[gdir.rgi_id, 'region'] = REGION
    exp_df.to_csv(os.path.join(cfg.PATHS['working_dir'], REGION + '_experiment.pkl'))

    ex_df = pd.DataFrame()
    model_df = pd.DataFrame()
    for dir in os.listdir(cfg.PATHS['working_dir']):
        if dir.startswith('temp'):
            p = os.path.join(cfg.PATHS['working_dir'],dir)
            df = pd.DataFrame()
            for file in os.listdir(p):
                if file.endswith('.pkl'):
                    rp = os.path.join(p,file)
                    df = df.append(pd.read_pickle(rp), ignore_index=False)
            ex_df.loc[:,dir] = df.mb_bias
            model_df.loc[:,dir] = df.ex_mod

    ex_df.loc[:, 'region'] = df.region
    ex_df = ex_df.sort_index()
    ex_df = ex_df[['region','temp_0', 'temp_-0.25','temp_-0.5','temp_-0.75','temp_-1']]

    model_df.loc[:, 'region'] = df.region
    model_df = model_df.sort_index()
    model_df.index.name='rgi_id'
    model_df = model_df[['region', 'temp_0', 'temp_-0.25', 'temp_-0.5', 'temp_-0.75','temp_-1']]

    model_df.to_pickle('/home/juliaeis/Dokumente/OGGM/work_dir/reconstruction/wgms/ex_mod.pkl')
    #ex_df.to_csv('/home/juliaeis/Dokumente/OGGM/work_dir/reconstruction/wgms/mb_bias.csv')


    wgms_df = pd.read_csv('/home/juliaeis/Dokumente/OGGM/work_dir/reconstruction/wgms/mb_bias.csv')
    lec_df = pd.read_csv('/home/juliaeis/Dokumente/OGGM/work_dir/reconstruction/leclercq_diff2/quick_experiment_df.csv')
    model_df = pd.read_pickle('/home/juliaeis/Dokumente/OGGM/work_dir/reconstruction/wgms/ex_mod.pkl')

    ##### calculate area diff for all wgms glaciers
    area_diff = pd.DataFrame()
    for region in model_df.region.unique():

        ids = model_df[model_df.region == region].drop('region', axis=1)

        # RGI file
        path = utils.get_rgi_region_file(str(region).zfill(2), version='61')
        rgidf = gpd.read_file(path).set_index('RGIId')
        rgidf.index.name='rgi_id'
        rgi_area = rgidf[rgidf.index.isin(ids.index)].Area


        area_df = ids.applymap(lambda x: x.area_km2_ts()[x.last_yr])
        area_diff = area_diff.append(area_df.sub(rgi_area, axis=0))


    ###### combine data about optimal mb bias from Leclercq and WGMS run ######
    wgms_df = wgms_df.sort_values(by='rgi_id')
    lec_df = lec_df.sort_values(by='rgi_id')

    wgms_df = wgms_df.set_index('rgi_id')
    lec_df = lec_df.set_index('rgi_id')

    # remove double entries
    lec_df = lec_df[~lec_df.index.isin(wgms_df.index)]

    df2 = pd.DataFrame()
    df3 = pd.DataFrame()
    for temp_bias in lec_df.temp_bias.unique():
        df = lec_df[lec_df.temp_bias==temp_bias]
        if temp_bias==0.0:
            temp_bias=0
        elif temp_bias==-1.0:
            temp_bias=-1
        col = 'temp_'+str(temp_bias)
        df2.loc[:, col] = df.bias
        df3.loc[:, col] = df.area_diff
    df2.loc[:,'RGIId'] = df2.index
    df2.loc[:,'region'] = df2.RGIId.apply(lambda x: x.split('RGI60-')[-1].split('.')[0])

    wgms_df.loc[:,'region'] = wgms_df.region.apply(lambda x: str(x).zfill(2))

    mb_df = wgms_df.append(df2)
    print(mb_df)
    area_diff = area_diff.append(df3)

    # sort by rgi_id
    mb_df.loc[:, 'RGIId'] = mb_df.index
    mb_df = mb_df.sort_values(by='RGIId')
    mb_df = mb_df.drop('RGIId',axis=1)

    mb_df = mb_df.replace([np.inf, -np.inf], np.nan)
    mb_df = mb_df.dropna()


    area_diff = area_diff.replace([np.inf, -np.inf], np.nan)
    area_diff = area_diff.dropna()

    mb_df = mb_df[['region', 'temp_-1', 'temp_-0.75', 'temp_-0.5', 'temp_-0.25', 'temp_0']]

    mb_df.to_csv('/home/juliaeis/Dokumente/OGGM/work_dir/reconstruction/wgms/all_mb_bias.csv')
    area_diff.to_csv('/home/juliaeis/Dokumente/OGGM/work_dir/reconstruction/wgms/all_area_diff.csv')

    '''
    mb_df = pd.read_csv('/home/juliaeis/Dokumente/OGGM/work_dir/reconstruction/wgms/all_mb_bias.csv').set_index('rgi_id')
    area_diff = pd.read_csv('/home/juliaeis/Dokumente/OGGM/work_dir/reconstruction/wgms/all_area_diff.csv').set_index('rgi_id')
    print(mb_df[mb_df==1000.0].loc[:,'temp_-0.75'].dropna())

    mb_df = mb_df[mb_df!=500.0]
    mb_df = mb_df[mb_df!=1000.0]
    mb_df = mb_df[mb_df!=250.0]
    mb_df = mb_df[mb_df!=-500.0]
    mb_df = mb_df[mb_df!=-1000.0]
    mb_df = mb_df[mb_df!=-250.0]
    mb_df = mb_df[mb_df!=0.0]
    area_diff = area_diff[~pd.isnull(mb_df)]

    mb_df = mb_df.drop('RGI60-05.00446')
    area_diff = area_diff.drop('RGI60-05.00446')


    plt.figure(figsize=(15, 15))
    grid = plt.GridSpec(2, 1, hspace=0.3, wspace=0.2)
    ax2 = plt.subplot(grid[0])
    ax1 = plt.subplot(grid[1])

    cmap = plt.get_cmap("tab10")
    diff = []
    for j,col in enumerate(mb_df.drop('region',axis=1).columns):
        sub = mb_df.loc[:, col].dropna()
        diff = np.append(diff,area_diff[area_diff.index.isin(sub.index)].loc[:,col].values)
        label = col.split('_')[-1]
        sns.distplot(sub, hist=False, kde=True,
                     kde_kws={'shade': True, 'linewidth': 3, 'alpha': 0.15},
                     bins=30,ax=ax2)
        ax2.plot(ax2.lines[-1].get_xydata()[:, 0],
                 ax2.lines[-1].get_xydata()[:, 1], color=cmap(j),
                 label=label)
        ax2.axvline(x=mb_df.loc[:,col].mean(), color=cmap(j), linestyle=':')
    ax1.text(-58,1.1e3,'n:'+str(len(diff)))
    ax1.hist(diff, bins=75)
    ax1.set_yscale('log')
    #ax1.set_xscale('symlog',linthreshx=1e-2)
    ax1.set_ylabel('Frequency', labelpad=30)
    ax1.set_xlabel(r'Area difference (km$^2$)')
    #ax1.set_xlim(-1.05, 1.05)
    ax2.set_ylim(0, None)
    ax2.set_yticks([0, 1e-3, 2e-3])
    ax2.set_xticks(np.arange(-1500, 2000, 500))
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
    plt.savefig(os.path.join(cfg.PATHS['plot_dir'], 'mb_bias.png'), dpi=300)
    plt.show()

    print(diff)
    print(len(diff))
    print(diff.mean())
    print(np.median(diff))
    print(np.std(diff))
    print(np.quantile(diff,0.95))
    print(np.quantile(diff,0.05))