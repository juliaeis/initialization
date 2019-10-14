
import os
import sys
from copy import deepcopy
sys.path.append('../')
from initialization.core import *
from plots_paper import *

import matplotlib.pyplot as plt
import salem
import geopandas as gpd
from oggm import cfg, workflow, utils, tasks
from oggm.core.flowline import FluxBasedModel
pd.options.mode.chained_assignment = None
import time
import matplotlib as mpl
import xarray as xr

mpl.rcParams['axes.linewidth'] = 3
mpl.rcParams['xtick.major.width'] = 3
mpl.rcParams['ytick.major.width'] = 3
mpl.rcParams['font.size'] =20
mpl.rcParams['font.weight'] = 'medium'
mpl.rcParams['axes.labelweight'] = 'medium'
mpl.rcParams['legend.fontsize'] = 20 #30
mpl.rcParams['lines.linewidth'] = 3

from scipy.stats import pearsonr

def calculate_pvalues(df):
    df = df._get_numeric_data()
    dfcols = pd.DataFrame(columns=df.columns)
    pvalues = dfcols.transpose().join(dfcols, how='outer')
    for r in df.columns:
        for c in df.columns:
            if c == r:
                df_corr = df[[r]].dropna()
            else:
                df_corr = df[[r,c]].dropna()
            pvalues[r][c] = pearsonr(df_corr[r], df_corr[c])[1]
    return pvalues

def corr_pvalue(df):
    numeric_df = df.dropna()._get_numeric_data()
    cols = numeric_df.columns
    mat = numeric_df.values

    arr = np.zeros((len(cols),len(cols)), dtype=object)

    for xi, x in enumerate(mat.T):
        for yi, y in enumerate(mat.T[xi:]):
            arr[xi, yi+xi] = map(lambda _: round(_,4), pearsonr(x,y))
            arr[yi+xi, xi] = arr[xi, yi+xi]

    return pd.DataFrame(arr, index=cols, columns=cols)

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


def read_results(gdirs):

    model_df = pd.DataFrame()
    pool = mp.Pool()
    list = pool.map(read_result_parallel, gdirs)
    pool.close()
    pool.join()
    model_df = model_df.append(list, ignore_index=True)
    return model_df


def plot_ratio_volume(df,acc_df,perc_df,ex_mod,gdir, ratio, t):

    fig, ax = plt.subplots(1,1, figsize=(10,7))
    norm = mpl.colors.LogNorm(vmin=0.01 / 125, vmax=10)
    cmap = mpl.cm.get_cmap('viridis')

    # df = df.sort_values('objective', ascending=False)
    df = df.sort_values('fitness', ascending=False)
    acc_min = acc_df.loc[acc_df.volume.idxmin(),'model'].volume_km3_ts()
    acc_max = acc_df.loc[acc_df.volume.idxmax(),'model'].volume_km3_ts()
    perc_min = perc_df.loc[perc_df.volume.idxmin(), 'model'].volume_km3_ts()
    perc_max = perc_df.loc[perc_df.volume.idxmax(), 'model'].volume_km3_ts()

    ax.fill_between(acc_max.index, acc_max.values,acc_min.values, color='grey', alpha=0.3)
    ax.fill_between(perc_max.index, perc_max.values, perc_min.values, color='C0', alpha=0.5)
    ex_mod.volume_km3_ts().plot(ax=ax, color='red', linestyle=':',
                                linewidth=3,
                                label='')
    ax.axvline(x=t, color='k')
    plt.title(gdir.rgi_id)
    plt.text(1950, 0.95 * acc_max.max(), 'measure: ' + str(ratio.round(4)))
    plt.text(1950, 0.9 * acc_max.max(), 't: ' + str(t))

    plt.ylabel('Volume (km$^3$)')
    plt.ylabel('Volume (km$^3$)')
    plt.xlabel('Time (years)')
    utils.mkdir(os.path.join(cfg.PATHS['plot_dir'],'ratio'))
    plt.savefig(os.path.join(cfg.PATHS['plot_dir'],'ratio',gdir.rgi_id+'.png'))
    plt.show()

def read_result_parallel(gdir):

    #try:
    rp = gdir.get_filepath('model_run', filesuffix='_experiment')
    ex_mod = FileModel(rp)

    df = pd.read_pickle(os.path.join(gdir.dir, 'result1850.pkl'),
                        compression='gzip')

    df.fitness = df.fitness / 125
    # replace all values smaller than 1e-4
    df.fitness[df.fitness < 1e-4] = 1e-4
    acc_df = df[df.fitness<1]
    perc_df = acc_df[acc_df.fitness <= round(acc_df.fitness.quantile(0.05),4)]

    v_acc = acc_df.model.apply(lambda x: x.volume_km3_ts())
    v_perc = perc_df.model.apply(lambda x: x.volume_km3_ts())

    try:
        ratio = 1-(v_perc.max()-v_perc.min())/(v_acc.max()[1850]-v_acc.min()[1850])
        t = int(ratio[ratio>0.98].index[0])
    except:
        t=2000

    plot_ratio_volume(df,acc_df, perc_df,ex_mod,gdir, ratio.loc[1850], t)

    inv_in = gdir.read_pickle('inversion_input')
    slope = []
    slope2 = []
    slope3 = []
    for i in range(len(inv_in)):
        n = int(len(inv_in[-1]['slope_angle'])/3)

        slope.extend(inv_in[i]['slope_angle'])
        slope2.extend(inv_in[i]['slope_angle'][-n:])
        slope3.extend(inv_in[i]['slope_angle'][-2*n:])


    ds = xr.open_dataset(gdir.get_filepath('model_diagnostics', filesuffix='_experiment'))
    diag = ds.to_dataframe()

    return pd.Series({'rgi':gdir.rgi_id, 'ratio':ratio, 'time':t,
                      'length':ex_mod.length_m_ts()[2000],
                      'area':ex_mod.area_km2_ts()[2000],
                      'volume':ex_mod.volume_km3_ts()[2000],
                      'slope_max':np.max(slope),
                      'slope_mean':np.mean(slope),
                      'slope_median':np.median(slope),
                      'slope_third': np.mean(slope2),
                      'slope_2thrid': np.mean(slope3),
                      'ela_2000':diag.ela_m[2000],
                      'ela_change':diag.ela_m[1850]/diag.ela_m[2000]})

    #except:
    #    return pd.Series({'rgi':gdir.rgi_id})

if __name__ == '__main__':
    cfg.initialize()

    ON_CLUSTER = False

    # Local paths
    if ON_CLUSTER:
        WORKING_DIR = os.environ.get("S_WORKDIR")
        WORKING_DIR = '/home/users/julia/initialization/out/paper_correction/paper_600'
        cfg.PATHS['working_dir'] = WORKING_DIR
    else:
        WORKING_DIR = '/home/juliaeis/Dokumente/OGGM/work_dir/reconstruction/600_paper_correction/'
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
    rgidf = rgidf[rgidf.RGIId.isin(['RGI60-11.00897','RGI60-11.00779', 'RGI60-11.00029', 'RGI60-11.00036', 'RGI60-11.00001','RGI60-11.00026','RGI60-11.00062'])]
    #rgidf = rgidf[rgidf.RGIId.isin(['RGI60-11.00074'])]
    # sort for efficient using
    rgidf = rgidf.sort_values('Area', ascending=False)
    gdirs = workflow.init_glacier_regions(rgidf)

    #df = read_results(gdirs).dropna()
    #df.to_pickle(os.path.join(cfg.PATHS['working_dir'],'ratio.pkl'),compression='gzip')


    df = pd.read_pickle(os.path.join(cfg.PATHS['working_dir'],'ratio.pkl'),compression='gzip').dropna()
    plt.figure(figsize=(15,8))
    grid = plt.GridSpec(1,2,top=0.75,bottom=0.25,wspace=0.4)
    ax1 = plt.subplot(grid[0])
    ax2 = plt.subplot(grid[1])
    df.ratio.hist(bins=20, ax=ax1)
    ax1.set_xlabel('Reconstructability measure')
    ax1.set_ylabel('Frequency')
    plt.suptitle('Reconstructability (Alps), n=2660')

    df = df[['rgi', 'ratio',  'length', 'area', 'volume', 'ela_2000',
             'ela_change', 'slope_max', 'slope_mean', 'slope_median', 'slope_2thrid',
             'slope_third']].dropna()
    print(df.corr()['ratio'])

        #for i, col in enumerate(df.drop(['rgi','ratio'],axis=1)):
        #fig, ax = plt.subplots(1, 1)
        #df.plot.scatter(x=col,y='ratio',ax=ax)
    #plt.show()


    p = calculate_pvalues(df)
    matrix = ax2.matshow(df.drop('rgi',axis=1).corr(),vmin=-1,vmax=1, cmap='RdBu_r')
    marks = ['Reconstr.', 'Length', 'Area','Volume', r'ELA$_{2000}$', r'$\Delta$ ELA ',
              r'Slope$_{max}$', r'Slope$_{mean}$', r'Slope$_{median}$',r'Slope$_{2/3}$',r'Slope$_{1/3}$']
    tick_marks = [i for i in range(len(df.drop('rgi',axis=1).columns))]
    plt.xticks(tick_marks, marks, rotation='vertical')
    plt.yticks(tick_marks, marks)
    plt.colorbar(matrix,fraction=0.046, pad=0.04, label='correlation coefficient')
    plt.tight_layout()
    plt.savefig('/home/juliaeis/Dropbox/Apps/Overleaf/reconstruction_paper/plots/measure.pdf')
    plt.show()


    '''

    from sklearn.linear_model import LinearRegression

    X = df.ratio.values.reshape(-1, 1)
    Y = df.slope_mean.values.reshape(-1, 1)
    linear_regressor = LinearRegression()  # create object for the class
    linear_regressor.fit(X, Y)  # perform linear regression
    Y_pred = linear_regressor.predict(X)


    fig,ax = plt.subplots(1,1)
    df.plot.scatter(x='ratio', y='slope_mean',color='C1',ax=ax,label='ratio')
    ax.plot(X, Y_pred, color='C1')
    plt.xlabel('reconstructability')
    plt.ylabel('Glacier Area (km²)')

    plt.xlabel('reconstructability measure')
    plt.ylabel('mean slope')

    plt.figure()
    df.time.plot.hist(bins=30)
    plt.show()


    print(gdir.rgi_id)
    print(len(df[df.fitness<1e-4]))
    print(len(df[(df.fitness > 1e-4) & (df.fitness< 1e-3)]))
    print(len(df[(df.fitness > 1e-3) & (df.fitness < 1e-2)]))
    print(len(df[(df.fitness > 1e-2) & (df.fitness < 1e-1)]))
    print(len(df[(df.fitness > 1e-1) & (df.fitness < 1e-0)]))
    print(len(df[(df.fitness > 1e0) & (df.fitness < 1e1)]))
    print(len(df[df.fitness > 1e1]))
    plt.show()

    df = pd.read_pickle(os.path.join(cfg.PATHS['working_dir'], 'ratio.pkl'),
                 compression='gzip')
    ratio = df.dropna().reconstructability
    ratio.plot.hist(bins=30)
    plt.yscale('log')
    plt.show()

    df = read_results(gdirs)
    if ON_CLUSTER:
        df.to_pickle(os.path.join(os.environ.get("S_WORKDIR"),'ratio.pkl'),compression='gzip')
    else:
        df.to_pickle(os.path.join(cfg.PATHS['working_dir'],'ratio.pkl'),compression='gzip')
    print(df)




                        med_mod, perc_min, perc_max = find_median(df)
                        min_mod = deepcopy(df.loc[df.fitness.idxmin(), 'model'])

                        df['fitness2'] = df.model.apply(lambda x: abs(x.area_km2_ts()[2000] - ex_mod.area_km2_ts()[2000]) ** 2)
                        df['fitness3'] = df.model.apply(lambda x: abs(x.length_m_ts()[2000] - ex_mod.length_m_ts()[2000]) ** 2)

                        # saves median state, minimum state and experiment model
                        model_df.loc[gdir.rgi_id, 'median'] = deepcopy(med_mod)
                        model_df.loc[gdir.rgi_id, 'minimum'] = deepcopy(min_mod)
                        model_df.loc[gdir.rgi_id, 'experiment'] = deepcopy(ex_mod)
                        model_df.loc[gdir.rgi_id, 'flowline'] = FluxBasedModel(flowlines=gdir.read_pickle('model_flowlines'))
                        model_df.loc[gdir.rgi_id, 'perc_min'] = deepcopy(perc_min)
                        model_df.loc[gdir.rgi_id, 'perc_max'] = deepcopy(perc_max)
                        model_df.loc[gdir.rgi_id, 'fit2'] = deepcopy(df.loc[df.fitness2.idxmin(), 'model'])
                        model_df.loc[gdir.rgi_id, 'fit3'] = deepcopy(df.loc[df.fitness3.idxmin(), 'model'])


                        # time_df
                        time_df.loc[gdir.rgi_id, 'time'] = time.time()-start

                        try:
                            # plots
                            plot_experiment(gdir, ex_mod, t_0, cfg.PATHS['plot_dir'])
                            plot_candidates(gdir, df, t_0, 'step3', cfg.PATHS['plot_dir'])
                            plot_fitness_values(gdir, df, ex_mod, t_0, cfg.PATHS['plot_dir'])
                            plot_median(gdir, df, epsilon, ex_mod, t_0, t_e, cfg.PATHS['plot_dir'])

                        except:
                            pass

            except:
                print(gdir.rgi_id+' failed')
    '''