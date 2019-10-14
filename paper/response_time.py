import os
import sys
import pandas as pd
import geopandas as gpd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from oggm import cfg, workflow, utils


sys.path.append('../')
from initialization.core import *
import matplotlib as mpl
from sklearn.linear_model import LinearRegression
from mpl_toolkits.axes_grid.anchored_artists import AnchoredText

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



def _response(gdir, model_df):

    fig,ax = plt.subplots(1,1)

    ex_mod2 = model_df.loc[gdir.rgi_id].experiment
    ax.plot(0,copy.deepcopy(ex_mod2.volume_km3),'o', label=r'experiment$_{1850}$')
    tasks.run_constant_climate(gdir, nyears=600, y0=1850, halfsize=15,
                               temperature_bias=-1,
                               store_monthly_step=False,
                               output_filesuffix='experiment_equilibrium',
                               init_model_fls=copy.deepcopy(ex_mod2.fls))

    rp = gdir.get_filepath('model_run', filesuffix='experiment_equilibrium')
    ex_mod = FileModel(rp)
    ex_mod.run_until(600)


    try:
        # scenario a
        tasks.run_constant_climate(gdir, nyears=600, y0=1865, halfsize=15,
                                   temperature_bias=-1.05,
                                   store_monthly_step=False,
                                   output_filesuffix='_-1.05',
                                   init_model_fls=copy.deepcopy(ex_mod.fls))
    except:
        pass

    try:

        # senario b
        tasks.run_constant_climate(gdir, nyears=600, y0=1865, halfsize=15,
                                   store_monthly_step=False, temperature_bias=-1.1,
                                   output_filesuffix='_-1.1',
                                   init_model_fls=copy.deepcopy(ex_mod.fls))

    except:
        pass

    try:

        # senario c
        tasks.run_constant_climate(gdir, nyears=600, y0=1865, halfsize=15,
                                   temperature_bias=-1,
                                   store_monthly_step=False,
                                   output_filesuffix='_-1',
                                   init_model_fls=copy.deepcopy(ex_mod.fls))
    except:
        pass

    try:
        # senario d
        tasks.run_constant_climate(gdir, nyears=600, y0=1865, halfsize=15,
                                   store_monthly_step=False, temperature_bias=-0.9,
                                   output_filesuffix='_-0.9',
                                   init_model_fls=copy.deepcopy(ex_mod.fls))
    except:
        pass

    try:
        # senario_e
        tasks.run_constant_climate(gdir, nyears=600, y0=1865, halfsize=15,temperature_bias=-0.95,
                                   store_monthly_step=False,
                                   output_filesuffix='_-0.95',
                                   init_model_fls=copy.deepcopy(ex_mod.fls))
    except:
        pass

    response = pd.DataFrame()

    colors = ['C0','C1','C2','C3','C4','C5',]
    for i, s in enumerate(['_-1.1', '_-1.05', '_-1', '_-0.95', '_-0.9']):

        #try:
        rp = gdir.get_filepath('model_run', filesuffix=s)
        mod2 = FileModel(rp)
        if mod2.volume_km3_ts()[600] > 0:

            diff = mod2.volume_km3_ts()[600] - ((
                                                mod2.volume_km3_ts()[600] -
                                                ex_mod.volume_km3_ts()[
                                                    600]) / np.exp(1))

            t = abs(mod2.volume_km3_ts() - diff).idxmin()
            response.at[gdir.rgi_id, s.split('_')[-1]] = t
            mod2.volume_km3_ts().plot(ax=ax, color=colors[i], label=s.split('_')[-1])
            ax.axvline(t,color=colors[i],linestyle=':')
        #except:
        #    pass
    plt.legend(loc='best')
    p = os.path.join(cfg.PATHS['plot_dir'],'response_time')
    utils.mkdir(p)
    #plt.savefig(os.path.join(p,str(gdir.rgi_id)+'.png'))
    plt.ylim((15,30))
    plt.show()
    return response

def value_to_color(val):
    n_colors = 15  # Use 256 colors for the diverging color palette
    palette = sns.color_palette("RdBu_r", 15)# sns.diverging_palette(20, 220, n=n_colors)  # Create the palette
    color_min, color_max = [-1,1]  # Range of values that will be mapped to the palette, i.e. min and max possible correlation
    val_position = float((val - color_min)) / (color_max - color_min) # position of value in the input range, relative to the length of the input range
    ind = int(val_position * (n_colors - 1)) # target index in the color palette
    return palette[ind]

def heatmap(x, y, size, color):
    plt.figure(figsize=(10,4))
    plot_grid = plt.GridSpec(5, 15, hspace=0.2, wspace=0.1)  # Setup a 1x15 grid
    ax = plt.subplot(plot_grid[:3, 2:])
    #fig, ax = plt.subplots(figsize=(5,1))

    # Mapping from column names to integer coordinates
    x_labels = x#[v for v in sorted(x.unique())]
    y_labels = y[::-1]#[v for v in sorted(y.unique())]
    x_to_num = {p[1]: p[0] for p in enumerate(x_labels)}
    y_to_num = {p[1]: p[0] for p in enumerate(y_labels)}
    import matplotlib

    size_scale = 500

    cmap = matplotlib.cm.get_cmap('RdBu_r',8)
    normalize = matplotlib.colors.Normalize(vmin=-1, vmax=1)
    colors = [cmap(normalize(value)) for value in color]

    sc = ax.scatter(
        x=x.map(x_to_num),  # Use mapping for x
        y=y.map(y_to_num),  # Use mapping for y
        s=size * size_scale,
        color=colors,
        # Vector of square sizes, proportional to size parameter
        marker='s'  # Use square as scatterplot marker
    )
    print(x.map(x_to_num),y.map(y_to_num))
    for i,col in enumerate(color):
        plt.text(x.map(x_to_num)[i]-.3,y.map(y_to_num)[i]-3.5 , round(col,2))
    # Show column labels on the axes
    ax.set_xticks([x_to_num[v] for v in x_labels])
    x_labels = ['Reconstructability', 'Response Time ', r'Length$_{2000}$', r'Area$_{2000}$', r'Volume$_{2000}$', r'ELA$_{2000}$', r'ELA$_{1850-2000}$',
              r'Slope$_{mean}^{2000}$',r'Slope$_{1/3}^{2000}$']
    ax.set_xticklabels(x_labels, rotation=45, horizontalalignment='right')
    ax.set_yticks([17,8])
    ax.set_yticklabels(['Reconstructability','Response Time'])
    ax.set_ylim(0,21)

    label = np.arange(-1,1.25,0.25)
    loc= copy.deepcopy(label)

    cax, _ = matplotlib.colorbar.make_axes(ax)
    cbar = matplotlib.colorbar.ColorbarBase(cax, cmap=cmap, norm=normalize,ticks=np.arange(-1,2))
    cbar.set_ticks(loc)
    cbar.set_ticklabels(label)
    plt.tight_layout()
    ax.set_title('Correlation (Alps), n='+str(len(df.dropna(subset=['0']))))
    plt.savefig(os.path.join(cfg.PATHS['working_dir'],'matrix.pdf'), dpi=300)

    plt.show()


def response_time(gdirs, model_df, job_nr):

    response_df = pd.DataFrame()

    pool = Pool()
    list = pool.map(partial(_response, model_df=model_df), gdirs)
    pool.close()
    pool.join()

    response_df = response_df.append(list, ignore_index=False)
    p = os.path.join(cfg.PATHS['working_dir'], 'response_time_'+str(job_nr)+'.pkl')
    response_df.to_pickle(p)
    return p

def slope(model):
    slope1850 = []
    slope2000 = []
    slope1850_third = []
    slope2000_third = []
    model.run_until(1850)
    for i in range(len(model.fls)):
        thick = model.fls[i].thick
        surface = model.fls[i].surface_h[thick>0]
        dx = model.fls[i].dx_meter
        if len(surface) > 1:
            gradient = -np.gradient(surface, dx)
            n = int(len(gradient)/3)
            slope1850.extend(gradient)
            slope1850_third.extend(gradient[-n:])

    model.run_until(2000)
    for i in range(len(model.fls)):
        thick = model.fls[i].thick
        surface = model.fls[i].surface_h[thick>0]
        dx = model.fls[i].dx_meter
        if len(surface) > 1:
            gradient = -np.gradient(surface, dx)
            n = int(len(gradient)/3)
            slope2000.extend(gradient)
            slope2000_third.extend(gradient[-n:])

    return np.mean(slope1850), np.max(slope1850),np.mean(slope1850_third), np.mean(slope2000), np.max(slope2000), np.mean(slope2000_third)
if __name__ == '__main__':
    cfg.initialize()

    ON_CLUSTER = False

    # Local paths
    if ON_CLUSTER:
        WORKING_DIR = os.environ.get("S_WORKDIR")
        cfg.PATHS['working_dir'] = WORKING_DIR
        job_nr = int(os.environ.get('I'))
    else:
        WORKING_DIR = '/home/juliaeis/Dokumente/OGGM/work_dir/reconstruction/response/'
        cfg.PATHS['working_dir'] = WORKING_DIR
        utils.mkdir(WORKING_DIR, reset=False)
        job_nr=0

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
    '''
    path = utils.get_rgi_region_file('11', version='61')
    rgidf = gpd.read_file(path)
    if ON_CLUSTER:
         p = os.path.join('/home/users/julia/initialization/out/paper_correction/paper_600', 'models_' + str(job_nr) + '.pkl')
         model_df = pd.read_pickle(p,compression='gzip')
    else:
        p = os.path.join(cfg.PATHS['working_dir'], 'models_0.pkl')
        model_df = pd.read_pickle(p)


    #rgidf = rgidf[rgidf.RGIId.isin(model_df.index)].sort_values(by='Area', ascending=False)
    #gdirs = workflow.init_glacier_regions(rgidf.head(1))

    #preprocessing(gdirs)
    #p = response_time(gdirs, model_df, job_nr)
    #print(pd.read_pickle(p))

    df = pd.DataFrame()
    for file in os.listdir(cfg.PATHS['working_dir']):

        if file.startswith('response_time') and not file.endswith('merge.pkl'):
            p = os.path.join(cfg.PATHS['working_dir'], file)
            df = df.append(pd.read_pickle(p))



    df.to_pickle(os.path.join(cfg.PATHS['working_dir'],'response_time_merge.pkl'))
    '''
    # read DataFrame with response times
    response_df = pd.read_pickle(os.path.join(cfg.PATHS['working_dir'],'response_time_merge.pkl'))
    # read stuff for correlation
    att_df = pd.read_pickle(os.path.join(cfg.PATHS['working_dir'], 'ratio.pkl'), compression='gzip')#.dropna().set_index('rgi')

    #att_df = att_df[att_df.index.isin(response_df.index)]
    #exp_df = pd.read_pickle(os.path.join(cfg.PATHS['working_dir'], 'experiment.pkl'))

    #exp_df = exp_df[att_df.index]
    #att_df['slope_mean1850'], att_df['slope_max1850'], att_df['slope_third1850'],att_df['slope_mean2000'], att_df['slope_max2000'], att_df['slope_third2000']=zip(*exp_df.map(slope))
    #att_df.to_pickle(os.path.join(cfg.PATHS['working_dir'],'ratio.pkl'), compression='gzip')
    #print(att_df[['slope_mean','slope_mean1850','slope_max1850','slope_third1850','slope_mean2000','slope_max2000', 'slope_third2000']])

    df =  response_df.join(att_df)

    #df1 = df[['ratio', '0', 'length', 'area', 'volume', 'ela_2000',  'slope_mean',
    #         'slope_third']]
    df1 = df[['ratio', '0', 'length', 'area', 'volume', 'ela_2000','ela_change', 'slope_mean2000', 'slope_third2000']]
    '''
    corr = df1.corr()
    corr = pd.melt(corr[['ratio','0']].reset_index(), id_vars='index') # Unpivot the dataframe, so we can get pair of arrays for x and y
    corr.columns = ['x', 'y', 'value']
    heatmap(
        x=corr['x'],
        y=corr['y'],
        size=corr['value'].abs(),
        color= corr['value']
    )



    df2 = df[['ratio', '-1.1', '-0.9', '-0.5', '0' ]]
    fig, ax = plt.subplots(1,1)
    print(df.isna().sum())
    print(df2.corr())

    matrix = ax.matshow(df2.corr(), vmin=-1, vmax=1, cmap='RdBu_r')
    marks = ['Reconstructability', '-1.1', '-0.9', '-0.5', '0']

    tick_marks = [i for i in range(len(df2.columns))]
    #marks = df2.columns
    plt.xticks(tick_marks, marks, rotation='vertical')
    plt.yticks(tick_marks, marks)
    plt.colorbar(matrix, fraction=0.046, pad=0.04,
                 label='correlation coefficient')
    plt.tight_layout()
    #plt.savefig(
    #    '/home/juliaeis/Dropbox/Apps/Overleaf/reconstruction_paper/plots/measure.pdf')
    plt.show()

    plt.figure(figsize=(15,10))
    df.ratio.plot.hist(bins=20)
    plt.xlabel('Reconstructability measure')
    plt.ylabel('Frequency')
    plt.title('Reconstructability (Alps), n=2660')
    plt.grid()
    plt.savefig('/home/juliaeis/Dropbox/Apps/Overleaf/reconstruction_paper/plots/reconstrucatibility_hist.pdf')
    plt.show()



    for col in df1.columns.drop('ratio'):

        X = df.dropna(subset=[col]).ratio.values.reshape(-1, 1)
        Y = df.loc[:, col].dropna().values.reshape(-1, 1)
        linear_regressor = LinearRegression()  # create object for the class
        linear_regressor.fit(X, Y)  # perform linear regression
        Y_pred = linear_regressor.predict(X)

        df1.plot.scatter(x='ratio', y=col,alpha=0.5)
        plt.plot(X, Y_pred, color='C0', linewidth=3)
        plt.xlabel('Reconstructability')
        if col == '0':
            col = 'response_time'
        plt.ylabel(col)
        #print(os.path.join(cfg.PATHS['plot_dir'],'test.png'))
        #plt.savefig(os.path.join(cfg.PATHS['plot_dir'],str(col)+'.png'), dpi=200)
        plt.show()

    plt.figure(figsize=(20,10))
    grid = plt.GridSpec(2, 3, hspace=0.3, wspace=0.4)
    ax0 = plt.subplot(grid[:,0])
    ax1 = plt.subplot(grid[0,1])
    ax2 = plt.subplot(grid[0, 2])
    ax3 = plt.subplot(grid[1,1])
    ax4 = plt.subplot(grid[1,2])

    df.ratio.plot.hist(bins=20, ax=ax0)

    # reconstr. vs. response time
    X = df.dropna(subset=['0']).ratio.values.reshape(-1, 1)
    Y = df.loc[:, '0'].dropna().values.reshape(-1, 1)
    linear_regressor = LinearRegression()  # create object for the class
    linear_regressor.fit(X, Y)  # perform linear regression
    Y_pred = linear_regressor.predict(X)
    ax1.plot( X, Y_pred, color='C0', linewidth=3)
    df1.plot.scatter(x='ratio',y='0', alpha=0.3, ax=ax1, color='C0')

    # reconst. vs ELA
    X = df.dropna(subset=['ela_2000']).ratio.values.reshape(-1, 1)
    Y = df.loc[:, 'ela_2000'].dropna().values.reshape(-1, 1)
    linear_regressor = LinearRegression()  # create object for the class
    linear_regressor.fit(X, Y)  # perform linear regression
    Y_pred = linear_regressor.predict(X)
    df1.plot.scatter(x='ratio', y='ela_2000', alpha=0.3, ax=ax2, color='C2')
    ax2.plot(X,Y_pred, color='C2', linewidth=3)

    # reconstr. vs. slope_mean
    X = df.dropna(subset=['slope_mean2000']).ratio.values.reshape(-1, 1)
    Y = df.loc[:, 'slope_mean2000'].dropna().values.reshape(-1, 1)
    linear_regressor = LinearRegression()  # create object for the class
    linear_regressor.fit(X, Y)  # perform linear regression
    Y_pred = linear_regressor.predict(X)
    df1.plot.scatter(x='ratio', y='slope_mean2000', alpha=0.3, ax=ax3, color='C1')
    ax3.plot(X,Y_pred,linewidth=3, color='C1')

    # reconstr. vs slope_third
    X = df.dropna(subset=['slope_third2000']).ratio.values.reshape(-1, 1)
    Y = df.loc[:, 'slope_third2000'].dropna().values.reshape(-1, 1)
    linear_regressor = LinearRegression()  # create object for the class
    linear_regressor.fit(X, Y)  # perform linear regression
    Y_pred = linear_regressor.predict(X)
    df1.plot.scatter(x='ratio', y='slope_third2000', alpha=0.3, ax=ax4, color='C3')
    ax4.plot(X,Y_pred,linewidth=3, color='C3')

    ax1.set_yticks([0,100,200,300,400])
    ax2.set_yticks([2000, 2500, 3000, 3500, 4000])
    ax3.set_yticks([0,0.5,1,1.5,2])
    ax4.set_yticks([0,0.5,1,1.5,2, 2.5])

    ax0.grid()
    ax1.grid()
    ax2.grid()
    ax3.grid()
    ax4.grid()

    ax0.set_xlabel('Reconstructability')
    ax1.set_xlabel('Reconstructability')
    ax2.set_xlabel('Reconstructability')
    ax3.set_xlabel('Reconstructability')
    ax4.set_xlabel('Reconstructability')

    ax1.set_ylabel('Response Time (years)')
    ax2.set_ylabel(r'ELA$_{2000}$ (m)')
    ax3.set_ylabel(r'Slope$_{mean}^{2000}$')
    ax4.set_ylabel(r'Slope$_{1/3}^{2000}$')

    add_at(ax0, r"a", loc=2)
    add_at(ax1, r"b", loc=2)
    add_at(ax2, r"c", loc=2)
    add_at(ax3, r"d", loc=2)
    add_at(ax4, r"e", loc=2)

    plt.savefig('/home/juliaeis/Dropbox/Apps/Overleaf/reconstruction_paper/plots/hist.pdf', dpi=150)
    '''
    df1 = df1.loc[['RGI60-11.01450','RGI60-11.00106','RGI60-11.00887','RGI60-11.00719','RGI60-11.00787','RGI60-11.00068','RGI60-11.00918','RGI60-11.00049', 'RGI60-11.00080','RGI60-11.00003','RGI60-11.00251','RGI60-11.00289','RGI60-11.00300','RGI60-11.00402','RGI60-11.00026', 'RGI60-11.02791', 'RGI60-11.02715',
                                    'RGI60-11.01246'],:]
    print(df1.sort_values(by='area', ascending=False).ratio)
    #print(df1[df1.ratio<0.1].sort_values(by=['ratio']).area)