
"""

This scipt calculates the errors from the model outputs, derived from the
PREVIOUSLY runs (e.g. run.py for the alps).

"""
import os
import pandas as pd
import geopandas as gpd
import numpy as np
import matplotlib.pyplot as plt
from oggm import cfg, workflow, utils
from plots_paper import *

if __name__ == '__main__':

    cfg.initialize()

    ON_CLUSTER = False

    # Local paths
    if ON_CLUSTER:
        WORKING_DIR = os.environ.get("S_WORKDIR")
        cfg.PATHS['working_dir'] = WORKING_DIR
        job_nr = int(os.environ.get('I'))
    else:
        WORKING_DIR = '/home/juliaeis/Dokumente/OGGM/work_dir/reconstruction/600_paper_correction/'
        cfg.PATHS['working_dir'] = WORKING_DIR
        utils.mkdir(WORKING_DIR, reset=False)

    cfg.PATHS['plot_dir'] = os.path.join(cfg.PATHS['working_dir'], 'plots')
    utils.mkdir(cfg.PATHS['plot_dir'], reset=False)
    '''
    model_df = pd.DataFrame()
    #time_df = pd.DataFrame()

    # merge all model outputs from different cluster jobs
    for file in os.listdir(cfg.PATHS['working_dir']):

        if file.startswith('models') and not file.endswith('merge.pkl'):
            p = os.path.join(cfg.PATHS['working_dir'], file)
            model_df = model_df.append(pd.read_pickle(p, compression='gzip'))
            print(p)
        if file.startswith('time') and not file.endswith('merge.pkl'):
            p = os.path.join(cfg.PATHS['working_dir'], file)
            #time_df = time_df.append(pd.read_pickle(p, compression='gzip'))

    # make sure that no entries are double
    model_df.loc[:, 'double'] = model_df.index.duplicated(keep='first')
    model_df = model_df[model_df['double'] == False].drop('double', axis=1)

    #time_df.loc[:, 'double'] = time_df.index.duplicated(keep='first')
    #time_df = time_df[time_df['double'] == False].drop('double', axis=1)

    model_df.to_pickle(os.path.join(cfg.PATHS['working_dir'], 'models_merge.pkl'))
    #time_df.to_pickle(os.path.join(cfg.PATHS['working_dir'], 'time_merge.pkl'))


    model_df = pd.read_pickle(os.path.join(cfg.PATHS['working_dir'], 'models_merge.pkl'))
    median_df = pd.read_pickle(os.path.join(cfg.PATHS['working_dir'], 'median.pkl'), compression='gzip').dropna()
    median_df = median_df.set_index('rgi').sort_index()
    model_df = model_df.sort_index()
    model_df['median']=median_df.loc[:,'median'].values



    #median_df.loc[:,'bias'] = median_df.loc[:,'median'].apply(lambda x: float(x.path.split('/')[-1].split('_')[-2]))
    #median_df.bias.hist(bins=40,figsize=(10,8))
    plt.ylabel('Frequecy')
    plt.xlabel('Temperature bias')
    plt.xticks(np.arange(-5,3,1))
    plt.savefig(os.path.join(cfg.PATHS['plot_dir'],'hist_median_bias.pdf'))
    plt.tight_layout()
    plt.show()

    v_1850 = model_df['experiment'].apply(lambda x: x.area_km2_ts()[1850])
    v_2000 = model_df['experiment'].apply(lambda x: x.area_km2_ts()[2000])


    print('glacierized area in 1850:', v_1850.sum())
    print('glacierized area in 2000:', v_2000.sum())

    path = utils.get_rgi_region_file('11', version='61')
    rgidf = gpd.read_file(path).set_index('RGIId')
    rgi_2000 = rgidf.loc[model_df.index].Area
    print('RGI area: ',rgi_2000.sum())

    (rgi_2000-v_2000).plot.hist(bins=50)
    plt.show()


    median = model_df['median'].apply(lambda x: x.volume_km3_ts())
    minimum = model_df['minimum'].apply(lambda x: x.volume_km3_ts())
    experiment = model_df['experiment'].apply(lambda x: x.volume_km3_ts())
    fit2 = model_df['fit2'].apply(lambda x: x.volume_km3_ts())
    fit3 = model_df['fit3'].apply(lambda x: x.volume_km3_ts())


    # absolute errors
    error1 = median - experiment
    error2 = minimum - experiment
    '''
    error1= pd.read_pickle(os.path.join(cfg.PATHS['working_dir'],'error1.pkl'))
    error2 = pd.read_pickle(os.path.join(cfg.PATHS['working_dir'], 'error2.pkl'))
    plot_absolute_error(error1, error2, cfg.PATHS['plot_dir'], all=True)
    '''
    # relative errors
    error3 = ((median - experiment) / experiment)*100
    error4 = ((minimum - experiment) / experiment)*100
    error3 = error3.replace([np.inf, -np.inf], np.nan).dropna()
    error4 = error4.replace([np.inf, -np.inf], np.nan).dropna()

    plot_relative_error_min_vs_med(error3, error4, cfg.PATHS['plot_dir'])


    # errors for different fitness functions
    error5 = ((fit2 - experiment) / experiment)*100
    error6 = ((fit3 - experiment) / experiment)*100

    plot_compare_fitness_functions(error4, error5, error6, cfg.PATHS['plot_dir'])


    #read measured time for each glacier in seconds
    time_df = pd.read_pickle(os.path.join(cfg.PATHS['working_dir'], 'time_merge.pkl'))
    print('Computational Times for the reconstructions in the Alps (n=)'+ str(len(time_df)))
    print('minimum time:', time_df.time.idxmin(), '\t', '{:10.2f}'.format(time_df.time.min()), ' seconds')
    print('maximum time:', time_df.time.idxmax(), '\t', '{:10.2f}'.format(time_df.time.max()/60), ' minutes')
    print('summarized time:', '\t \t \t \t', '{:10.2f}'.format(time_df.time.sum() / 60/60/24), ' days')


    print(error3.loc[:,1850].mean())
    print(error4.loc[:,1850].mean())

    error = error4.loc[:,1850]
    print(len(error[error<0]),len(error[error>0]) )
    '''