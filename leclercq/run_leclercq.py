
import os
import sys
import pandas as pd
import geopandas as gpd
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from oggm import cfg, utils, workflow, tasks
from oggm.utils._downloads import get_demo_file
from oggm.core.flowline import FluxBasedModel, FileModel
import oggm
from leclercq.leclercq_plots import *
sys.path.append('../')
sys.path.append('../../')
from  initialization.core import *
import advanced_experiments
import copy

mpl.rcParams['axes.linewidth'] = 3
mpl.rcParams['xtick.major.width'] = 3
mpl.rcParams['ytick.major.width'] = 3
mpl.rcParams['font.size'] =25
mpl.rcParams['font.weight'] = 'medium'
mpl.rcParams['axes.labelweight'] = 'medium'
mpl.rcParams['legend.fontsize'] = 20 #30
mpl.rcParams['lines.linewidth'] = 3


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

        return deepcopy(df.loc[df.fitness.idxmin()].model), None, None


def _run_experiment(gdir, temp_bias, bias, ys,ye):
    """
    Creates the synthetic experiment for one glacier. model_run_experiment.nc
    will be saved in working directory.
    """

    # check, if this experiment already exists
    try:
        rp = gdir.get_filepath('model_run', filesuffix='_advanced_experiment_' +str(temp_bias)+'_'+ str(bias))
        model = FileModel(rp)

    # otherwise create experiment
    except:

        fls = gdir.read_pickle('model_flowlines')

        try:
            model = tasks.run_random_climate(gdir, nyears=600, y0=ys, bias=bias, seed=1,
                                             temperature_bias=temp_bias,
                                             init_model_fls=fls,output_filesuffix='_random_experiment_'+str(temp_bias)+'_'+str(bias) )

            # construct observed glacier, previous glacier will be run forward from
            # 1917 - rgi_date with past climate file

            fls = deepcopy(model.fls)
            tasks.run_from_climate_data(gdir, ys=ys, ye=ye, init_model_fls=fls,bias=bias,
                                        output_filesuffix='_advanced_experiment_'+str(temp_bias)+'_'+str(bias))
            # to return FileModel
            rp = gdir.get_filepath('model_run',filesuffix='_advanced_experiment_' + str(
                                       temp_bias) + '_' + str(bias))
            model = FileModel(rp)


        except:
            pass

    return model


def find_residual(gdir, temp_bias_list, ys, a=-2000, b=2000):
    best_df = pd.DataFrame()

    fls = gdir.read_pickle('model_flowlines')
    mod = FluxBasedModel(flowlines=fls)

    for temp_bias in temp_bias_list:

        try:
            ye = gdir.rgi_date
            max_it = 15
            i = 0
            bounds = [a, b]

            df = pd.DataFrame()

            while i < max_it:
                bias = round((bounds[0] + bounds[1]) / 2, 1)

                ex_mod2 = _run_experiment(gdir, temp_bias, bias, ys, ye)

                diff = mod.area_km2 - ex_mod2.area_km2_ts()[ye]

                df = df.append(pd.Series({'bias': bias, 'area_diff': diff}),
                               ignore_index=True)

                if (abs(diff) < 1e-4) or bounds[1] - bounds[0] <= 1:
                    break

                elif ex_mod2.area_km2_ts()[ye] > mod.area_km2:
                    bounds[0] = bias
                else:
                    bounds[1] = bias
                i += 1

            # best bias found
            bias = df.iloc[df.area_diff.abs().idxmin()].bias

            rp = gdir.get_filepath('model_run',
                                   filesuffix='_advanced_experiment_' + str(
                                       temp_bias) + '_' + str(bias))
            model = FileModel(rp)

            diff = gdir.rgi_area_km2 - model.area_km2_ts()[gdir.rgi_date]

            series = pd.Series(
                {'rgi_id': gdir.rgi_id, 'bias': bias, 'iterations': i,
                 'area_diff': diff, 'model': model, 'temp_bias': temp_bias})

        except:
            series = pd.Series({'rgi_id': gdir.rgi_id, 'temp_bias': temp_bias})
        best_df = best_df.append(series, ignore_index=True)

    return best_df


def get_ref_length_data(gdir):
    """Get the glacier lenght data from P. Leclercq's data base.

     https://folk.uio.no/paulwl/data.php

     For some glaciers only!
     """

    df = pd.read_csv('rgi_leclercq_links_2014_RGIV6.csv')
    df = df.loc[df.RGI_ID == gdir.rgi_id]
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
    return df


def advanced_experiments(gdirs, temp_bias_list ,ys , region):

    exp_df = pd.DataFrame()

    pool = Pool()
    list = pool.map(partial(find_residual,temp_bias_list=temp_bias_list,ys=ys),gdirs)
    pool.close()
    pool.join()

    exp_df = exp_df.append(list, ignore_index=True)
    p = os.path.join(cfg.PATHS['working_dir'], str(region)+'_advanced_experiments.pkl')
    exp_df.to_pickle(p, compression='gzip')
    return p

if __name__ == '__main__':

    cfg.initialize()

    ON_CLUSTER = True
    REGION = '11'

    # Local paths
    if ON_CLUSTER:
        WORKING_DIR = os.environ.get("S_WORKDIR")
        cfg.PATHS['working_dir'] = WORKING_DIR
        TEMP_BIAS = float(os.environ.get('TEMP_BIAS'))
        REGION = str(os.environ.get('I')).zfill(2)
        if float(REGION)<=19:
            REGION_NAME = REGION+'a'
        elif float(REGION)<=38:
            REGION = str(int(REGION)-19).zfill(2)
            REGION_NAME = REGION+'b'
        elif float(REGION)<=57:
            REGION = str(int(REGION)-38).zfill(2)
            REGION_NAME = REGION+'c'
        else:
            REGION = str(int(REGION)-57).zfill(2)
            REGION_NAME = REGION+'d'

    else:
	WORKING_DIR = '/home/juliaeis/Dokumente/OGGM/work_dir/reconstruction/leclercq'
        cfg.PATHS['working_dir'] = WORKING_DIR
        utils.mkdir(WORKING_DIR, reset=False)
    print(REGION_NAME)
    print(TEMP_BIAS)

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

    # read leclercq links
    lec = pd.read_csv('rgi_leclercq_links_2014_RGIV6.csv')
    lec['REGION'] = lec.RGI_ID.apply(lambda x: x.split('-')[-1].split('.')[0])

    # We use intersects
    db = utils.get_rgi_intersects_region_file(version='61', region=REGION)
    cfg.set_intersects_db(db)

    # RGI file
    path = utils.get_rgi_region_file(REGION, version='61')
    rgidf = gpd.read_file(path)

    # only the ones with leclercq observation
    rgidf = rgidf[rgidf.RGIId.isin(lec[lec.REGION==REGION].RGI_ID.values)]

     # exclude non-landterminating glaciers
    rgidf = rgidf[rgidf.TermType==0]
    rgidf = rgidf[rgidf.Connect !=2]

    # sort for efficient using
    rgidf = rgidf.sort_values('Area', ascending=True)

    gdirs = workflow.init_glacier_regions(rgidf)
    if REGION_NAME.endswith('a'):
        gdirs = gdirs[0:len(gdirs):4]
    elif REGION_NAME.endswith('b'):
        gdirs = gdirs[1:len(gdirs):4]
    elif REGION_NAME.endswith('c'):
        gdirs = gdirs[2:len(gdirs):4]
    else:
	gdirs = gdirs[3:len(gdirs):4]

    preprocessing(gdirs)
    p =advanced_experiments(gdirs,[TEMP_BIAS] ,1917 ,REGION_NAME)
    df = pd.read_pickle(p, compression='gzip')
    df = df.set_index('rgi_id')

    diff = pd.DataFrame()
    for gdir in gdirs:
        try:
            ye = gdir.rgi_date

            bias = df.loc[gdir.rgi_id].bias
            ex_mod = df.loc[gdir.rgi_id].model

            lec = get_ref_length_data(gdir).dL
            lec = lec[lec.index <= ye]
            lec = (lec - lec.iloc[-1]) + ex_mod.length_m_ts(rollmin=5)[
                lec.index[-1]]

            ini_df = find_possible_glaciers(gdir, 1917,ye, 200, ex_mod, bias)
            ini_df.fitness = pd.to_numeric(ini_df.fitness/125)
            ini_df = ini_df.dropna(subset=['fitness'])

            med_mod, perc_min, perc_max = find_median(ini_df)

            lec.loc[1917]= np.nan
            lec = lec.sort_index().interpolate()[lec.index >= 1917]

            rmse = np.sqrt(((lec - med_mod.length_m_ts(rollmin=5)[lec.index]) ** 2).mean())
            rmspe = np.sqrt((((lec - med_mod.length_m_ts(rollmin=5)[lec.index])/lec) ** 2).mean())*100
            error = (lec - med_mod.length_m_ts(rollmin=5)[lec.index]).mean()
            perc_error = ((lec - med_mod.length_m_ts(rollmin=5)[lec.index])/lec).mean()

            max = (lec - med_mod.length_m_ts(rollmin=5)[lec.index]).max()
            min = (lec - med_mod.length_m_ts(rollmin=5)[lec.index]).min()
            if abs(max)>abs(min):
                max_diff = max
            else:
                max_diff = min

            temp_bias = cfg.PATHS['working_dir'].split('_')[-1]

            diff.loc[gdir.rgi_id, 'region'] = REGION
            diff.loc[gdir.rgi_id,'rmse']=rmse
            diff.loc[gdir.rgi_id, 'rmspe'] = rmspe
            diff.loc[gdir.rgi_id, 'error'] = error
            diff.loc[gdir.rgi_id, 'perc_error'] = perc_error
            diff.loc[gdir.rgi_id,'max_diff'] = max_diff


            plot_fitness_values(gdir, ini_df, ex_mod, 1917, ye, lec, cfg.PATHS['plot_dir'])
            plot_median(gdir, ini_df, 125, ex_mod, 1917, ye, lec, cfg.PATHS['plot_dir'])

        except:
            pass

    diff.to_pickle(os.path.join(cfg.PATHS['working_dir'], REGION_NAME+'_leclercq_difference.pkl'))


'''


    #except:
    #    pass


    diff = pd.read_pickle(os.path.join(cfg.PATHS['working_dir'],'leclercq_difference.pkl'))
    diff.sort_values(by='region')
    diff.rmse = diff.rmse/1000
    diff.max_diff = diff.max_diff/1000

    print(diff[diff.region=='01'].sort_values(by='rmse').rmse)



    fig,ax = plt.subplots(1,1)
    bp = diff.boxplot(column='rmse', by='region',return_type='dict',ax=ax)
    # Calculate number of obs per group & median to position labels

    nobs = diff['region'].value_counts().sort_index().values
    nobs = [str(x) for x in nobs.tolist()]
    nobs = ["n: " + i for i in nobs]

    # Add it to the plot
    pos = range(len(nobs))
    for tick, label in zip(pos, ax.get_xticklabels()):

        ax.text(pos[tick]+1, -0.15, nobs[tick],
                horizontalalignment='center', size='x-small', color='k',
                weight='semibold')

    [[item.set_linewidth(2) for item in key[item]] for key in bp for item in bp['rmse'].keys()]
    [[item.set_color('red') for item in key['medians']] for key in bp]


    #plt.yticks(np.arange(0, 15, step=1))
    plt.ylabel('Root Mean Square Error (km)')
    plt.xlabel('Region')
    plt.title('Comparison with Leclercq')
    plt.suptitle('')

    # maximum difference
    fig, ax = plt.subplots(1, 1)
    bp = diff.boxplot(column='max_diff', by='region', return_type='dict', ax=ax)
    # Calculate number of obs per group & median to position labels

    nobs = diff['region'].value_counts().sort_index().values
    nobs = [str(x) for x in nobs.tolist()]
    nobs = ["n: " + i for i in nobs]

    # Add it to the plot
    pos = range(len(nobs))
    for tick, label in zip(pos, ax.get_xticklabels()):
        ax.text(pos[tick] + 1, -0.4, nobs[tick],
                horizontalalignment='center', size='x-small', color='k',
                weight='semibold')

    [[item.set_linewidth(2) for item in key[item]] for key in bp for item in
     bp['max_diff'].keys()]
    [[item.set_color('red') for item in key['medians']] for key in bp]

    # plt.yticks(np.arange(0, 15, step=1))
    plt.ylabel('Root Mean Square Error (km)')
    plt.xlabel('Region')
    plt.title('Comparison with Leclercq')
    plt.suptitle('')


    plt.show()


    p = os.path.join(cfg.PATHS['working_dir'],'11_advanced_experiments.pkl')
    df = pd.read_pickle(p,compression='gzip')
    df = df.set_index('rgi_id')
    df.fitness = df.fitness/125
    print(len(df.dropna()))

    df.bias.plot.hist()
    plt.title('mb bias, Alps (n=85)')
    plt.show()


    df = pd.read_pickle(os.path.join(cfg.PATHS['working_dir'],'09_leclercq_difference.pkl'))
    print(df)
'''