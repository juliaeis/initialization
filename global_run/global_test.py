import sys
sys.path.append('../')
from initialization.core import *
from paper.plots_paper import *



import geopandas as gpd
from oggm import cfg, utils
pd.options.mode.chained_assignment = None


if __name__ == '__main__':
    cfg.initialize()

    ON_CLUSTER = False

    # Local paths
    if ON_CLUSTER:
        WORKING_DIR = os.environ.get("S_WORKDIR")
        cfg.PATHS['working_dir'] = WORKING_DIR
        OUT_DIR = os.environ.get("OUTDIR")
        REGION = str(os.environ.get('REGION')).zfill(2)
        ID =  os.environ.get("I")
    else:
        WORKING_DIR = '/home/juliaeis/Dokumente/OGGM/work_dir/reconstruction/global/'
        OUT_DIR = WORKING_DIR
        cfg.PATHS['working_dir'] = WORKING_DIR
        utils.mkdir(WORKING_DIR, reset=False)
        REGION='01'

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
    rgidf = rgidf[rgidf.RGIId == 'RGI60-01.23037']
    #select HEF
    gdirs = workflow.init_glacier_regions(rgidf)


    t_0 = 1917
    epsilon = 125
    #preprocessing(gdirs)
    #advanced_experiments(gdirs, [0], 1917, REGION)

    for gdir in gdirs:

        #dir = os.path.join(OUT_DIR, gdir.dir.split('/global/')[-1])
        '''
        dir = os.path.join(OUT_DIR, gdir.dir.split('/global/')[-1])
        os.system('cp -rf '+dir+'/* '+ gdir.dir)

        t_e = gdir.rgi_date
        ex = [f for f in os.listdir(gdir.dir) if f.startswith('model_run_ad')]
        if len(ex)==1 :
            dst = os.path.join(gdir.dir,ex[0])
            ex_mod = FileModel(dst)
            bias = float(ex[0].split('_')[-1].split('.nc')[0])
            try:
                df = find_possible_glaciers(gdir, t_0, t_e, 200, ex_mod, bias, delete=True)
                print(df)
                #save = pickle.load(open(os.path.join(gdir.dir, 'initialization_output.pkl'),'rb'))
                print(save)
            except Exception as e:
                print(e)
        '''
        save = {}
        results = pd.read_pickle(os.path.join(gdir.dir,'result1917.pkl'), compression='gzip')


        # minimum
        save.update({'minimum_exp': results.loc[results.fitness.idxmin(), 'model'].split('/')[-1]})
        save.update({'minimum_fls': results.loc[results.fitness_fls.idxmin(), 'model'].split('/')[-1]})

        # acceptable
        results_exp = results[results.fitness <=1]
        if len(results_exp)>0:
            # acceptable
            save.update({'acc_min_exp':results_exp.loc[results_exp.length.idxmin(),'model'].split('/')[-1]})
            save.update({'acc_max_exp': results_exp.loc[results_exp.length.idxmax(), 'model'].split('/')[-1]})

            # 5th percentile
            results_exp = results_exp[results_exp.fitness <= results_exp.fitness.quantile(0.05)]
            save.update({'perc_min_exp': results_exp.loc[results_exp.length.idxmin(), 'model'].split('/')[-1]})
            save.update({'perc_max_exp': results_exp.loc[results_exp.length.idxmax(), 'model'].split('/')[-1]})

            # median
            results_exp = results_exp.sort_values(by='length')
            l1 = len(results_exp)
            if l1 % 2:
                index_exp = int((l1 - 1) / 2)
            else:
                index_exp = int(l1 / 2)
            save.update({'median_exp': results_exp.iloc[index_exp].model.split('/')[-1]})

        results_fls = results[results.fitness_fls <= 1]
        if len(results_fls)>0:
            save.update({'acc_min_fls': results_fls.loc[results_fls.length.idxmin(), 'model'].split('/')[-1]})
            save.update({'acc_max_fls': results_fls.loc[results_fls.length.idxmax(), 'model'].split('/')[-1]})

            results_fls = results_fls[results_fls.fitness_fls <= results_fls.fitness_fls.quantile(0.05)]
            save.update({'perc_min_fls': results_fls.loc[results_fls.length.idxmin(), 'model'].split('/')[-1]})
            save.update({'perc_max_fls': results_fls.loc[results_fls.length.idxmax(), 'model'].split('/')[-1]})


            results_fls = results_fls.sort_values(by='length')
            l2 = len(results_fls)
            if l2 % 2:
                index_fls = int((l2 - 1) / 2)
            else:
                index_fls = int(l2 / 2)

            save.update({'median_fls': results_exp.iloc[index_fls].model.split('/')[-1]})




