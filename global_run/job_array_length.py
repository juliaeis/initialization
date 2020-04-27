from oggm import cfg, utils
import geopandas as gpd
import pandas as pd
import math


df = pd.DataFrame()

for region in range(1,19):
    region = str(region).zfill(2)
    print(region)
    # RGI file
    path = utils.get_rgi_region_file(region, version='61')
    rgidf = gpd.read_file(path)
    #rgidf = rgidf.sort_values('Area', ascending=False)

    # exclude non-landterminating glaciers
    #rgidf = rgidf[rgidf.TermType == 0]
    #rgidf = rgidf[rgidf.Connect != 2]

    rgidf = rgidf.sort_values('Area', ascending=False)

    l = len(rgidf[rgidf.Area>=1])+(len(rgidf[rgidf.Area<1])/30)
    df.loc[region, 'n_array1'] = int(math.ceil((len(rgidf)/50)-1))
    df.loc[region, 'n_array2'] = int(math.ceil(l))
    df.loc[region, 'n_glacier'] = len(rgidf)
print(df.n_array2)


    #print(rgidf[int(6401):len(rgidf):math.ceil(len(rgidf[rgidf.Area < 1]) / 30)])


'''
    if I <= len(rgidf[rgidf.Area >= 1]):
        rgidf = rgidf[I:I + 1]
    else:
        rgidf = rgidf[I:len(rgidf):math.ceil(len(rgidf[rgidf.Area < 1]) / 30)]
    print(rgidf)




df.index.name = 'REGION'
df = df.astype({'n_array1': 'int32','n_array2':'int32'})
print(df.n_glacier.sum())
#df.to_csv('number_job_arrays.csv')
'''