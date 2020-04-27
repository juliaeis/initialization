import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

def tranform_date(yr):
    if np.isnan(yr):
        return yr
    else:
        return int(str(yr)[:4])

rgi='60'

dir = '/home/juliaeis/PycharmProjects/initialization/advanced_experiments/DOI-WGMS-FoG-2018-11'

if rgi=='60':
    df_links = pd.read_csv(os.path.join(dir, '00_rgi60_links.csv'))
    df_links = df_links.rename(columns={"RGIId": "RGI_ID", "FoGId": "WGMS_ID"})

else:
    df_links = pd.read_csv(os.path.join(dir, 'WGMS-FoG-2018-11-AA-GLACIER-ID-LUT.csv'),encoding='iso8859_15')

df_len = pd.read_csv(os.path.join(dir, 'WGMS-FoG-2018-11-C-FRONT-VARIATION.csv'), encoding='iso8859_15')

df_links = df_links.dropna(subset=['RGI_ID'])  # keep the ones with a valid RGI ID
print('Total number of RGI links: {}'.format(len(df_links)))

#df_links.loc[~df_links.duplicated('RGI_ID', keep='first')]

len(df_len[df_len.WGMS_ID.isin(df_links.WGMS_ID)])

'''
for i in df_links.WGMS_ID:
    rgi_id = df_links[df_links.WGMS_ID == i].RGI_ID.iloc[-1]

    for i in df_len[df_len.WGMS_ID ==i].index:
        df_len.loc[i,'RGI_ID'] = rgi_id

df_len = df_len.dropna(subset=['RGI_ID'])
if rgi=='60':
    df_len.to_pickle(os.path.join(dir,'length_obs_60.pkl'))
else:
    df_len.to_pickle(os.path.join(dir, 'length_obs.pkl'))

'''
if rgi=='60':
    df_len = pd.read_pickle(os.path.join(dir,'length_obs_60.pkl'))
else:
    df_len = pd.read_pickle(os.path.join(dir, 'length_obs.pkl'))

counts = df_len.RGI_ID.value_counts()
miss_df = pd.DataFrame()
'''

for rgi in counts.index:

    try:
        df = df_len[df_len.RGI_ID.isin([rgi])]
        #df = df_len[df_len.RGI_ID.isin(['RGI50-11.00897'])]
        #df = df_len[df_len.RGI_ID.isin(['RGI50-13.43207'])]

        df = df[df.WGMS_ID.isin([df.WGMS_ID.value_counts().index[0]])].dropna(subset=['FRONT_VARIATION'])
        df.REFERENCE_DATE = df.REFERENCE_DATE.apply(tranform_date)
        print(df.WGMS_ID.iloc[0], df.RGI_ID.iloc[0])

        if np.isnan(df.iloc[0].REFERENCE_DATE):

            sum_df = pd.DataFrame({'Year': df.iloc[0].Year,
                                   'dL': float(0),'remark':'reference date unknown'}, index=[0])
        else:
            sum_df = pd.DataFrame({'Year': int(str(df.iloc[0].REFERENCE_DATE)[:4]),
                                   'dL':float(0)}, index=[0])


        for i in df.index:
            ref_yr = df.loc[i,'REFERENCE_DATE']

            if not np.isnan(ref_yr):

                # get index of reference year
                ref_yr = int(str(ref_yr)[:4])
                j = sum_df[sum_df['Year']==ref_yr].index

                if len(j) !=0:
                    dl = (sum_df.loc[j, 'dL'].values + float(df.loc[i,'FRONT_VARIATION']))[0]
                else:

                    lf = sum_df[sum_df.Year < ref_yr]
                    if len(lf) !=0:

                        prev = df[df.Year < ref_yr].index[-1]
                        x1 = [df.loc[prev,'Year'],df.loc[i,'Year']]
                        x2 = [df.loc[prev,'REFERENCE_DATE'],df.loc[i,'REFERENCE_DATE']]
                        if np.isnan(x2[0]):
                            x2[0]=df.iloc[df.index.get_loc(prev)-1].Year
                        y = np.divide([df.loc[prev,'FRONT_VARIATION'],df.loc[i,'FRONT_VARIATION']],np.subtract(x1, x2))

                        f = interpolate.interp1d(x1, y)
                        prev_dL = sum_df.loc[sum_df[sum_df.Year==x1[0]].index,'dL'].values[0]
                        sum_df = sum_df.append({'Year':int(ref_yr),'dL':prev_dL+f(ref_yr)*(ref_yr-x1[0]), 'remark':'reference date not included'}, ignore_index=True)
                    else:
                        # append first reference year
                        sum_df = sum_df.append({'Year': int(ref_yr), 'dL': 0, 'remark':np.nan}, ignore_index=True)
                    j = sum_df[sum_df['Year'] == ref_yr].index
                    dl = (sum_df.loc[j, 'dL'].values + float(df.loc[i, 'FRONT_VARIATION']))[0]
                sum_df = sum_df.append({'Year':int(df.loc[i,'Year']),'dL':dl,'remark':np.nan},ignore_index=True)

            # if reference year unknown, take row before as reference
            else:

                dl = sum_df.iloc[-1].dL + float(df.loc[i, 'FRONT_VARIATION'])
                sum_df = sum_df.append({'Year': int(df.loc[i, 'Year']), 'dL': dl,'remark':'reference date unknown'}, ignore_index=True)


        sum_df = sum_df.sort_values('Year').reset_index(drop=True)

        miss = sum_df[sum_df.remark=='reference date not included']

        miss_df = miss_df.append({'WGMS_ID': int(df.WGMS_ID.iloc[0]), 'RGI_ID':df.RGI_ID.iloc[0],'missing dL': int(len(miss)),'Years': ",".join(str(int(x)) for x in miss.Year.values)}, ignore_index=True)


        fig, ax = plt.subplots(1,1)
        # all points
        sum_df.plot.scatter(x='Year',y='dL',ax=ax)
        # reference date not included

        try:
            sum_df[sum_df.remark=='reference date not included'].plot.scatter(x='Year',y='dL',ax=ax,color='r')
        except:
            pass
        try:
            sum_df[sum_df.remark == 'reference date unknown'].plot.scatter(x='Year', y='dL', ax=ax, color='orange')
        except:
            pass

        # plot observations
        sum_df.plot(x='Year', y='dL', zorder=0, ax=ax, color='C0', legend=None)

        try:
            sub = sum_df[sum_df.remark == 'reference date not included']
            sub.plot.scatter(x='Year', y='dL', ax=ax, color='C3', zorder=2,label='reference date not included')
            for i in sub.index:
                sum_df.loc[i-1:i].plot(x='Year',y='dL',ax=ax,color='C3', zorder=0, legend=None)
        except:
            pass


        try:
            sub = sum_df[sum_df.remark == 'reference date unknown']
            sub.plot.scatter(x='Year', y='dL', ax=ax, color='orange', zorder=2, label='reference unknown')
            for i in sub.index:
                sum_df.loc[i - 1:i].plot(x='Year', y='dL', ax=ax, color='orange',
                                           zorder=0, legend=None)

        except:
            pass
        plt.title(df.iloc[-1].NAME)
        plt.savefig(os.path.join(dir, 'plots', df.RGI_ID.iloc[0] + '.png'))
        plt.close()
        #plt.show()

    except:
       pass

miss_df = miss_df.sort_values(by='missing dL', ascending=False)
miss_df = miss_df[['WGMS_ID','RGI_ID','missing dL', 'Years']]

miss_df = miss_df.astype({"WGMS_ID": int, "missing dL":int})

miss_df = miss_df.set_index('WGMS_ID')
print(miss_df)
if rgi=='60':
    miss_df.to_csv(os.path.join(dir,'60_missing_reference_years.csv'))
else:
    miss_df.to_csv(os.path.join(dir, 'missing_reference_years.csv'))



    no_ref = list()
    for ix,ref_yr in zip(hf.REFERENCE_DATE.dropna().index,hf.REFERENCE_DATE.dropna()):
        ref_yr = int(str(ref_yr)[:4])

        # if ref_yr not in df
        if hf[hf.Year.isin([int(ref_yr)])].empty:

            index = int(hf.index[-1])+1
            data = pd.DataFrame({'Year':ref_yr,'FRONT_VARIATION':0,
                                 'RGI_ID':hf.loc[ix,'RGI_ID'],
                                 'WGMS_ID':hf.loc[ix,'WGMS_ID']},
                                index=[index])
            hf = hf.append(data,sort=True)
            no_ref.append(index)
    hf = hf.sort_values('Year')

    plt.figure()
    for ix in hf[['Year','REFERENCE_DATE']].dropna().index:
        yr = hf.loc[ix,'Year']
        cl = hf.FRONT_VARIATION.cumsum()
        ref_yr = int(str(hf.loc[ix, 'REFERENCE_DATE'])[:4])
        ref_ix = hf[hf.Year.isin([int(ref_yr)])].index
        for ref_x in ref_ix:
            plt.plot([ref_yr,yr],[cl.loc[ref_x],cl.loc[ix]],color='C0')

    #filter "REFERENCE_YEAR unknown"
    unkw = hf.REFERENCE.str.find('REFERENCE_YEAR unknown').dropna()
    unkw = unkw[unkw!=-1].index


    plt.plot(hf.Year.values, hf.FRONT_VARIATION.cumsum().values,'o', color='C0')
    plt.plot(hf.loc[no_ref].Year.values, hf.FRONT_VARIATION.cumsum().loc[no_ref].values,'o', color='C1', label='reference year not included')
    plt.plot(hf.loc[unkw].Year.values, hf.FRONT_VARIATION.cumsum().loc[unkw].values,'o', color='C3', label='reference year unknown')
    plt.legend(loc='best')
    plt.title(hf.iloc[-1].NAME)
    plt.savefig(os.path.join(dir,'plots',hf.loc[ix,'RGI_ID']+'.png'))

print(hf[['RGI_ID','Year','REFERENCE_DATE','FRONT_VARIATION']])

'''
