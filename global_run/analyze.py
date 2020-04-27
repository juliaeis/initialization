import os
import pandas as pd
import numpy
import matplotlib.pyplot as plt



WORKING_DIR = '/home/juliaeis/Dokumente/OGGM/work_dir/reconstruction/global'

df = pd.DataFrame()
for i in range(1,19):
    df_i = pd.read_pickle(os.path.join(WORKING_DIR,'volume'+str(i).zfill(2)+'.pkl'))
    df_i = df_i.set_index('rgi_id')
    df = df.append(df_i)
    df_i = df_i.dropna(subset=[1917.0])
    df_i = df_i.dropna(axis=1)
    #df_i.sum().plot()
    #plt.show()
print(len(df)-len(df.dropna()))

volume1917 = df.loc[:,1917.0].sum()
mass1917 = volume1917*0.9167
sll1917 = mass1917*(1/361.8)

print(volume1917)
volume2009 = df.loc[:,2016.0].sum()
mass2009 = volume2009*0.9167
sll2009 = mass2009*(1/361.8)
print(volume2009)

print(170*(1/361.8))

df = df.dropna()
df.sum().plot()
plt.xlabel('Time (years)')
plt.ylabel(r'Volume (km$^3$)')
plt.title('Global Glacier Volume, n='+str(len(df)))
plt.show()