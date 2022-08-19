#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 18 09:44:28 2022

@author: thorbjornostenbymoe

"""
import time
start = time.time()

import scipy.io
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

def filtering(data,sigma):
    return(gaussian_filter(data, sigma=sigma))

mypath = '/Users/thorbjornostenbymoe/Desktop/Initialization_of_model/MODEL_SETUP/'
mat = scipy.io.loadmat(mypath+'InputProfileCentralIcelandSea.mat') #load data

###############################################################################
################################## - Sizes - ##################################
###############################################################################

n_profiles = 300
max_depth_data = 2180
max_depth_domain = 1152


###############################################################################
############ - Calculate mean Termperature and Salinity profiles - ############
###############################################################################

temp = mat['IS']['temp'] #retrieve temp, sal and press from dataset
press = mat['IS']['press']
sal = mat['IS']['sal']

#arbitrary depth/zero vectors of maximum depth in CTD's
h = np.linspace(0,max_depth_data-1,max_depth_data)
h1 = np.linspace(0,0,max_depth_data)


#create empty dataset with 0-vec for temp, sal and index to calculate mean
df_temp = pd.DataFrame(data=h1,index=h)
for i in range(n_profiles):
    df_temp[str(i)]=h1
df_temp = df_temp.drop(0,axis=1)

df_salt = pd.DataFrame(data=h1,index=h)
for i in range(n_profiles):
    df_salt[str(i)]=h1
df_salt = df_salt.drop(0,axis=1)

df_press = pd.DataFrame(data=h1,index=h)
for i in range(n_profiles):
    df_press[str(i)]=h1
df_press = df_press.drop(0,axis=1)

#simpler averaging loop heavily dependent on indexing
for i,j,k,ind in zip(press[0],temp[0],sal[0],range(len(press[0]))):
    b_t = np.nan * np.ones(shape=(max_depth_data))
    b_t[i.astype(int)] = j
    df_temp.iloc[:,ind] = b_t

    b_s = np.nan * np.ones(shape=(max_depth_data))
    b_s[i.astype(int)] = k
    df_salt.iloc[:,ind] = b_s
    
    b_p = np.nan * np.ones(shape=(max_depth_data))
    b_p[i.astype(int)] = i
    df_press.iloc[:,ind] = b_p
    
#get number of datapoints
index = df_press.iloc[0:max_depth_domain,:].sum(axis=1)
index[:] = index[:]/np.arange(0,max_depth_domain,1) 

#mean
mean_df_s1 = df_salt.iloc[0:max_depth_domain,:].mean(axis=1)

mean_df_t1 = df_temp.iloc[0:max_depth_domain,:].mean(axis=1)

mean_df_p = df_press.iloc[0:max_depth_domain,:].mean(axis=1)

#filter
mean_df_s = filtering(filtering(mean_df_s1,[3]),[3])
mean_df_t = filtering(filtering(mean_df_t1,[3]),[3])

"""
plt.figure(figsize=[5,7])
plt.plot(mean_df_t,-mean_df_p,'.',markersize=1)
plt.ylim([-max_depth_domain,0])

plt.figure(figsize=[5,7])
plt.plot(mean_df_s,-mean_df_p,'.',markersize=1)
plt.ylim([-max_depth_domain,0])
"""

###############################################################################
############################# - Save to NetCDF - ##############################
###############################################################################



fn = 'Final_mean_CTD.nc'
ds = nc.Dataset(mypath+fn, 'w', format='NETCDF4')

pressure = ds.createDimension('pressure', max_depth_domain)

pressure = ds.createVariable('pressure', 'f4', ('pressure',))
Temperature = ds.createVariable('Temperature', 'f4', ('pressure'))
Temperature.units = 'Celsius'
Salinity = ds.createVariable('Salinity', 'f4', ('pressure'))
Salinity.units = 'PSU'

pressure[:] = mean_df_p.to_numpy()

Temperature[:] = mean_df_t

Salinity[:] = mean_df_s

ds.close()

#timing loop
end = time.time()
total_time = end - start
print("\n"+ str(total_time))




