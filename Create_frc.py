#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 19 11:39:22 2022

@author: thorbjornostenbymoe
"""
""" 

This script creates a forcing file test for ROMS idealized model runs. (basically a copy of the input script)

Currently, it simply uses the uppermost levels of temp and salt from the init file (created in Create_input.py) and 
appropriatly dimensioned zero matrices for all other variables.

For later use it should be sufficient to change the variables, as the writing to netCDF is 
made in accordance to the ROMS standards.

"""
###############################################################################
################################## - Packages - ###############################
###############################################################################

import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt

###############################################################################
############################### - User options - ##############################
###############################################################################

#File directory/names
directory       = '/Users/thorbjornostenbymoe/Desktop/Initialization_of_model/MODEL_SETUP/' #User directory/path
frc_file        = 'ROMS_in/frc_file.nc'
grid_file       = 'ROMS_in/grid_file.nc'
Input_file      = 'ROMS_in/init_file.nc'

ds_grid = nc.Dataset(directory+grid_file)
ds_init = nc.Dataset(directory+Input_file)


#file description
Descrip_grd   = 'Iceland slope - ROMS idealized model run - frc file'  
Author        = 'Thorbjoern Oestenby Moe'

#logical switches
create_frc_file = 1                                                   # create input (1 = true)
plot = 1

#Dims
x = 268
y = 401
z = 30                                
#%%  

###############################################################################
################################# - func(s) - #################################
###############################################################################

def double_trouble(old_var):             # add arbitrary dimension (dim = 2)
    return np.array([old_var,old_var])

###############################################################################
################################ - load input - ###############################
###############################################################################

salt = ds_init['salt'][0,29,:,:]
temp = ds_init['temp'][0,29,:,:]

###############################################################################
############################ - create variables - #############################
###############################################################################
""" all variables except temp and salt are set to zero """

# salt/temp
SSS_var = np.array([salt,salt])
SST_var = np.array([temp,temp])

#time
sms_time_var = np.array([0, 180])
sss_time_var = np.array([0, 180])
sst_time_var = np.array([0, 180])
#momentum
sustr_var = np.zeros((2, 268, 401))
svstr_var = np.zeros((2, 267, 402))
# heat flux
dQdSST_var = np.ones((2, 268, 402))*-40


#%%

if plot == 1:
    X,Y = np.meshgrid(np.linspace(0,402*1.5,402),np.linspace(0,268*1.5,268))
    fig,axs = plt.subplots(1, 2, figsize=(7,2))
    
    c1 = axs[0].contourf(X,Y,SST_var[0,:,:])
    plt.colorbar(c1,ax=axs[0])
    
    c2 = axs[1].contourf(X,Y,SSS_var[0,:,:])
    plt.colorbar(c2,ax=axs[1])
    
#%%
if create_frc_file == 1:
    ###############################################################################
    ################################### - file - ##################################
    ###############################################################################
    
    path_string = directory+frc_file
    ncid=nc.Dataset(path_string,mode='w',format='NETCDF4')
    ncid.description=Descrip_grd
    ncid.author=Author
    
    ###############################################################################
    ############################### - dimentions - ################################
    ###############################################################################
    
    # Meridional: eta
    eta_rho=ncid.createDimension('eta_rho',x)
    eta_u=ncid.createDimension('eta_u',x)
    eta_v=ncid.createDimension('eta_v',x-1)
    eta_psi=ncid.createDimension('eta_psi',x-1)
    
    # Zonal: xi
    xi_rho=ncid.createDimension('xi_rho',y+1)
    xi_u=ncid.createDimension('xi_u',y)
    xi_v=ncid.createDimension('xi_v',y+1)
    xi_psi=ncid.createDimension('xi_psi',y)
    
    #time
    sms_time=ncid.createDimension('sms_time',2)
    sss_time=ncid.createDimension('sss_time',2)
    sst_time=ncid.createDimension('sst_time',2)
    
    ###############################################################################
    ############################### - Variables - #################################
    ###############################################################################
    
    #time
    sms_time=ncid.createVariable('sms_time','f8',('sms_time'))
    sms_time.units = 'day'
    sms_time.cycle_length = 360.
    sms_time.field = 'sms_time, scalar, series'
    sms_time.long_name = 'surface momentum stress time'
    
    sss_time=ncid.createVariable('sss_time','f8',('sss_time'))
    sss_time.units = 'day'
    sss_time.cycle_length = 360.
    sss_time.field = 'sms_time, scalar, series'
    sss_time.long_name = 'surface salinity time'

    sst_time=ncid.createVariable('sst_time','f8',('sst_time'))
    sst_time.units = 'day'
    sst_time.cycle_length = 360.
    sst_time.field = 'sms_time, scalar, series'
    sst_time.long_name = 'surface temperature time'
    
    #momentum
    sustr=ncid.createVariable('sustr','f8',('sms_time','eta_u','xi_u',))
    sustr.long_name='surface u-momentum stress'
    sustr.time = 'sms_time'
    sustr.units='N/m2'
    sustr.field='surface u-momentum stress, scalar, series'
    
    svstr=ncid.createVariable('svstr','f8',('sms_time','eta_v','xi_v',))
    svstr.long_name='surface V-momentum stress'
    svstr.time = 'sms_time'
    svstr.units='N/m2'
    svstr.field='surface v-momentum stress, scalar, series'
    
    #temp/salt
    SSS=ncid.createVariable('SSS','f8',('sss_time','eta_rho','xi_rho',)) 
    SSS.long_name='sea surface salinity climatology'
    SSS.time = 'sss_time'
    SSS.units='PSU'
    SSS.field='SSS,scalar,series'    
    
    SST=ncid.createVariable('SST','f8',('sst_time','eta_rho','xi_rho',))
    SST.long_name='sea surface temperature climatology'
    SST.units='Celsius'
    SST.time = 'sst_time'
    SST.field='SST,scalar,series' 
    
    #heat flux
    dQdSST=ncid.createVariable('dQdSST','f8',('sst_time','eta_rho','xi_rho',))
    dQdSST.long_name='surface net heat flux sensitivity to SST'
    dQdSST.time='sst_time'
    dQdSST.units='Watts/m2/Celsius'
    dQdSST.field='dQdSST, scalar, series'
    
    ###############################################################################
    ############################## - Assign values - ##############################
    ###############################################################################
    
    sms_time[:] = sms_time_var
    sss_time[:] = sss_time_var
    sst_time[:] = sst_time_var
    
    sustr[:,:,:] = sustr_var
    svstr[:,:,:] = svstr_var
    
    SSS[:,:,:] = SSS_var
    SST[:,:,:] = SST_var
    dQdSST[:,:,:] = dQdSST_var
    
    ncid.close()










