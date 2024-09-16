#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 14 15:10:36 2023

@author: thorbjorn
"""

""" 

This script creates a forcing file test for ROMS idealized model runs. (basically a copy of the input script)

Currently, it simply uses the uppermost levels of temp and salt from the init file (created in Make_input.py) and 
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
from GlobParam import *                # User input comes from here

print('###############################################################################')
print('############################# - Making frc file - #############################')
print('###############################################################################')

###############################################################################
########################### - Paths and settings - ############################
###############################################################################

#File directory/names
common_dir          = GlobParam.common_dir
directory           = GlobParam.EXP_dir

frc_file            = directory + GlobParam.frc_file
input_file          = directory + GlobParam.init_file
grid_file           = directory + GlobParam.grid_file


grid_ds             = nc.Dataset(grid_file)
init_ds             = nc.Dataset(input_file)

Descrip_grd         = GlobParam.General_desc + GlobParam.spec_desc_frc
Author              = GlobParam.Author
create_frc_file     = GlobParam.Create_files
plotting            = GlobParam.Gen_plot
z                   = GlobParam.z
res                 = GlobParam.Res

#Dims
x   = np.shape(grid_ds['h'])[0]
y   = np.shape(grid_ds['h'])[1]

###############################################################################
################################ - load input - ###############################
###############################################################################

# Extract surface fields from initial file and add arbitrary time dimension of size 2

SSS_var = init_ds['salt'][:,-1,:,:].repeat(2,0)
SST_var = init_ds['temp'][:,-1,:,:].repeat(2,0)

###############################################################################
############################ - create variables - #############################
###############################################################################
""" Variables (except dQ) are set to zero """

#time
sms_time_var    = np.array([0, 180])
sss_time_var    = np.array([0, 180])
sst_time_var    = np.array([0, 180])

#momentum
sustr_var       = np.zeros((2, x, y-1))
svstr_var       = np.zeros((2, x-1, y))

# heat flux
dQdSST_var      = np.ones((2, x, y))*-40


#%%
if create_frc_file == 1:
    ###############################################################################
    ################################### - file - ##################################
    ###############################################################################
    
    ncid                    = nc.Dataset(frc_file,mode='w',format='NETCDF4')
    ncid.description        = Descrip_grd
    ncid.author             = Author
    
    ###############################################################################
    ############################### - dimentions - ################################
    ###############################################################################
    
    # Meridional: eta
    eta_rho                 = ncid.createDimension('eta_rho',  x)
    eta_u                   = ncid.createDimension('eta_u',    x)
    eta_v                   = ncid.createDimension('eta_v',    x-1)
    eta_psi                 = ncid.createDimension('eta_psi',  x-1)
    
    # Zonal: xi
    xi_rho                  = ncid.createDimension('xi_rho',   y)
    xi_u                    = ncid.createDimension('xi_u',     y-1)
    xi_v                    = ncid.createDimension('xi_v',     y)
    xi_psi                  = ncid.createDimension('xi_psi',   y-1)
    
    #time
    sms_time                = ncid.createDimension('sms_time', 2)
    sss_time                = ncid.createDimension('sss_time', 2)
    sst_time                = ncid.createDimension('sst_time', 2)
    
    ###############################################################################
    ############################### - Variables - #################################
    ###############################################################################
    
    #time
    sms_time              = ncid.createVariable('sms_time','f8',('sms_time'))
    sms_time.units          = 'day'
    sms_time.cycle_length   = 360.
    sms_time.field          = 'sms_time, scalar, series'
    sms_time.long_name      = 'surface momentum stress time'
    
    sss_time              = ncid.createVariable('sss_time','f8',('sss_time'))
    sss_time.units          = 'day'
    sss_time.cycle_length   = 360.
    sss_time.field          = 'sms_time, scalar, series'
    sss_time.long_name      = 'surface salinity time'

    sst_time              = ncid.createVariable('sst_time','f8',('sst_time'))
    sst_time.units          = 'day'
    sst_time.cycle_length   = 360.
    sst_time.field          = 'sms_time, scalar, series'
    sst_time.long_name      = 'surface temperature time'
    
    #momentum
    sustr                 = ncid.createVariable('sustr','f8',('sms_time','eta_u','xi_u',))
    sustr.long_name         = 'surface u-momentum stress'
    sustr.time              = 'sms_time'
    sustr.units             = 'N/m2'
    sustr.field             = 'surface u-momentum stress, scalar, series'
    
    svstr                 = ncid.createVariable('svstr','f8',('sms_time','eta_v','xi_v',))
    svstr.long_name         = 'surface V-momentum stress'
    svstr.time              = 'sms_time'
    svstr.units             = 'N/m2'
    svstr.field             = 'surface v-momentum stress, scalar, series'
    
    #temp/salt
    SSS                   = ncid.createVariable('SSS','f8',('sss_time','eta_rho','xi_rho',)) 
    SSS.long_name           = 'sea surface salinity climatology'
    SSS.time                = 'sss_time'
    SSS.units               = 'PSU'
    SSS.field               = 'SSS,scalar,series'    
    
    SST                   = ncid.createVariable('SST','f8',('sst_time','eta_rho','xi_rho',))
    SST.long_name           = 'sea surface temperature climatology'
    SST.units               = 'Celsius'
    SST.time                = 'sst_time'
    SST.field               = 'SST,scalar,series' 
    
    #heat flux
    dQdSST                = ncid.createVariable('dQdSST','f8',('sst_time','eta_rho','xi_rho',))
    dQdSST.long_name        = 'surface net heat flux sensitivity to SST'
    dQdSST.time             = 'sst_time'
    dQdSST.units            = 'Watts/m2/Celsius'
    dQdSST.field            = 'dQdSST, scalar, series'
    
    ###############################################################################
    ############################## - Assign values - ##############################
    ###############################################################################
    
    sms_time[:]     = sms_time_var
    sss_time[:]     = sss_time_var
    sst_time[:]     = sst_time_var
    
    sustr[:,:,:]    = sustr_var
    svstr[:,:,:]    = svstr_var
    
    SSS[:,:,:]      = SSS_var
    SST[:,:,:]      = SST_var
    
    dQdSST[:,:,:]   = dQdSST_var
    
    ncid.close()



