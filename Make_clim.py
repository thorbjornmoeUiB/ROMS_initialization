#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 14 15:12:03 2023

@author: thorbjorn
"""
""" 
This script creates a climatology file for ROMS idealized model runs. (basically a copy of input script)

In this test file, all variables except T & S are set as appropriately dimensioned 0 matrices.

The transformation to sigma coordinates are made from the stretching vector (retrieved from the model) 
and the depth coordinates from the grid file (made in Make_grid.py).

For later use, the "Loop to create 3D temp/sali domains" secton will have to be significantly modified along with all the 
sections retriving data used in the loop. Additionally, the "write empty zeta, u, v, ubar, vbar data" section must be 
exchanged with actual observations. 

For later use it should be sufficient to change the variables, as the writing to netCDF is 
made in accordance to the ROMS standards.

"""
###############################################################################
################################## - Packages - ###############################
###############################################################################

import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
from GlobParam import *                # This is where I change user outputs now

print('###############################################################################')
print('############################# - Making clm file - #############################')
print('###############################################################################')
###############################################################################
############################### - User options - ##############################
###############################################################################

#File directory/names
common_dir          = GlobParam.common_dir
directory           = GlobParam.EXP_dir

clm_file            = directory + GlobParam.clm_file
grid_file           = directory + GlobParam.grid_file
Input_file          = directory + GlobParam.init_file
bry_file            = directory + GlobParam.bry_file
stretching_file     = common_dir + GlobParam.Streching


grid_ds             = nc.Dataset(grid_file)
init_ds             = nc.Dataset(Input_file)
bnry_ds             = nc.Dataset(bry_file)
stre_ds             = nc.Dataset(stretching_file)

Descrip_grd         = GlobParam.General_desc + GlobParam.spec_desc_clm
Author              = GlobParam.Author
create_clm_file     = GlobParam.Create_files
plotting            = GlobParam.Gen_plot
Iplottin            = GlobParam.Imp_plot
z                   = GlobParam.z
res                 = GlobParam.Res


pNIIC_w         = GlobParam.pNIIC_w

#Dims
x               = np.shape(grid_ds['h'])[0]                             # dim meridional ~ 400 km   (1 km resolution)  
y               = np.shape(grid_ds['h'])[1]                             # dim zonal      ~ 650 km  

#nudge zone 
nud             = 17   
spe             = 0.10

    
###############################################################################
###################### - load temp/salt from init file - ######################
###############################################################################

salt_from_init = init_ds['salt'][:,:,:,:].repeat(2,0)
temp_from_init = init_ds['temp'][:,:,:,:].repeat(2,0)

###############################################################################
########## - Empty arrays for velocities and sea surface elevation - ##########
###############################################################################
#%%
ds_u_3d         = np.zeros((z, x, y-1))
ds_u_4d         = np.zeros((2, z, x, y-1))
ds_ubar_3d      = np.zeros((2, x, y-1))

ds_v_3d         = np.zeros((z, x-1, y))
ds_v_4d         = np.zeros((2, z, x-1, y))
ds_vbar_3d      = np.zeros((2, x-1, y))

ds_zeta_3d      = np.shape((2, x, y))

###############################################################################
######################## - Load properties from bry - #########################
###############################################################################

# test boudnary file with temp and velocity
# bnry_ds = nc.Dataset('.../Iceland_Sea_Mod/NewConfig/INPUT/VarSlo_WINTER_2xDomain/pNIIC_6C_NIJ/bry_file.nc')

u_west   ,u_east    = bnry_ds['u_west'][:],    bnry_ds['u_east'][:]
ubar_west,ubar_east = bnry_ds['ubar_west'][:], bnry_ds['ubar_east'][:]
temp_west,temp_east = bnry_ds['temp_west'][:], bnry_ds['temp_east'][:]

###############################################################################
################## - Apply boundary fields in nudging zone - ##################
###############################################################################

ds_u_4d[:,:,:,:nud]         = u_west[:,:,:,np.newaxis].repeat(nud,3)
ds_u_4d[:,:,:,(y-nud-1):]   = u_east[:,:,:,np.newaxis].repeat(nud,3)

ds_ubar_3d[:,:,:nud]        = ubar_west[:,:,np.newaxis].repeat(nud,2)
ds_ubar_3d[:,:,(y-nud-1):]  = ubar_east[:,:,np.newaxis].repeat(nud,2)

temp_from_init[:,:,:,:nud]      = temp_west[:,:,:,np.newaxis].repeat(nud,3)
temp_from_init[:,:,:,(y-nud):]  = temp_east[:,:,:,np.newaxis].repeat(nud,3)


if create_clm_file == 1:
    ###############################################################################
    ################################### - file - ##################################
    ###############################################################################
    
    ncid                = nc.Dataset(clm_file,mode='w',format='NETCDF4')
    ncid.description    = Descrip_grd
    ncid.author         = Author
    
    ###############################################################################
    ############################### - dimensions - ################################
    ###############################################################################
    
    # Meridional: eta
    eta_rho             = ncid.createDimension('eta_rho',  x)
    eta_u               = ncid.createDimension('eta_u',    x)
    eta_v               = ncid.createDimension('eta_v',    x-1)
    eta_psi             = ncid.createDimension('eta_psi',  x-1)
    
    # Zonal: xi
    xi_rho              = ncid.createDimension('xi_rho',   y)
    xi_u                = ncid.createDimension('xi_u',     y-1)
    xi_v                = ncid.createDimension('xi_v',     y)
    xi_psi              = ncid.createDimension('xi_psi',   y-1)
    
    s_rho               = ncid.createDimension('s_rho',    z)
    time                = ncid.createDimension('time',     2)
    one                 = ncid.createDimension('one',      1)
    
    ###############################################################################
    ############################### - Variables - #################################
    ###############################################################################
    
    clm_time              = ncid.createVariable('clm_time','f8',('time'))
    clm_time.units          = 'days'
    clm_time.time           = 'clm_time'
    clm_time.cycle_length   = 360.
    clm_time.field          = 'time, scalar, series'


    salt                  = ncid.createVariable('salt','f8',('time','s_rho','eta_rho','xi_rho',))
    salt.long_name          = 'salinity'
    salt.time='clm_time'
    
    temp                  = ncid.createVariable('temp','f8',('time','s_rho','eta_rho','xi_rho',))
    temp.long_name          = 'Temperature'
    temp.time               = 'clm_time'
    temp.units              = 'Celsius'
    
    
    u                     = ncid.createVariable('u','f8',('time','s_rho','eta_u','xi_u',))
    u.long_name             = '3D U-momentum'
    u.time                  = 'clm_time'
    u.units                 = 'm/s'
    
    v                     = ncid.createVariable('v','f8',('time','s_rho','eta_v','xi_v',))
    v.long_name             = '3D V-momentum'
    v.time                  = 'clm_time'
    v.units                 = 'm/s'
    
    
    ubar                  = ncid.createVariable('ubar','f8',('time','eta_u','xi_u',))
    ubar.long_name          = '2D U-momentum'
    ubar.time               = 'clm_time'
    ubar.units              = 'm/s'
    
    vbar                  = ncid.createVariable('vbar','f8',('time','eta_v','xi_v',))
    vbar.long_name          = '2D v-momentum'
    vbar.time               = 'clm_time'
    vbar.units              = 'm/s'
    
    
    zeta                  = ncid.createVariable('zeta','f8',('time','eta_rho','xi_rho',))
    zeta.long_name          = 'free surface'
    zeta.time               = 'clm_time'
    zeta.units              = 'm'
    
    ###############################################################################
    ############################## - Assign values - ##############################
    ###############################################################################
    
    clm_time[:]             = np.array([0, 180])
    
    salt[:,:,:,:]           = salt_from_init
    temp[:,:,:,:]           = temp_from_init
    
    u[:,:,:,:]              = ds_u_4d
    v[:,:,:,:]              = ds_v_4d
    
    ubar[:,:,:]             = ds_ubar_3d
    vbar[:,:,:]             = ds_vbar_3d
    
    zeta[:,:,:]             = ds_zeta_3d
    
    ncid.close()
    
