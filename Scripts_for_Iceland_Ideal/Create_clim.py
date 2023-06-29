#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 18 16:33:25 2022

@author: thorbjornostenbymoe
"""
""" 
This script creates a climatology file test for ROMS idealized model runs. (basically a copy of input script)

In this test file, all variables except T & S are set as appropriately dimensioned 0 matrices.

The transformation to sigma coordinates are made from the stretching vector (made in stretching.py), the depth coordinates 
from the grid file (made in Create_grid.py), and the averaged CTD values are from the final ctd file (made in CTD_mean.py).

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

###############################################################################
############################### - User options - ##############################
###############################################################################

#File directory/names
directory       = '/Users/thorbjornostenbymoe/Desktop/Initialization_of_model/MODEL_SETUP/' #User directory/path
clm_file        = 'ROMS_in/clm_file.nc'
grid_file       = 'ROMS_in/grid_file.nc'
Input_file      = 'ROMS_in/init_file.nc'
Stretching_file = 'ROMS_in/Stretching.nc'

ds_grid         = nc.Dataset(directory+grid_file)
ds_init         = nc.Dataset(directory+Input_file)
ds_stre         = nc.Dataset(directory+Stretching_file)

#file description
Descrip_grd     = 'Iceland slope - ROMS idealized model run - clm file'  
Author          = 'Thorbjoern Oestenby Moe'

#logical switches
create_clm_file = 1                                       # create input (1 = true)
plotting        = 1

#Experiments in "idealized modeling of the North Icelandic Jet" master thesis
pNIJ            = 1                                       # EXP 1 == prescribed NIJ outflow in west
pNIIC           = 1                                       # EXP 2 == prescribed NIIC inflow in west
pNIIC2          = 0                                       # EXP 2.5 == altered prescribed NIIC inflow in west
pNIJ_md         = 0                                       # EXP 3 == mid-depth confined prescribed NIJ outflow in west
pIFSJ           = 0                                       # EXP 4 == prescribed IFSJ outflow in east
onslope         = 0                                       # Sets degree to which pIFSJ is prescribed on the boundary, very onslope,1=onslope, 0=offslope
NIIC_EB         = 0                                       # EXP 5 == prescribed NIIC outflow in east

#Dims
x    = 268                                                # dim meridional ~ 400 km   (1.5 km resolution)  
y    = 401                                                # dim zonal      ~ 600 km  
z    = 30                                                 # dim vertical   = varying 
  
#nudge zone 
nud  = 17                                                 # NB! Has to be the same as in 'Create_clim_nudge.py'!!! 
#speed of prescribed currents
spe  = 0.10

###############################################################################
############################### - Functions - #################################
###############################################################################

def double_trouble(old_var):             # add arbitrary dimension (dim = 2)
    return np.array([old_var,old_var])
    
###############################################################################
###################### - load temp/salt from init file - ######################
###############################################################################

salt_from_init = double_trouble(ds_init['salt'][0,:,:,:])
temp_from_init = double_trouble(ds_init['temp'][0,:,:,:])

###############################################################################
########## - Empty arrays for velocities and sea surface elevation - ##########
###############################################################################

ds_u_3d         = np.zeros((z, x, y))
ds_u_4d         = np.zeros((2, z, x, y))
ds_ubar_3d      = np.zeros((2, x, y))

ds_v_3d         = np.zeros((z, x-1, y+1))
ds_v_4d         = np.zeros((2, z, x-1, y+1))
ds_vbar_3d      = np.zeros((2, x-1, y+1))

ds_zeta_3d      = np.shape((2, x, y+1))

###############################################################################
############################ - prescribe currents - ###########################
###############################################################################
"""
From now on, the code goes into specific configurations for previous configurations of
the model, skip/customize/adapt if set-up is different

The code under simply defines a box of positive/negative velocities on the boundary (and neighbouring nudging region), 
meant to prescribe idealized versions of currents

The extents are identical to the boundary file (Create_bry.py)


-- Should be expanded with physically consistent sea surface elevation... --
"""

#3D fields!!

if pNIJ == 1:
    ds_u_4d[:,0:13,100:118,0:nud] = double_trouble(np.ones([13,18,nud])*-spe)
    ds_u_4d_UBARCALC = ds_u_4d
    
    if plotting == 1:
        plt.figure()
        plt.title('horizontal extent of "box" in clm file (at 5th sigma level)')
        X,Y = np.meshgrid(np.linspace(0,y*1.5,y),np.linspace(0,x*1.5,x))
        plt.contourf(X,Y,ds_u_4d[0,5,:,:])
        cbarax = plt.colorbar()
        cbarax.set_label('velocity [m/s]')
        plt.xlabel('distance east [km]')
        plt.ylabel('distance north [km]')  

if pNIIC == 1:
    ds_u_4d[:,0:z,50:95,0:nud] = double_trouble(np.ones([z,45,nud])*spe)
    ds_u_4d_UBARCALC = ds_u_4d
    
    if plotting == 1:
        plt.figure()
        plt.title('horizontal extent of "box" in clm file (at 5th sigma level)')
        X,Y = np.meshgrid(np.linspace(0,y*1.5,y),np.linspace(0,x*1.5,x))
        plt.contourf(X,Y,ds_u_4d[0,5,:,:])
        cbarax = plt.colorbar()
        cbarax.set_label('velocity [m/s]')
        plt.xlabel('distance east [km]')
        plt.ylabel('distance north [km]')
    
if pNIIC2 == 1:
    ds_u_4d[:,0:z,40:85,(y-nud):] = double_trouble(np.ones([z,45,nud])*spe)
    ds_u_4d_UBARCALC = ds_u_4d
    
    if plotting == 1:
        plt.figure()
        plt.title('horizontal extent of "box" in clm file (at 5th sigma level)')
        X,Y = np.meshgrid(np.linspace(0,y*1.5,y),np.linspace(0,x*1.5,x))
        plt.contourf(X,Y,ds_u_4d[0,5,:,:])
        cbarax = plt.colorbar()
        cbarax.set_label('velocity [m/s]')
        plt.xlabel('distance east [km]')
        plt.ylabel('distance north [km]')
  
if pNIJ_md == 1:
    ds_u_4d[:,7:15,100:128,0:nud] = double_trouble(np.ones([8,28,nud])*-spe)
    ds_u_4d_UBARCALC = ds_u_4d
    
    if plotting == 1:
        plt.figure()
        plt.title('horizontal extent of "box" in clm file (at 5th sigma level)')
        X,Y = np.meshgrid(np.linspace(0,y*1.5,y),np.linspace(0,x*1.5,x))
        plt.contourf(X,Y,ds_u_4d[0,5,:,:])
        cbarax = plt.colorbar()
        cbarax.set_label('velocity [m/s]')
        plt.xlabel('distance east [km]')
        plt.ylabel('distance north [km]')  
        
if pIFSJ == 1:
    if onslope == 1:
        ds_u_4d[:,0:8,105:120,(y-nud):] = double_trouble(np.ones([8,15,nud])*spe)
        ds_u_4d_UBARCALC = ds_u_4d
    
    elif onslope == 2:
        ds_u_4d[:,0:12,99:118,(y-nud):] = double_trouble(np.ones([12,118-99,nud])*spe)  
        ds_u_4d_UBARCALC = ds_u_4d
    
    else:
        ds_u_4d[:,0:8,115:128,(y-nud):] = double_trouble(np.ones([8,128-115,nud])*spe)
        ds_u_4d_UBARCALC = ds_u_4d
        
    if plotting == 1:
        plt.figure()
        plt.title('horizontal extent of "box" in clm file (at 5th sigma level)')
        X,Y = np.meshgrid(np.linspace(0,y*1.5,y),np.linspace(0,x*1.5,x))
        plt.contourf(X,Y,ds_u_4d[0,5,:,:])
        cbarax = plt.colorbar()
        cbarax.set_label('velocity [m/s]')
        plt.xlabel('distance east [km]')
        plt.ylabel('distance north [km]')
        

#2D fields!! - (ubar)

if pNIJ == 1 or pNIIC == 1 or pNIJ_md == 1:
    Cs_w = ds_stre['Cs_r'][:]
    diff_Cs_w = np.diff(Cs_w)
    
    for k in range(0,2):
        for j in range(0,nud,1):
            for i in range(0,x,1):
                un = ds_u_4d_UBARCALC[0,:,:,j]
                weighted = un[:-1,i]*diff_Cs_w
                weighted_mean = np.nansum(weighted)
                ds_ubar_3d[k,i,j] = weighted_mean
    plt.figure()
    plt.title('Ubar at western boundary')
    plt.plot(ds_ubar_3d[0,:,0])
                
if pIFSJ == 1:  #could also be an 'if-statement' in last part, but I think this is more clear (on eastern boudnary that is)
    Cs_w = ds_stre['Cs_r'][:]
    diff_Cs_w = np.diff(Cs_w)
    
    for k in range(0,2):
        for j in range(y-nud,y,1):
            for i in range(0,x,1):
                un = ds_u_4d_UBARCALC[0,:,:,j]
                weighted = un[:-1,i]*diff_Cs_w
                weighted_mean = np.nansum(weighted)
                ds_ubar_3d[k,i,j] = weighted_mean
    
    plt.figure()
    plt.title('Ubar at eastern boundary')
    plt.plot(ds_ubar_3d[0,:,-1])
    


###############################################################################
################################ - make file - ################################
###############################################################################

if create_clm_file == 1:
    ###############################################################################
    ################################### - file - ##################################
    ###############################################################################
    
    path_string         = directory+clm_file
    ncid                = nc.Dataset(path_string,mode='w',format='NETCDF4')
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
    xi_rho              = ncid.createDimension('xi_rho',   y+1)
    xi_u                = ncid.createDimension('xi_u',     y)
    xi_v                = ncid.createDimension('xi_v',     y+1)
    xi_psi              = ncid.createDimension('xi_psi',   y)
    
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
    
