#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 14 15:12:03 2023

@author: thorbjorn
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
from GlobalParameters import *                # This is where I change user outputs now

print('###############################################################################')
print('############################# - Making clm file - #############################')
print('###############################################################################')

###############################################################################
############################### - User options - ##############################
###############################################################################

#File directory/names
common_dir          = GlobalParameters.common_dir
directory           = GlobalParameters.exp_dir

clm_file            = directory + GlobalParameters.clm_file
grid_file           = directory + GlobalParameters.grid_file
Input_file          = directory + GlobalParameters.init_file
if GlobalParameters.Use_stretch == 1:
    stretching_file = directory + GlobalParameters.stretching
else:
    stretching_file = common_dir + GlobalParameters.stretching_alt

grid_ds             = nc.Dataset(grid_file)
init_ds             = nc.Dataset(Input_file)
stre_ds             = nc.Dataset(stretching_file)

Descrip_grd         = GlobalParameters.General_desc + GlobalParameters.spec_desc_clm
Author              = GlobalParameters.Author
create_clm_file     = GlobalParameters.Create_files
plotting            = GlobalParameters.Gen_plot
Iplottin            = GlobalParameters.Imp_plot
z                   = GlobalParameters.z
res                 = GlobalParameters.Res

#Experiments in "idealized modeling of the North Icelandic Jet" master thesis
pNIJ            = GlobalParameters.pNIJ
pNIIC           = GlobalParameters.pNIIC
pNIIC_hydro     = GlobalParameters.pNIIC_hydro
pNIIC2          = GlobalParameters.pNIIC2
pNIJ_md         = GlobalParameters.pNIJ_md
pIFSJ           = GlobalParameters.pIFSJ
onslope         = GlobalParameters.onslope
NIIC_EB         = GlobalParameters.NIIC_EB

#Dims
x               = np.shape(grid_ds['h'])[0]                             # dim meridional ~ 400 km   (1 km resolution)  
y               = np.shape(grid_ds['h'])[1]                             # dim zonal      ~ 650 km  

#nudge zone 
if res   == 1.5:
    nud         = 25
elif res == 1:
    nud         = 17   
elif res == 3:
    nud         = 8

spe             = 0.10

###############################################################################
############################### - Functions - #################################
###############################################################################

def double_trouble(old_var):             # add arbitrary dimension (dim = 2)
    return np.array([old_var,old_var])
    
###############################################################################
###################### - load temp/salt from init file - ######################
###############################################################################

salt_from_init = double_trouble(init_ds['salt'][0,:,:,:])
temp_from_init = double_trouble(init_ds['temp'][0,:,:,:])

###############################################################################
########## - Empty arrays for velocities and sea surface elevation - ##########
###############################################################################

ds_u_3d         = np.zeros((z, x, y-1))
ds_u_4d         = np.zeros((2, z, x, y-1))
ds_ubar_3d      = np.zeros((2, x, y-1))

ds_v_3d         = np.zeros((z, x-1, y))
ds_v_4d         = np.zeros((2, z, x-1, y))
ds_vbar_3d      = np.zeros((2, x-1, y))

ds_zeta_3d      = np.shape((2, x, y))

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
    if GlobalParameters.steepnes_factor == 7:
        if res == 3:
            ds_u_4d[:,0:13,50:59,0:nud] = double_trouble(np.ones([13,9,nud])*-0.0955)
        elif res == 1.5:
            if GlobalParameters.Use_stretch == 1:
                if GlobalParameters.theta_b == 2:
                    ds_u_4d[:,0:13,100:116,0:nud] = double_trouble(np.ones([13,16,nud])*-0.1025)
                elif GlobalParameters.theta_b == 14:
                    ds_u_4d[:,0:13,100:118,0:nud] = double_trouble(np.ones([13,18,nud])*-0.102)
            else:
                if GlobalParameters.New_S_exag == 0:
                    ds_u_4d[:,0:13,100:119,0:nud] = double_trouble(np.ones([13,19,nud])*-0.10)
                elif GlobalParameters.New_S_exag == 1:
                    ds_u_4d[:,0:15,100:120,0:nud] = double_trouble(np.ones([15,20,nud])*-0.101)
                elif GlobalParameters.New_S_exag == 2:
                    ds_u_4d[:,0:14,100:119,0:nud] = double_trouble(np.ones([14,19,nud])*-0.101)
        elif res == 1:
            if GlobalParameters.Use_stretch == 1:
                if GlobalParameters.theta_b == 2:
                    ds_u_4d[:,0:13,150:174,0:nud] = double_trouble(np.ones([13,24,nud])*-spe)
                elif GlobalParameters.theta_b == 14:
                    ds_u_4d[:,0:13,150:177,0:nud] = double_trouble(np.ones([13,27,nud])*-spe)
            else:
                if GlobalParameters.New_S_exag == 0:
                    ds_u_4d[:,0:13,150:178,0:nud] = double_trouble(np.ones([13,28,nud])*-0.105)
                elif GlobalParameters.New_S_exag == 1:
                    ds_u_4d[:,0:15,150:179,0:nud] = double_trouble(np.ones([15,29,nud])*-0.103)
                elif GlobalParameters.New_S_exag == 2:
                    ds_u_4d[:,0:15,150:179,0:nud] = double_trouble(np.ones([15,29,nud])*-0.103)
    elif GlobalParameters.steepnes_factor == 15:
        if res == 3:
            ds_u_4d[:,0:13,50:59,0:nud] = double_trouble(np.ones([13,9,nud])*-0.0955)
        elif res == 1.5:
            if GlobalParameters.Use_stretch == 1:
                if GlobalParameters.theta_b == 2:
                    ds_u_4d[:,0:13,100:116,0:nud] = double_trouble(np.ones([13,16,nud])*-0.1025)
                elif GlobalParameters.theta_b == 14:
                    ds_u_4d[:,0:13,100:118,0:nud] = double_trouble(np.ones([13,18,nud])*-0.102)
            else:
                if GlobalParameters.New_S_exag == 0:
                    ds_u_4d[:,0:14,100:119,0:nud] = double_trouble(np.ones([14,19,nud])*-0.101)
                elif GlobalParameters.New_S_exag == 1:
                    ds_u_4d[:,0:15,100:120,0:nud] = double_trouble(np.ones([15,20,nud])*-0.101)
                elif GlobalParameters.New_S_exag == 2:
                    ds_u_4d[:,0:14,100:119,0:nud] = double_trouble(np.ones([14,19,nud])*-0.101)
        elif res == 1:
            if GlobalParameters.Use_stretch == 1:
                if GlobalParameters.theta_b == 2:
                    ds_u_4d[:,0:13,150:174,0:nud] = double_trouble(np.ones([13,24,nud])*-spe)
                elif GlobalParameters.theta_b == 14:
                    ds_u_4d[:,0:13,150:177,0:nud] = double_trouble(np.ones([13,27,nud])*-spe)
            else:
                if GlobalParameters.New_S_exag == 0:
                    ds_u_4d[:,0:13,150:178,0:nud] = double_trouble(np.ones([13,28,nud])*-0.105)
                elif GlobalParameters.New_S_exag == 1:
                    ds_u_4d[:,0:15,150:179,0:nud] = double_trouble(np.ones([15,29,nud])*-0.103)
                elif GlobalParameters.New_S_exag == 2:
                    ds_u_4d[:,0:15,150:179,0:nud] = double_trouble(np.ones([15,29,nud])*-0.103)
        
    ds_u_4d_UBARCALC = ds_u_4d.copy()
    
    if Iplottin == 1:
        plt.figure()
        plt.suptitle('FROM Make_clm.py')
        plt.title('horizontal extent of "box" in clm file (at 5th sigma level)')
        X,Y = np.meshgrid(np.linspace(0,(y-1)*res,y-1),np.linspace(0,x*res,x))
        plt.contourf(X,Y,ds_u_4d[0,5,:,:])
        cbarax = plt.colorbar()
        cbarax.set_label('velocity [m/s]')
        plt.xlabel('distance east [km]')
        plt.ylabel('distance north [km]')  


if pNIIC == 1:
    ds_u_4d[:,0:z,50:95,0:nud] = double_trouble(np.ones([z,45,nud])*spe)
    ds_u_4d_UBARCALC = ds_u_4d
    
    if Iplottin == 1:
        plt.figure()
        plt.suptitle('FROM Make_clm.py')
        plt.title('horizontal extent of "box" in clm file (at 5th sigma level)')
        X,Y = np.meshgrid(np.linspace(0,y*res,y),np.linspace(0,x*res,x))
        plt.contourf(X,Y,ds_u_4d[0,5,:,:])
        cbarax = plt.colorbar()
        cbarax.set_label('velocity [m/s]')
        plt.xlabel('distance east [km]')
        plt.ylabel('distance north [km]')
        
if pNIIC_hydro == 1:
    z1,z2,y1,y2,vel,temp,salt = 0,30,70,142, 0.116, 3 , 35.2
    
    # Velocity
    ds_u_4d[:,z1:z2,y1:y2,0:nud] = double_trouble(np.ones([z2-z1,y2-y1,nud])*vel)
    ds_u_4d_UBARCALC = ds_u_4d
    
    z1,z2,y1,y2,vel,temp,salt = 10,30,70,135, 0.107, 3 , 35.2
    
    # temperataure
    temp_from_init[:,z1:z2,y1:y2,0:5] = double_trouble(np.ones([z2-z1,y2-y1,5])*temp)
    
    # salinity
    salt_from_init[:,z1:z2,y1:y2,0:5] = double_trouble(np.ones([z2-z1,y2-y1,5])*salt)
    
    if Iplottin == 1:
        
        fig = plt.figure(figsize=[7,10])
        plt.suptitle('FROM Make_clm.py',y=0.94)
        grid = plt.GridSpec(3,1,hspace=0.4)
        
        ax = fig.add_subplot(grid[0, 0])
        
        ax.set_title('Velocity at sigma lvl 11')
        X,Y = np.meshgrid(np.linspace(0,y*res,y-1),np.linspace(0,x*res,x))
        colorplot = ax.contourf(X,Y,ds_u_4d[0,11,:,:])
        cbarax = plt.colorbar(colorplot,ax=ax)
        cbarax.set_label('velocity [m/s]')
        ax.set_ylabel('distance north [km]')
        
        ax = fig.add_subplot(grid[1, 0])
        ax.set_title('Temp at sigma lvl 11')
        X,Y = np.meshgrid(np.linspace(0,y*res,y),np.linspace(0,x*res,x))
        colorplot = ax.contourf(X,Y,temp_from_init[0,11,:,:])
        cbarax = plt.colorbar(colorplot,ax=ax)
        cbarax.set_label('temp')
        ax.set_ylabel('distance north [km]')
        
        ax = fig.add_subplot(grid[2, 0])
        ax.set_title('Salinity at sigma lvl 11')
        X,Y = np.meshgrid(np.linspace(0,y*res,y),np.linspace(0,x*res,x))
        colorplot = ax.contourf(X,Y,salt_from_init[0,11,:,:])
        cbarax = plt.colorbar(colorplot,ax=ax)
        cbarax.set_label('salt')
        ax.set_xlabel('distance east [km]')
        ax.set_ylabel('distance north [km]')
    
    
if pNIIC2 == 1:
    ds_u_4d[:,0:z,40:85,(y-nud):] = double_trouble(np.ones([z,45,nud])*spe)
    ds_u_4d_UBARCALC = ds_u_4d
    
    if Iplottin == 1:
        plt.figure()
        plt.suptitle('FROM Make_clm.py')
        plt.title('horizontal extent of "box" in clm file (at 5th sigma level)')
        X,Y = np.meshgrid(np.linspace(0,y*res,y),np.linspace(0,x*res,x))
        plt.contourf(X,Y,ds_u_4d[0,5,:,:])
        cbarax = plt.colorbar()
        cbarax.set_label('velocity [m/s]')
        plt.xlabel('distance east [km]')
        plt.ylabel('distance north [km]')
  
if pNIJ_md == 1:
    ds_u_4d[:,7:15,100:128,0:nud] = double_trouble(np.ones([8,28,nud])*-spe)
    ds_u_4d_UBARCALC = ds_u_4d
    
    if Iplottin == 1:
        plt.figure()
        plt.suptitle('FROM Make_clm.py')
        plt.title('horizontal extent of "box" in clm file (at 5th sigma level)')
        X,Y = np.meshgrid(np.linspace(0,y*res,y),np.linspace(0,x*res,x))
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
        
    if Iplottin == 1:
        plt.figure()
        plt.suptitle('FROM Make_clm.py')
        plt.title('horizontal extent of "box" in clm file (at 5th sigma level)')
        X,Y = np.meshgrid(np.linspace(0,y*res,y),np.linspace(0,x*res,x))
        plt.contourf(X,Y,ds_u_4d[0,5,:,:])
        cbarax = plt.colorbar()
        cbarax.set_label('velocity [m/s]')
        plt.xlabel('distance east [km]')
        plt.ylabel('distance north [km]')
        

#2D fields!! - (ubar)

if pNIJ == 1 or pNIIC == 1 or pNIIC_hydro == 1 or pNIJ_md == 1:
    Cs_w = stre_ds['Cs_r'][:]
    diff_Cs_w = np.diff(Cs_w)
    
    for k in range(0,2):
        for j in range(0,nud,1):
            for i in range(0,x,1):
                un = ds_u_4d_UBARCALC[0,:,:,j]
                weighted = un[:-1,i]*diff_Cs_w
                weighted_mean = np.nansum(weighted)
                ds_ubar_3d[k,i,j] = weighted_mean
    if plotting == 1:
        plt.figure()
        plt.suptitle('FROM Make_clm.py')
        plt.title('Ubar at western boundary')
        plt.plot(np.linspace(0,x*res,x),ds_ubar_3d[0,:,0])
                
if pIFSJ == 1:  #could also be an 'if-statement' in last part, but I think this is more clear (on eastern boudnary that is)
    Cs_w = stre_ds['Cs_r'][:]
    diff_Cs_w = np.diff(Cs_w)
    
    for k in range(0,2):
        for j in range(y-nud,y,1):
            for i in range(0,x,1):
                un = ds_u_4d_UBARCALC[0,:,:,j]
                weighted = un[:-1,i]*diff_Cs_w
                weighted_mean = np.nansum(weighted)
                ds_ubar_3d[k,i,j] = weighted_mean
    
    if plotting == 1:
        plt.figure()
        plt.suptitle('FROM Make_clm.py')
        plt.title('Ubar at eastern boundary')
        plt.plot(np.linspace(0,x*res,x),ds_ubar_3d[0,:,-1])
    


###############################################################################
################################ - make file - ################################
###############################################################################

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
    
