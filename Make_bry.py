#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 14 15:11:13 2023

@author: thorbjorn
"""

""" 
This script creates a boundary file for idealized model runs using ROMS. 


This script currently simply uses the input fields as the basis for the hydrographic 
boundary conditions. This can be changed by simply setting the desired variables equal to some other,
predetermined value (i.e. temp_east[:,:] = temp_from_file or otherwise)

The script currently only specifies boundary conditions for the western and eastern edges, as the 
meridional boundaries are walls.

The writing to netCDF is made in accordance to the ROMS standards (e.g. name of variables)
Completely different set-ups can still use a large part of this script as a template
"""

###############################################################################
################################ - Packages - #################################
###############################################################################

import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy.matlib as repmaty
from GlobParam import *                # This is where I change user outputs now
from scipy.ndimage import gaussian_filter

print('###############################################################################')
print('############################# - Making bry file - #############################')
print('###############################################################################')

###############################################################################
################################ - Functions - ################################
###############################################################################

def make_XY(ds,zd,xd,vs=False):                                    # get coordinates for plotting
    stretch = stretch_ds['Cs_r'][:]                       # stretching
    heighta = ds['h'][:,:]                                # depth array
    
    bath = heighta[:,:]*stretch[:,np.newaxis,np.newaxis]
    X = repmaty.repmat(np.linspace(0,x,xd),zd,1)

    return(bath,X)

def get_sizes(mnz,mxz,mnh,mxh,velo): # function for making "box" on boundary == outflow/inflow
    z   = mxz - mnz
    h   = mxh - mnh
    box = np.ones([2,z,h])*velo
    return(mnz,mxz,mnh,mxh,box)


###############################################################################
############################### - User options - ##############################
###############################################################################

#File directory/names
common_dir          = GlobParam.common_dir
directory           = GlobParam.EXP_dir

bry_file            = directory + GlobParam.bry_file
grid_file           = directory + GlobParam.grid_file
Input_file          = directory + GlobParam.init_file
stretching_file     = common_dir + GlobParam.Streching

grid_ds             = nc.Dataset(grid_file)
init_ds             = nc.Dataset(Input_file)
stretch_ds          = nc.Dataset(stretching_file)

Descrip_grd         = GlobParam.General_desc + GlobParam.spec_desc_bry
Author              = GlobParam.Author
create_bry_file     = GlobParam.Create_files
plotting            = GlobParam.Gen_plot
Iplottin            = GlobParam.Imp_plot
z                   = GlobParam.z
res                 = GlobParam.Res

#Experiments in "idealized modeling of the North Icelandic Jet" master thesis
pNIJ            = GlobParam.pNIJ
pNIIC           = GlobParam.pNIIC
pNIIC_f         = GlobParam.pNIIC_f
pNIIC_w         = GlobParam.pNIIC_w
pIFSJ           = GlobParam.pIFSJ

InflowTemp      = GlobParam.temp_of_wNIIC

#Dims
x   = np.shape(grid_ds['h'])[0]                           # dim meridional ~ 400 km   (1 km resolution)  
y   = np.shape(grid_ds['h'])[1]                           # dim zonal      ~ 650 km 

#%%

###############################################################################
############################ - load input fields - ############################
###############################################################################

salt    = init_ds['salt']                                      # salinity
temp    = init_ds['temp']                                      # temperature
u       = init_ds['u']                                         # zonal velocity
v       = init_ds['v']                                         # meridional velocity
ubar    = init_ds['ubar']                                      # depth integrated u
vbar    = init_ds['vbar']                                      # depth integrated v
zeta    = init_ds['zeta']                                      # sea surface elevation

###############################################################################
################################ - load grid - ################################
###############################################################################

lon_rho = grid_ds['lon_rho']                                  # coordinates
lat_rho = grid_ds['lat_rho']                                  # ...
lon_u   = grid_ds['lon_u']                                    #
lat_u   = grid_ds['lat_u']                                    #
lon_v   = grid_ds['lon_v']                                    #
lat_v   = grid_ds['lat_v']                                    #

###############################################################################
########### - Strictly define variables according to input fields - ###########
###############################################################################

eas = -1                                       # coordinate of eastern boundary
wes = 1                                        # coordinate of western boundary

bry_time_var = np.array([0, 180])                    # arbitrary time dimension

"""
Here we start defining the boundary fields directly from the input fields.

One variable is created for each boundary (east/west in this case, can be 
expanded to north/south if these are open)

Modifications to these can be made underway or, as done here, afterwards. 
Depends on how complex the experiments are

"""

salt_east_var    = salt[:,:,:,eas].repeat(2,0)
salt_west_var    = salt[:,:,:,wes].repeat(2,0)

temp_east_var    = temp[:,:,:,eas].repeat(2,0)
temp_west_var    = temp[:,:,:,wes].repeat(2,0)


u_east_var       = u[:,:,:,eas].repeat(2,0)
u_west_var       = u[:,:,:,wes].repeat(2,0)

v_east_var       = v[:,:,:,eas].repeat(2,0)
v_west_var       = v[:,:,:,wes].repeat(2,0)


ubar_east_var    = ubar[:,:,eas].repeat(2,0)
ubar_west_var    = ubar[:,:,wes].repeat(2,0)

vbar_east_var    = vbar[:,:,eas].repeat(2,0)
vbar_west_var    = vbar[:,:,wes].repeat(2,0)


zeta_east_var    = zeta[:,:,eas].repeat(2,0)
zeta_west_var    = zeta[:,:,wes].repeat(2,0)


lon_rho_east_var = lon_rho[:,eas]
lon_rho_west_var = lon_rho[:,wes]

lat_rho_east_var = lat_rho[:,eas]
lat_rho_west_var = lat_rho[:,wes]


lon_u_east_var   = lon_u[:,eas]
lon_u_west_var   = lon_u[:,wes]

lat_u_east_var   = lat_u[:,eas]
lat_u_west_var   = lat_u[:,wes]


lon_v_east_var   = lon_v[:,eas]
lon_v_west_var   = lon_v[:,wes]

lat_v_east_var   = lat_v[:,eas]
lat_v_west_var   = lat_v[:,wes]


###############################################################################
############################ - prescribe currents - ###########################
###############################################################################
"""
From now on, the code goes into specific configurations for previous configurations of
the model, skip/customize/adapt if set-up is different

The code under simply defines a box of positive/negative velocities on the boundary, 
meant to prescribe idealized versions of currents
"""

u_west_var_for_ubar_calc = np.zeros([2,30,x])
u_east_var_for_ubar_calc = np.zeros([2,30,x])

    
#prescribing the current(s) - mutliple are allowed.
""" 
For now I have simply ran this script mulitple times to get the balance between 
velocity, area, and volume transport. If this was figured out in a better way 
these options could be moved to global params

"""

if pNIJ == 1:
    # STEEP
    if GlobParam.steepnes_factor == 7 and GlobParam.VarSlope == 0:
        z1,z2,h1,h2,box = get_sizes(0,15,150,176,-0.101) #(0,15,150,176,-0.101) #(0,15,150,168,-0.103) = 1.2
    
    # MODERATE or VarSlo
    elif GlobParam.steepnes_factor == 15 or GlobParam.VarSlope == 1:
        z1,z2,h1,h2,box = get_sizes(0,15,150,178,-0.104)

    # Apply BC
    u_west_var[:,z1:z2,h1:h2]               = box
    u_west_var_for_ubar_calc[:,z1:z2,h1:h2] = box


        
if pNIIC == 1:
    # STEEP
    if GlobParam.steepnes_factor == 7 and GlobParam.VarSlope == 0:
        z1,z2,h1,h2,box = get_sizes(0,30,75,144,0.103) 
        
    # MODERATE or VarSlo
    elif GlobParam.steepnes_factor == 15 or GlobParam.VarSlope == 1:
        z1,z2,h1,h2,box = get_sizes(0,30,82,144,0.10) 
    
    # Apply BC
    u_west_var[:,z1:z2,h1:h2]               = box
    u_west_var_for_ubar_calc[:,z1:z2,h1:h2] = box

    
if pNIIC_f == 1 and pNIIC_w == 0:
    # STEEP
    if GlobParam.steepnes_factor == 7 and GlobParam.VarSlope == 0:
        # 2xInflow
        if GlobParam.speed_of_wNIIC == '2':
            z1,z2,h1,h2,box = get_sizes(0,30,110,146,0.198)
        # 3xInflow
        elif GlobParam.speed_of_wNIIC == '3':
            z1,z2,h1,h2,box = get_sizes(0,30,129,148,0.305)
        else:
            print('pNIIC_f not applied due to lacking speed specification. Look into "Make_bry.py" and make appropriate changes :/')
    
    # MODERATE or VarSlo
    elif GlobParam.steepnes_factor == 15 or GlobParam.VarSlope == 1:
        # 2xInflow
        if GlobParam.speed_of_wNIIC == '2':
            z1,z2,h1,h2,box = get_sizes(0,30,124,146,0.2) 
        # 3xInflow
        elif GlobParam.speed_of_wNIIC == '3':
            z1,z2,h1,h2,box = get_sizes(0,30,129,148,0.305)
        else:
            print('pNIIC_f not applied due to lacking speed specification. Look into "Make_bry.py" and make appropriate changes :/')
            
    # Apply BC
    u_west_var[:,z1:z2,h1:h2]               = box
    u_west_var_for_ubar_calc[:,z1:z2,h1:h2] = box
    
    
if pNIIC_w == 1 and pNIIC_f == 0:
    # STEEP
    if GlobParam.steepnes_factor == 7:
        z1,z2,y1,y2,vel = 0,30,115,148,0.11
        z1,z2,h1,h2,box_u = get_sizes(z1,z2,y1,y2,vel) 
        z1,z2,h1,h2,box_t = get_sizes(z1,z2,y1,y2,temp) 
        
    # MODERATE
    elif GlobParam.steepnes_factor == 15:
        z1,z2,y1,y2,vel = 0,30,90,145,0.11
        z1,z2,h1,h2,box_u = get_sizes(z1,z2,y1,y2,vel) 
        z1,z2,h1,h2,box_t = get_sizes(z1,z2,y1,y2,InflowTemp) 
    
    # Apply BC
    u_west_var[:,z1:z2,h1:h2]               = box_u     # Velocity
    u_west_var_for_ubar_calc[:,z1:z2,h1:h2] = box_u

    temp_west_var[:,z1:z2,h1:h2]            = box_t     # temperataure


if pNIIC_f == 1 and pNIIC_w == 1:
    # STEEP
    if GlobParam.steepnes_factor == 7 and GlobParam.VarSlope == 0:
        # 2xInflow
        if GlobParam.speed_of_wNIIC == '2':
            z1,z2,y1,y2,vel = 0,30,110,146,0.198
            
            z1,z2,h1,h2,box_u = get_sizes(z1,z2,y1,y2,vel)
            z1,z2,h1,h2,box_t = get_sizes(z1,z2,y1,y2,InflowTemp)
        # 3xInflow
        elif GlobParam.speed_of_wNIIC == '3':
            z1,z2,y1,y2,vel = 0,30,129,148,0.305
            
            z1,z2,h1,h2,box_u = get_sizes(z1,z2,y1,y2,vel)
            z1,z2,h1,h2,box_t = get_sizes(z1,z2,y1,y2,InflowTemp)
            
        else:
            print('pNIIC_f not applied due to lacking speed specification. Look into "Make_bry.py" and make appropriate changes :/')
    
    # MODERATE or VarSlo
    elif GlobParam.steepnes_factor == 15 or GlobParam.VarSlope == 1:
        # 2xInflow
        if GlobParam.speed_of_wNIIC == '2':
            z1,z2,y1,y2,vel = 0,30,124,146,0.2
            
            z1,z2,h1,h2,box_u = get_sizes(z1,z2,y1,y2,vel)
            z1,z2,h1,h2,box_t = get_sizes(z1,z2,y1,y2,InflowTemp)
        # 3xInflow
        elif GlobParam.speed_of_wNIIC == '3':
            z1,z2,y1,y2,vel = 0,30,129,148,0.305
            
            z1,z2,h1,h2,box_u = get_sizes(z1,z2,y1,y2,vel)
            z1,z2,h1,h2,box_t = get_sizes(z1,z2,y1,y2,InflowTemp)
        else:
            print('pNIIC_f not applied due to lacking speed specification. Look into "Make_bry.py" and make appropriate changes :/')
            
    # Apply BC
    u_west_var[:,z1:z2,h1:h2]               = box_u     # Velocity
    u_west_var_for_ubar_calc[:,z1:z2,h1:h2] = box_u

    temp_west_var[:,z1:z2,h1:h2]            = box_t     # temperataure


if pIFSJ == 1:
    #z1,z2,h1,h2,box = get_sizes(0,8,105,120,0.095)
    #z1,z2,h1,h2,box = get_sizes(0,12,99,118,0.104)
    z1,z2,h1,h2,box = get_sizes(0,9,152,175,0.13)
   
    u_east_var[:,z1:z2,h1:h2]               = np.array([box,box])
    u_east_var_for_ubar_calc[:,z1:z2,h1:h2] = np.array([box,box])
    

""" FROM OLD TEST ON PRESCRIBING THE NIJ THROUGH THE EASTERN BOUNDARY """
#if pNIJ_East == 1:
#    z1,z2,h1,h2,box = get_sizes(0,15,150,176,-0.101)
#    
#    u_east_var[:,z1:z2,h1:h2]               = np.array([box,box])
#    u_east_var_for_ubar_calc[:,z1:z2,h1:h2] = np.array([box,box])

if Iplottin == 1:
###############################################################################
########################## - plot volume transport - ##########################
###############################################################################
    
    if pNIJ == 1 or pNIIC == 1 or pNIIC_w == 1 or pNIIC_f == 1  or pIFSJ == 1:
        
        # make area
        h = grid_ds['h'][:,:]
        c  = stretch_ds['Cs_r'][:]
        c2  = stretch_ds['Cs_w'][:]
        
        Y = h[np.newaxis,:,0]*c[:,np.newaxis]
        X = repmaty.repmat(np.linspace(0,401,401),30,1)

        # Volume transport
        area = np.diff(h[np.newaxis,:,0]*c2[:,np.newaxis],axis=0)*1000
        
        
        
        u_W_spec = u_west_var[0,:,:]
        u_E_spec = u_west_var[0,:,:]
        u_E_EB   = u_east_var[0,:,:]

        
        #westward and eastward seperately
        u_W_spec = np.where(u_W_spec>0,u_W_spec,np.zeros(np.shape(u_W_spec)))
        u_E_spec = np.where(u_E_spec<0,u_E_spec,np.zeros(np.shape(u_E_spec)))
        
        # calc Volume Transport (Cross sectional area * current speed)
        vol_trans_W = area * u_W_spec
        vol_trans_E = area * u_E_spec
        vol_trans_EB= area * u_E_EB

        X1,Y1,vol_trans1 = X[:,:],Y[:,:],vol_trans_W[:,:]
        
        #plot 
        if pIFSJ == 1: 
            fig = plt.figure(figsize=[14,8])
            plt.suptitle('FROM Make_bry.py')
            grid = plt.GridSpec(3,3,hspace=0.5)
        else:
            fig = plt.figure(figsize=[14,7])
            plt.suptitle('FROM Make_bry.py')
            grid = plt.GridSpec(3,2,hspace=0.4)
        
        
        ax1 = fig.add_subplot(grid[0, 0])
        ax1.set_title('Eastward (NIIC): \n area of each grid cell')
        z1 = ax1.contourf(X,Y,area)
        plt.colorbar(z1,ax=ax1)
        
        ax1 = fig.add_subplot(grid[1, 0])
        ax1.set_title('Eastward (NIIC): current')
        z1 = ax1.contourf(X1,Y1,u_W_spec[:,:],cmap='Reds')
        plt.colorbar(z1,ax=ax1)
        
        ax3 = fig.add_subplot(grid[2, 0])
        ax3.set_title('volume transport \n total: '+str(round(np.nansum(vol_trans1)/1000000,3))+' Sv')
        z3 = ax3.contourf(X1,Y1,vol_trans1,cmap='Reds')
        plt.colorbar(z3,ax=ax3)
        
        


        ax1 = fig.add_subplot(grid[0, 1])
        ax1.set_title('Westward (NIJ): \narea of each grid cell')
        z1 = ax1.contourf(X,Y,area)
        plt.colorbar(z1,ax=ax1)
        
        ax1 = fig.add_subplot(grid[1, 1])
        ax1.set_title('Westward (NIJ): current')
        z1 = ax1.contourf(X1,Y1,u_E_spec[:,:],cmap='Blues_r')
        plt.colorbar(z1,ax=ax1)

        X1,Y1,vol_trans1 = X[:,:],Y[:,:],vol_trans_E[:,:]
        ax3 = fig.add_subplot(grid[2, 1])
        ax3.set_title('volume transport \n total: '+str(round(np.nansum(vol_trans1)/1000000,3))+' Sv')
        z3 = ax3.contourf(X1,Y1,vol_trans1,cmap='Blues_r')
        plt.colorbar(z3,ax=ax3)
        
        
        if pIFSJ == 1: 
            clr = 'Reds'
            ax1 = fig.add_subplot(grid[0, 2])
            ax1.set_title('Velocity at Eastern boundary: \n area of each grid cell')
            z1 = ax1.contourf(X,Y,area)
            plt.colorbar(z1,ax=ax1)
            
            ax1 = fig.add_subplot(grid[1, 2])
            ax1.set_title('Velocity')
            z1 = ax1.contourf(X1,Y1,u_E_EB[:,:],cmap=clr)
            plt.colorbar(z1,ax=ax1)
            
            X1,Y1,vol_trans1 = X[:,:],Y[:,:],vol_trans_E[:,:]
            
            ax3 = fig.add_subplot(grid[2, 2])
            ax3.set_title('volume transport \n total: '+str(round(np.nansum(vol_trans_EB)/1000000,3))+' Sv')
            z3 = ax3.contourf(X1,Y1,vol_trans_EB,cmap=clr)
            plt.colorbar(z3,ax=ax3)
        
    
###############################################################################
####################### - plot depth averaged velocity - ######################
###############################################################################

if pNIJ == 1 or pNIIC == 1 or pNIIC_w == 1 or pNIIC_w == 1:
    Cs_w = stretch_ds['Cs_r'][:]
    diff_Cs_w= np.diff(Cs_w)
    
    ubarc = np.zeros([x])
    
    for i in range(x):
        un = u_west_var_for_ubar_calc[0,:-1,i]
        weighted_mean = np.nansum(un*diff_Cs_w) / np.nansum(diff_Cs_w)
        ubarc[i] = weighted_mean

    ubar_west_var = np.array([ubarc,ubarc])
    
    if plotting == 1 or Iplottin == 1:
        plt.figure()
        plt.suptitle('FROM Make_bry.py')
        plt.plot(np.linspace(0,x*res,x),ubarc[:])
        plt.title('western Ubar BC from bry file')
        plt.ylabel('velocity [m/s]')
        plt.xlabel('distance north [km]')
    
if pNIIC_w == 1:
    if plotting == 1 or Iplottin == 1:
        z10 = np.zeros([z,x])
        for i in range(x):
            z10[:,i] = deptha*depth_a[i]
        
        pp,X,Y = make_XY(grid_ds,z,x-1)
        X,Y    = X[:,:],Y[:,1:,0]
        
        fig = plt.figure(figsize=[7,6])
        plt.suptitle('FROM Make_bry.py \nonly pNIIC_(NIJ)_hydro')
        grid = plt.GridSpec(3,1,hspace=0.2)
        
        ax = fig.add_subplot(grid[0, 0])
        z1 = ax.contourf(X,Y,u_west_var[0,:,1:])
        plt.colorbar(z1,ax=ax)
        
        ax = fig.add_subplot(grid[1, 0])
        z1 = ax.contourf(X,Y,temp_west_var[0,:,1:],levels=np.arange(-1,6.11,0.1))
        plt.colorbar(z1,ax=ax)
        
        ax = fig.add_subplot(grid[2, 0])
        z1 = ax.contourf(X,Y,salt_west_var[0,:,1:],levels=np.arange(34.74,35.02,0.01))
        plt.colorbar(z1,ax=ax)

if pIFSJ == 1:
    Cs_w = stretch_ds['Cs_w'][:]
    
    diff_Cs_w=np.diff(Cs_w)
    
    ubarc = np.zeros([x])
    for i in range(0,x,1):
        un = u_east_var_for_ubar_calc[0,:,i]
        weighted = un*diff_Cs_w
        weighted_mean = np.nansum(weighted)
        ubarc[i] = weighted_mean
        
    ubar_east_var = np.array([ubarc,ubarc])
        
    if plotting == 1 or Iplottin == 1:
        plt.figure()
        plt.suptitle('FROM Make_bry.py')
        plt.plot(np.linspace(0,x,x),ubar_east_var[0,:])
        plt.title('eastern Ubar BC from bry file')
        plt.ylabel('velocity [m/s]')
        plt.xlabel('distance north [km]')
    

#%%

if create_bry_file == 1:
    
    ###############################################################################
    ################################### - file - ##################################
    ###############################################################################
    
    ncid                = nc.Dataset(bry_file,mode='w',format='NETCDF4')
    ncid.description    = Descrip_grd
    ncid.author         = Author
    
    ###############################################################################
    ############################### - Dimensions - ################################
    ###############################################################################
    
    # Meridional: eta
    eta_rho             = ncid.createDimension('eta_rho',x)
    eta_u               = ncid.createDimension('eta_u',  x)
    eta_v               = ncid.createDimension('eta_v',  x-1)
    
    # Zonal: xi
    xi_rho              = ncid.createDimension('xi_rho', y)
    xi_u                = ncid.createDimension('xi_u',   y-1)
    xi_v                = ncid.createDimension('xi_v',   y)
    
    s_rho               = ncid.createDimension('s_rho',  z)
    s_w                 = ncid.createDimension('s_w',    z+1)
    
    time                = ncid.createDimension('time',   2)
    tracer              = ncid.createDimension('tracer', 2)
    one                 = ncid.createDimension('one',    1)
    
    ###############################################################################
    ############################### - Variables - #################################
    ###############################################################################
    
    bry_time              = ncid.createVariable('bry_time','f8',('time'))
    bry_time.units          = 'days'
    bry_time.long_name      = 'time since initialization'
    bry_time.cycle_length   = 360.
    bry_time.time           = 'bry_time'
    bry_time.field          = 'time, scalar, series'
    
    #######################################
    salt_east             = ncid.createVariable('salt_east','f8',('time','s_rho','eta_rho',))
    salt_east.long_name     = 'salinity_east'
    salt_east.time          = 'bry_time'
    
    salt_west             = ncid.createVariable('salt_west','f8',('time','s_rho','eta_rho',))
    salt_west.long_name     = 'salinity_east'
    salt_west.time          = 'bry_time'
    
    #######################################
    temp_east             = ncid.createVariable('temp_east','f8',('time','s_rho','eta_rho',))
    temp_east.long_name     = 'Temperature'
    temp_east.units         = 'Celsius'
    temp_east.time          = 'bry_time'
    
    temp_west             = ncid.createVariable('temp_west','f8',('time','s_rho','eta_rho',))
    temp_west.long_name     = 'Temperature'
    temp_west.units         = 'Celsius'
    temp_west.time          = 'bry_time'
    
    #######################################
    u_east                = ncid.createVariable('u_east','f8',('time','s_rho','eta_u',))
    u_east.long_name        = '3D U-momentum'
    u_east.units            = 'm/s'
    u_east.time             = 'bry_time'
    
    u_west                = ncid.createVariable('u_west','f8',('time','s_rho','eta_u',))
    u_west.long_name        = '3D U-momentum'
    u_west.units            = 'm/s'
    u_west.time             = 'bry_time'
    
    #######################################
    v_east                = ncid.createVariable('v_east','f8',('time','s_rho','eta_v',))
    v_east.long_name        = '3D V-momentum'
    v_east.units            = 'm/s'
    v_east.time             = 'bry_time'
    
    v_west                = ncid.createVariable('v_west','f8',('time','s_rho','eta_v',))
    v_west.long_name        = '3D V-momentum'
    v_west.units            = 'm/s'
    v_west.time             = 'bry_time'
    
    #######################################
    ubar_east             = ncid.createVariable('ubar_east','f8',('time','eta_u',))
    ubar_east.long_name     = '2D U-momentum'
    ubar_east.units         = 'm/s'
    ubar_east.time          = 'bry_time'
    
    ubar_west             = ncid.createVariable('ubar_west','f8',('time','eta_u',))
    ubar_west.long_name     = '2D U-momentum'
    ubar_west.units         = 'm/s'
    ubar_west.time          = 'bry_time'
    
    #######################################
    vbar_east             = ncid.createVariable('vbar_east','f8',('time','eta_v',))
    vbar_east.long_name     = '2D v-momentum'
    vbar_east.units         = 'm/s'
    vbar_east.time          = 'bry_time'
    
    vbar_west             = ncid.createVariable('vbar_west','f8',('time','eta_v',))
    vbar_west.long_name     = '2D v-momentum'
    vbar_west.units         = 'm/s'
    vbar_west.time          = 'bry_time'
    
    #######################################
    zeta_east             = ncid.createVariable('zeta_east','f8',('time','eta_rho',))
    zeta_east.long_name     = 'free surface'
    zeta_east.time          = 'bry_time'
    
    zeta_west             = ncid.createVariable('zeta_west','f8',('time','eta_rho',))
    zeta_west.long_name     = 'free surface'
    zeta_west.time          = 'bry_time'
    
    #######################################
    lon_rho_west          = ncid.createVariable('lon_rho_west','f8',('eta_rho',))
    lon_rho_west.long_name  = 'longitude of RHO-points, western boundary condition'
    lon_rho_west.units      = 'degree_east'
    
    lat_rho_west          = ncid.createVariable('lat_rho_west','f8',('eta_rho',))
    lat_rho_west.long_name  = 'latitude of RHO-points, western boundary condition'
    lat_rho_west.units      = 'degree_north'
    
    #######################################
    lon_rho_east          = ncid.createVariable('lon_rho_east','f8',('eta_rho',))
    lon_rho_east.long_name  = 'longitude of RHO-points, eastern boundary condition'
    lon_rho_east.units      = 'degree_east'
    
    lat_rho_east          = ncid.createVariable('lat_rho_east','f8',('eta_rho',))
    lat_rho_east.long_name  = 'latitude of RHO-points, eastern boundary condition'
    lat_rho_east.units      = 'degree_north'
    
    
    #######################################
    lon_u_west            = ncid.createVariable('lon_u_west','f8',('eta_u',))
    lon_u_west.long_name    = 'longitude of U-points, western boundary condition'
    lon_u_west.units        = 'degree_east'
    
    lat_u_west            = ncid.createVariable('lat_u_west','f8',('eta_u',))
    lat_u_west.long_name    = 'latitude of U-points, western boundary condition'
    lat_u_west.units        = 'degree_north'
    
    #######################################
    lon_u_east            = ncid.createVariable('lon_u_east','f8',('eta_u',))
    lon_u_east.long_name    = 'longitude of U-points, eastern boundary condition'
    lon_u_east.units        = 'degree_east'
    
    lat_u_east            = ncid.createVariable('lat_u_east','f8',('eta_u',))
    lat_u_east.long_name    = 'latitude of U-points, eastern boundary condition'
    lat_u_east.units        = 'degree_north'
    
    
    #######################################
    lon_v_west            = ncid.createVariable('lon_v_west','f8',('eta_v',))
    lon_v_west.long_name    = 'longitude of V-points, western boundary condition'
    lon_v_west.units        = 'degree_east'
    
    lat_v_west            = ncid.createVariable('lat_v_west','f8',('eta_v',))
    lat_v_west.long_name    = 'latitude of V-points, western boundary condition'
    lat_v_west.units        = 'degree_north'
    
    #######################################
    lon_v_east            = ncid.createVariable('lon_v_east','f8',('eta_v',))
    lon_v_east.long_name    = 'longitude of V-points, eastern boundary condition'
    lon_v_east.units        = 'degree_east'
    
    lat_v_east            = ncid.createVariable('lat_v_east','f8',('eta_v',))
    lat_v_east.long_name    = 'latitude of V-points, eastern boundary condition'
    lat_v_east.units        = 'degree_north'
    
    ###############################################################################
    ############################## - Assign values - ##############################
    ###############################################################################
    
    bry_time[:]         = bry_time_var
    
    
    salt_east[:,:,:]    = salt_east_var
    salt_west[:,:,:]    = salt_west_var
    
    temp_east[:,:,:]    = temp_east_var
    temp_west[:,:,:]    = temp_west_var
    
    
    u_east[:,:,:]       = u_east_var
    u_west[:,:,:]       = u_west_var
        
    v_east[:,:,:]       = v_east_var
    v_west[:,:,:]       = v_west_var
    
    
    
    ubar_east[:,:]      = ubar_east_var
    ubar_west[:,:]      = ubar_west_var
    
    vbar_east[:,:]      = vbar_east_var
    vbar_west[:,:]      = vbar_west_var
    
    zeta_east[:,:]      = zeta_east_var
    zeta_west[:,:]      = zeta_west_var
    


    lon_rho_east[:]     = lon_rho_east_var
    lon_rho_west[:]     = lon_rho_west_var
    
    lat_rho_east[:]     = lat_rho_east_var
    lat_rho_west[:]     = lat_rho_west_var
    
    lon_u_east[:]       = lon_u_east_var
    lon_u_west[:]       = lon_u_west_var
    
    lat_u_east[:]       = lat_u_east_var
    lat_u_west[:]       = lat_u_west_var
    
    lon_v_east[:]       = lon_v_east_var
    lon_v_west[:]       = lon_v_west_var
    
    lat_v_east[:]       = lat_v_east_var
    lat_v_west[:]       = lat_v_west_var
    
    ncid.close()
