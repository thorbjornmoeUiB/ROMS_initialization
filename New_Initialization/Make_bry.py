#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 14 15:11:13 2023

@author: thorbjorn
"""

""" 
This script creates a boundary file for idealized model runs using ROMS. 


This script currently simply uses the input fields as boundary conditions. This 
can be changes by simply setting the desired variables equal to some other,
predetermined value (i.e. temp_east[:,:] = temp_from_file or otherwise)


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
from GlobalParameters import *                # This is where I change user outputs now

print('###############################################################################')
print('############################# - Making bry file - #############################')
print('###############################################################################')

###############################################################################
################################ - Functions - ################################
###############################################################################

def make_XY(ds,zd,xd):                                    # get coordinates for plotting
    stretch = stretch_ds['Cs_r'][:]                       # stretching
    heighta = ds['h'][:,:]                                # depth array
    
    bath = np.zeros((30,x,y))
    for i in range(x):
        for j in range(y):
            bath[:,i,j] = heighta[i,j]*stretch            # stretching * depth = coordinate

    X = repmaty.repmat(np.linspace(0,x,xd),zd,1)
    Y = bath

    return(bath,X,Y)

def get_sizes(mnz,mxz,mnh,mxh,velo): # function for making "box" on boundary == outflow/inflow
    z   = mxz - mnz
    h   = mxh - mnh
    box = np.ones([z,h])*velo
    return(mnz,mxz,mnh,mxh,box)

def double_trouble(old_var):             # add arbitrary dimension (dim = 2)
    return np.array([old_var,old_var])

###############################################################################
############################### - User options - ##############################
###############################################################################

#File directory/names
common_dir          = GlobalParameters.common_dir
directory           = GlobalParameters.exp_dir

bry_file            = directory + GlobalParameters.bry_file
grid_file           = directory + GlobalParameters.grid_file
Input_file          = directory + GlobalParameters.init_file
if GlobalParameters.Use_stretch == 1:
    stretching_file = directory + GlobalParameters.stretching
else:
    stretching_file = common_dir + GlobalParameters.stretching_alt

grid_ds             = nc.Dataset(grid_file)
init_ds             = nc.Dataset(Input_file)
stretch_ds          = nc.Dataset(stretching_file)

Descrip_grd         = GlobalParameters.General_desc + GlobalParameters.spec_desc_bry
Author              = GlobalParameters.Author
create_bry_file     = GlobalParameters.Create_files
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
x   = np.shape(grid_ds['h'])[0]                           # dim meridional ~ 400 km   (1 km resolution)  
y   = np.shape(grid_ds['h'])[1]                           # dim zonal      ~ 650 km 

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

lon_rho = grid_ds['lon_rho']                                # coordinates
lat_rho = grid_ds['lat_rho']                                # ...
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

salt_east_var    = double_trouble(salt[0,:,:,eas])
salt_west_var    = double_trouble(salt[0,:,:,wes])

temp_east_var    = double_trouble(temp[0,:,:,eas])
temp_west_var    = double_trouble(temp[0,:,:,wes])


u_east_var       = double_trouble(u[0,:,:,eas])
u_west_var       = double_trouble(u[0,:,:,wes])

v_east_var       = double_trouble(v[0,:,:,eas])
v_west_var       = double_trouble(v[0,:,:,wes])


ubar_east_var    = double_trouble(ubar[:,:,eas])
ubar_west_var    = double_trouble(ubar[:,:,wes])

vbar_east_var    = double_trouble(vbar[:,:,eas])
vbar_west_var    = double_trouble(vbar[:,:,wes])


zeta_east_var    = double_trouble(zeta[:,:,eas])
zeta_west_var    = double_trouble(zeta[:,:,wes])


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
if pNIJ == 1:
    if GlobalParameters.steepnes_factor == 7:
        if res == 3:
            z1,z2,h1,h2,box = get_sizes(0,13,50,59,-0.0955)  
        elif res == 1.5:
            if GlobalParameters.Use_stretch == 1:
                if GlobalParameters.theta_b == 2:
                    z1,z2,h1,h2,box = get_sizes(0,13,100,116,-0.1025)
                elif GlobalParameters.theta_b == 14:
                    z1,z2,h1,h2,box = get_sizes(0,13,100,118,-0.102)
            else:
                if GlobalParameters.New_S_exag == 0:
                    z1,z2,h1,h2,box = get_sizes(0,13,100,119,-0.1)
                elif GlobalParameters.New_S_exag == 1:
                    z1,z2,h1,h2,box = get_sizes(0,15,100,120,-0.101)
                elif GlobalParameters.New_S_exag == 2:
                    z1,z2,h1,h2,box = get_sizes(0,14,100,119,-0.101)
        elif res == 1:
            if GlobalParameters.Use_stretch == 1:
                if GlobalParameters.theta_b == 2:
                    z1,z2,h1,h2,box = get_sizes(0,13,150,174,-0.10)
                elif GlobalParameters.theta_b == 14:
                    z1,z2,h1,h2,box = get_sizes(0,13,150,177,-0.10)
            else:
                if GlobalParameters.New_S_exag == 0:
                    z1,z2,h1,h2,box = get_sizes(0,13,150,178,-0.105)
                elif GlobalParameters.New_S_exag == 1:
                    z1,z2,h1,h2,box = get_sizes(0,15,150,179,-0.103)
                elif GlobalParameters.New_S_exag == 2:
                    z1,z2,h1,h2,box = get_sizes(0,15,150,178,-0.10)
    elif GlobalParameters.steepnes_factor == 15:
        if res == 3:
            z1,z2,h1,h2,box = get_sizes(0,13,50,59,-0.0955)  
        elif res == 1.5:
            if GlobalParameters.Use_stretch == 1:
                if GlobalParameters.theta_b == 2:
                    z1,z2,h1,h2,box = get_sizes(0,13,100,116,-0.1025)
                elif GlobalParameters.theta_b == 14:
                    z1,z2,h1,h2,box = get_sizes(0,13,100,118,-0.102)
            else:
                if GlobalParameters.New_S_exag == 0:
                    z1,z2,h1,h2,box = get_sizes(0,14,100,119,-0.101)
                elif GlobalParameters.New_S_exag == 1:
                    z1,z2,h1,h2,box = get_sizes(0,15,100,120,-0.101)
                elif GlobalParameters.New_S_exag == 2:
                    z1,z2,h1,h2,box = get_sizes(0,14,100,119,-0.101)
        elif res == 1:
            if GlobalParameters.Use_stretch == 1:
                if GlobalParameters.theta_b == 2:
                    z1,z2,h1,h2,box = get_sizes(0,13,150,174,-0.10)
                elif GlobalParameters.theta_b == 14:
                    z1,z2,h1,h2,box = get_sizes(0,13,150,177,-0.10)
            else:
                if GlobalParameters.New_S_exag == 0:
                    z1,z2,h1,h2,box = get_sizes(0,13,150,178,-0.105)
                elif GlobalParameters.New_S_exag == 1:
                    z1,z2,h1,h2,box = get_sizes(0,15,150,179,-0.103)
                elif GlobalParameters.New_S_exag == 2:
                    z1,z2,h1,h2,box = get_sizes(0,15,150,178,-0.10)

    u_west_var[:,z1:z2,h1:h2]               = np.array([box,box])
    u_west_var_for_ubar_calc[:,z1:z2,h1:h2] = np.array([box,box])
        
if pNIIC == 1:
    z1,z2,h1,h2,box = get_sizes(0,30,75,142,0.10) 
    
    u_west_var[:,z1:z2,h1:h2]               = np.array([box,box])
    u_west_var_for_ubar_calc[:,z1:z2,h1:h2] = np.array([box,box])

if pNIIC_hydro == 1:
    z1,z2,y1,y2,vel,temp,salt = 0,30,70,142, 0.107, 3 , 35.2 #0,30,75,142, 0.116, 3 , 35.2
    
    # Velocity
    z1,z2,h1,h2,box = get_sizes(z1,z2,y1,y2,vel) 
    u_west_var[:,z1:z2,h1:h2]               = np.array([box,box])
    u_west_var_for_ubar_calc[:,z1:z2,h1:h2] = np.array([box,box])
    
    
    z1,z2,y1,y2,vel,temp,salt = 10,30,70,135, 0.107, 3 , 35.2
    
    # temperataure
    z1,z2,h1,h2,box = get_sizes(z1,z2,y1,y2,temp) 
    temp_west_var[:,z1:z2,h1:h2]            = np.array([box,box])
    
    # salinity
    z1,z2,h1,h2,box = get_sizes(z1,z2,y1,y2,salt) 
    salt_west_var[:,z1:z2,h1:h2]            = np.array([box,box])

if pNIIC2 == 1:
    z1,z2,h1,h2,box = get_sizes(0,30,60,127,0.08)
    
    u_east_var[:,z1:z2,h1:h2]               = np.array([box,box])
    u_east_var_for_ubar_calc[:,z1:z2,h1:h2] = np.array([box,box])
        
if pNIJ_md == 1:
    
    z1,z2,h1,h2,box = get_sizes(7,15,150,192,-0.10) 
    
    u_west_var[:,z1:z2,h1:h2]               = np.array([box,box])
    u_west_var_for_ubar_calc[:,z1:z2,h1:h2] = np.array([box,box])

if pIFSJ == 1:
    if onslope == 1:
        z1,z2,h1,h2,box = get_sizes(0,8,105,120,0.095)
    elif onslope == 2:
        z1,z2,h1,h2,box = get_sizes(0,12,99,118,0.104)
    else:
        z1,z2,h1,h2,box = get_sizes(0,8,115,128,0.10)
   
    u_east_var[:,z1:z2,h1:h2]               = np.array([box,box])
    u_east_var_for_ubar_calc[:,z1:z2,h1:h2] = np.array([box,box])


if Iplottin == 1:
###############################################################################
########################## - plot volume transport - ##########################
###############################################################################
    
    if pNIJ == 1 or pNIIC == 1 or pNIIC_hydro == 1 or pNIIC2 == 1 or pNIJ_md == 1 or pIFSJ == 1:
        
        # make area
        depth_a = grid_ds['h'][:,0]
        deptha  = stretch_ds['Cs_r'][:]
        
        z10 = np.zeros([z,x])
        for i in range(x):
            z10[:,i] = deptha*depth_a[i]
        
        depth = (z10*-1)[::-1]
        
        pp,X,Y = make_XY(grid_ds,30,x-1)
        X,Y    = X[:-1,4:],Y[:-1,5:,0]
        
        u_W_spec = u_west_var[0,:,5:]
        u_E_spec = u_west_var[0,:,5:]
        
        diffs = np.zeros([29,x-1])
        
        for i in range(x-1):
            diffs[:,i] = np.diff(depth[:,i])
        
        #westward and eastward seperately
        u_W_spec = np.where(u_W_spec>0,u_W_spec,np.zeros(np.shape(u_W_spec)))[:-1,:]
        u_E_spec = np.where(u_E_spec<0,u_E_spec,np.zeros(np.shape(u_E_spec)))[:-1,:]
        
        # calc Volume Transport (Cross sectional area * current speed)
        area = diffs*1000*res
        area = area[:,4:]
        area = area[::-1]
        vol_trans_W = area * u_W_spec
        vol_trans_E = area * u_E_spec
        
        
        
        
        
        
        
        X1,Y1,vol_trans1 = X[:,:],Y[:,:],vol_trans_W[:,:]
        
        #plot 
        fig = plt.figure(figsize=[14,12])
        plt.suptitle('FROM Make_bry.py')
        grid = plt.GridSpec(4,2,hspace=0.4)
        
        ax1 = fig.add_subplot(grid[0, 0])
        ax1.set_title('Eastward (NIIC): current')
        z1 = ax1.contourf(X1,Y1,u_W_spec[:,:],cmap='Reds')
        plt.colorbar(z1,ax=ax1)
        
        ax1 = fig.add_subplot(grid[0, 1])
        ax1.set_title('Westward (NIJ): current')
        z1 = ax1.contourf(X1,Y1,u_E_spec[:,:],cmap='Blues_r')
        plt.colorbar(z1,ax=ax1)
        
        
        
        ax1 = fig.add_subplot(grid[1, 0])
        ax1.set_title('Eastward (NIIC): \n area of each grid cell')
        z1 = ax1.contourf(X,Y,area)
        plt.colorbar(z1,ax=ax1)
        
        
        ax2 = fig.add_subplot(grid[2, 0])
        ax2.set_title('u')
        z2 = ax2.contourf(X,Y,u_W_spec,cmap='Reds')
        plt.colorbar(z2,ax=ax2)
        
        
        
        
        ax3 = fig.add_subplot(grid[3, 0])
        ax3.set_title('volume transport \n total: '+str(round(np.nansum(vol_trans1)/1000000,3))+' Sv')
        z3 = ax3.contourf(X1,Y1,vol_trans1,cmap='Reds')
        plt.colorbar(z3,ax=ax3)
        
        ax1 = fig.add_subplot(grid[1, 1])
        ax1.set_title('Westward (NIJ): \narea of each grid cell')
        z1 = ax1.contourf(X,Y,area)
        plt.colorbar(z1,ax=ax1)
        
        
        ax2 = fig.add_subplot(grid[2, 1])
        ax2.set_title('u')
        z2 = ax2.contourf(X,Y,u_E_spec,cmap='Blues_r')
        plt.colorbar(z2,ax=ax2)
        
        
        X1,Y1,vol_trans1 = X[:,:],Y[:,:],vol_trans_E[:,:]
        ax3 = fig.add_subplot(grid[3, 1])
        ax3.set_title('volume transport \n total: '+str(round(np.nansum(vol_trans1)/1000000,3))+' Sv')
        z3 = ax3.contourf(X1,Y1,vol_trans1,cmap='Blues_r')
        plt.colorbar(z3,ax=ax3)
    
    
###############################################################################
####################### - plot depth averaged velocity - ######################
###############################################################################

if pNIJ == 1 or pNIIC == 1 or pNIIC_hydro == 1 or pNIJ_md == 1:
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
    
    

if pIFSJ == 1:
    Cs_w = stretch_ds['Cs_r'][:]
    
    diff_Cs_w=[]
    for i in range(1,31,1):
        diff_Cs_w.append(Cs_w[i]-Cs_w[i-1])
    
    ubarc = np.zeros([y])
    for i in range(0,y,1):
        un = u_east_var_for_ubar_calc[0,:,i]
        weighted = un*diff_Cs_w
        weighted_mean = np.nansum(weighted)
        ubarc[i] = weighted_mean
        
    ubar_east_var = np.array([ubarc,ubarc])
        
    if plotting == 1 or Iplottin == 1:
        plt.figure()
        plt.suptitle('FROM Make_bry.py')
        plt.plot(np.linspace(0,y,y),ubarc)
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
