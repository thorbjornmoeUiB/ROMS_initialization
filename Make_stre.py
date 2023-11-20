#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 14 15:08:37 2023

@author: thorbjorn
"""

###############################################################################
################################# - Packages - ################################
###############################################################################

import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import xarray as xr
from GlobalParameters import *                # This is where I change user outputs now

print('###############################################################################')
print('############################ - Making stre file - #############################')
print('###############################################################################')

###############################################################################
##################### - Options from "GlobalParameters" - #####################
###############################################################################

directory            = GlobalParameters.exp_dir
grid_file_name       = GlobalParameters.grid_file
Stretching_file_name = GlobalParameters.stretching

Use_stretch          = GlobalParameters.Use_stretch
Create_file          = GlobalParameters.Create_files
plotting             = GlobalParameters.Gen_plot
Iplottin             = GlobalParameters.Imp_plot

N                    = GlobalParameters.z       
theta_s              = GlobalParameters.theta_s       
theta_b              = GlobalParameters.theta_b         

###############################################################################
############################### - User options - ##############################
###############################################################################

# Parameters for the stretching vector, see romsWIKI or ocean.in 
vstretch             = 2         # which stretching func to use, vstretch 1 is funky, but we only use 2 
hc                   = 40        # critical depth
kgrid                = 1         # What grid to produce (middle or top/bottom of grid)

###############################################################################
############################# - files and shapes - ############################
###############################################################################

ds                   = nc.Dataset(directory + grid_file_name)

x = np.shape(ds['h'])[0]
y = np.shape(ds['h'])[1]

###############################################################################
############################ - stretching vector - ############################
###############################################################################

# function for grid spacing vector 
# adapted from matlab script, possibly poorly translated...
def streching(Vstretching,theta_s,theta_b,hc,N,kgrid):
    s = []
    c = []
    #parameters can be added here as in matlab script, but do not se how that's necessary for this
   
    Np = N + 1
    
    if Vstretching == 1: #Song and Haidvogel (1994) Version of stretching
        
        ds=1.0/N 
	  
        if kgrid == 1:  #r or w?
            Nlev    = Np 
            lev     = np.arange(0,N,1)
            s       = (lev-N)*ds
        else: #the other one (r or w)
            Nlev    = N 
            lev     = np.arange(0,N,1) - 0.5 
            s       = (lev-N) * ds 
        if theta_s > 0:     #func for more grid cells at surface
            Ptheta  = (np.sinh(theta_s*s))/np.sinh(theta_s)
            Rtheta  = np.tanh(theta_s*(s+0.5))/(2*np.tanh(0.5*theta_s))-0.5
            c       = ((1-theta_b)*Ptheta)+(theta_b*Rtheta)
        else:
            c = s
        
    # the one we use """ 
    elif Vstretching == 2: #A. Shchepetkin (2005) version
        alfa    = 1.0
        beta    = 1.0
        ds      = 1/N
        if kgrid == 1: #r or w? 
            Nlev    = Np 
            lev     = np.arange(0,N,1)
            s       = (lev-N)*ds
        else: #the other one (r or w)
            Nlev    = N 
            lev     = np.arange(0,N,1)-0.5 
            s       = (lev-N)*ds 
        if theta_s > 0:
            Csur=(1.0-np.cosh(theta_s*s))/(np.cosh(theta_s)-1.0)
            if theta_b > 0:
                Cbot    = -1.0+np.sinh(theta_b*(s+1))/np.sinh(theta_b)
                weigth  = ((s+1)**alfa) * (1+(alfa/beta)*((1-(s+1)**beta)))
                c       = weigth*Csur+(1-weigth)*Cbot
            else:
                c = Csur
        else:
            c = s
    return(s,c)

s,c = streching(vstretch,theta_s,theta_b,hc,N,kgrid)

###############################################################################
####################### - Plotting layer distribution - #######################
###############################################################################

#Compute new layers for every grid point
def Compute_depth_layers(ds, hmin=0.1):
    #compute depths of ROMS vertical levels (Vtransform = 2)
   
    # compute vertical transformation functional
    S_rho = (ds.hc * ds.s_rho + ds.Cs_r * ds.h) / (ds.hc + ds.h)
   
    # compute depth of rho (layers)
    z_rho = ds.zeta + (ds.zeta + ds.h) * S_rho
   
    # transpose arrays and fill NaNs with a minimal depth
    ds['z_rho'] = z_rho.transpose(*('s_rho','eta_rho','xi_rho','Cs_r'),
                                  transpose_coords=False).fillna(hmin)
    
    # add z_rho to xarray coordinates
    ds = ds.set_coords(['z_rho'])
    return ds


# this part works better with xarray
ds          = xr.open_dataset(directory + grid_file_name)
ds['hc']    = hc
ds['s_rho'] = s
ds['Cs_r']  = c
ds['h1']    = 600
ds['zeta']  = 0
ds1         = Compute_depth_layers(ds, hmin=0.1)


#In order to plot with pre-made rutine 
air = ds.z_rho
air.attrs = ds.z_rho.attrs
if plotting == 1:
    # plot line for every vertical grid-cell (arbitrary south-north transect)
    fig = plt.figure(figsize=[9,6])
    liste = ['#af1446', '#c1274a', '#d23a4e', '#de4c4b', '#e85b48', '#f26944', '#f67c4a', '#f99355', '#fca85e', '#fdb96a', '#fec877', '#feda86', '#fee695', '#fff0a6', '#fffab6', '#fbfdb8', '#f3faac', '#ebf7a0', '#dff299', '#caea9e', '#b8e2a1', '#a4daa4', '#8fd2a4', '#76c8a5', '#62bda7', '#52abae', '#4199b6', '#3585bb', '#4273b3', '#5061aa']
    num = 0
    num1=0
    for i in range(0,30):
        air1d1 = air.isel(s_rho=29,xi_rho=15,Cs_r=num)
        air1d1.plot(color=liste[num1],alpha=0.8)
        num = num + 1
        num1 = num1 + 1
    plt.suptitle('FROM Make_stre.py')
    #writing info on plot
    text1 = x
    for i,j,k in zip(range(1,31),c*1151,liste):
        if i < 19:
            plt.text(text1,j,i,fontsize=8,color = 'k')
        elif (i % 2) != 0 and i < 24:
            plt.text(text1,j,i,fontsize=8,color = 'k')
        if i == 30:
            plt.text(text1,j-20,i,fontsize=8,color = 'k')
    plt.xticks(ticks=np.linspace(0,268,268)[1::int(100/1.5)], labels=[0,100,200,300,400])
    plt.xlim([1,280])
    plt.xlabel('Distance North [km]')
    plt.ylabel('Depth [m]')
    plt.ylim([-1200,15])
#%%
if plotting == 1:
    
    fig = plt.figure(figsize=[8,6])
    grid = plt.GridSpec(2,1, hspace=0.3,height_ratios=[2,1])
    plt.suptitle('FROM Make_stre.py')
    
    ax1 = fig.add_subplot(grid[0, 0])
    ax1.plot(c*1151,'*',color='Teal')
    ax1.set_ylabel('Depth [m]')
    ax1.set_xlabel('# layer')
    ax1.set_xlim([-1,30])
    ax1.set_ylim([-1200,50])
    ax1.grid(alpha=0.4)
    ax1.text(23,-1150,'Layer depth',fontsize=15)
    ax1.text(-0.5,-300,'Theta_b = '+str(theta_b)+'\nTheta_s = '+str(theta_s),fontsize=15)
    
    ax1 = fig.add_subplot(grid[1, 0])
    ax1.plot(np.diff(c*1151),'*',color='Teal')
    ax1.set_ylabel('Depth [m]')
    ax1.set_xlabel('# layer')
    ax1.set_xlim([-1,30])
    ax1.set_ylim([0,80])
    ax1.grid(alpha=0.4)
    ax1.text(20,60,'diff between layer',fontsize=15)
    
    depth = 1151

#%%
if Use_stretch == 1 and Create_file == 1:
    ###############################################################################
    ################################ - Transect - #################################
    ###############################################################################
    
    H = nc.Dataset(directory+grid_file_name)['h'][:,0]
    
    Ht = np.zeros([N,x])
    for i in range(x):
        Ht[:,i] = c*H[i]
        
    ###############################################################################
    ############################## - Create .nc file - ############################
    ###############################################################################
    
    fn = directory+Stretching_file_name
    ds = nc.Dataset(fn, 'w', format='NETCDF4')
    
    ###############################################################################
    ############################ - Create dimensions - ############################
    ###############################################################################
    
    z_l             = ds.createDimension('z_l', N)
    x_l             = ds.createDimension('x_l', x)
    
    ###############################################################################
    ############################ - Define variables - #############################
    ###############################################################################
    
    Cs_r            = ds.createVariable('Cs_r', 'f4', ('z_l',))
    Cs_r.units      = 'dimless (0-1)'
    
    Z_lvl_max       = ds.createVariable('z_max', 'f4', ('z_l',))
    Z_lvl_max.units = 'meters'
    
    cross_sec       = ds.createVariable('cross_sec', 'f4', ('z_l','x_l',))
    cross_sec.units = 'meters'
    
    ###############################################################################
    ############################ - Assign variables - #############################
    ###############################################################################
    
    Cs_r[:]         = c
    Z_lvl_max[:]    = c*1151
    cross_sec[:,:]  = Ht
    
    ds.close()
    
    
    
    
    