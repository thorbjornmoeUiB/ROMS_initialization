#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 12 11:50:35 2022

@author: thorbjornostenbymoe
"""

""" 
This script creates a input file test for ROMS idealized model runs.

In this test file, all variables except T & S are set as appropriately dimensioned 0 matrices.

The T&S fields are made horizontally homogeneous from the average hydrography from 1980-2010 in the Iceland Sea.

The transformation to sigma coordinates are made from the stretching vector (made in stretching.py), the depth coordinates 
from the grid file (made in Create_grid.py), and the averaged CTD values are from the final ctd file (made in CTD_mean.py).

For later use it should be sufficient to change the variables, as the writing to netCDF is 
made in accordance to the ROMS standards.
"""
###############################################################################
################################## - Packages - ###############################
###############################################################################

import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import cmocean
import gsw.conversions as gswc
import gsw.density as gswd
import scipy.io as sci

###############################################################################
################################ - Function - #################################
###############################################################################

def rho(T,S,z):
    
    press = gswc.p_from_z(z,68)
    SA    = gswc.SA_from_SP(S, press,-10,68)
    
    return(gswd.sigma0(SA,T))

###############################################################################
############################### - User options - ##############################
###############################################################################

#File directory/names
directory           = '/Users/thorbjornostenbymoe/Desktop/Initialization_of_model/MODEL_SETUP/ROMS_in/' #User directory/path
input_file          = 'init_file.nc'
grid_file           = 'grid_file.nc'
stretching_file     = 'Stretching.nc'

#data paths
variable_mean_file  = 'New_depth_profiles_for_model.nc'
Raw_profiles        = '/Users/thorbjornostenbymoe/Downloads/InputProfileCentralIS_Thorbjørn/InputProfileCentralIcelandSea.mat'

grid_ds             = nc.Dataset(directory + grid_file)
stre_ds             = nc.Dataset(directory + stretching_file)
prof_ds             = nc.Dataset(directory + variable_mean_file)
mat_rp              = sci.loadmat(Raw_profiles)


#file description
Descrip_grd         = 'Iceland slope - ROMS idealized model run - Input file'  
Author              = 'Thorbjoern Oestenby Moe'

#logical switches
create_input_file   = 1                                                   # create input (1 = true)
plotting            = 1                                                   # plot input   (1 = true)

#Dimensions
x  = 268
y  = 401
z  = 30

################################################################
##################### - Load mean profile - ####################
################################################################

temp_d  = prof_ds['CT'][:] 
sal_d   = prof_ds['SA'][:] 
depth_d = prof_ds['z'][:]

#convert to practical salinity
pres    = gswc.p_from_z(depth_d,68)
salt_d  = gswc.SP_from_SA(sal_d,pres,-10,68)

temp    = np.array(temp_d)
salt    = np.array(salt_d)

################################################################
################## - make sigma depth array - ##################
################################################################

depths  = grid_ds['h'][:,0]  #depths
stretc  = stre_ds['Cs_r'][:] #stretching vec
H       = np.zeros([30,268])

for i in range(268):
    H[:,i] = - depths[i]*stretc

################################################################
########## - Apply profiles homogeneously to domain - ##########
################################################################

# apply profiles to correct depth along meridional cross sect (in sigma coords)
temp_2d = np.zeros([z, x])
salt_2d = np.zeros([z, x])
for i in range(0,z):
    for j in range(x):
        temp_2d[i,j] = temp[int(H[i,j])]
        salt_2d[i,j] = salt[int(H[i,j])]

# spread across zonal domain
temp_3d = np.zeros([z,x,y+1])
salt_3d = np.zeros([z,x,y+1])
for i in range(y+1):
    temp_3d[:,:,i] = temp_2d
    salt_3d[:,:,i] = salt_2d
    
# add arbitrary time dimension
temperature = np.array([temp_3d])
salinity    = np.array([salt_3d])

###############################################################################
################# - Make empty zeta, u, v, ubar, vbar data - ##################
###############################################################################

#U
ds_u_3d    = np.zeros((30, 268, 401))
ds_u_4d    = np.zeros((1, 30, 268, 401))
ds_ubar_3d = np.zeros((1, 268, 401))

#V
ds_v_3d    = np.zeros((30, 267, 402))
ds_v_4d    = np.zeros((1, 30, 267, 402))
ds_vbar_3d = np.zeros((1, 267, 402))

#Z
ds_zeta_3d = np.zeros((1, 268, 402))

################################################################
##################### - Plotting density - #####################
################################################################

if plotting == 1:
    X = np.matlib.repmat(np.linspace(0,x*1.5,x),30,1)
    Y = - H
    r = rho(temp_2d,salt_2d,Y)
    
    lvlt = np.linspace(-0.75,0.15,13)
    lvls = np.linspace(34.76,34.92,17)
    
    
    fig  = plt.figure(figsize=[8,6])
    grid = plt.GridSpec(1,1, hspace=0.15)
    
    ax1 = fig.add_subplot(grid[0, 0])
    ax1.set_ylim([-1200,0])
    ax1.set_facecolor('lightgray')
    ax1.plot(X[0,:],Y[0,:],linewidth=2,color='gray')
    z1 = ax1.contourf(X,Y,r,np.linspace(27.9,28.075,33),cmap=cmocean.cm.dense) 
    cbar = fig.colorbar(z1,ax=ax1)
    
    c = ax1.contour(X,Y,r,[28.04,28.06,28.07,28.08],colors='w',alpha=0.4,linewidths=0.75)
    plt.clabel(c,inline=1)

################################################################
############# - Plotting temperature and salinity - ############
################################################################
if plotting == 1:
    
    # load raw profiles 
    IcSe = mat_rp['IS']
    IST  = IcSe['temp'][0]
    ISP  = IcSe['press'][0]
    ISS  = IcSe['sal'][0]
    
    # levels
    r_lvl = [27.94,28.,28.03,28.05,28.06]
    lvlt  = np.linspace(-0.75,0.15,13)
    lvls  = np.linspace(34.76,34.92,17)
    
    #fig
    fig = plt.figure(figsize=[9,1.54])
    grid = plt.GridSpec(1,5, hspace=0.15,width_ratios=(2,0.6,0.15,2,0.6),wspace=0.03)
    
    
    ax1 = fig.add_subplot(grid[0, 0])
    ax1.set_facecolor('lightgray')
    ax1.plot(X[0,:]-2,Y[0,:]-12,linewidth=2,color='gray')
    z1 = ax1.contourf(X,Y,temp_2d,lvlt,cmap=cmocean.cm.thermal) 
    
    xloc = 360
    manloc=[(360,-50),(360,-150),(360,-300),(360,-600),(360,-800)]
    c = ax1.contour(X,Y,r,levels=r_lvl,colors='k',linewidths=0.75)
    plt.clabel(c,inline=1,manual=manloc)
    
    ax1.text(1,-1150,'(a)')
    ax1.set_ylim([-1200,0])
    ax1.set_yticks([-1200,-1000,-800,-600,-400,-200,0])
    ax1.yaxis.set_ticklabels(['1200','','800','','400','','0'])
    ax1.set_ylabel('Depth [m]')
    ax1.set_xticks([0,100,200,300,401])
    ax1.set_xlabel('Distance North [km]')
    ax1.xaxis.set_ticklabels([0,100,200,300,''])
    ax1.set_xlim([0,402])
    
    cbar_ax = fig.add_axes([0.125, -0.35, 0.33, 0.1])
    c = fig.colorbar(z1,cax=cbar_ax,orientation='horizontal')
    c.set_label('Conservative Temperature [C'+r'$^{\circ}$'+']')
    
    
    
    
    ax11 = fig.add_subplot(grid[0, 1])
    ax11.plot([0,0],[1200,0],'--',linewidth=1,color='dimgray')
    ax11.plot([-1,-1],[1200,0],'--',linewidth=1,color='dimgray')
    ax11.plot([1,1],[1200,0],'--',linewidth=1,color='dimgray')
    
    for i in range(300):
        ax11.plot(IST[i],ISP[i],'*',color='lightgray',alpha=1,markersize=0.5)
    ax11.set_ylim(1200,0)
    ax11.yaxis.set_ticklabels([])
    ax11.set_yticks([])
    ax11.set_xlim([-2,2])
    ax11.set_xticks([-2,-1,0,1,2])
    ax11.xaxis.set_ticklabels(['$-2$','',0,'','2'])
    ax11.plot(temp,pres,color='k')
    ax11.text(1,1150,'(b)')
    
    
    
    
    ax2 = fig.add_subplot(grid[0, 3])
    ax2.set_facecolor('lightgray')
    ax2.plot(X[0,:]-2,Y[0,:]-12,linewidth=2,color='gray')
    z2 = ax2.contourf(X,Y,salt_2d,lvls,cmap=cmocean.cm.haline)
    
    c = ax2.contour(X,Y,r,levels=r_lvl,colors='k',linewidths=0.75)
    plt.clabel(c,inline=1,manual=manloc)
    
    tikky = np.arange(34.76,34.92,0.02)
    cbar.set_label('Practical Salinity')
    ax2.text(1,-1150,'(c)')
    ax2.set_yticks([-1200,-800,-400,0])
    ax2.yaxis.set_ticklabels([])
    ax2.set_ylim([-1200,0])
    ax2.set_xlabel('Distance North [km]')
    ax2.set_xticks([0,100,200,300])
    ax2.xaxis.set_ticklabels([0,100,200,300])
    ax2.set_yticks([-1200,-1000,-800,-600,-400,-200,0])
    ax2.set_xlim([0,402])
    
    cbar_ax = fig.add_axes([0.528, -0.35, 0.33, 0.1])
    c = fig.colorbar(z2,cax=cbar_ax,orientation='horizontal')
    c.set_ticks(lvls[::3])
    c.set_label('Practical Salinity')
    
    
    
    
    
    ax22 = fig.add_subplot(grid[0, 4])
    ax22.plot([34.625,34.625],[1200,0],'--',linewidth=1,color='dimgray')
    ax22.plot([34.75,34.75],[1200,0],'--',linewidth=1,color='dimgray')
    ax22.plot([34.875,34.875],[1200,0],'--',linewidth=1,color='dimgray')
    
    for i in range(300):
        ax22.plot(ISS[i],ISP[i],'*',color='lightgray',alpha=1,markersize=0.5)
    ax22.set_ylim(1200,0)
    ax22.yaxis.set_ticklabels([])
    ax22.set_yticks([])
    ax22.set_xlim([34.5,35])
    ax22.set_xticks([34.5,34.625,34.75,34.875,35])
    ax22.xaxis.set_ticklabels([34.5,'','','',35])
    ax22.plot(salt,pres,color='k')
    ax22.text(34.495,1150,'(d)')
    
    #fig.savefig('/Users/thorbjornostenbymoe/Desktop/PUSH_init',dpi=300, facecolor='w', edgecolor='w', orientation='landscape',format=None,transparent=False,bbox_inches='tight',pad_inches=0.25)

###############################################################################
############################# - write to NetCDF - #############################
###############################################################################

if create_input_file == 1:
    ###############################################################################
    ################################### - file - ##################################
    ###############################################################################
    
    path_string               = directory+input_file
    ncid                      = nc.Dataset(path_string,mode='w',format='NETCDF4')
    ncid.description          = Descrip_grd
    ncid.author               = Author
    
    ###############################################################################
    ############################### - Dimensions - ################################
    ###############################################################################
    
    # Meridional: eta
    eta_rho                   = ncid.createDimension('eta_rho',x)
    eta_u                     = ncid.createDimension('eta_u',  x)
    eta_v                     = ncid.createDimension('eta_v',  x-1)
    eta_psi                   = ncid.createDimension('eta_psi',x-1)
    
    # Zonal: xi
    xi_rho                    = ncid.createDimension('xi_rho', y+1)
    xi_u                      = ncid.createDimension('xi_u',   y)
    xi_v                      = ncid.createDimension('xi_v',   y+1)
    xi_psi                    = ncid.createDimension('xi_psi', y)

    # z and time
    s_rho                     = ncid.createDimension('s_rho',  z)
    time                      = ncid.createDimension('time',   1)
    one                       = ncid.createDimension('one',    1)
    
    ###############################################################################
    ############################### - Variables - #################################
    ###############################################################################
    
    ocean_time              = ncid.createVariable('ocean_time','f8',('time'))
    ocean_time.units          = 'day'
    ocean_time.time           = 'ocean_time'
    ocean_time.cycle_length   = 360. 
    ocean_time.field          = 'time, scalar, series'
    ocean_time.long_name      = 'day since initialization'
    
    salt                    = ncid.createVariable('salt','f8',('time','s_rho','eta_rho','xi_rho',))
    salt.long_name            = 'salinity'
    salt.time                 = 'ocean_time'

    temp                    = ncid.createVariable('temp','f8',('time','s_rho','eta_rho','xi_rho',))
    temp.long_name            = 'temperature'
    temp.units                = 'degree Celsius'
    temp.time                 = 'ocean_time'
    
    u                       = ncid.createVariable('u','f8',('time','s_rho','eta_u','xi_u',))
    u.long_name               = '3D U-momentum'
    u.units                   = 'm/s'
    u.time                    = 'ocean_time'
    
    v                       = ncid.createVariable('v','f8',('time','s_rho','eta_v','xi_v',))
    v.long_name               = '3D V-momentum'
    v.units                   = 'm/s'
    v.time                    = 'ocean_time'
    
    ubar                    = ncid.createVariable('ubar','f8',('time','eta_u','xi_u',))
    ubar.long_name            = '2D U-momentum'
    ubar.time                 = 'ocean_time'
    ubar.units                = 'm/s'
    
    vbar                    = ncid.createVariable('vbar','f8',('time','eta_v','xi_v',))
    vbar.long_name            = '2D v-momentum'
    vbar.time                 = 'ocean_time'
    vbar.units                = 'm/s'
    
    zeta                    = ncid.createVariable('zeta','f8',('time','eta_rho','xi_rho',))
    zeta.long_name            = 'free surface'
    zeta.time                 = 'ocean_time'
    
    ###############################################################################
    ############################## - Assign values - ##############################
    ###############################################################################
    
    ocean_time[:] = [0.0]
    
    salt[:,:,:,:] = salinity
    temp[:,:,:,:] = temperature
    
    u[:,:,:,:]    = ds_u_4d
    v[:,:,:,:]    = ds_v_4d
    
    ubar[:,:,:]   = ds_ubar_3d
    vbar[:,:,:]   = ds_vbar_3d
    
    zeta[:,:,:]   = ds_zeta_3d
    
    
    ncid.close()

    
    
