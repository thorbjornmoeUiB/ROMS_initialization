#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 14 15:09:59 2023

@author: thorbjorn
"""

""" 
This script creates a input file for ROMS idealized model runs.

In this file, all variables except T & S are set as appropriately dimensioned 0 matrices.

The T&S fields are made horizontally homogeneous from the average hydrography from 1980-2010 in the 
Iceland Sea (Brakstad et al., 2023).

The transformation to sigma coordinates are made from the stretching vector (retrieved from the model) 
and the depth coordinates from the grid file (made in Make_grid.py).

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
import numpy.matlib as repmaty
from GlobParam import *                # This is where I change user outputs now

print('###############################################################################')
print('############################ - Making init file - #############################')
print('###############################################################################')

###############################################################################
################################ - Functions - ################################
###############################################################################

def get_sizes(mnz,mxz,mnh,mxh,velo): # function for making "box" on boundary == outflow/inflow
    z   = mxz - mnz
    h   = mxh - mnh
    box = np.ones([z,h])*velo
    return(mnz,mxz,mnh,mxh,box)

###############################################################################
############################### - User options - ##############################
###############################################################################

# Paths 
common_dir          = GlobParam.common_dir
directory           = GlobParam.EXP_dir
input_file          = directory + GlobParam.init_file
grid_file           = directory + GlobParam.grid_file
stre_file           = common_dir + GlobParam.Streching

grid_ds             = nc.Dataset(grid_file)
stre_ds             = nc.Dataset(stre_file)

# hydrography
hydrography = GlobParam.hydrography
if hydrography == 'w':
    variable_mean_file  = common_dir    + GlobParam.Mean_profiles_w
    Raw_profiles        = common_dir + GlobParam.Raw_profiles_w
    prof_ds             = nc.Dataset(variable_mean_file)
elif hydrography == 's':
    variable_mean_file  = common_dir    + GlobParam.Mean_profiles_s
    Raw_profiles        = common_dir + GlobParam.Raw_profiles_s
    prof_ds             = sci.loadmat(variable_mean_file)
mat_rp              = sci.loadmat(Raw_profiles)

Descrip_grd         = GlobParam.General_desc + GlobParam.spec_desc_ini
Author              = GlobParam.Author
create_input_file   = GlobParam.Create_files
plotting            = GlobParam.Gen_plot
Iplottin            = GlobParam.Imp_plot
z                   = GlobParam.z
res                 = GlobParam.Res
VS                  = GlobParam.VarSlope

#Dimensions
x   = np.shape(grid_ds['h'])[0]
y   = np.shape(grid_ds['h'])[1]
tdim = 1

#%%
################################################################
##################### - Load mean profile - ####################
################################################################

if hydrography == 'w':
    ConsTemp  = prof_ds['CT'][:] 
    AbsoSalt   = prof_ds['SA'][:] 
    depth_d = prof_ds['z'][:]
elif hydrography == 's':
    ConsTemp  = prof_ds['CT_av']
    AbsoSalt   = prof_ds['SA_av']
    depth_d = -1*prof_ds['dpt'].T

# Convert away from TEOS-10
pressure    = gswc.p_from_z(depth_d,68)
salt        = gswc.SP_from_SA(AbsoSalt,pressure,-10,68)
temp        = gswc.t_from_CT(AbsoSalt,ConsTemp,68)
sigma       = gswd.sigma0(salt,temp)

#%%

if VS == 0:
    ################################################################
    ################## - make sigma depth array - ##################
    ################################################################
    
    depths  = grid_ds['h'][:,0]  #depths
    stretc  = stre_ds['Cs_r'][:] #stretching vec
    
    H = -depths[:]*stretc[:,np.newaxis]


    ################################################################
    ########## - Apply profiles homogeneously to domain - ##########
    ################################################################
    
    # apply profiles to correct depth along meridional cross sect (in sigma coords)
    temp_2d = np.zeros([tdim, z, x])
    salt_2d = np.zeros([tdim, z, x])
    for i in range(0,z):
        for j in range(x):
            temp_2d[0,i,j] = temp[int(H[i,j])]
            salt_2d[0,i,j] = salt[int(H[i,j])]
    
    temperature = temp_2d[:,:,:,np.newaxis].repeat(y,3)
    salinity    = salt_2d[:,:,:,np.newaxis].repeat(y,3)

else:
    ################################################################
    ################## - make sigma depth array - ##################
    ################################################################
    
    depths  = grid_ds['h'][:,:]  #depths
    stretc  = stre_ds['Cs_r'][:] #stretching vec
    
    H = -depths[:,:]*stretc[:,np.newaxis,np.newaxis]
        
    ################################################################
    ########## - Apply profiles homogeneously to domain - ##########
    ################################################################
    
    temperature = np.zeros((tdim,z,x,y))
    salinity    = np.zeros((tdim,z,x,y))
    for i in range(z):
        for j in range(x):
            for k in range(y):
                temperature[0,i,j,k] = temp[int(H[i,j,k])]
                salinity[0,i,j,k]    = salt[int(H[i,j,k])]
               

###############################################################################
################# - Make empty zeta, u, v, ubar, vbar data - ##################
###############################################################################

#U
ds_u_4d    = np.zeros((tdim, z, x, y-1))
ds_ubar_3d = np.zeros((tdim, x, y-1))

#V
ds_v_4d    = np.zeros((tdim, z, x-1, y))
ds_vbar_3d = np.zeros((tdim, x-1, y))

#Z
ds_zeta_3d = np.zeros((tdim, x, y))

################################################################
##################### - Plotting density - #####################
################################################################

if VS == 1:
    Y_east = - H[:,:,-1]
    Y_west = - H[:,:,0]
else:
    Y_east = - H[:,:]
    Y_west = - H[:,:]
    
#convert back to TEOS-10 for plotting

pressure_west    = gswc.p_from_z(Y_west,68)
pressure_east    = gswc.p_from_z(Y_east,68)
ConsTemp_3d_east,AbsoSalt_3d_east = gswc.CT_from_t(salinity[0,:,:,-1],temperature[0,:,:,-1],68),gswc.SA_from_SP(salinity[0,:,:,-1],pressure_east,-10,68)
ConsTemp_3d_west,AbsoSalt_3d_west = gswc.CT_from_t(salinity[0,:,:,0], temperature[0,:,:,0],68), gswc.SA_from_SP(salinity[0,:,:,0], pressure_west,-10,68)

sigma_east     = gswd.sigma0(AbsoSalt_3d_east,ConsTemp_3d_east)
sigma_west     = gswd.sigma0(AbsoSalt_3d_west,ConsTemp_3d_west)

X = repmaty.repmat(np.linspace(0,x*res,x),z,1)

EAST = ['EAST',Y_east,sigma_east,ConsTemp_3d_east,AbsoSalt_3d_east]
WEST = ['WEST',Y_west,sigma_west,ConsTemp_3d_west,AbsoSalt_3d_west]

things = [EAST,WEST]

#%%
################################################################
############# - Plotting temperature and salinity - ############
################################################################
if Iplottin == 1:
    
    # load raw profiles 
    IcSe = mat_rp['IS']
    IST  = IcSe['temp'][0]
    ISP  = IcSe['press'][0]
    ISS  = IcSe['sal'][0]
    ISt = IcSe['time'][0]


    r_lvl = [27.94,28.,28.03,28.05,28.06,28.07]
    lvlt  = np.arange(-0.75,0.3,.075)
    lvls  = np.arange(34.93,35.1,0.01)
    for side in things:
        name,Y,rho,temp1,salt1 = side
        #fig
        

        fig = plt.figure(figsize=[10,1.54])
        grid = plt.GridSpec(1,2, hspace=0.15,wspace=0.2)
        #
        plt.suptitle('FROM Make_init.py',y=1.05)
        
        ax1 = fig.add_subplot(grid[0, 0])
        ax1.fill_between(X[0,:],Y[0,:],-1200,linewidth=1.5,color='lightgray')
        ax1.plot(X[0,:],Y[0,:]-12,linewidth=2,color='gray')
        z1 = ax1.contourf(X,Y,temp1,lvlt,cmap=cmocean.cm.thermal,extend='max') 
    
        
        ax1.text(1,-1150,'(a)')
        ax1.set_ylim([-1200,0])
        ax1.set_yticks([-1200,-1000,-800,-600,-400,-200,0])
        ax1.yaxis.set_ticklabels(['1200','','800','','400','','0'])
        ax1.set_ylabel('Depth [m]')
        ax1.set_xticks([0,100,200,300,401])
        ax1.set_xlabel('Distance North [km]')
        ax1.xaxis.set_ticklabels([0,100,200,300,''])
        ax1.set_xlim([0,402])
        

        
        ax2 = fig.add_subplot(grid[0, 1])
        ax2.fill_between(X[0,:],Y[0,:],-1200,linewidth=1.5,color='lightgray')
        ax2.plot(X[0,:],Y[0,:]-12,linewidth=2,color='gray')
        z2 = ax2.contourf(X,Y,salt1,lvls,cmap=cmocean.cm.haline)
        
        tikky = np.arange(34.76,34.92,0.02)
        ax2.text(1,-1150,'(b)')
        ax2.set_yticks([-1200,-800,-400,0])
        ax2.yaxis.set_ticklabels([])
        ax2.set_ylim([-1200,0])
        ax2.set_xlabel('Distance North [km]')
        ax2.set_xticks([0,100,200,300])
        ax2.xaxis.set_ticklabels([0,100,200,300])
        ax2.set_yticks([-1200,-1000,-800,-600,-400,-200,0])
        ax2.set_xlim([0,402])
        
        
        cbar_ax = fig.add_axes([0.125, -0.35, 0.35, 0.1])
        c = fig.colorbar(z1,cax=cbar_ax,orientation='horizontal')
        c.set_label('Conservative Temperature'+' [C'+r'$^{\circ}$'+']')
        
        cbar_ax = fig.add_axes([0.55, -0.35, 0.35, 0.1])
        c = fig.colorbar(z2,cax=cbar_ax,orientation='horizontal')
        c.set_ticks(lvls[::5])
        c.set_label('Absolute Salinity')
        
    
        
#%%
###############################################################################
############################# - write to NetCDF - #############################
###############################################################################

if create_input_file == 1:
    ###############################################################################
    ################################### - file - ##################################
    ###############################################################################

    ncid                      = nc.Dataset(input_file,mode='w',format='NETCDF4')
    ncid.description          = Descrip_grd
    ncid.author               = Author
    
    ###############################################################################
    ############################### - dimentions - ################################
    ###############################################################################
    
    # Meridional: eta
    eta_rho                   = ncid.createDimension('eta_rho',x)
    eta_u                     = ncid.createDimension('eta_u',  x)
    eta_v                     = ncid.createDimension('eta_v',  x-1)
    eta_psi                   = ncid.createDimension('eta_psi',x-1)
    
    # Zonal: xi
    xi_rho                    = ncid.createDimension('xi_rho', y)
    xi_u                      = ncid.createDimension('xi_u',   y-1)
    xi_v                      = ncid.createDimension('xi_v',   y)
    xi_psi                    = ncid.createDimension('xi_psi', y-1)
    
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
