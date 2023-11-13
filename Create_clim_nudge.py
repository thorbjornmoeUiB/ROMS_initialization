#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 25 10:43:04 2022

@author: thorbjornostenbymoe
"""

""" 
This script creates a nudging file for ROMS idealized model runs.

The nudging coefficients are set to be vertically uniform and constant along the boundary
for M2 (not vetical), M3, temp, and salt. Script should be modifiable enough from the user options section.

For later use it should be sufficient to change the variables, as the writing to netCDF is 
made in accordance to the ROMS standards.
"""
###############################################################################
################################## - Packages - ###############################
###############################################################################

import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy.matlib as repmaty

###############################################################################
############################### - User options - ##############################
###############################################################################

#File directory/names
directory           = '/Users/thorbjornostenbymoe/Desktop/PhD_Iceland_sea_model/Initialization/TEST6_ALL/' #User directory/path

clm_nudge_file_name = 'clm_nudge_file.nc.nosync'
grid_file           = 'grid_file.nc.nosync'
Input_file          = 'init_file.nc.nosync'

grid_ds             = nc.Dataset(directory+grid_file)

#logical switches                               
plot                = 1                                   # Logical switch 1/0
create_file         = 1                                   # Logical switch 1/0

#Experiments in "idealized modeling of the North Icelandic Jet" master thesis
pNIJ                = 1                                   # EXP 1 == prescribed NIJ outflow in west
pNIIC               = 0                                   # EXP 2 == prescribed NIIC inflow in west
pNIIC2              = 0                                   # EXP 2.5 == altered prescribed NIIC inflow in west
pNIJ_md             = 0                                   # EXP 3 == mid-depth confined prescribed NIJ outflow in west
pIFSJ               = 0                                   # EXP 4 == prescribed IFSJ outflow in east
onslope             = 0                                   # Sets degree to which pIFSJ is prescribed on the boundary, very onslope,1=onslope, 0=offslope
NIIC_EB             = 0                                   # EXP 5 == prescribed NIIC outflow in east

#Dim and plot switch
x   = np.shape(grid_ds['h'])[0]
y   = np.shape(grid_ds['h'])[1]
z   = 30
res = 1

#Trigonometric coefficients for the nudging function
nudgemax        = 0.2                                      # Max nudging coefficient, taken from previous nudging files (other values could be used)
amplitude       = nudgemax/2
pi_factor       = 5
nudge_cells     = 25 #(1.5 km = 17, 1 km = 25, 3 km = 8)    # Number of cells with nudging (away from boundary)

#.nc file info                   
Descrip_grd   = 'Iceland slope - ROMS idealized model run - nudging file, currently using 17 cell padding with nudging coeff approximated from Vår'  #Description for grid .nc file
Author        = 'Thorbjoern Oestenby Moe'

###############################################################################
############################## - Load lat/lon - ###############################
###############################################################################

lat_rho_var1 = grid_ds['lat_rho'][:,:] #using standard ROMS format
lon_rho_var1 = grid_ds['lon_rho'][:,:]

###############################################################################
########################### - Create M2_NudgeCoef - ###########################
###############################################################################

#Lon_rho y dim 
dim_y = y

#Create zero df of dim (x,y+1)
empty_row   = np.zeros((x, y))
empty_col   = np.zeros((x, y))
empty_fin   = np.zeros((x, y))

#Create zero vector of dim y+1
empty_array = np.zeros(dim_y)

#Create slope of decreasing nudging away from boundaries
x_lis       = np.linspace(0,nudgemax,nudge_cells)
nudge_array = amplitude * np.cos(x_lis * pi_factor * np.pi) + amplitude

#Apply to start/end of zero vec
empty_array[0:nudge_cells]              = nudge_array
empty_array[(dim_y-nudge_cells):dim_y]  = np.flip(nudge_array)

#Fill two df's with nudging vector (they are currently just correct in x & y direction) 
for row in range(0,x):
    if row < nudge_cells:
        empty_row[row,:] =  empty_array
        
        dummy = np.zeros(dim_y)
        dummy[:] = nudge_array[row]
        empty_col[row,:] = dummy
        empty_col[(x-row-1),:] = dummy
    
    else:
        empty_row[row,:] =  empty_array

#Combine the two arrays to get one with nudging coeffs on all 4 boundaries
for row in range(x):    #must be better way to combine two dataframes and chose the highest values
    for column in range(y):
        if empty_row[row,column]  > empty_col[row,column]:
            empty_fin[row,column] = empty_row[row,column]
        else:
            empty_fin[row,column] = empty_col[row,column]


#empty_fin[251:,11:391] = empty_fin[251:,11:391]*0.75  #reduce nudging in the north

if plot == 1: #plot
    plt.figure()
    plt.contourf(empty_fin,30)
    plt.colorbar()
    plt.show()  
    
#Create vertically unifrom nudging, in this preliminary script this is used for M3 coeff as well as temp & salt
empty_final_3d = np.array([empty_fin]*30)

#%%
                
if pNIJ == 1:
    empty_final_3d_momentum = empty_final_3d
    if res == 3:
        empty_final_3d_momentum[:13,50:59,0:nudge_cells] = empty_final_3d_momentum[:13,50:59,0:nudge_cells]*2
    elif res == 1.5:
        #empty_final_3d_momentum[:13,100:116,0:nudge_cells] = empty_final_3d_momentum[:13,100:116,0:nudge_cells]*2
        empty_final_3d_momentum[:13,100:118,0:nudge_cells] = empty_final_3d_momentum[:13,100:118,0:nudge_cells]*2
    elif res == 1:
        #empty_final_3d_momentum[:13,150:174,0:nudge_cells] = empty_final_3d_momentum[:13,150:174,0:nudge_cells]*2
        empty_final_3d_momentum[:13,150:177,0:nudge_cells] = empty_final_3d_momentum[:13,150:177,0:nudge_cells]*2
    
    
    X,Y = np.meshgrid(np.linspace(0,x*res,x),np.linspace(0,z,z))
    plt.figure()
    plt.contourf(X,Y,empty_final_3d_momentum[:,:,0])
    plt.colorbar()
    
    
    X,Y = np.meshgrid(np.linspace(0,y*res,y),np.linspace(0,x*res,x))
    plt.figure()
    plt.title('Velocity nudging coefficient in clm_nudge file')
    c1 = plt.contourf(X,Y,empty_final_3d_momentum[5,:,:])
    plt.xlabel('distance east [km]')
    plt.ylabel('distance north [km]')
    cbarax = plt.colorbar(c1)
    cbarax.set_label('nudging coefficient')
    
if pNIIC == 1:
    empty_final_3d_momentum = empty_final_3d
    empty_final_3d_momentum[0:30,50:95,0:nudge_cells] = empty_final_3d_momentum[0:30,50:95,0:nudge_cells]*2
    
    X,Y = np.meshgrid(np.linspace(0,y*res,y),np.linspace(0,z,z))
    plt.figure()
    plt.contourf(X,Y,empty_final_3d_momentum[:,:,0])
    plt.colorbar()
    
    
    X,Y = np.meshgrid(np.linspace(0,x*res,x),np.linspace(0,y*res,y))
    plt.figure()
    plt.title('Velocity nudging coefficient in clm_nudge file')
    c1 = plt.contourf(X,Y,empty_final_3d_momentum[5,:,:])
    plt.xlabel('distance east [km]')
    plt.ylabel('distance north [km]')
    cbarax = plt.colorbar(c1)
    cbarax.set_label('nudging coefficient')
    
if pNIIC2 == 1:
    empty_final_3d_momentum = empty_final_3d
    empty_final_3d_momentum[0:30,40:85,(y-nudge_cells):] = empty_final_3d_momentum[0:30,40:85,(y-nudge_cells):]*2
    
    X,Y = np.meshgrid(np.linspace(0,y*res,y),np.linspace(0,z,z))
    plt.figure()
    plt.contourf(X,Y,empty_final_3d_momentum[:,:,0])
    plt.colorbar()
    
    
    X,Y = np.meshgrid(np.linspace(0,x*res,x),np.linspace(0,y*res,y))
    plt.figure()
    plt.title('Velocity nudging coefficient in clm_nudge file')
    c1 = plt.contourf(X,Y,empty_final_3d_momentum[5,:,:])
    plt.xlabel('distance east [km]')
    plt.ylabel('distance north [km]')
    cbarax = plt.colorbar(c1)
    cbarax.set_label('nudging coefficient')
    
if pNIJ_md == 1:
    
    empty_final_3d_momentum = empty_final_3d 
    empty_final_3d_momentum[7:15,100:128,0:nudge_cells] = empty_final_3d_momentum[7:15,100:128,0:nudge_cells]*2
    
    X,Y = np.meshgrid(np.linspace(0,y*res,y),np.linspace(0,z,z))
    plt.figure()
    plt.contourf(X,Y,empty_final_3d_momentum[:,:,0])
    plt.colorbar()
    
    X,Y = np.meshgrid(np.linspace(0,x*res,x),np.linspace(0,y*res,y))
    plt.figure()
    plt.title('Velocity nudging coefficient in clm_nudge file')
    c1 = plt.contourf(X,Y,empty_final_3d_momentum[10,:,:])
    plt.xlabel('distance east [km]')
    plt.ylabel('distance north [km]')
    cbarax = plt.colorbar(c1)
    cbarax.set_label('nudging coefficient')
    
if pIFSJ == 1:
    empty_final_3d_momentum = empty_final_3d
    
    if onslope == 1:
        empty_final_3d_momentum[0:8,105:120,(y-nudge_cells):] = empty_final_3d_momentum[0:8,105:120,(y-nudge_cells):]*2
    elif onslope == 2:
        empty_final_3d_momentum[0:12,99:118,(y-nudge_cells):] = empty_final_3d_momentum[0:12,99:118,(y-nudge_cells):]*2
    else:
        empty_final_3d_momentum[0:8,115:128,(y-nudge_cells):] = empty_final_3d_momentum[0:8,115:128,(y-nudge_cells):]*2
        
        
    X,Y = np.meshgrid(np.linspace(0,y*res,y),np.linspace(0,z,z))
    plt.figure()
    plt.contourf(X,Y,empty_final_3d_momentum[:,:,-1])
    plt.colorbar()
    
    X,Y = np.meshgrid(np.linspace(0,x*res,x),np.linspace(0,y*res,y))
    plt.figure()
    plt.title('Velocity nudging coefficient in clm_nudge file')
    c1 = plt.contourf(X,Y,empty_final_3d_momentum[5,:,:])
    plt.xlabel('distance east [km]')
    plt.ylabel('distance north [km]')
    cbarax = plt.colorbar(c1)
    cbarax.set_label('nudging coefficient')
    
    
#%%
if create_file == 1:

    ###############################################################################
    ############################# - Create .nc file - #############################
    ###############################################################################
    
    path_string                     = directory+clm_nudge_file_name
    ncid                            = nc.Dataset(path_string,mode='w',format='NETCDF4')
    ncid.description                = Descrip_grd
    ncid.author                     = Author
    
    ###############################################################################
    ########################### - Create dimensions - #############################
    ###############################################################################
    
    eta_rho                         = ncid.createDimension('eta_rho',x)
    xi_rho                          = ncid.createDimension('xi_rho', y)
    s_rho                           = ncid.createDimension('s_rho',  z)
    one                             = ncid.createDimension('one',    1) 
    
    ###############################################################################
    ########################### - Create variables - ##############################
    ###############################################################################
    
    spherical_var                 = ncid.createVariable('spherical','i')
    spherical_var.long_name         = 'dgrid type logical switch'
    spherical_var.flag_values       = '[0 1]'
    spherical_var.flag_meaning      = 'Cartesian Spherical'
    
    lon_rho_var                   = ncid.createVariable('lon_rho','f8',('eta_rho', 'xi_rho',))
    lon_rho_var.long_name           = 'longitude of RHO-points'
    lon_rho_var.units               ='degree_east'
    lon_rho_var.standard_name       = 'longitude'
    
    lat_rho_var                   = ncid.createVariable('lat_rho','f8',('eta_rho', 'xi_rho',))
    lat_rho_var.long_name           = 'latitude of RHO-points'
    lat_rho_var.units               = 'degree_north'
    lat_rho_var.standard_name       = 'latitude'
    
    M2_NudgeCoef_var              = ncid.createVariable('M2_NudgeCoef','f8',('eta_rho', 'xi_rho',))
    M2_NudgeCoef_var.long_name      = '2D momentum inverse nudging coefficients'
    M2_NudgeCoef_var.units          = 'day-1'
    
    M3_NudgeCoef_var              = ncid.createVariable('M3_NudgeCoef','f8',('s_rho','eta_rho', 'xi_rho',))
    M3_NudgeCoef_var.long_name      = '3D momentum inverse nudging coefficients'
    M3_NudgeCoef_var.units          = 'day-1'
    
    temp_NudgeCoef_var            = ncid.createVariable('temp_NudgeCoef','f8',('s_rho','eta_rho', 'xi_rho',))
    temp_NudgeCoef_var.long_name    = 'temp inverse nudging coefficients'
    temp_NudgeCoef_var.units        = 'day-1'
    
    salt_NudgeCoef_var            = ncid.createVariable('salt_NudgeCoef','f8',('s_rho','eta_rho', 'xi_rho',))
    salt_NudgeCoef_var.long_name    = 'temp inverse nudging coefficients'
    salt_NudgeCoef_var.units        ='day-1'
    
    ###############################################################################
    ############################## - Assign values - ##############################
    ###############################################################################
    
    spherical_var[:]                    = 1
    lon_rho_var[:,:]                    = lon_rho_var1
    lat_rho_var[:,:]                    = lat_rho_var1
    
    M2_NudgeCoef_var[:,:]               = empty_fin
    
    M3_NudgeCoef_var[:,:,:]             = empty_final_3d
    temp_NudgeCoef_var[:,:,:]           = empty_final_3d
    salt_NudgeCoef_var[:,:,:]           = empty_final_3d
    
    ncid.close()
    
    

















