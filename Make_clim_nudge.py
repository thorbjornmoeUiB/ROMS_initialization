#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 14 15:05:07 2023

@author: thorbjorn
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
from GlobParam import *                # This is where I change user outputs now

print('###############################################################################')
print('########################## - Making clm nudge file - ##########################')
print('###############################################################################')

###############################################################################
############################### - User options - ##############################
###############################################################################

#File directory/names
common_dir          = GlobParam.common_dir
directory           = GlobParam.EXP_dir

clm_nudge_file_name = directory + GlobParam.clm_nudge_file
grid_file           = directory + GlobParam.grid_file
bry_file            = directory + GlobParam.bry_file

grid_ds             = nc.Dataset(grid_file)
bnry_ds             = nc.Dataset(bry_file)

Descrip_grd         = GlobParam.General_desc + GlobParam.spec_desc_nud
Author              = GlobParam.Author
create_file         = GlobParam.Create_files
plotting            = GlobParam.Gen_plot
Iplottin            = GlobParam.Imp_plot
z                   = GlobParam.z
res                 = GlobParam.Res

pIFSJ               = GlobParam.pIFSJ
pNIIC_w             = GlobParam.pNIIC_w

#Dim and plot switch
x                   = np.shape(grid_ds['h'])[0]
y                   = np.shape(grid_ds['h'])[1]

#Trigonometric coefficients for the nudging function
nudgemax            = 0.2                                      # Max nudging coefficient, taken from previous nudging files (other values could be used)
amplitude           = nudgemax/2
pi_factor           = 5

#nudge zone 
nudge_cells     = 17  
#%% 

###############################################################################
############################## - Load lat/lon - ###############################
###############################################################################

lat_rho_var1 = grid_ds['lat_rho'][:,:] #using standard ROMS format
lon_rho_var1 = grid_ds['lon_rho'][:,:]

###############################################################################
########################## - Create M2/M3_NudgeCoef - #########################
###############################################################################

#Lon_rho y dim 
dim_y = y

#Create zero df of dim (x,y+1)
empty_row   = np.zeros((x, y))
empty_col   = np.zeros((x, y))
M2_NC   = np.zeros((x, y))

#Create zero vector of dim y+1
empty_array = np.zeros(dim_y)

#Create slope of decreasing nudging away from boundaries
x_lis       = np.linspace(0,nudgemax,nudge_cells)
nudge_array = amplitude * np.cos(x_lis * pi_factor * np.pi) + amplitude

#Apply to start/end of zero vec
empty_array[0:nudge_cells]              = nudge_array
empty_array[(dim_y-nudge_cells):]       = nudge_array[::-1]

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
for row in range(x):    #must be better way to combine two sets and chose the highest values
    for column in range(y):
        if empty_row[row,column]  > empty_col[row,column]:
            M2_NC[row,column] = empty_row[row,column]
        else:
            M2_NC[row,column] = empty_col[row,column]

#Create vertically unifrom nudging, in this preliminary script this is used for M3 coeff as well as temp & salt
M3_NC = M2_NC[np.newaxis,:,:].repeat(30,0)
M3_NC_nonedit = M3_NC.copy()

if plotting == 1: #plot
    plt.figure(figsize=[10,3])
    plt.suptitle('FROM Make_clm_nudge.py')
    plt.contourf(M3_NC[0,:,:],30)
    plt.colorbar()
    plt.show()  
    
#%%
# test boudnary file with temp and velocity
# bnry_ds = nc.Dataset('.../Iceland_Sea_Mod/NewConfig/INPUT/VarSlo_WINTER_2xDomain/pNIIC_6C_NIJ/bry_file.nc')

###############################################################################
############### - Increase nudging where current is prescribed - ##############
###############################################################################

# Load properties from bry 
u_west, u_east = bnry_ds['u_west'][:], bnry_ds['u_east'][:]

# Make mask for where there are modifications to velocity
# 1. It only doubles the nudging where a current is prescribed
# 2. Doesnt matter with temperature because the two regions overlap

bottom_mask_west = np.where(u_west[0,0,:] < 0,2,1)
middep_mask_west = np.where(u_west[0,16,:] != 0,2,1)

bottom_mask_east = np.where(u_east[0,0,:] < 0,2,1)
middep_mask_east = np.where(u_east[0,16,:] != 0,2,1)

# Apply mask to nudging
# only to level 15 for bottom mask because this is the shallowest level the NIJ is prescribed

M3_NC[:,:,:nudge_cells] = M3_NC[:,:,:nudge_cells]*middep_mask_west[np.newaxis,:,np.newaxis]
M3_NC[:15,:,:nudge_cells] = M3_NC[:15,:,:nudge_cells]*bottom_mask_west[np.newaxis,:,np.newaxis]


M3_NC[:,:,(y-nudge_cells):] = M3_NC[:,:,(y-nudge_cells):]*middep_mask_east[np.newaxis,:,np.newaxis]
M3_NC[:15,:,(y-nudge_cells):] = M3_NC[:15,:,(y-nudge_cells):]*bottom_mask_east[np.newaxis,:,np.newaxis]

# Plot for safety

plt.figure()
plt.contourf(M3_NC[0,:,:])
plt.colorbar()

###############################################################################
########################### - same for 3d nudging - ###########################
###############################################################################

mask_west = np.where(u_west[0,0,:] != 0,2,1)
mask_east = np.where(u_east[0,0,:] != 0,2,1)

M2_NC[:,:nudge_cells]     = M2_NC[:,:nudge_cells]*mask_west[:,np.newaxis]
M2_NC[:,(y-nudge_cells):] = M2_NC[:,(y-nudge_cells):]*mask_east[:,np.newaxis]

plt.figure()
plt.contourf(M2_NC[:,:])
plt.colorbar()
#%%
if create_file == 1:

    ###############################################################################
    ############################# - Create .nc file - #############################
    ###############################################################################
    
    ncid                            = nc.Dataset(clm_nudge_file_name,mode='w',format='NETCDF4')
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
    
    M2_NudgeCoef              = ncid.createVariable('M2_NudgeCoef','f8',('eta_rho', 'xi_rho',))
    M2_NudgeCoef.long_name      = '2D momentum inverse nudging coefficients'
    M2_NudgeCoef.units          = 'day-1'
    
    M3_NudgeCoef              = ncid.createVariable('M3_NudgeCoef','f8',('s_rho','eta_rho', 'xi_rho',))
    M3_NudgeCoef.long_name      = '3D momentum inverse nudging coefficients'
    M3_NudgeCoef.units          = 'day-1'
    
    temp_NudgeCoef            = ncid.createVariable('temp_NudgeCoef','f8',('s_rho','eta_rho', 'xi_rho',))
    temp_NudgeCoef.long_name    = 'temp inverse nudging coefficients'
    temp_NudgeCoef.units        = 'day-1'
    
    salt_NudgeCoef            = ncid.createVariable('salt_NudgeCoef','f8',('s_rho','eta_rho', 'xi_rho',))
    salt_NudgeCoef.long_name    = 'temp inverse nudging coefficients'
    salt_NudgeCoef.units        ='day-1'
    
    ###############################################################################
    ############################## - Assign values - ##############################
    ###############################################################################
    
    spherical_var[:]                    = 1
    lon_rho_var[:,:]                    = lon_rho_var1
    lat_rho_var[:,:]                    = lat_rho_var1
    
    M2_NudgeCoef[:,:]               = M2_NC
    M3_NudgeCoef[:,:,:]             = M3_NC
    salt_NudgeCoef[:,:,:]           = M3_NC_nonedit
    
    if pNIIC_w == 1:
        temp_NudgeCoef[:,:,:]       = M3_NC
    else:
        temp_NudgeCoef[:,:,:]       = M3_NC_nonedit
    
    ncid.close()
    
    


