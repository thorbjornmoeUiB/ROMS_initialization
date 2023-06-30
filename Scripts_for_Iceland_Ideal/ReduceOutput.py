
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 18 16:56:38 2023

@author: thorbjornostenbymoe
"""
"""
Script for reducing size of output file, can surely be solved in the supercomputer 
with nc.dump, but this gives a lot of choices....

Interpolation from u and v points to rho points should probably be done in the horizontal
region part of the script    :)
"""

###############################################################################
################################# - Packages - ################################
###############################################################################

import numpy as np
import netCDF4 as nc

###############################################################################
############################### - User options - ##############################
###############################################################################
""" NB - nothing outside of user options should have to be changed unless a new variable is added"""

# file path - split up to fit with writing of new file
Main_path   = '/Users/thorbjornostenbymoe/Desktop/Initialization_of_model/MODEL_SETUP/NEW_OUTPUT_HHH'
Experiment  = '/MOD/pNIIC_NIJ/'
Filename    = 'roms_avg_0004.nc.nosync'

ds          = nc.Dataset(Main_path+Experiment+Filename)

#original dimensions
Orig_t = np.shape(ds['temp'])[0]
Orig_z = np.shape(ds['temp'])[1]
Orig_x = np.shape(ds['temp'])[2]
Orig_y = np.shape(ds['temp'])[3]

# Print dimensions along the way?
PrintDim = 1

#Variables to include (1 = include)
incl_Temp       = 1       # Temperature
incl_Salt       = 1       # Salinity
incl_Uvel       = 1       # Zonal velocity
incl_Vvel       = 1       # Meridional Velocity
incl_Ubar       = 1       # Depth-integrated zonal velocity
incl_Vbar       = 1       # Depth-integrated meridional velocity
incl_Zeta       = 1       # Sea surface height
incl_Wvel       = 1       # Vertical Velocity
incl_Omeg       = 1       # Vertical momentum
incl_shfx       = 1       # surface heat flux
incl_ssfx       = 1       # surface salt flux
only_time_3D    = 0       # Only do horizontal averages for 4D variables, 3D variables are only treated in time 
# if new variables are desired, simply add here and below according to dimension of variable

#recomended to always include
incl_Cs_r       = 1       # rho-coordinate stretching
incl_Cs_w       = 1       # W-coordinate stretching
incl_dept       = 1       # depth array

# Time averages?
Time_avg        = 0       # 1 = time average
Time_a_start    = 6       # When does the time average start? (month) 
Time_a_end      = 12      # When does the time average end? (month)

# Specific period? 
# (If only one time-step is desired, simply set Time_p_start = Time_p_end)
Time_per        = 1       # 1 = specific period
Time_p_start    = 6       # When does the time period start? (month) 
Time_p_end      = 12      # When does the time period end? (month)
time_p_dim      = Time_p_end - Time_p_start #dimension of newly defined period


# Horizontal average? 
Horiz_avg       = 0       # 1 = Horizontal average
Horiz_avg_d     = 'x'     # 'x' = zonal average, 'y' = meridional average
Horiz_avg_s     = 125     # Where does the horizontal average start?
Horiz_avg_e     = 175     # Where does the horizontal average end?
Horiz_avg_dim   = Horiz_avg_e - Horiz_avg_s

# Horizontal region? (retains only a specified region)
Horiz_reg       = 1       # 1 = Horizontal average
Horiz_reg_both  = 1       # 1 = average over both dims
Horiz_reg_dim   = 'y'     # 'x' = zonal average, 'y' = meridional average
#  1 is for zonal 
Horiz_reg_s1    = 125     # Where does the horizontal region start?
Horiz_reg_e1    = 175     # Where does the horizontal region end?
Horiz_reg_dim1  = Horiz_reg_e1 - Horiz_reg_s1
# 2 is for meridional 
Horiz_reg_s2    = 100     # Where does the horizontal region start?
Horiz_reg_e2    = 200     # Where does the horizontal region end?
Horiz_reg_dim2  = Horiz_reg_e2 - Horiz_reg_s2

###############################################################################
################################ - Functions - ################################
###############################################################################

#function to get desired data, should work for all dimensions and output (as long as the dimension is 3D or 4D)
def get_avgs_4D(VARIABLE,name,counter):
    orig_z = np.shape(VARIABLE)[1] #original dimensions
    orig_x = np.shape(VARIABLE)[2]
    orig_y = np.shape(VARIABLE)[3]
    
    ###################### - TIME - ######################
    if Time_avg == 1: # if user wants a time average, start and end is defined in 'User Options'
        VARIABLE = np.array([np.nanmean(VARIABLE[Time_a_start:Time_a_end,:,:,:],axis=0)]) # np.array to fit NetCDF file
        tdim = 1
        
    elif Time_per == 1: # if user wants a certain time period, start and end is defined in 'User Options'
        if Time_p_start == Time_p_end: # if start = end, only one timestep is retrieved
            VARIABLE = np.array([VARIABLE[Time_p_start,:,:,:]])
            tdim = 1
            
        else: # else a time period is chosen
            VARIABLE = VARIABLE[Time_p_start:Time_p_end,:,:,:]
            tdim = time_p_dim
           
    else: #get full data period
        VARIABLE = VARIABLE
        tdim = Orig_t
    
    ###################### - SPACE - ######################
    if Horiz_avg == 1: # retrieve horizontal transect averaged over chosen extent, zonal or meridional
        if Horiz_avg_d == 'x': # zonal average
            VARIABLE = np.reshape(np.array([np.nanmean(VARIABLE[:,:,Horiz_avg_s:Horiz_avg_e,:],axis=2)]),newshape=[tdim,orig_z,1,orig_y]) # Reshape to fit NetCDF file
        elif Horiz_avg_d == 'y': # meridional average 
            VARIABLE = np.reshape(np.array([np.nanmean(VARIABLE[:,:,:,Horiz_avg_s:Horiz_avg_e],axis=3)]),newshape=[tdim,orig_z,orig_x,1]) # Reshape to fit NetCDF file

    elif Horiz_reg == 1: # retrieve only a chosen region of the domain, could be used to for instance exclude nudging regions
        if Horiz_reg_both == 1: # both zonal and meridional selected region
            VARIABLE = VARIABLE[:,:,Horiz_reg_s1:Horiz_reg_e1,Horiz_reg_s2:Horiz_reg_e2]
        else:
            if Horiz_reg_dim == 'x': # only zonal region
                VARIABLE = VARIABLE[:,:,Horiz_reg_s1:Horiz_reg_e1,:]
            if Horiz_reg_dim == 'y': # only meridional region
                VARIABLE = VARIABLE[:,:,:,Horiz_reg_s2:Horiz_reg_e2]
                
    if PrintDim == 1: #print dimension
        print(str(counter)+'. '+name+'(t,z,x,y): '+str(np.shape(VARIABLE)))
        counter = counter + 1
    return VARIABLE,counter

# same as last function for 3D variables
def get_avgs_3D(VARIABLE,name,counter):
    orig_x = np.shape(VARIABLE)[1]
    orig_y = np.shape(VARIABLE)[2]
    if   Time_avg == 1:
        VARIABLE = np.array([np.nanmean(VARIABLE[Time_a_start:Time_a_end,:,:],axis=0)])
        tdim = 1
    elif Time_per == 1:
        if Time_p_start == Time_p_end: 
            VARIABLE =  np.array([VARIABLE[Time_p_start,:,:]])
            tdim = 1
        else:
            VARIABLE = VARIABLE[Time_p_start:Time_p_end,:,:]
            tdim = time_p_dim
    else:
        VARIABLE = VARIABLE
        tdim = Orig_t

    #space
    if Horiz_avg == 1:
        if Horiz_avg_d == 'x':
            VARIABLE = np.reshape(np.array([np.nanmean(VARIABLE[:,Horiz_avg_s:Horiz_avg_e,:],axis=1)]),newshape=[tdim,1,orig_y])
        elif Horiz_avg_d == 'y':
            VARIABLE = np.reshape(np.array([np.nanmean(VARIABLE[:,:,Horiz_avg_s:Horiz_avg_e],axis=2)]),newshape=[tdim,orig_x,1])

    elif Horiz_reg == 1:
        if Horiz_reg_both == 1:
            VARIABLE = VARIABLE[:,Horiz_reg_s1:Horiz_reg_e1,Horiz_reg_s2:Horiz_reg_e2]
        else:
            if Horiz_reg_dim == 'x':
                VARIABLE = VARIABLE[:,Horiz_reg_s1:Horiz_reg_e1,:]
            if Horiz_reg_dim == 'y':
                VARIABLE = VARIABLE[:,:,Horiz_reg_s2:Horiz_reg_e2]
                
    if PrintDim == 1:
        print(str(counter)+'. '+name+'(t,x,y):   '+str(np.shape(VARIABLE)))
        counter = counter + 1
    return VARIABLE,counter


def get_avgs_3D_time(VARIABLE,name,counter): # same as previous but only for time
    if   Time_avg == 1:
        VARIABLE = np.array([np.nanmean(VARIABLE[Time_a_start:Time_a_end,:,:],axis=0)])
    elif Time_per == 1:
        if Time_p_start == Time_p_end: 
            VARIABLE =  np.array([VARIABLE[Time_p_start,:,:]])
        else:
            VARIABLE = VARIABLE[Time_p_start:Time_p_end,:,:]
    else:
        VARIABLE = VARIABLE
                
    if PrintDim == 1:
        print(str(counter)+'. '+name+'(t,x,y):   '+str(np.shape(VARIABLE)))
        counter = counter + 1
    return VARIABLE,counter


###############################################################################
################################### - file - ##################################
###############################################################################

path_string         = Main_path+Experiment+'Reduced_'+Filename
ncid                = nc.Dataset(path_string,mode='w',format='NETCDF4')
ncid.description    = 'Reduced output of file: '+Main_path+Experiment+Filename
ncid.author         = 'Thorbjoern Oestenby Moe'

###############################################################################
############################# - Set dimensions - ##############################
###############################################################################

# Original dimensions 
# Time
tdim = Orig_t
# Depth
z1   = ds['Cs_r'][:].shape[0]
z2   = ds['Cs_w'][:].shape[0]
# Horizontal
x1   = ds['v'][:].shape[2]
x2   = ds['u'][:].shape[2]
y1   = ds['u'][:].shape[3]
y2   = ds['v'][:].shape[3]
# Dim to use if a dimension is averaged over (to keep the shape constant)
dim0 = 1

# Setting dimensions of NetCDF file according to the desired data structure
# This means keeping the original dimensions for variables or dimensions not changed
# and reducing dimensions when neccessary

if Time_avg == 1 or Time_p_start == Time_p_end:
    tdim = dim0
elif Time_per == 1:
    tdim = time_p_dim
else:
    tdim = Orig_t

if Horiz_avg == 1 and Horiz_avg_d == 'y':
    x1 = x1
    x2 = x2
    y1 = dim0
    y2 = dim0
    
elif Horiz_avg == 1 and Horiz_avg_d == 'x':
    x1 = dim0
    x2 = dim0
    y1 = y1
    y2 = y2   

elif Horiz_reg == 1 and Horiz_reg_both == 0 and Horiz_reg_dim == 'y':
    x1 = x1
    x2 = x2
    y1 = Horiz_reg_dim2
    y2 = Horiz_reg_dim2
    
elif Horiz_reg == 1 and Horiz_reg_both == 0 and Horiz_reg_dim == 'x':
    x1 = Horiz_reg_dim1
    x2 = Horiz_reg_dim1
    y1 = y1
    y2 = y2

elif Horiz_reg == 1 and Horiz_reg_both == 1:
    x1 = Horiz_reg_dim1
    x2 = Horiz_reg_dim1
    y1 = Horiz_reg_dim2
    y2 = Horiz_reg_dim2
    
###############################################################################
############################# - Make dimensions - #############################
###############################################################################
    
dt          = ncid.createDimension('dt',  tdim)
dz1         = ncid.createDimension('dz1', z1)
dz2         = ncid.createDimension('dz2', z2)
dx1         = ncid.createDimension('dx1', x1)
dx2         = ncid.createDimension('dx2', x2)
dy1         = ncid.createDimension('dy1', y1)
dy2         = ncid.createDimension('dy2', y2)

hdx         = ncid.createDimension('hdx', 268)
hdy         = ncid.createDimension('hdy', 402) 

if only_time_3D == 1:
    x3   = ds['v'][:].shape[2]
    x4   = ds['u'][:].shape[2]
    y3   = ds['u'][:].shape[3]
    y4   = ds['v'][:].shape[3]
    
    dx3         = ncid.createDimension('dx3', x3)
    dx4         = ncid.createDimension('dx4', x4)
    dy3         = ncid.createDimension('dy3', y3)
    dy4         = ncid.createDimension('dy4', y4)

###############################################################################
######### - Get data, average, define variable, and assign variable - #########
###############################################################################

#Template to add another variable:
"""
4D:
if incl_NEW_VAR == 1:
    # Get data
    NEW_VAR_raw = ds['new_var'][:,:,:,:]
    # Make fields
    NEW_VAR,counter    = get_avgs_4D(NEW_VAR_raw,'NEW_VAR short name',counter)
    
    # Create variable:
    #'dt': should be constant
    #'dz1': change to dz2 if variable is located on vertical w-point (e.g. W, omega,..) 
    #'dx2': change to dx1 if variable is located on eta_vi-point
    #'dy2': change to dy1 if variable is located on xi_u-point 
    
    new_var            = ncid.createVariable('new_var','f8',('dt','dz1','dx2','dy2',))
    new_var.long_name  = 'NEW_VAR long name'
    # Assign variable
    new_var[:,:,:,:]   = NEW_VAR
  
3D:
if incl_NEW_VAR == 1:
    NEW_VAR_raw = ds['new_var'][:,:,:]
    
    if only_time_3D == 0:
        NEW_VAR,counter    = get_avgs_3D(NEW_VAR_raw,NEW_VAR short name,counter) 
        
        new_var            = ncid.createVariable('new_var','f8',('dt','dx2','dy1',))
        new_var.long_name  = 'NEW_VAR long name'
        new_var[:,:,:]     = NEW_VAR
    
    elif only_time_3D == 1:
        NEW_VAR,counter    = get_avgs_3D_time(NEW_VAR_raw,NEW_VAR short name,counter)
        
        new_var            = ncid.createVariable('new_var','f8',('dt','dx4','dy3',))
        new_var.long_name  = 'NEW_VAR long name'
        new_var[:,:,:]     = NEW_VAR
"""

if PrintDim == 1:
    print('#####################################################################')
    print('################ - Variables included & dimension - #################')
    print('##################################################################### \n')
counter = 1
if incl_Temp == 1:
    # Get data
    TEMP_raw = ds['temp'][:,:,:,:]
    # Make fields
    TEMP,counter    = get_avgs_4D(TEMP_raw,'T',counter)
    
    # Create variable
    temp            = ncid.createVariable('temp','f8',('dt','dz1','dx2','dy2',))
    temp.long_name  = 'Temperature'
    # Assign variable
    temp[:,:,:,:]   = TEMP
       
if incl_Salt == 1: #same comments for all others
    SALT_raw = ds['salt'][:,:,:,:]
    SALT,counter    = get_avgs_4D(SALT_raw,'S',counter)
    
    salt            = ncid.createVariable('salt','f8',('dt','dz1','dx2','dy2',))
    salt.long_name  = 'Salinity'
    salt[:,:,:,:]   = SALT
  
    
if incl_Uvel == 1:
    UVEL_raw = ds['u'][:,:,:,:]
    UVEL,counter    = get_avgs_4D(UVEL_raw,'u',counter)
    
    u               = ncid.createVariable('u','f8',('dt','dz1','dx2','dy1',))
    u.long_name     = 'Zonal Velocity'
    u[:,:,:,:]      = UVEL

if incl_Vvel == 1:
    VVEL_raw = ds['v'][:,:,:,:]
    VVEL,counter    = get_avgs_4D(VVEL_raw,'v',counter)
    
    v               = ncid.createVariable('v','f8',('dt','dz1','dx1','dy2',))
    v.long_name     = 'Meridional Velocity'
    v[:,:,:,:]      = VVEL
    
if incl_Wvel == 1:
    WVEL_raw = ds['w'][:,:,:,:]
    WVEL,counter    = get_avgs_4D(WVEL_raw,'w',counter)
    
    w               = ncid.createVariable('w','f8',('dt','dz2','dx2','dy2',))
    w.long_name     = 'Vertical Velocity'
    w[:,:,:,:]      = WVEL
    
if incl_Omeg == 1:
    OMEG_raw = ds['omega'][:,:,:,:]
    OMEG,counter    = get_avgs_4D(OMEG_raw,'ω',counter)
    
    omega           = ncid.createVariable('omega','f8',('dt','dz2','dx2','dy2',))
    omega.long_name = 'Vertical Momentum'
    omega[:,:,:,:]  = OMEG

if incl_Ubar == 1:
    UBAR_raw = ds['ubar'][:,:,:]
    
    if only_time_3D == 0:
        UBAR,counter    = get_avgs_3D(UBAR_raw,u'u\u0305',counter) 
        
        ubar            = ncid.createVariable('ubar','f8',('dt','dx2','dy1',))
        ubar.long_name  = 'Depth-Integrated Zonal Velocity'
        ubar[:,:,:]     = UBAR
    
    elif only_time_3D == 1:
        UBAR,counter    = get_avgs_3D_time(UBAR_raw,u'u\u0305',counter)
        
        ubar            = ncid.createVariable('ubar','f8',('dt','dx4','dy3',))
        ubar.long_name  = 'Depth-Integrated Zonal Velocity'
        ubar[:,:,:]     = UBAR
    
if incl_Vbar == 1:
    VBAR_raw = ds['vbar'][:,:,:]
    
    if only_time_3D == 0:
        VBAR,counter    = get_avgs_3D(VBAR_raw,u'v\u0305',counter)  
        
        vbar            = ncid.createVariable('vbar','f8',('dt','dx1','dy2',))
        vbar.long_name  = 'Depth-Integrated Meridional Velocity'
        vbar[:,:,:]     = VBAR
        
    elif only_time_3D == 1:  
        VBAR,counter    = get_avgs_3D_time(VBAR_raw,u'v\u0305',counter)  
        
        vbar            = ncid.createVariable('vbar','f8',('dt','dx3','dy4',))
        vbar.long_name  = 'Depth-Integrated Meridional Velocity'
        vbar[:,:,:]     = VBAR
    
if incl_Zeta == 1:
    ZETA_raw = ds['zeta'][:,:,:]
    
    if only_time_3D == 0:
        ZETA,counter    = get_avgs_3D(ZETA_raw,'η',counter) 
        
        zeta            = ncid.createVariable('zeta','f8',('dt','dx2','dy2',))
        zeta.long_name  = 'Sea Surface Elevation'
        zeta[:,:,:]     = ZETA
        
    elif only_time_3D == 1:  
        ZETA,counter    = get_avgs_3D_time(ZETA_raw,'η',counter) 
        
        zeta            = ncid.createVariable('zeta','f8',('dt','dx4','dy4',))
        zeta.long_name  = 'Sea Surface Elevation'
        zeta[:,:,:]     = ZETA
        
if incl_shfx == 1:
    SHFX_raw = ds['shflux'][:,:,:]
    
    if only_time_3D == 0:
        SHFX,counter      = get_avgs_3D(SHFX_raw,'dQ',counter)
        
        shflux            = ncid.createVariable('shflux','f8',('dt','dx2','dy2',))
        shflux.long_name  = 'heat flux'
        shflux[:,:,:]     = SHFX
        
    elif only_time_3D == 1:
        SHFX,counter      = get_avgs_3D_time(SHFX_raw,'dQ',counter)
        
        shflux            = ncid.createVariable('shflux','f8',('dt','dx4','dy4',))
        shflux.long_name  = 'heat flux'
        shflux[:,:,:]     = SHFX
    
if incl_ssfx == 1:
    SSFX_raw = ds['ssflux'][:,:,:]
    if only_time_3D == 0:
            SSFX,counter      = get_avgs_3D(SSFX_raw,'dS',counter) 
            
            ssflux            = ncid.createVariable('ssflux','f8',('dt','dx2','dy2',))
            ssflux.long_name  = 'salt flux'
            ssflux[:,:,:]     = SSFX
            
    elif only_time_3D == 1:       
            SSFX,counter      = get_avgs_3D_time(SSFX_raw,'dS',counter) 
            
            ssflux            = ncid.createVariable('ssflux','f8',('dt','dx4','dy4',))
            ssflux.long_name  = 'salt flux'
            ssflux[:,:,:]     = SSFX
    
if incl_Cs_r == 1:
    Cs_r            = ncid.createVariable('Cs_r','f8',('dz1',))
    Cs_r.long_name  = 'stretching on rho points'
    Cs_r[:]         = ds['Cs_r'][:]
    
    
if incl_Cs_w == 1:
    Cs_w            = ncid.createVariable('Cs_w','f8',('dz2',))
    Cs_w.long_name  = 'stretching on w points'
    Cs_w[:]         = ds['Cs_w'][:]

if incl_dept == 1:
    h            = ncid.createVariable('h','f8',('hdx','hdy',))
    h.long_name  = 'salt flux'
    h[:,:]       = ds['h'][:,:]


ncid.close()














