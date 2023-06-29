#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 19 16:37:49 2022

@author: thorbjornostenbymoe
"""

""" 
This script creates a grid file for an idealized ROMS test run. 

The script should be modifiable enough from the user options section. 

Large modifcations are needed to the bathymetry and masking sections if y dimension is not constant. 

Many grid parameters not needed for the preliminary test run are set to illogical values or copied 
from a more complicated ROMS simulation. 

For later use it should be sufficient to change the variables, as the writing to netCDF is 
made in accordance to the ROMS standards.
"""

###############################################################################
################################ - packages - #################################
###############################################################################

import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np
import gsw as gsw

###############################################################################
############################### - User options - ##############################
###############################################################################

#File directory/names/details
directory           = '/Users/thorbjornostenbymoe/Desktop/Initialization_of_model/MODEL_SETUP/mach4/' #User directory/path
grid_file_name      = 'grid_file.nc' 

#file description
Descrip_grd         = 'Iceland slope - ROMS idealized model run - Grid file'  
Author              = 'Thorbjoern Oestenby Moe'

#logical switches. 1 = True, 0 = False
create_grid         = 1           # Create GRID. Turn OFF to work with a previously created grid (i.e. grid variables existing on Workspace)
plot_grid           = 1           # Plot grid
save_grid           = 1           # Save GRID in a NetCDF file
plot_mask           = 0           # plot masks
plotting            = 1

# Geographical and Grid parameters --------
lat         = 66                      # Latitude  (degrees) of the bottom-left corner of the grid.
lon         = -16                     # Longitude (degrees) of the bottom-left corner of the grid. 

rotangle    = 0                       # Angle (degrees) to rotate the grid conterclock-wise
res         = 1.5                     # km per grid cell width in east and north direction 
resol       = res * 1000              # Horizontal cell width (i.e. Resolution)in meters. Grid cells are forced to be (almost) square.
N           = 30                      # Number of vertical levels
X           = 600000+resol            # Width of domain (meters)
Y           = 400000+resol            # Length of domain (meters)    
x           = np.arange(0,X+1,resol) 
y           = np.arange(0,Y,resol)

hh          = 1500                   # Analytical Depth (meters) used to create a uniform-depth grid. If using a bathymetry file, leave hh = nan;
minh        = 4                      # Arbitrary depth offset in meters (see above). minh should be a little more than the maximum expected tidal variation.

# Trigonometric coefficients for creation of bathymetry
amplitude       = -500
offset          = -651
phase           = -150 
steepnes_factor = 15  #(7 er steep slope, 15 er modslope, 25 er mild slope) (11 er meep (steep+mod)/2)


#Trigonometric coefficients for sponge region
nudgemax        = 3.2                                                    #Max nudging coefficient, taken from previous nudging files (other values could be used)
amplitude       = nudgemax/2
pi_factor       = 0.3
nudge_cells     = 10                                                     #Amount of cells with nudging (away from boundary)
nudge_cells_s   = 0

#plotting
plot_depth      = -1200

#masking
mask_n          = 5 #number of cells to mask at southern boundary of domain

###############################################################################
############################# - Make bathymetry - #############################
###############################################################################

h  = np.zeros([y.shape[0],x.shape[0]])
i  = np.arange(0,y.shape[0],1)

#actual (meridional) bathymetry function (fit to bathymetry northeast of Iceland)
hy = offset+amplitude*((np.tanh((i*res+phase)/steepnes_factor)))  #works for 400X200 km grid

#spread bathymetry across zonal domain
for ind in np.arange(0,x.shape[0],1):
    h[:,ind]=-hy


if plot_grid==1:
    plotting_y = np.arange(0, y.shape[0],1)*res
    fig, ax= plt.subplots()
    ax.plot(plotting_y,hy,color='gray')
    ax.fill_between(plotting_y,np.linspace(plot_depth,plot_depth,len(plotting_y)),hy,color='gray')
    ax.fill_between(plotting_y,np.linspace(0,0,len(plotting_y)),hy,color='lightblue',alpha=0.6)
    ax.set_xlabel('distance north [km]')
    ax.set_ylabel('Depth [m]')
    ax.set_xlim([0,400])
    ax.set_ylim([plot_depth,5])
    ax.set_title('Transect')
    ax.grid(alpha=0.3)
    
    xcord = np.arange(0,x.shape[0],1)*res
    ycord = np.arange(0,y.shape[0],1)*res
    [Xcord,Ycord] = np.meshgrid(xcord,ycord)
    
    fig, ax= plt.subplots(figsize=[7,4])
    cs = ax.pcolormesh(Xcord,Ycord,h)
    cbar=fig.colorbar(cs)
    cbar.set_label('Depth [m]')
    ax.set_title('Depth')
    ax.set_ylabel('Latitude [km]')
    ax.set_xlabel('Longitude [km]')

###############################################################################
################################ - Make masks - ###############################
###############################################################################

x_dim = np.shape(x)[0]
y_dim = np.shape(y)[0]

mask_rho = np.ones([y_dim, x_dim])
mask_rho[0:mask_n,:] = 0

mask_psi = np.ones([y_dim-1, x_dim-1])
mask_psi[0:mask_n,:] = 0

mask_u = np.ones([y_dim, x_dim-1])
mask_u[0:mask_n,:] = 0

mask_v = np.ones([y_dim-1, x_dim])
mask_v[0:mask_n,:] = 0

if plot_mask==1:
    plt.figure()
    plt.contourf(mask_v,1,vmax=1)
    plt.colorbar()
    plt.ylabel('Latitude [km]')
    plt.xlabel('Longitude [km]')
    plt.title('masks (1=ocean)')
    plt.show()

###############################################################################
############################ - Spone zone - ###################################
###############################################################################

#Lon_rho y dim 
dim_y = y.shape[0]

#Create zero array of dim (x,y+1)
empty_row = np.zeros((x.shape[0], y.shape[0]))
empty_col = np.zeros((x.shape[0], y.shape[0]))
empty_fin = np.zeros((x.shape[0], y.shape[0]))

#Create zero vector of dim y+1
empty_array = np.zeros(dim_y)

#Create slope of decreasing nudging away from boundaries
x_lis = np.linspace(0,nudgemax,nudge_cells)
nudge_array = amplitude * np.cos(x_lis * pi_factor * np.pi) + amplitude

x_lis = np.linspace(0,nudgemax,nudge_cells_s)
nudge_array_s = amplitude * np.cos(x_lis * pi_factor * np.pi) + amplitude

if plotting == 1:
    plt.figure()
    plt.plot(nudge_array)

#Apply to start/end of zero vec
empty_array[:nudge_cells_s] = nudgemax
empty_array[nudge_cells_s:nudge_cells_s+nudge_cells] = nudge_array
empty_array[(dim_y-nudge_cells):dim_y] = np.flip(nudge_array)

if plotting == 1:
    plt.figure()
    plt.plot(empty_array)

#Fill two df's with nudging vector (they are currently just correct in x & y direction) 
for row in range(x.shape[0]):
    if row < nudge_cells:
        empty_row[row,:] =  empty_array[:]
        
        dummy = np.zeros(dim_y)
        dummy[:] = nudge_array[row]
        empty_col[row,:] = dummy
        empty_col[(x.shape[0]-row-1),:] = dummy
    
    else:
        empty_row[row,:] =  empty_array[:]
        
#Combine the two df's to get one with nudging coeffs on all 4 boundaries
for row in range(x.shape[0]):    #must be better way to combine two dataframes and chose the highest values
    for column in range(y.shape[0]):
        if empty_row[row,column] > empty_col[row,column]:
            empty_fin[row,column] = empty_row[row,column]
        else:
            empty_fin[row,column] = empty_col[row,column]

if plotting == 1:
    plt.figure()
    plt.contourf(empty_fin.T,30)
    plt.colorbar()
    plt.show()

###############################################################################
############################### - Conversions - ###############################
###############################################################################

rotangle = rotangle/180*np.pi                 # Convert Angle for grid rotation from degrees to radians
latdist  = gsw.distance([lon,lon],[lat,lat+1]) # Length (in meters) of 1 degree of latitude

###############################################################################
############################### - Create grid - ###############################
###############################################################################

if create_grid == 1: # Only create grid if switch (in USER SETTINGS) is turned ON

    Lm   = x.shape[0]-2
    Lp   = Lm +2
    L    = Lp-1
    Mm   = y.shape[0]-2
    Mp   = Mm +2
    M    = Mp-1
    
    ################## ----- RHO GRID ----- ##################
    #Create non-georeferenced grid in meters (origin = 0,0)
    x_rho = np.tile(x,[y.shape[0],1])
    y_rho = np.transpose(np.tile(y,[x.shape[0],1]))
    
    #Rotate grid ...See: http://en.wikipedia.org/wiki/Rotation_(mathematics)
    Rx_rho = x_rho * np.cos(rotangle) - y_rho * np.sin(rotangle)
    Ry_rho = x_rho * np.sin(rotangle) + y_rho * np.cos(rotangle)
    
    #Estimate Latitude and Longitude of each Grid point
    lat_rho = np.zeros(x_rho.shape) #initialize
    lon_rho = np.zeros(x_rho.shape) #initialize
    lat_rho[0,0] = lat + (Ry_rho[0,0] / latdist)
    lon_rho[0,0] = lon + (Rx_rho[0,0] / gsw.distance([lon,lon+1],[lat_rho[0,0],lat_rho[0,0]]))
    
    for i in np.arange(1,Ry_rho.shape[0],1):
        lat_rho[i,0] = lat + (Ry_rho[i,0] / latdist)
        lon_rho[i,0] = lon_rho[i-1,0] + ((Rx_rho[i,0]-Rx_rho[i-1,0]) / gsw.distance([lon_rho[i-1,0],lon_rho[i-1,0]+1],[lat_rho[i,0],lat_rho[i,0]]))
        
    for i in np.arange(0,Ry_rho.shape[0],1):
        for j in np.arange(1,Ry_rho.shape[1],1):
            lat_rho[i,j] = lat + (Ry_rho[i,j] / latdist);
            lon_rho[i,j] = lon_rho[i,j-1] + ((Rx_rho[i,j]-Rx_rho[i,j-1]) / gsw.distance([lon_rho[i,j-1],lon_rho[i,j-1]+1],[lat_rho[i,j],lat_rho[i,j]]))
    del Rx_rho, Ry_rho
    
    
    
    ################## ----- U GRID ----- ##################
    #Create non-georeferenced grid in meters (origin = 0,0)
    x_u = (x_rho[:,0:L] + x_rho[:,1:Lp])/2
    y_u = (y_rho[:,0:L] + y_rho[:,1:Lp])/2
    
    #Rotate grid ...See: http://en.wikipedia.org/wiki/Rotation_(mathematics)
    Rx_u = x_u * np.cos(rotangle) - y_u * np.sin(rotangle)
    Ry_u = x_u * np.sin(rotangle) + y_u * np.cos(rotangle)
    
    #Estimate Latitude and Longitude of each Grid point
    lat_u = np.zeros(np.shape(Ry_u)) #initialize
    lon_u = np.zeros(np.shape(Ry_u)) #initialize
    lat_u[0,0] = lat + (Ry_u[0,0] / latdist)
    lon_u[0,0] = lon + (Rx_u[0,0] / gsw.distance([lon,lon+1],[lat_u[0,0],lat_u[0,0]]))
    for i in np.arange(1,Ry_u.shape[0],1):
        lat_u[i,0] = lat + (Ry_u[i,0] / latdist)
        lon_u[i,0] = lon_u[i-1,0] + ((Rx_u[i,0]-Rx_u[i-1,0]) / gsw.distance([lon_u[i-1,0],lon_u[i-1,0]+1],[lat_u[i,1],lat_u[i,0]]))
        
    for i in np.arange(0,Ry_u.shape[0],1):
        for j in np.arange(1,Ry_u.shape[1],1):
            lat_u[i,j] = lat + (Ry_u[i,j] / latdist)
            lon_u[i,j] = lon_u[i,j-1] + ((Rx_u[i,j]-Rx_u[i,j-1]) / gsw.distance([lon_u[i,j-1],lon_u[i,j-1]+1],[lat_u[i,j],lat_u[i,j]]));   
    del x_u, y_u, Rx_u, Ry_u
    
    ################## ----- V GRID ----- ##################
    #Create non-georeferenced grid in meters (origin = 0,0)
    x_v   = (x_rho[0:M,:]   + x_rho[1:Mp,:])/2
    y_v   = (y_rho[0:M,:]   + y_rho[1:Mp,:])/2

    #Rotate grid ...See: http://en.wikipedia.org/wiki/Rotation_(mathematics)
    Rx_v = x_v * np.cos(rotangle) - y_v * np.sin(rotangle);
    Ry_v = x_v * np.sin(rotangle) + y_v * np.cos(rotangle);

    #Estimate Latitude and Longitude of each Grid point
    lat_v = np.zeros(np.shape(Ry_v)) #initialize
    lon_v = np.zeros(np.shape(Ry_v)) #initialize
    lat_v[0,0] = lat + (Ry_v[0,0] / latdist)
    lon_v[0,0] = lon + (Rx_v[0,0] / gsw.distance([lon,lon+1],[lat_v[0,0],lat_v[0,0]]))
    for i in np.arange(1,Ry_v.shape[0],1):
        lat_v[i,0] = lat + (Ry_v[i,0] / latdist)
        lon_v[i,0] = lon_v[i-1,0] + ((Rx_v[i,0]-Rx_v[i-1,0]) / gsw.distance([lon_v[i-1,0],lon_v[i-1,0]+1],[lat_v[i,0],lat_v[i,0]]))

    for i in np.arange(0,Ry_v.shape[0]):
        for j in np.arange(1,Ry_v.shape[1]):
            lat_v[i,j] = lat + (Ry_v[i,j] / latdist);
            lon_v[i,j] = lon_v[i,j-1] + ((Rx_v[i,j]-Rx_v[i,j-1]) / gsw.distance([lon_v[i,j-1],lon_v[i,j-1]+1],[lat_v[i,j],lat_v[i,j]]))
    del x_v, y_v, Rx_v, Ry_v
    
    
    
    ################## ----- PSI GRID ----- ##################
    #Create non-georeferenced grid in meters (origin = 0,0)
    x_psi = (x_rho[0:M,0:L] + x_rho[1:Mp,1:Lp])/2
    y_psi = (y_rho[0:M,0:L] + y_rho[1:Mp,1:Lp])/2

    #Rotate grid ...See: http://en.wikipedia.org/wiki/Rotation_(mathematics)
    Rx_psi = x_psi * np.cos(rotangle) - y_psi * np.sin(rotangle)
    Ry_psi = x_psi * np.sin(rotangle) + y_psi * np.cos(rotangle)

    #Estimate Latitude and Longitude of each Grid point
    lat_psi = np.zeros(np.shape(Ry_psi)) #initialize
    lon_psi = np.zeros(np.shape(Ry_psi)) #initialize
    lat_psi[0,0] = lat + (Ry_psi[0,0] / latdist)
    lon_psi[0,0] = lon + (Rx_psi[0,0] / gsw.distance([lon,lon+1],[lat_psi[0,0],lat_psi[0,0]]))
    for i in np.arange(1,Ry_psi.shape[0]):
        lat_psi[i,0] = lat + (Ry_psi[i,0] / latdist)
        lon_psi[i,0] = lon_psi[i-1,0] + ((Rx_psi[i,0]-Rx_psi[i-1,0]) / gsw.distance([lon_psi[i-1,0],lon_psi[i-1,0]+1],[lat_psi[i,0],lat_psi[i,0]]))

    for i in np.arange(0,Ry_psi.shape[0]):
        for j in np.arange(1,Ry_psi.shape[1]):
            lat_psi[i,j] = lat + (Ry_psi[i,j] / latdist)
            lon_psi[i,j] = lon_psi[i,j-1] + ((Rx_psi[i,j]-Rx_psi[i,j-1]) / gsw.distance([lon_psi[i,j-1],lon_psi[i,j-1]+1],[lat_psi[i,j],lat_psi[i,j]]))
    del x_rho, y_rho, x_psi, y_psi, Rx_psi, Ry_psi, i, j

    
    ################## ----- spacing and parameters ----- ##################
    el      = lat_u[-1,0] - lat_u[0,0]
    xl      = lon_v[1,-1] - lon_v[1,1]

    dx = np.zeros([Mp,Lp])
    dy = np.zeros([Mp,Lp])

    for i in np.arange(0,lon_u.shape[0],1):
        dx[i,1:L-1]= gsw.distance(lon_u[i,1:L],lat_u[i,1:L]) #sperical distance calculation

    dx[:,0]  = dx[:,1]  
    dx[:,-2] = dx[:,-3]  
    dx[:,-1] = dx[:,-2] 

    for j in np.arange(0,lon_v.shape[1],1):
        dy[1:M-1,j]= gsw.distance(lon_v[1:M,j],lat_v[1:M,j]) #sperical distance calculation

    dy[0,:]  = dy[1,:]
    dy[-2,:] = dy[-3,:]
    dy[-1,:] = dy[-2,:]
    pm  = 1/dx
    pn  = 1/dy

    dndx = np.zeros(pm.shape)
    dmde = dndx
    dndx[1:M,1:L] = 0.5*(1/pn[1:M,2:Lp] - 1/pn[1:M,0:Lm])
    dmde[1:M,1:L] = 0.5*(1/pm[2:Mp,1:L] - 1/pm[0:Mm,1:L])

    angle = np.ones(pm.shape)*rotangle
    
    # Coriolis 
    # f = 2 * 7.29E-5 * np.sin(lat_rho * (np.pi/180)) #Estimation of Coriolis over the grid domain. OMEGA=7.29E-5
    #More info: http://en.wikipedia.org/wiki/Coriolis_effect#Formula
    #set f constant
    f = 1e-4*np.ones(lat_rho.shape)
    
    
###############################################################################
############################## - Create .nc file - ############################
###############################################################################
    
if save_grid == 1: 
    path_string         = directory+'mach4/grid_file.nc' 
    ncid                = nc.Dataset(path_string,mode='w',format='NETCDF4')
    ncid.description    = Descrip_grd
    ncid.author         = Author
        
    ###############################################################################
    ############################ - Create dimensions - ############################
    ###############################################################################
    
    #dimesions for variables, see notebook.
    dim_x2              = x.shape[0]
    dim_x1              = dim_x2 - 1 
    dim_y2              = y.shape[0]
    dim_y1              = dim_y2 - 1
    dim1                = 1
    
    # Meridional: eta
    eta_rho             = ncid.createDimension('eta_rho',dim_y2)
    eta_u               = ncid.createDimension('eta_u',dim_y2)
    eta_v               = ncid.createDimension('eta_v',dim_y1)
    eta_psi             = ncid.createDimension('eta_psi',dim_y1)
    
    # Zonal: xi
    xi_rho              = ncid.createDimension('xi_rho',dim_x2)
    xi_u                = ncid.createDimension('xi_u',dim_x1)
    xi_v                = ncid.createDimension('xi_v',dim_x2)
    xi_psi              = ncid.createDimension('xi_psi',dim_x1)
    
    one                 = ncid.createDimension('one',dim1)
    
    ###############################################################################
    ############################ - Define variables - #############################
    ###############################################################################
    
    #preliminary: factors, domain length & bathymetry
    dmde_var                      = ncid.createVariable('dmde','f8',('eta_rho','xi_rho'))
    dmde_var.long_name              = 'eta derivative of inverse metric factor pm'
    dmde_var.units                  = 'meter'
        
    dndx_var                      = ncid.createVariable('dndx','f8',('eta_rho','xi_rho'))
    dndx_var.long_name              = 'xi derivative of inverse metric factor pn'
    dndx_var.units                  = 'meter'
                
    el_var                        = ncid.createVariable('el','f8',('one'))
    el_var.long_name                = 'domain length in the ETA-direction'
    el_var.units                    = 'degrees'
                
    f_var                         = ncid.createVariable('f','f8',('eta_rho','xi_rho'))
    f_var.long_name                 = 'Coriolis parameter at RHO-points'
    f_var.units                     = 'second-1'
            
    h_var                         = ncid.createVariable('h','f8',('eta_rho','xi_rho'))
    h_var.long_name                 = 'bathymetry at RHO-points'
    h_var.units                     = 'meter'
        
    zice_var                      = ncid.createVariable('zice','f8',('eta_rho','xi_rho'))
    zice_var.long_name              = 'bathymetry at RHO-points'
    zice_var.units                  = 'meter'
        
    
    #lats
    lat_rho_var                   = ncid.createVariable('lat_rho','f8',('eta_rho','xi_rho'))
    lat_rho_var.long_name           = 'latitude of RHO-points'
    lat_rho_var.units               = 'degree_north'
        
    lat_psi_var                   = ncid.createVariable('lat_psi','f8',('eta_psi','xi_psi'))
    lat_psi_var.long_name           = 'latitude of PSI-points'
    lat_psi_var.units               = 'degree_north'
        
    lat_u_var                     = ncid.createVariable('lat_u','f8',('eta_u','xi_u'))
    lat_u_var.long_name             = 'latitude of U-points'
    lat_u_var.units                 = 'degree_north'
        
    lat_v_var                     = ncid.createVariable('lat_v','f8',('eta_v','xi_v'))
    lat_v_var.long_name             = 'latitude of V-points'
    lat_v_var.units                 = 'degree_north'
      
    
    #lons
    lon_rho_var                   = ncid.createVariable('lon_rho','f8',('eta_rho','xi_rho'))
    lon_rho_var.long_name           = 'longitude of RHO-points'
    lon_rho_var.units               = 'degree_east'
        
    lon_psi_var                   = ncid.createVariable('lon_psi','f8',('eta_psi','xi_psi'))
    lon_psi_var.long_name           = 'longitude of PSI-points'
    lon_psi_var.units               = 'degree_east'
        
    lon_u_var                     = ncid.createVariable('lon_u','f8',('eta_u','xi_u'))
    lon_u_var.long_name             = 'longitude of U-points'
    lon_u_var.units                 = 'degree_east'
        
    lon_v_var                     = ncid.createVariable('lon_v','f8',('eta_v','xi_v'))
    lon_v_var.long_name             = 'longitude of V-points'
    lon_v_var.units                 = 'degree_east'
        
    
    #masks
    mask_rho_var                  = ncid.createVariable('mask_rho','f8',('eta_rho','xi_rho'))
    mask_rho_var.long_name          = 'mask on RHO-points'
    mask_rho_var.units              = '(option_0:land,option_1:water)'
        
    mask_psi_var                  = ncid.createVariable('mask_psi','f8',('eta_psi','xi_psi'))
    mask_psi_var.long_name          = 'mask on PSI-points'
    mask_psi_var.units              = '(option_0:land,option_1:water)'
            
    mask_u_var                    = ncid.createVariable('mask_u','f8',('eta_u','xi_u'))
    mask_u_var.long_name            = 'mask on U-points'
    mask_u_var.units                = '(option_0:land,option_1:water)'
        
    mask_v_var                    = ncid.createVariable('mask_v','f8',('eta_v','xi_v'))
    mask_v_var.long_name            = 'mask on V-points'
    mask_v_var.units                = '(option_0:land,option_1:water)'
        
    
    #div
    pm_var                        = ncid.createVariable('pm','f8',('eta_rho','xi_rho'))
    pm_var.long_name                = 'curvilinear coordinate metric in XI'
    pm_var.units                    = 'meter-1'
        
    pn_var                        = ncid.createVariable('pn','f8',('eta_rho','xi_rho'))
    pn_var.long_name                = 'curvilinear coordinate metric in ETA'
    pn_var.units                    = 'meter-1'
        
    spherical_var                 = ncid.createVariable('spherical','f8',('one'))
    spherical_var.long_name         ='Grid type logical switch'
    spherical_var.option            = 'option_t: spherical'
        
    xl_var                        = ncid.createVariable('xl','f8',('one'))
    xl_var.long_name                = 'domain length in the XI-direction'
    xl_var.units                    = 'degrees'
        
    X_var                         = ncid.createVariable('X','f8',('one'))
    X_var.long_name                 = 'width of domain (degrees)'
        
    Y_var                         = ncid.createVariable('Y','f8',('one'))
    Y_var.long_name                 = 'width of domain (degrees)'
        
    dx_var                        = ncid.createVariable('dx','f8',('one'))
    dx_var.long_name                = 'resolution in x (degrees)'
            
    d_var                         = ncid.createVariable('d','f8',('one'))
    d_var.long_name                 = 'resolution in y (degrees)'
        
    angle_var                     = ncid.createVariable('angle','f8',('eta_rho','xi_rho'))
    angle_var.long_name             = 'Grid rotation angle'
    angle_var.units                 = 'rad'
        
    visc_factor_var               = ncid.createVariable('visc_factor','f8',('eta_rho','xi_rho'))
    visc_factor_var.long_name       ='horizontal viscosity sponge factor'
    visc_factor_var.coordinates     = 'lon_rho lat_rho'
            
    diff_factor_var               = ncid.createVariable('diff_factor','f8',('eta_rho','xi_rho'))
    diff_factor_var.long_name       = 'horizontal diffusivity sponge factor'
    diff_factor_var.coordinates     = 'lon_rho lat_rho'
    
    
    ###############################################################################
    ############################ - Assign variables - #############################
    ###############################################################################
    
    zeroes              = np.zeros([dim_y2,dim_x2])
    

    dmde_var[:]         = dmde
    dndx_var[:]         = dndx
    el_var[:]           = 4.5157011 #(value from Kjersti)
    f_var[:]            = f
    h_var[:]            = h              
    zice_var[:]         = np.zeros([dim_y2,dim_x2])             
    
    #lats
    lat_rho_var[:]      = lat_rho
    lat_psi_var[:]      = lat_psi
    lat_u_var[:]        = lat_u
    lat_v_var[:]        = lat_v
    
    #lons
    lon_rho_var[:]      = lon_rho
    lon_psi_var[:]      = lon_psi
    lon_u_var[:]        = lon_u
    lon_v_var[:]        = lon_v
    
    #masks
    mask_rho_var[:]     = mask_rho
    mask_psi_var[:]     = mask_psi      
    mask_u_var[:]       = mask_u          
    mask_v_var[:]       = mask_v          
    
    #div
    pm_var[:]           = pm
    pn_var[:]           = pn
    xl_var[:]           = xl
    X_var[:]            = X
    Y_var[:]            = Y
    angle_var[:]        = angle
    
    visc_factor_var[:]  = empty_fin.T
    diff_factor_var[:]  = empty_fin.T
    spherical_var[:]    = 5 
    dx_var[:]           = 437.86664147
    d_var[:]            = 0.01349892
    
    ncid.close()
