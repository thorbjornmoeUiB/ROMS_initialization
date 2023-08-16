#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  8 12:45:19 2023

@author: thorbjornostenbymoe
"""

"""
This script includes a short template for plotting output from ROMS. 
The script goes through:
    - Making and plotting the depth array (in rho coordinates)
    - Plotting a transect of temperature and density from the middle of the domain
    - Plotting horizontal depth-integrated circulation with arrows and sea surface height/speed
        - Having to interpolate u and v coordinates to rho
    - Calculate and plot volume transport
    - Making GIFS
"""

###############################################################################
################################ - PACKAGES - #################################
###############################################################################
from scipy.interpolate import interp1d
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import cmocean
import numpy.matlib as repmaty
import gsw.conversions as gswc
import gsw.density as gswd
from matplotlib.colors import TwoSlopeNorm
import imageio
cT=cmocean.cm.thermal


###############################################################################
################################# - DATASET - #################################
###############################################################################
path = '/Users/thorbjornostenbymoe/Desktop/Initialization_of_model/MODEL_SETUP/NEW_OUTPUT_HHH/roms_avg_0004.nc'
ds = nc.Dataset(path)

""" 
The dataset includes all variables purposfully included in the ocean.in/roms.in
script used for running the model. The variables and shapes can be seen by 
printing ds
"""
#print(ds)
"""
In this script only ubar, vbar, temp and salt is used, as these include the most 
central aspects of plotting ROMS output. Additionally, Z and Cs_r are used 
for plotting
"""

###############################################################################
######################### - Preliminary def and var - #########################
###############################################################################

# Resolution - needed to plot in km instead of # grid boxes. 
# Change to resolution chosen in the grid file
res = 1.5

# needed variables
# Variables are given in shape T = T(time,depth,y/meridional,x/zonal)
T = ds['temp'][:,:,:,:] # temperature
S = ds['salt'][:,:,:,:] # salt
U = ds['ubar'][:,:,:]   # depth integrated zonal velocity
V = ds['vbar'][:,:,:]   # depth integrated meridional velocity
h = ds['h'][:,:]        # depth 
Cs_r = ds['Cs_r'][:]    # stretching vector in rho coordinate
Cs_w = ds['Cs_w'][:]    # stretching vector in w coordinate

# rho dimensions (ROMS uses a staggered grid, read more at: https://www.myroms.org/wiki/Numerical_Solution_Technique#Horizontal_Discretization)
t_dim = np.shape(T)[0] # time
z_rho = np.shape(T)[1] # depth
y_rho = np.shape(T)[2] # y
x_rho = np.shape(T)[3] # x

# u dimensions (t and z are the same for u and v (not same for w!))
y_U = np.shape(U)[1]
x_U = np.shape(U)[2]

# v dimensions
y_V = np.shape(V)[1]
x_V = np.shape(V)[2]

#%%
###############################################################################
############################ - Make depth array - #############################
###############################################################################
""" 
if domain is constant in one dimension, this dimension can be neglected. 
The depth can be arbitrarily chosen to be negative or positive (*-1)
"""
H = np.zeros([z_rho,y_rho,x_rho]) # depth array

for y in range(y_rho): # loop over y
    for x in range(x_rho): # loop over x
        H[:,y,x] = h[y,x] * Cs_r #multiple depth at certain (x,y) point by stretching vector

#%%
###############################################################################
############################ - Plot depth array - #############################
###############################################################################

""" Plotting the depth array """
Xh,Yh = np.meshgrid(np.linspace(0,x_rho*res,x_rho),np.linspace(0,y_rho*res,y_rho)) #making custom dimensions for plotting

#make fig
fig  = plt.figure(figsize=[(x_rho*res)/100,(y_rho*res)/100])
grid = plt.GridSpec(1,1)
ax   = fig.add_subplot(grid[0,0])


z = ax.contourf(Xh,Yh,H[0,:,:],cmap='Blues_r') #plot depth
cbar = plt.colorbar(z,ax=ax) #cbar
cbar.set_label('Depth [m]')

ax.set_ylabel('Distance "North" from \nsouthwestern corner in km')
ax.set_xlabel('Distance "East" from southwestern \ncorner in km')

#%%
###############################################################################
############################## - Plot transect - ##############################
###############################################################################
"""
Here we're going to make a Transect of temperature and density from an arbitrary 
distance west in the domain (here x = middle)
"""
# position and timestep for transect
tra_x = int(x_rho/2)
timestep = 0

# transect dimensions
Xt = repmaty.repmat(np.linspace(0,y_rho*res,y_rho),z_rho,1)
Yt = H[:,:,tra_x]

temp_tra = T[timestep,:,:,tra_x]
salt_tra = S[timestep,:,:,tra_x]
# levels for plot, better to customize or make better automatic, but this works 
# if you don't know the range
levels_T = np.linspace(np.round(np.min(temp_tra),3),np.round(np.max(temp_tra),3),10)



# fig
fig  = plt.figure(figsize=[7,4])
grid = plt.GridSpec(1,1)
ax   = fig.add_subplot(grid[0,0])

# add gray fill and gray line under bottom
ax.set_facecolor('lightgray')
ax.plot(Xt[0,:],Yt[0,:],linewidth=2,color='dimgray')

# plot T
z = ax.contourf(Xt,Yt,temp_tra,levels=levels_T,cmap=cT)
cbar = plt.colorbar(z,ax=ax) #cbar
cbar.set_label(r'$\theta$'+' [C'+r'$^{\circ}$'+']')

# calc density
press = gswc.p_from_z(Yt,68)
SA = gswc.SA_from_SP(salt_tra, press,-10,68)
r = gswd.sigma0(SA,temp_tra)

# plotting density as contour lines
cont = ax.contour(Xt,Yt,r,colors='k',linewidths=.75)
plt.clabel(cont,inline=1, fontsize=8,fmt=' {:.2f} '.format)
#cont.collections[2].set_linewidth(1) #for if you want to highlight a specific contour level

# misc
ax.set_ylim([-1200,0])
ax.set_ylabel('Depth [m]')
ax.set_xlabel('Distance North [km]')
#%%
###############################################################################
############### - Plot Horizontal depth-integrated circulation - ##############
###############################################################################

# Optional function for removing small data (if not used, these are plotted at small dots)
def exclude_small(U,V,lim):     # input is u and v velocities and the limit for plotting
    s = np.ones(np.shape(U))    
    masku = np.ones(np.shape(U))
    UV = np.sqrt(U**2 + V**2)   # Absolute speed used as threshold
    masku[UV<lim] = np.nan      # remove values under threshold
    U = U*masku
    V = V*masku
    s = s*masku
    return(U,V,s)


T_dim_hor   = T[timestep,0,:,:]
V_hor       = V[timestep,:,:]
U_hor       = U[timestep,:,:]

"""
Interpolating both U and V to rho coordinates
"""
old = np.linspace(0,1,np.shape(U_hor)[1]) # U coordinate we want to change, in this case dimension 1
new = np.linspace(0,1,np.shape(T_dim_hor)[1]) # rho coordinate we want to change it to
u = interp1d(old, U_hor, axis=1)(new) #int old to new coord

X = np.linspace(0,1,np.shape(V_hor)[0]) # V coordinate we want to change, in this case dimension 0
Y = np.linspace(0,1,np.shape(T_dim_hor)[0]) # rho coordinate we want to change it to
v = interp1d(X, V_hor, axis=0)(Y) #int old to new coord

# If we plot all datapoints, the whole plot is going to be black, this therefore
# includes only the nth data-point. "You" can experiment with different sizes...
freq = 10 # distance between each plotted arrow
offs = 6  # offset of merridional variable due to mask in southern end
Xh_quiv,Yh_quiv,u_quiv,v_quiv = Xh[offs::freq,::freq],Yh[offs::freq,::freq],u[offs::freq,::freq],v[offs::freq,::freq]

# using the function in the top of the section to remove the smaller values
u_quiv1,v_quiv1,s = exclude_small(u_quiv,v_quiv,0.01) # remove small values


# Fig
fig  = plt.figure(figsize=[(x_rho*res)/100,(y_rho*res)/100])
grid = plt.GridSpec(1,1)
ax   = fig.add_subplot(grid[0,0])

# add colormap underr arrows, could chose sea surface elevation ('zeta'), zonal/meridional 
# velocity, etc. Here the absolute speed is plotted
c = ax.contourf(Xh[offs:,:],Yh[offs:,:],np.sqrt(u[offs:,:]**2 + v[offs:,:]**2),cmap='Reds')
cb = plt.colorbar(c,ax=ax)
cb.set_label('Abs. Speed [m s$^{-1}$]')

# plot arrows, many option could be chosen in terms of arrow size and shape etc.
quiv = ax.quiver(Xh_quiv,Yh_quiv,u_quiv1,v_quiv1)
plt.quiverkey(quiv, 0.9, 0.9, 0.1, label='10 cm s$^{-1}$',) # add scale for arrows

# misc
ax.set_ylim([0,y_rho*res])
ax.set_xlim([0,x_rho*res])
ax.set_xlabel('Distance east [km]')
ax.set_ylabel('Distance north [km]')


#%%
###############################################################################
############################# - Volume Transport - ############################
###############################################################################

# new depth array for w levels (staggered)
H1 = np.zeros([z_rho+1,y_rho,x_rho]) # depth array

for y in range(y_rho): # loop over y
    for x in range(x_rho): # loop over x
        H1[:,y,x] = h[y,x] * Cs_w # multiple depth at certain (x,y) point by stretching vector


# calc diffs to find depth of each level and multiply by horizontal dim to get area
area = np.zeros([z_rho,y_rho])
for i in range(y_rho):
    area[:,i] = np.diff(H1[:,i,0])*res*1000 # times resolution and 1000m to convert to meter
# area is now the are of each grid box
#%%
# levels for area plot
lvlA = np.arange(0,125000,5000)

# plot to see that the are makes sense (in terms of the stretching vector)
# if the are is largest where the depth of each grid cell is largest, then the
# area's are correct
fig  = plt.figure(figsize=[7,4])
grid = plt.GridSpec(1,2,width_ratios=(1,0.2),wspace=0.05)
ax   = fig.add_subplot(grid[0,0])

# plot transect of areas
z = ax.contourf(Xt,Yt,area,levels=lvlA)

# misc
ax.set_ylabel('Depth [m]')
ax.set_ylim([-1200,0])
ax.set_xlabel('Distance north [km]')

# plot distance between levels
ax   = fig.add_subplot(grid[0,1])
ax.plot(np.diff(Cs_w)*-np.min(H1),np.linspace(np.min(H1),0,30),'.')

# colorbar for area plot
cb = plt.colorbar(z,ax=ax)
cb.set_label('Area of each grid cell [m^2]')

# misc
ax.tick_params(labelbottom=True, labeltop=False, labelleft=False, labelright=False,
                     bottom=True, top=False, left=False, right=False)
ax.set_xlabel('depth of each \ngrid cell [m]')
ax.set_ylim([-1200,0])

#%%
""" Plotting VT and speed to get a sense of size etc."""

# Calculate volume transport (in Sv)
Ut = ds['u'][timestep,:,:,tra_x] # get transect of u from middle :)
VT = Ut*area                     # calc volume transport area * speed [m^2 * m/s = m^3/s]
VTs = VT/1000000                 # VT in sverdrups


# normalize colorbar
norm = TwoSlopeNorm(vmin=np.min(VTs), vcenter=0, vmax=-np.min(VTs))

# fig
fig  = plt.figure(figsize=[7,6])
grid = plt.GridSpec(2,1,hspace=0.1)
ax   = fig.add_subplot(grid[0,0])

# plot VT
z = ax.contourf(Xt,Yt,VTs,cmap='RdBu_r',norm=norm)
cbar1 = plt.colorbar(z,ax=ax)
cbar1.set_label('Volume transport in each grid box [Sv]')

# misc
ax.set_ylabel('Depth [m]')
ax.set_ylim([-1200,0])
ax.tick_params(labelbottom=False, labeltop=False, labelleft=True, labelright=False,
                     bottom=False, top=False, left=True, right=False)



# norm and levels for zonal current speed
norm = TwoSlopeNorm(vmin=np.min(Ut*100), vcenter=0, vmax=-np.min(Ut*100))
levels = np.arange(-0.10,0.061,0.01)*100 # remove or customize
ax   = fig.add_subplot(grid[1,0])

#plot speed
z = ax.contourf(Xt,Yt,Ut*100,cmap='RdBu_r',norm=norm,levels=levels)
cbar2 = plt.colorbar(z,ax=ax)
cbar2.set_label('Along-slope velocity [cm/s]')

#misc
ax.set_ylabel('Depth [m]')
ax.set_ylim([-1200,0])
ax.set_xlabel('Distance north [km]')

#%%
""" 
choose speed threshold for volume transport within currents (could be u = 0 but 
this is often misrepresenting the VT within the current ....)
"""

thresh = -0.03 #m/s
mask = np.zeros(np.shape(Ut))
Vt_curr = np.where(Ut < thresh,VTs,mask)


# check that the sign is right (does it retain only the region with the deisred current?)
fig  = plt.figure(figsize=[7,4])
grid = plt.GridSpec(1,1)
ax   = fig.add_subplot(grid[0,0])
z = ax.contourf(Xt,Yt,Vt_curr,cmap='Reds_r')
plt.colorbar(z,ax=ax)
ax.set_ylabel('Depth [m]')
ax.set_ylim([-1200,0])
ax.set_xlabel('Distance north [km]')

# Finally we have the total volume transport within the current over the speed threshold
print('Volume transport of current (u < '+str(thresh)+' m/s) is: '+str(np.sum(Vt_curr))+' Sv')

#%%
###############################################################################
################################### - GIF - ###################################
###############################################################################

""" chose first what to make a gif of, then create the figures and save a 
series of figures in a logical way """

U_GIF = ds['ubar'][:,:,:]
V_GIF = ds['vbar'][:,:,:]
Z_GIF = ds['zeta'][:,:,:]

lvlZ = np.arange(-0.048,-0.012,0.002)


"""
Doing the same steps as for the horizontal fields but keeping the time dimension throughout
"""
old = np.linspace(0,1,np.shape(U_GIF)[2]) # U coordinate we want to change, in this case dimension 1
new = np.linspace(0,1,np.shape(Z_GIF)[2]) # rho coordinate we want to change it to
u = interp1d(old, U_GIF, axis=2)(new) #int old to new coord

X = np.linspace(0,1,np.shape(V_GIF)[1]) # V coordinate we want to change, in this case dimension 0
Y = np.linspace(0,1,np.shape(Z_GIF)[1]) # rho coordinate we want to change it to
v = interp1d(X, V_GIF, axis=1)(Y) #int old to new coord

freq = 10 # frequency of plotted arrows
offs = 6  # offset of merridional variable due to mask in southern end
Xh_quiv,Yh_quiv,u_quiv,v_quiv = Xh[offs::freq,::freq],Yh[offs::freq,::freq],u[:,offs::freq,::freq],v[:,offs::freq,::freq]

u_quiv1,v_quiv1,s = exclude_small(u_quiv,v_quiv,0.01)

#loop to make and save figs
for i in range(t_dim):
    fig  = plt.figure(figsize=[7,4])
    grid = plt.GridSpec(1,1)
    ax   = fig.add_subplot(grid[0,0])
    
    z = ax.contourf(Xh,Yh,Z_GIF[i,:,:],cmap='Blues_r',levels=lvlZ)
    cb = plt.colorbar(z,ax=ax)
    cb.set_label('Free surface [m]')
    
    quiv,ax.quiver(Xh_quiv,Yh_quiv,u_quiv1[i,:,:],v_quiv1[i,:,:])
    plt.quiverkey(quiv, 0.05, 0.05, 0.1, label='10 cm s$^{-1}$',) # add scale for arrows
    
    ax.set_ylabel('Distance north [km]')
    ax.set_xlabel('Distance east [km]')
    
    fig.savefig('/Users/thorbjornostenbymoe/Desktop/GIFS_test/pic'+str(i)+'.png',dpi=300, facecolor='w', edgecolor='w', orientation='landscape',format=None,transparent=False,bbox_inches='tight',pad_inches=0.25)

#%%
""" Then use this short code to make a GIF of the saved figures """


with imageio.get_writer('/Users/thorbjornostenbymoe/Desktop/GIFS_test/movie.gif', mode='I',duration=0.5) as writer:
    for filename in range(t_dim):
        image = imageio.imread('/Users/thorbjornostenbymoe/Desktop/GIFS_test/pic' + str(filename)+'.png')
        writer.append_data(image)

#%%


