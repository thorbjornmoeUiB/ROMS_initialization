#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  6 12:27:53 2024

@author: thorbjorn
"""

###############################################################################
################################## - Packages - ###############################
###############################################################################

import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy.matlib as repmaty
from GlobParam import *                

print('###############################################################################')
print('###################### - Double_checking_input_fields - #######################')
print('###############################################################################')

###############################################################################
############################### - User options - ##############################
###############################################################################

#File directory/names
common_dir          = GlobParam.common_dir
directory           = GlobParam.EXP_dir

grd     = directory + GlobParam.grid_file
bry     = directory + GlobParam.bry_file
clm     = directory + GlobParam.clm_file
clm_nud = directory + GlobParam.clm_nudge_file
stre    = common_dir + GlobParam.Streching

grd_ds = nc.Dataset(grd)
bry_ds = nc.Dataset(bry)
clm_ds = nc.Dataset(clm)
nud_ds = nc.Dataset(clm_nud)
str_ds = nc.Dataset(stre)

pIFSJ               = GlobParam.pIFSJ
pNIIC_w             = GlobParam.pNIIC_w
pIFSJ=1
pNIIC_w=1
###############################################################################
########################### - Make bathymetry array - #########################
###############################################################################

h  = grd_ds['h'][:]
c  = str_ds['Cs_r'][:]
c2 = str_ds['Cs_w'][:]

#plotting
Y_west = h[np.newaxis,:,0]*c[:,np.newaxis]
Y_east = h[np.newaxis,:,-1]*c[:,np.newaxis]
X = repmaty.repmat(np.linspace(0,401,401),30,1)

# Volume transport
VT_west = np.diff(h[np.newaxis,:,0]*c2[:,np.newaxis],axis=0)*1000
VT_east = np.diff(h[np.newaxis,:,-1]*c2[:,np.newaxis],axis=0)*1000

###############################################################################
############################ - PLOTTING BRY FILES - ###########################
###############################################################################
if pIFSJ == 0 and pNIIC_w == 0:
    fig = plt.figure(figsize=[7,4])
    grid = plt.GridSpec(2,1,hspace=0.2)
    
elif pIFSJ == 1 and pNIIC_w == 0:
    fig = plt.figure(figsize=[14,4])
    grid = plt.GridSpec(2,2,hspace=0.2)
    
elif pIFSJ == 0 and pNIIC_w == 1:
    fig = plt.figure(figsize=[7,6])
    grid = plt.GridSpec(3,1,hspace=0.2)

elif pIFSJ == 1 and pNIIC_w == 1: 
    fig = plt.figure(figsize=[14,6])
    grid = plt.GridSpec(3,2,hspace=0.2)

plt.suptitle('Testing the Boundary file (bry_file.nc)')

ax = fig.add_subplot(grid[0, 0])
ax.set_title('West')
ax.plot(X[0,:],Y_west[0,:],color='dimgray')
ax.fill_between(X[0,:],Y_west[0,:],-1200,color='lightgray')
z = ax.contourf(X,Y_west,bry_ds['u_west'][0,:,:],cmap='RdBu_r',levels=np.arange(-0.15,0.16,0.01),extend='both')
cbar = plt.colorbar(z,ax=ax)
cbar.set_label('Along-slope \nvelocity [m s$^{-1}$]')
ax.set(ylim=[-1200,0],yticks=[-1200,-800,-400,0],xticklabels=[],ylabel='Depth [m]')

ax = fig.add_subplot(grid[1, 0])
ax.plot(X[0,:],Y_west[0,:],color='dimgray')
ax.fill_between(X[0,:],Y_west[0,:],-1200,color='lightgray')
z = ax.contourf(X,Y_west,bry_ds['u_west'][0,:,:]*VT_west,cmap='PuOr',levels=np.arange(-7000,7010,1000),extend='both')
cbar = plt.colorbar(z,ax=ax)
cbar.set_label('Volume transport [m$^{3}$ s$^{-1}$]')
ax.text(10,-900,'NIIC: '+str(np.round(np.nansum(np.where(bry_ds['u_west'][0,:,:]*VT_west>0,bry_ds['u_west'][0,:,:]*VT_west,0))/1000000,3))+' Sv')
ax.text(10,-1050,'NIJ: '+str(np.round(np.nansum(np.where(bry_ds['u_west'][0,:,:]*VT_west<0,bry_ds['u_west'][0,:,:]*VT_west,0))/1000000,3))+' Sv')

if pNIIC_w == 0:
    ax.set(ylim=[-1200,0],yticks=[-1200,-800,-400,0],ylabel='Depth [m]',xlabel='Distance North [km]')

elif pNIIC_w == 1:
    ax.set(ylim=[-1200,0],yticks=[-1200,-800,-400,0],xticklabels=[],ylabel='Depth [m]')
    
    ax = fig.add_subplot(grid[2, 0])
    ax.plot(X[0,:],Y_west[0,:],color='dimgray')
    ax.fill_between(X[0,:],Y_west[0,:],-1200,color='lightgray')
    z = ax.contourf(X,Y_west,bry_ds['temp_west'][0,:,:],cmap='cmo.thermal',levels=np.arange(-0.8,1,0.1),extend='max')
    cbar = plt.colorbar(z,ax=ax)
    cbar.set_label('Temperature [C]')
    ax.set(ylim=[-1200,0],yticks=[-1200,-800,-400,0],ylabel='Depth [m]',xlabel='Distance North [km]')
    
if pIFSJ == 1:
    ax = fig.add_subplot(grid[0, 1])
    ax.set_title('East')
    ax.plot(X[0,:],Y_east[0,:],color='dimgray')
    ax.fill_between(X[0,:],Y_east[0,:],-1200,color='lightgray')
    z = ax.contourf(X,Y_east,bry_ds['u_east'][0,:,:],cmap='RdBu_r',levels=np.arange(-0.15,0.16,0.01),extend='both')
    cbar = plt.colorbar(z,ax=ax)
    cbar.set_label('Along-slope \nvelocity [m s$^{-1}$]')
    ax.set(ylim=[-1200,0],yticks=[-1200,-800,-400,0],xticklabels=[])
    
    ax = fig.add_subplot(grid[1, 1])
    ax.plot(X[0,:],Y_east[0,:],color='dimgray')
    ax.fill_between(X[0,:],Y_east[0,:],-1200,color='lightgray')
    z = ax.contourf(X,Y_east,bry_ds['u_east'][0,:,:]*VT_east,cmap='PuOr',levels=np.arange(-7000,7010,1000),extend='both')
    cbar = plt.colorbar(z,ax=ax)
    cbar.set_label('Volume transport [m$^{3}$ s$^{-1}$]')
    ax.text(10,-900,'IFSJ: '+str(np.round(np.nansum(np.where(bry_ds['u_east'][0,:,:]*VT_east>0,bry_ds['u_east'][0,:,:]*VT_east,0))/1000000,3))+' Sv')
    ax.text(10,-1050,'NIJ: '+str(np.round(np.nansum(np.where(bry_ds['u_east'][0,:,:]*VT_east<0,bry_ds['u_east'][0,:,:]*VT_east,0))/1000000,3))+' Sv')
    ax.set(ylim=[-1200,0],yticks=[-1200,-800,-400,0])
    
    
#%%
if pIFSJ == 0:
    fig = plt.figure(figsize=[7,3])
    grid = plt.GridSpec(1,1,hspace=0.2)

elif pIFSJ == 1: 
    fig = plt.figure(figsize=[7,6])
    grid = plt.GridSpec(2,1,hspace=0.2)
    
ax = fig.add_subplot(grid[0, 0])
ax.set_title('Ubar (if clm, nud, and bry file line up then we are good)')
ax.plot(clm_ds['ubar'][0,:,0],'k',marker='o',label='clm file')
ax.plot(bry_ds['ubar_west'][0],'Teal',label='bry file',linewidth=4)
ax.plot(nud_ds['M2_NudgeCoef'][:,0]/2-0.1,'red',label='clim nudge$^{[1]}$')

ax.text(190,-0.07,'$^{[1]}$clim nudge should always be positive')

ax.set(xlim=[0,401],ylabel='depth-averaged velocity [m s$^{-1}$] \n (and "normalized" nudging)')
ax.legend()    
if pIFSJ == 0:
    ax.set_xlabel('Distance north [km]')

if pIFSJ == 1:
    ax.set_xticks([])
    ax = fig.add_subplot(grid[1, 0])
    ax.set_title('Ubar east')
    ax.plot(clm_ds['ubar'][0,:,-1],'k',marker='o',label='clm file')
    ax.plot(bry_ds['ubar_east'][0],'Teal',label='bry file',linewidth=4)
    ax.plot(nud_ds['M2_NudgeCoef'][:,-1]/2-0.1,'red',label='clim nudge')

    ax.set(xlim=[0,401],ylabel='depth-averaged velocity [m s$^{-1}$]',xlabel='Distance north [km]')
    ax.legend()  
