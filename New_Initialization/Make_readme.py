#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 21 15:41:55 2023

@author: thorbjorn
"""
import netCDF4 as nc
from GlobalParameters import * 

print('###############################################################################')
print('########################### - Making readme file - ############################')
print('###############################################################################')

path = GlobalParameters.exp_dir + 'readme.txt'

string1 = 'Experiment with parameters: \n'
string2 = 'Slope:      '+str(GlobalParameters.slope)+'\n'
string3 = 'Forcing:    '+str(GlobalParameters.exp)+'\n'
string4 = 'Resolution: '+str(GlobalParameters.Res)+' km \n'
string5 = 'Dimensions: '+str(int(GlobalParameters.X/1000))+'x'+str(int(GlobalParameters.Y/1000))+' km \n'
string6 = GlobalParameters.statement

string  = string1 + string2 + string3 + string4 + string5 + string6

with open(path, 'w') as f:
    f.writelines(string)
    f.close()
    
    


###############################################################################
############################# - Create .nc file - #############################
###############################################################################

ncid            = nc.Dataset(GlobalParameters.exp_dir + 'readme_file.nc',mode='w',format='NETCDF4')
one             = ncid.createDimension('one',    1) 


xpoint_      = ncid.createVariable('xpoint_','f8')
ypoint_      = ncid.createVariable('ypoint_','f8') 

tileX_       = ncid.createVariable('tileX_','f8')
tileY_       = ncid.createVariable('tileY_','f8') 

Vstretch_    = ncid.createVariable('Vstretch_','f8')
thetaS_      = ncid.createVariable('thetaS_','f8')
thetaB_      = ncid.createVariable('thetaB_','f8') 
Tcline_      = ncid.createVariable('Tcline_','f8')


if GlobalParameters.Res == 1.5:
    xpoint_[:] = 400
    ypoint_[:] = 266
else:
    xpoint_[:] = 600
    ypoint_[:] = 399

tileX_[:] = 4
tileY_[:] = 4

if GlobalParameters.Use_stretch == 1:
    Vstretch_[:] = 2
    thetaS_[:]   = float(GlobalParameters.theta_s)
    thetaB_[:]   = float(GlobalParameters.theta_b)
    Tcline_[:]   = 40.0
elif GlobalParameters.Use_stretch == 0:
    if GlobalParameters.New_S_exag == 0:
        Vstretch_[:] = 3
        thetaS_[:]   = 4.0
        thetaB_[:]   = 0.9
        Tcline_[:]   = 20.0
    elif GlobalParameters.New_S_exag == 1:
        Vstretch_[:] = 3
        thetaS_[:]   = 4.0
        thetaB_[:]   = 2.0
        Tcline_[:]   = 20.0
    if GlobalParameters.New_S_exag == 2:
        Vstretch_[:] = 3
        thetaS_[:]   = 4.0
        thetaB_[:]   = 1.2
        Tcline_[:]   = 20.0



ncid.close()
