#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 21 15:41:55 2023

@author: thorbjorn
"""
import netCDF4 as nc
from GlobParam import * 

print('###############################################################################')
print('########################### - Making readme file - ############################')
print('###############################################################################')

path = GlobParam.EXP_dir + GlobParam.SubEXP_dir + 'readme.txt'

stringp1  = 'You have now made an experiment or stumbled upon a readme file for the IcelandIdeal model in ROMS. '
stringp2  = 'This file belongs to a host of experiments designed to study the mechanisms behind '
stringp3  = 'the emergence of the NIJ. Good luck with experimenting! :)\n'
string1 = 'Experiment with parameters: \n'
string2 = 'Slope:      '+str(GlobParam.slope)[:-1]+'\n'
string3 = 'Forcing:    '+str(GlobParam.exp)+'\n'
string4 = 'Resolution: '+str(GlobParam.Res)+' km \n'
string5 = 'Dimensions: '+str(int(GlobParam.X/1000))+'x'+str(int(GlobParam.Y/1000))+' km \n'
string6 = GlobParam.statement

string  = stringp1+stringp2+stringp3+string1 + string2 + string3 + string4 + string5 + string6

with open(path, 'w') as f:
    f.writelines(string)
    f.close()
    
    


###############################################################################
############################# - Create .nc file - #############################
###############################################################################

ncid            = nc.Dataset(GlobParam.EXP_dir + GlobParam.SubEXP_dir + 'ExpParams_file.nc',mode='w',format='NETCDF4')
one             = ncid.createDimension('one',    1) 


xpoint_      = ncid.createVariable('xpoint_','f8')
ypoint_      = ncid.createVariable('ypoint_','f8') 

tileX_       = ncid.createVariable('tileX_','f8')
tileY_       = ncid.createVariable('tileY_','f8') 

Vstretch_    = ncid.createVariable('Vstretch_','f8')
thetaS_      = ncid.createVariable('thetaS_','f8')
thetaB_      = ncid.createVariable('thetaB_','f8') 
Tcline_      = ncid.createVariable('Tcline_','f8')


if GlobParam.Res == 1.5:
    xpoint_[:] = 400
    ypoint_[:] = 266
elif GlobParam.x2_domain == 1:
    xpoint_[:] = 1200
    ypoint_[:] = 399
else:
    xpoint_[:] = 600
    ypoint_[:] = 399


tileX_[:] = 4
tileY_[:] = 4



Vstretch_[:] = 3
thetaS_[:]   = 4.0
thetaB_[:]   = 1.2
Tcline_[:]   = 20.0



ncid.close()
