#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 13 08:53:14 2023

@author: thorbjorn
"""
class GlobParam(object):
    import getpass
    import os
    import sys
    import netCDF4 as nc
    
    work_p      = '/cluster/work/users/huk006/metroms_ideal_run'
    proj_p      = '/cluster/projects/nn9608k/Iceland_Ideal'
    home_p      = '/cluster/home/huk006/metroms_ideal/apps'
    gen_dir     = '/PHD_ii_HH_V2'
    exp_dir     = '/TEST1_new_params/TEST_FRAMEWORK'
    readF       = '/readme_file.nc'
    ds          = nc.Dataset(proj_p+gen_dir+exp_dir+readF)

    xpoint      = str(int(ds['xpoint_'][0]))
    ypoint      = str(int(ds['ypoint_'][0]))
    
    tileX       = str(int(ds['tileX_'][0]))
    tileY       = str(int(ds['tileY_'][0]))
    
    Vstretch    = str(int(ds['Vstretch_'][0]))
    thetaS      = str(ds['thetaS_'][0])+'d0'
    thetaB      = str(ds['thetaB_'][0])+'d0'
    Tcline      = str(ds['Tcline_'][0])+'d0'
    
    grdFN       = '/grid_file.nc'
    iniFN       = '/init_file.nc'
    bryFN       = '/bry_file.nc'
    clmFN       = '/clm_file.nc'
    nudFN       = '/clm_nudge_file.nc'
    frcFN       = '/frc_file.nc'
    

    
    
    
    
