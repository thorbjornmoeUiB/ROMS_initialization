#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 14 15:12:32 2023

@author: thorbjorn
"""

class GlobalParameters(object):
    import getpass
    import os
    import sys
    
    ###############################################################################
    ############################# - Directory paths - #############################
    ###############################################################################
    
    old_dir         = '/Users/thorbjorn/Library/CloudStorage/OneDrive-UniversityofBergen/Documents/Master_mod/MODEL_SETUP/'
    common_dir      = '/Users/thorbjorn/Library/CloudStorage/OneDrive-UniversityofBergen/Documents/PhD/Iceland_Sea_Mod/'
    exp_dir         = common_dir + 'New_test_of_params/INPUT/Test_pNIIC&NIJ_warm/'
    
    ###############################################################################
    ########################## - Name of all used files - #########################
    ###############################################################################
    
    grid_file       = 'grid_file.nc'
    stretching      = 'stet_file.nc'
    init_file       = 'init_file.nc'
    frc_file        = 'frc_file.nc'
    bry_file        = 'bry_file.nc'
    clm_file        = 'clm_file.nc'
    clm_nudge_file  = 'clm_nudge_file.nc'
    stretching_alt  = 'Testing_params_OUTPUT/Stretch_testing/Vs3_Vt1_s4_b09_hc20.nc'
    Mean_profiles   = 'New_depth_profiles_for_model.nc'
    Raw_profiles    = 'InputProfileCentralIcelandSea.mat'
    
    ###############################################################################
    ############################# - File descriptors - ############################
    ###############################################################################
    
    General_desc    = 'Iceland Ideal - ROMS idealized model run - '
    spec_desc_grd   = 'Grid file'
    spec_desc_ini   = 'Initial file'
    spec_desc_frc   = 'Forcing file'
    spec_desc_bry   = 'Boundary file'
    spec_desc_clm   = 'Climatology file'
    spec_desc_nud   = 'Climatology nudging file'
    Author          = 'Thorbjoern Oestenby Moe'
    
    ###############################################################################
    ############################# - Logical switches - ############################
    ###############################################################################
    
    Use_stretch     = 1         # Alternatives:           1 (=yes)      0 (=no) --> (use stretching from ROMS output)
    Gen_plot        = 0         # plotting relatively unimportant figures from all scripts
    Imp_plot        = 1         # plotting the important figures from all scripts
    Create_files    = 1         # Create all the files
    create_grid     = 1         # only for grid file (making the grid takes a lot of time) Alternatives: 1 (=yes)      0 (=no)
    
    ###############################################################################
    ############################# - Major parameters - ############################
    ###############################################################################
    
    z               = 30        # Depth layers                  - Always 30
    Res             = 1.0       # Resolution                    - Alternatives = 1             1.5             (3)
    X               = 600000    # Zonal dimension               - Alternatives = 650000        600000
    Y               = 400000    # Meridional dimension          - Always 400000
    steepnes_factor = 7         # Steepness of slope            - Alternatives = 7 (steep)     15 (moderate)   (25 (mild slope))
    theta_s         = 7         # Layer distribution at surface - Always 7
    theta_b         = 14        # Layer distribution at bottom  - Alternatives = 2             14
    
    ###############################################################################
    ########################## - Experiment composition - #########################
    ###############################################################################
    
    # if all of these are = 0 then its a ref!
    pNIJ            = 1                                       # EXP 1 == prescribed NIJ outflow in west
    pNIIC           = 0                                       # EXP 2 == prescribed NIIC inflow in west
    pNIIC_hydro     = 1                                       # EXP 2.33 == Same pNIIC but with warm water flowing in
    pNIIC2          = 0                                       # EXP 2.67 == altered prescribed NIIC inflow in west
    pNIJ_md         = 0                                       # EXP 3 == mid-depth confined prescribed NIJ outflow in west
    pIFSJ           = 0                                       # EXP 4 == prescribed IFSJ outflow in east
    onslope         = 0                                       # Sets degree to which pIFSJ is prescribed on the boundary, very onslope,1=onslope, 0=offslope
    NIIC_EB         = 0                                       # EXP 5 == prescribed NIIC outflow in east 
    
    if steepnes_factor == 7:
        slope = 'Steep_'
    else:
        slope =  'Moderate_'
    
    if pNIJ == 0 and pNIIC == 0 and pNIIC2 == 0 and pNIJ_md == 0 and pIFSJ == 0 and NIIC_EB == 0:
        exp = 'REF'
    elif pNIJ == 1 and pNIIC == 0:
        exp = 'pNIJ'
    elif pNIJ == 1 and pNIIC == 1:
        exp = 'pNIIC&NIJ_warm'
    elif pNIJ == 1 and pNIIC_hydro == 1:
        exp = 'pNIIC&NIJ'
    elif pNIJ_md == 1 and pNIIC == 0:
        exp = 'pNIJ_md'
    elif pNIJ_md == 1 and pNIIC == 1:
        exp = 'pNIIC&NIJ_md'
    
    if Use_stretch == 1:
        statement = 'With OLD STRETCHING (Vstretchin = 2, theta_s = '+str(theta_s)+', and theta_b = '+str(theta_b)
    else:
        statement = 'With NEW STRETCHING (Vstretchin = 3, theta_s = 4, and theta_b = 0.9)'
    
    print('\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
    print('Making run for: \n\nResolution  = '+str(Res)+' km \n\nDomain size = '+str(int(X/1000))+'x'+str(int(Y/1000))+' km \n\nExperiment  = '+slope+exp+'\n\n'+statement+'\n\nIn the directory: \n'+exp_dir+'<--')
    print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
    
    
    
    
    
    
    
    
    
    
