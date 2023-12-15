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
    exp_dir         = common_dir + 'New_test_of_params/INPUT/REF_MOD_OON3/'
    
    ###############################################################################
    ########################## - Name of all used files - #########################
    ###############################################################################
    
    grid_file       = 'grid_file.nc'
    stretching      = 'Stetching.nc'
    init_file       = 'init_file.nc'
    frc_file        = 'frc_file.nc'
    bry_file        = 'bry_file.nc'
    clm_file        = 'clm_file.nc'
    clm_nudge_file  = 'clm_nudge_file.nc'
    stretching_alt1 = 'Testing_params_OUTPUT/Stretch_testing/Vs3_Vt1_s4_b09_hc20.nc'
    stretching_alt2 = 'Testing_params_OUTPUT/Stretch_testing/Vs3_Vt1_s4_s2_hc20.nc'
    stretching_alt3 = 'Testing_params_OUTPUT/Stretch_testing/Vs3_Vt1_s4_b12_hc20.nc'
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
    
    Use_stretch     = 0         # Alternatives:           1 (=yes)      0 (=no) --> (use stretching from ROMS output)
    New_S_exag      = 2         # Which alternative stretching to use: 0 = higher bottom res, 1 = very much higher bottom res, 2 = alternative 3
    Gen_plot        = 1         # plotting relatively unimportant figures from all scripts
    Imp_plot        = 1         # plotting the important figures from all scripts
    Create_files    = 1         # Create all the files
    create_grid     = 1         # only for grid file (making the grid takes a lot of time) Alternatives: 1 (=yes)      0 (=no)
    VarSlope        = 0         # do we use VarSlo grid? 1 = yes, 0 = no
    
    ###############################################################################
    ###################### - Alternative stretching option - ######################
    ###############################################################################
    
    if New_S_exag == 0:
        stretching_alt = stretching_alt1
    elif New_S_exag == 1:
        stretching_alt = stretching_alt2
    elif New_S_exag == 2:
        stretching_alt = stretching_alt3
    
    ###############################################################################
    ############################# - Major parameters - ############################
    ###############################################################################
    
    z               = 30        # Depth layers                  - Always 30
    Res             = 1.5       # Resolution                    - Alternatives = 1             1.5             (3)
    X               = 600000    # Zonal dimension               - Alternatives = 650000        600000
    Y               = 400000    # Meridional dimension          - Always 400000
    steepnes_factor = 7         # Steepness of slope            - Alternatives = 7 (steep)     15 (moderate)   (25 (mild slope))
    theta_s         = 7         # Layer distribution at surface - Always 7
    theta_b         = 14        # Layer distribution at bottom  - Alternatives = 2             14
    
    ###############################################################################
    ########################### - Varying slope option - ##########################
    ###############################################################################
    
    same_top        = 1                     # same top or same middle (1 = top, 0 = middle)
    transition      = 90                    # how wide is the transition zone (only options are: 30, 40, 50, 60, 70, 80, 90, and 100)
    West_slope      = 15                    # west slope  #(7 er steep slope, 15 er modslope, 25 er mild slope)
    East_slope      = 7                     # east slope
    zonal_center    = int(X/(1000*Res*2))   # Zonal center, not an option
    
    if same_top == 0:                       # meridional shift of the (western?) slope 
        merid_shift = 0
    else: 
        merid_shift = 15.955
    
    ###############################################################################
    ########################## - Experiment composition - #########################
    ###############################################################################
    
    # if all of these are = 0 then its a ref!
    pNIJ            = 0                                       # EXP 1 == prescribed NIJ outflow in west
    pNIIC           = 0                                       # EXP 2 == prescribed NIIC inflow in west
    pNIIC_hydro     = 0                                       # EXP 2.33 == Same pNIIC but with warm water flowing in
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
    elif pNIJ == 1 and pNIIC == 0 and pNIIC_hydro == 0:
        exp = 'pNIJ'
    elif pNIJ == 1 and pNIIC == 1:
        exp = 'pNIIC&NIJ'
    elif pNIJ == 1 and pNIIC_hydro == 1:
        exp = 'pNIIC&NIJ_warm'
    elif pNIJ_md == 1 and pNIIC == 0:
        exp = 'pNIJ_md'
    elif pNIJ_md == 1 and pNIIC == 1:
        exp = 'pNIIC&NIJ_md'
    
    if Use_stretch == 1:
        statement = 'With OLD STRETCHING (Vstretchin = 2, theta_s = '+str(theta_s)+', and theta_b = '+str(theta_b)
    else:
        if New_S_exag == 0:
            statement = 'With NEW STRETCHING (Vstretchin = 3, theta_s = 4, and theta_b = 0.9)'
        elif New_S_exag == 1:
            statement = 'With NEW STRETCHING 2 (Vstretchin = 3, theta_s = 4, and theta_b = 2)'
        elif New_S_exag == 2:
            statement = 'With NEW STRETCHING 3 (Vstretchin = 3, theta_s = 4, and theta_b = 1.2)'
    
    print('\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
    print('Making run for: \n\nResolution  = '+str(Res)+' km \n\nDomain size = '+str(int(X/1000))+'x'+str(int(Y/1000))+' km \n\nExperiment  = '+slope+exp+'\n\n'+statement+'\n\nIn the directory: \n'+exp_dir+'<--')
    print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
    
    
    
    
    
    
    
    
    
    
