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

    common_dir      = '.../Iceland_Sea_Mod/'                    # Home directory for the model
    EXP_dir         = common_dir + 'CONFIG/INPUT/'+'MODERATE/'    # Change for the slope configuration
    SubEXP_dir      = 'REF/'                                    # Change for experiments (e.g., IN_OUTFLOW)
    ref_dir         = 'REF/'                                    # Never change
    
    ###############################################################################
    ########################WINT## - Name of all used files - #########################
    ###############################################################################
    
    grid_file       =  ref_dir    + 'grid_file.nc'              # Files strictly for SPINUPS
    stretching      =  ref_dir    + 'Stetching.nc'              # Files strictly for SPINUPS
    init_file       =  ref_dir    + 'init_file.nc'              # Files strictly for SPINUPS
    frc_file        =  ref_dir    + 'frc_file.nc'               # Files strictly for SPINUPS
    bry_file        =  SubEXP_dir + 'bry_file.nc'               # Files for every individual experiment (and spuinup)
    clm_file        =  SubEXP_dir + 'clm_file.nc'               # Files for every individual experiment (and spuinup)
    clm_nudge_file  =  SubEXP_dir + 'clm_nudge_file.nc'         # Files for every individual experiment (and spuinup)
    
    
    Streching       = 'STRETCHING.nc'                   # Stretching file
    
    Mean_profiles_w = 'ProfilesMean_WINTER.nc'          # Hydrography: mean winter
    Raw_profiles_w  = 'ProfilesAll_WINTER.mat'          # Hydrography: all winter
    
    Mean_profiles_s = 'ProfilesMean_SUMMER.mat'         # Hydrography: mean summer
    Raw_profiles_s  = 'ProfilesAll_SUMMER.mat'          # Hydrography: all summer

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
    Author          = 'TOeM'
    
    ###############################################################################
    ############################# - Logical switches - ############################
    ###############################################################################
    
    hydrography     ='w'        # Alternatives: 'w' == winter hydrography, 's' == summer hydrography
    Gen_plot        = 1         # Plotting relatively unimportant figures from "Make_*" scripts
    Imp_plot        = 1         # Plotting the important figures from "Make_*" scripts
    Create_files    = 1         # Create all the files (turn off when testing)
    create_grid     = 1         # Only for grid file (making the grid takes a lot of time) Alternatives: 1 (=yes)      0 (=no)
    VarSlope        = 0         # Do we use VarSlo grid? 1 = yes, 0 = no
    Tilted_shelf    = 0         # Include an inclined plane on the shelf. Alternatives: 1 (=yes)      0 (=no) (not included in paper)
    x2_domain       = 1         # Is the domain doubled relative to original configuration of 600x400. if 1 (=yes) the domain is 1200x400
    
    ###############################################################################
    ############################# - NIIC properties - #############################
    ###############################################################################
    
    # must be applied along with the pNIIC_f and/or pNIIC_w below
    # NB! the temperature of the Inflow is simply set as a float/int and is directly used for temperature,
    #     while the speed of the Inflow must be given as a str ('2' or '3') due to volume transport conserns
    
    temp_of_wNIIC   = 0         # Temperature of "Atlantic Water" in Inflow (i.e., 3cInflow, 6cInflow)
    speed_of_wNIIC  = '2'       # Speed of inflow: Alternatives '2' (=0.2 m/s) or '3' (=0.3 m/s) (i.e., 2xInflow or 3xInflow)
    
    # ^ MAKE SURE TO SEE IF VOLUME TRANSPORT OF INFLOW IS CORRECT
    ###############################################################################
    ############################### - STEEPNESS - #################################
    ###############################################################################
    
    steepnes_factor = 15        # Steepness of slope. Alternatives = 7 (steep)     15 (moderate and VarSlo)
    
    ###############################################################################
    ########################## - Experiment composition - #########################
    ###############################################################################
    
    # if all of these are = 0 then its a REF!
    
    #variable       #switch     # experiment name               # Explanation
    pNIJ            = 0         # = Outflow                     # prescribed NIJ outflow in west
    pNIIC           = 0         # = In- & Outflow               # prescribed NIIC inflow in west
    pNIIC_f         = 0         # = 2/3xInflow                  # prescribed NIIC inflow in west that is 2x faster than normal (pNIIC)
    pNIIC_w         = 0         # = 3/6cInflow                  # prescribed NIIC inflow in west that is warmer than normal (pNIIC)
    pIFSJ           = 0         # = not included in paper       # prescribed IFSJ outflow in EAST NB!

    ###############################################################################
    ##################### - Major parameters (Not changed) - ######################
    ###############################################################################
    
    z               = 30        # Depth layers                  - Always 30
    Res             = 1         # Resolution                    - Alternatives = 1             1.5             (3)
    X               = 600000    # Zonal dimension
    if x2_domain == 1:
        X           = X*2
    Y               = 400000    # Meridional dimension          - Always 400000
    theta_s         = 7         # Layer distribution at surface - Always 7
    theta_b         = 14        # Layer distribution at bottom  - Alternatives = 2             14
    
    ###############################################################################
    ########################### - Varying slope option - ##########################
    ###############################################################################
    
    same_top        = 1                     # same top or same middle (1 = top, 0 = middle)
    transition      = 90                    # how wide is the transition zone (only options are: 30, 40, 50, 60, 70, 80, 90, and 100)
    West_slope      = 15                    # west slope  #(7 er steep slope, 15 er modslope, 25 er mild slope)
    East_slope      = 7                     # east slope
    zonal_center    = int(X/(1000*Res*2))   # Zonal center, switch to int(X/(1000*Res*4)) for 2x domain
    
    if same_top == 0:                       # meridional shift of the (western?) slope 
        merid_shift = 0
    else: 
        merid_shift = 15.955
        
    ###############################################################################
    ############## - Shelf tilting options (not included in paper) - ##############
    ###############################################################################
    
    tilt = 100          # min = 100, med = 115, max = 130
    tilt_loc = 66       # min = 50,  med = 66,  max = 83

    
    if Tilted_shelf == 1:
        ts = 'TS_'
    else:
        ts = ''
    
    
    if VarSlope == 1:
        slope = 'VarSlo_'
    else:
        if steepnes_factor == 7:
            slope = 'Steep_'
        else:
            slope =  'Moderate_'
    
    if hydrography == 'w':
        hyd = 'Winter_'
    else:
        hyd = 'Summer_'
        
    exp = '' 
    if pNIJ == 0 and pNIIC == 0 and pNIIC_f == 0 and pNIIC_w == 0 and pIFSJ == 0:
        exp = 'REF'
    if pNIJ == 1 and pNIIC == 0 and pNIIC_f == 0 and pNIIC_w == 0 and pIFSJ == 0:
        exp = 'OUTFLOW'
    if pNIJ == 0 and pNIIC == 1 and pNIIC_f == 0 and pNIIC_w == 0 and pIFSJ == 0:
        exp = 'INFLOW'
    if pNIJ == 1 and pNIIC == 1 and pNIIC_f == 0 and pNIIC_w == 0 and pIFSJ == 0:
        exp = 'IN_OUTFLOW'
    if pNIJ == 1 and pNIIC == 0 and pNIIC_f == 1 and pNIIC_w == 0 and pIFSJ == 0:
        if speed_of_wNIIC == 0.2:
            exp = '2xINFLOW'
        elif speed_of_wNIIC == 0.3:
            exp = '3xINFLOW'
        else:
            exp = str(speed_of_wNIIC*10)+'xINFLOW'
    if pNIJ == 1 and pNIIC == 0 and pNIIC_f == 0 and pNIIC_w == 1 and pIFSJ == 0:
        if temp_of_wNIIC == 3:
            exp = '3cINFLOW'
        elif temp_of_wNIIC == 6:
            exp = '6cINFLOW'
        else:
            exp = str(temp_of_wNIIC)+'cINFLOW'
    if pIFSJ == 1:
        exp += '_EASTOUTFLOW'
    
    

    statement = 'With STRETCHING: Vstretchin = 3, theta_s = 4, and theta_b = 1.2'
    
    print('\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
    print('Making run for: \n\nResolution  = '+str(Res)+' km \n\nDomain size = '+str(int(X/1000))+'x'+str(int(Y/1000))+' km \n\nExperiment  = '+ts+slope+hyd+exp+'\n\n'+statement+'\n\nIn the directory: \n'+EXP_dir+SubEXP_dir+'<--')
    print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
    
    
    
    
    
    
    
    
    
