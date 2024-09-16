#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 14 15:12:47 2023

@author: thorbjorn
"""

from GlobParam import * 

fn = '/Users/thorbjorn/Library/CloudStorage/OneDrive-UniversityofBergen/Documents/PhD/Iceland_Sea_Mod/New_Input_scripts/GITHUB/'
#%%

# Only run for Spinup

exec(open(fn+"Make_grid.py").read())
exec(open(fn+"Make_stre.py").read())
exec(open(fn+"Make_init.py").read())
exec(open(fn+"Make_frc.py").read())

#%%

# Run for Restart

exec(open(fn+"Make_bry.py").read())
exec(open(fn+"Make_clim.py").read())
exec(open(fn+"Make_clim_nudge.py").read())
exec(open(fn+"Make_DoubleCheck.py").read())

if GlobParam.Create_files == 1:
    exec(open(fn+"Make_readme.py").read())


