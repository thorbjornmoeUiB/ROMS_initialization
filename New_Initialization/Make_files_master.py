#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 14 15:12:47 2023

@author: thorbjorn
"""

from GlobalParameters import * 

exec(open("Make_grid.py").read())
exec(open("Make_stre.py").read())
exec(open("Make_init.py").read())
exec(open("Make_frc.py").read())

#%% # in case modifications are needed to bry file

exec(open("Make_bry.py").read())
exec(open("Make_clim.py").read())
exec(open("Make_clim_nudge.py").read())
exec(open("Make_readme.py").read())

