#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 14 15:12:47 2023

@author: thorbjorn
"""

from GlobalParameters import * 

#  No need to run in right order, it will stop and rerun after doing the right order :)

exec(open("Make_grid.py").read())
exec(open("Make_stre.py").read())
exec(open("Make_init.py").read())
exec(open("Make_frc.py").read())
exec(open("Make_bry.py").read())
exec(open("Make_clim.py").read())
exec(open("Make_clim_nudge.py").read())
