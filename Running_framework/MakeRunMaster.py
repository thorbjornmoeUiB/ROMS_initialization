#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 13 10:43:11 2023

@author: thorbjorn
"""

from GlobParam import * 
import os
path = os.getcwd()

exec(open(path+"/Make_oceanin.py").read())
exec(open(path+"/Make_runon_fram.py").read())
