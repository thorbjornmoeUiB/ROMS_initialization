#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 11 16:57:55 2023

@author: thorbjorn
"""

import os
import shutil
from GlobParam import *

#VARNAME = '/cluster/projects/nn9608k/Iceland_Ideal/varinfo.dat'

KEYWORDLIST=[['RUNNAME',"'PHD_ii_HH_V2'"], 
            ['RUNID',"'TEST1_new_params/TEST_FRAMEWORK'"]]


p  = GlobParam.home_p + GlobParam.gen_dir + GlobParam.exp_dir
p1 = GlobParam.proj_p + GlobParam.gen_dir  


dest_file = p + "/runon_fram_testing.sh"
shutil.copy(p1+'/MakeRun/runon_fram_testing.sh',dest_file)

with open(p+'/runon_fram_testing.sh', 'r') as f:
    newlines = f.read()
for key,value in KEYWORDLIST:
    newlines = newlines.replace(key,value)
with open(dest_file, 'w') as f:
    f.write(newlines)

dest_file = p + "/runon_fram.sh"
shutil.copy(p1+'/MakeRun/runon_fram.sh',dest_file)

with open(p+'/runon_fram.sh', 'r') as f:
    newlines = f.read()
for key,value in KEYWORDLIST:
    newlines = newlines.replace(key,value)
with open(dest_file, 'w') as f:
    f.write(newlines)
