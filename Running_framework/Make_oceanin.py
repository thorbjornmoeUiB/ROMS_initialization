#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 11 15:50:12 2023

@author: th
"""

import os
import shutil
from GlobParam import *

proj_dir = GlobParam.proj_p + GlobParam.gen_dir + GlobParam.exp_dir

KEYWORDLIST=[['APPTITLE',"ROMS 3.8 IceSlope, CATS2008 k1o1m2s2, initializing preliminary ROMS run"], 
            ['MYAPPCPPNAME',"Iceland_Ideal"],
            ['VARFILE',"/cluster/projects/nn9608k/Iceland_Ideal/varinfo.dat"],
            
            ['XPOINTS',GlobParam.xpoint],  
            ['YPOINTS',GlobParam.ypoint],
            ['NLEVELS',"30"], 
            
            ['XCPU',GlobParam.tileX],
            ['YCPU',GlobParam.tileY],
            
            ['TSTEPS',"1728000.0"],
            ['DELTAT',"180"],
            ['RATIO',"25"],
            
            ['IRESTART',"0"],
            ['LcycleRST1',"T"],
            ['RSTSTEP',"172800.0"],
            ['STASTEP',"20.0"],
            ['NFLT1',"1"],
            ['INFOSTEP',"10.0"],
            
            ['_TNU2_',"1.0d0 1.0d0"],
            ['_TNU4_',"1.0d0 1.0d0"],
            ['_TNU2_ad',"1.0d0 1.0d0"], #1 1
            ['_TNU4_ad',"1.0d0 1.0d0"], #1 1
            ['_VISC2_',"1.0d0"],
            ['_VISC4_',"1.0d0"],
            ['_VISC2_ad',"1.0d0"], #1 1
            ['_VISC4_ad',"1.0d0"], #1 1
            
            ['V_TRANS',"2"],
            ['V_STRETCH',GlobParam.Vstretch],
            ['GRDTHETAS',GlobParam.thetaS],
            ['GRDTHETAB',GlobParam.thetaB],
            ['GRDTCLINE',GlobParam.Tcline],
            
            ['STARTTIME',"0"],
            ['TIDEREFd0',"1d0"],
            ['TIMEREF',"-1"],
            
            ['_TNUDG_',"2*9.0d0"],
            ['_ZNUDG_',"180.0d0"],
            ['_M2NUDG_',"2*9.0d0"],
            ['_M3NUDG_',"2*9.0d0"],
            
            ['_GRDNAME_',proj_dir + GlobParam.grdFN],
            ['_ININAME_',proj_dir + GlobParam.iniFN],
            ['_BRYNAME_',proj_dir + GlobParam.bryFN],
            ['_CLMNAME_',proj_dir + GlobParam.clmFN],
            ['_NUDNAME_',proj_dir + GlobParam.nudFN],
            ['_FRCNAME_',proj_dir + GlobParam.frcFN],          
]


p_old = GlobParam.proj_p + GlobParam.gen_dir + '/MakeRun'
p_new = GlobParam.work_p + GlobParam.gen_dir + GlobParam.exp_dir


dest_file = p_new + "/ocean.in"
shutil.copy(p_old+'/ocean.in',dest_file)

with open(p_old+'/ocean.in', 'r') as f:
    newlines = f.read()
for key,value in KEYWORDLIST:
    newlines = newlines.replace(key,value)
with open(dest_file, 'w') as f:
    f.write(newlines)
    
    
shutil.copy(p_old+'/romsM',p_new+'/romsM')

        
