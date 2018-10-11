#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 27 09:42:03 2017

This script 

(1) loads the maximum U and V wind components calculated from all
NARCliM output files over 1979-2013 in: 
    /home/z3131380/hdrive/PhD/PartB/Scripts/NARCliM_max_winds.sh

(2) calculates the max U and V components over all loaded files

(3) calculates the max resulting wind speed component

(4) calculates the time an air parcel would take, at this max windspeed, to 
travel across a 50km grid cell.

@author: z3131380
"""

import numpy as np
from netCDF4 import Dataset 

fdir=r'/srv/ccrc/data03/z3131380/PartB/Model_testing/Max_windspeed/'

U_vals = []; V_vals = [];

for year in np.arange(1979,2013):
    for month in np.arange(1,13):
        m = '%02d' % month
        file = fdir+'wrfout_d01_'+str(year)+'-'+m+'_UV_max.nc'
        fh = Dataset(file, mode='r')
        U = fh.variables['U'][:]
        V = fh.variables['V'][:]
        fh.close()  
        U_vals = np.append(U_vals,U)
        V_vals = np.append(V_vals,V)

max_U = U_vals.max()
max_V = V_vals.max()

WS_max = np.sqrt((max_U**2) + (max_V**2))

TS = (50000/WS_max)/60 #[mins]