# -*- coding: utf-8 -*-
"""
Created on Tue Apr  9 13:03:59 2019

@author: z3131380
"""

import pandas as pd
import numpy as np

dir_out = r'/srv/ccrc/data19/z3131380/PartB/Output/Australia/100parcels/TS10min/exp01/Processed/Rainfall_recycling/Yearly/'
region='Australia'
timeblock='annual'
df = pd.read_csv(dir_out+region+'_'+timeblock+'_rainfall_recycling_1979-2013.csv')
df['Incoming_mm'] = df['outregion_ocean_mm'] + df['region_land_mm']
lowpctl = np.percentile(df['Incoming_mm'], 20)
highpctl = np.percentile(df['Incoming_mm'], 80)
df['pctl'] = "in between"
df['pctl'] = np.where(df['Incoming_mm']<=lowpctl,"<=lowpctl",df['pctl'])
df['pctl'] = np.where(df['Incoming_mm']>=highpctl,">=highpctl",df['pctl'])
df['Type'] = "average"
df['Type'] = np.where(df['P_anom']<=-0.1,"Dry",df['Type'])
df['Type'] = np.where(df['P_anom']>=0.1,"Wet",df['Type'])       
        
        
