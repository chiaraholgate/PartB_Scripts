# -*- coding: utf-8 -*-
"""
Created on Thu Apr  4 14:13:04 2019

This script creates: Table 1. Mean annual contribution of oceanic and terrestrial moisture to precipitation 
by region (1979-2013).

@author: z3131380
"""
import pandas as pd
import os.path

#==============================================================================
# Definitions
#==============================================================================
dir_in = r'/srv/ccrc/data19/z3131380/PartB/Output/Australia/100parcels/TS10min/exp01/Processed/Rainfall_recycling/Yearly/'
dir_out = r'/home/z3131380/hdrive/PhD/PartB/Reporting/Paper_1/Figures/'

regions = ['Australia','CarpentariaCoast','LakeEyreBasin','NorthEastCoast','NorthWesternPlateau','MurrayDarlingBasin',\
                'SouthWestCoast','PilbaraGascoyne','SouthAustralianGulf','SouthEastCoastNSW','SouthEastCoastVictoria',\
                'SouthWesternPlateau','TanamiTimorSeaCoast','Tasmania']
                
#==============================================================================
# Create table
#==============================================================================
df = pd.DataFrame()

for region in regions:
    fname = dir_in+region+'_annual_rainfall_recycling_1979-2013.csv'
    if os.path.isfile(fname) == True:
        print region
        df_in = pd.read_csv(fname)    
        df = df.append(pd.DataFrame({'Basin':region,'RR':[df_in['RR'].mean()],\
                    'outregion_ocean_pct':[df_in['outregion_ocean_pct'].mean()],\
                    'outregion_land_Aus_pct':[df_in['outregion_land_Aus_pct'].mean()],\
                    'terrestrial_cont':[df_in['RR'].mean() + df_in['outregion_land_Aus_pct'].mean()],\
                    'check_sum':[df_in['RR'].mean() + df_in['outregion_ocean_pct'].mean() + \
                        df_in['outregion_land_Aus_pct'].mean()]}))
outname = dir_out+'Table1_mean_annual_land_ocean_pct.csv'
df.to_csv(outname)                              