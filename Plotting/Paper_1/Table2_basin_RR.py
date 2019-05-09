# -*- coding: utf-8 -*-
"""
Created on Thu Apr  4 14:13:04 2019

This script creates: Table 2. Annual and seasonal mean recycling ratio.

The power law scaling by basin area is taken from :
Dirmeyer, P.A. and K.L. Brubaker, 2007: Characterization of the Global Hydrologic Cycle from a 
Back-Trajectory Analysis of Atmospheric Water Vapor. J. Hydrometeor., 8, 20–37, https://doi.org/10.1175/JHM557.1

**TO DO:
-Check!


@author: z3131380
"""
import pandas as pd
import os.path
import numpy as np

#==============================================================================
# Definitions
#==============================================================================
dir_in = r'/srv/ccrc/data19/z3131380/PartB/Output/Australia/100parcels/TS10min/exp01/Processed/Rainfall_recycling/'
dir_out = r'/home/z3131380/hdrive/PhD/PartB/Reporting/Paper_1/Figures/'

regions = ['CarpentariaCoast','LakeEyreBasin','NorthEastCoast','NorthWesternPlateau','MurrayDarlingBasin',\
                'SouthWestCoast','PilbaraGascoyne','SouthAustralianGulf','SouthEastCoastNSW','SouthEastCoastVictoria',\
                'SouthWesternPlateau','TanamiTimorSeaCoast','Tasmania']

# Basin area calculated in qgis from shapefile:
#/srv/ccrc/data03/z3131380/PartB/Masks/Aus_Drainage_Divisions/geofabric/Divisions/HR_Regions_river_region_Divisions.shp
dir_basin_area = r'/srv/ccrc/data03/z3131380/PartB/Masks/Aus_Drainage_Divisions/geofabric/Divisions/'
df_area = pd.read_csv(dir_basin_area+'Divisions_area_10^10_m2.csv')   

                
#==============================================================================
# Functions
#==============================================================================
# Returns rainfall recycling ratio (%) based on area (km2) and Dirmeyer & Brubaker 2007 scaling.
b = 0.457 # From Dirmeyer & Brubaker 2007
common_area = 10**5 #(km2)
def powerlaw(a,area):    
    RR = a*area**b
    return RR


#==============================================================================
# Create table
#==============================================================================
df = pd.DataFrame()
df['Basin'] = regions
df['Area_km2'] = np.nan
df['Summer'] = np.nan
df['Fall'] = np.nan
df['Winter'] = np.nan
df['Spring'] = np.nan
df['Annual'] = np.nan
df['Annual_DBscale'] = np.nan

for region in regions:
    
    # Basin area    
    area_km2 = df_area.loc[df_area['Division']==region,['Area_10^10_m2']].values[0]*10**10/10**6
    df.ix[df['Basin']==region,['Area_km2']] = area_km2
    
    fname = dir_in+'Seasonal/'+region+'_seasonal_rainfall_recycling_1979-2013.csv'
    if os.path.isfile(fname) == True:
        df_in = pd.read_csv(fname)    
        groups = df_in.groupby('Season')
        df.ix[df['Basin']==region,['Summer']] = groups.get_group('DJF')['RR'].mean()
        df.ix[df['Basin']==region,['Fall']] = groups.get_group('MAM')['RR'].mean()
        df.ix[df['Basin']==region,['Winter']] = groups.get_group('JJA')['RR'].mean()
        df.ix[df['Basin']==region,['Spring']] = groups.get_group('SON')['RR'].mean()
        
    fname = dir_in+'Yearly/'+region+'_annual_rainfall_recycling_1979-2013.csv'
    if os.path.isfile(fname) == True:
        df_in = pd.read_csv(fname)    
        df.ix[df['Basin']==region,['Annual']] = df_in['RR'].mean()
    
    df.ix[df['Basin']==region,['Annual_DBscale']] = powerlaw(df.loc[df['Basin']==region,['Annual']].values[0] / area_km2**b, common_area)
    
outname = dir_out+'Table2_basin_RR.csv'
df.to_csv(outname)                              