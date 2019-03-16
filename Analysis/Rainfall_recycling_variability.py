# -*- coding: utf-8 -*-
"""
Created on Thu Dec 13 15:12:51 2018

@author: z3131380
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

dir_out = '/srv/ccrc/data19/z3131380/PartB/Output/Australia/100parcels/TS10min/Processed/Rainfall_recycling/'
domain = 'Australia'
region = 'Australia'
timeblock = 'annual'

#==============================================================================
# Load rainfall recycling data frame
#==============================================================================
dataf = pd.read_csv(dir_out+region+'_'+timeblock+'_rainfall_recycling_1979-2013.csv')

#==============================================================================
# Q: What proportion of the rainfall interannual variability comes from the land and ocean?
# A: Australia: 25% land, 75% ocean. Error ~0.3%
#     MDB: 13% land, 74% ocean, 14% outregion land (Aus). Error ~ -1.4% 
#     SWWA: 3% land, 92% ocean, 9% outregion land (Aus). Error -4%
#==============================================================================
df = pd.DataFrame()
df['Year'] = dataf['Year']
P_mm_mean = dataf['P_total'].mean()
region_land_mm_mean = dataf['region_land_mm'].mean()
outregion_land_Aus_mm_mean = dataf['outregion_land_Aus_mm'].mean()
outregion_ocean_mm_mean = dataf['outregion_ocean_mm'].mean()
df['P_mm_anom'] = dataf['P_total'] - P_mm_mean
df['region_land_mm_anom'] = dataf['region_land_mm'] - region_land_mm_mean
df['outregion_land_Aus_mm_anom'] = dataf['outregion_land_Aus_mm'] - outregion_land_Aus_mm_mean
df['outregion_ocean_mm_anom'] = dataf['outregion_ocean_mm'] - outregion_ocean_mm_mean
df['zero_line'] = 0

# Find area under curves - assume width of rectangle is 1?
P_mm_anom_area =sum(np.abs(df['P_mm_anom']) - df['zero_line']) # *1 etc
region_land_mm_anom_area = sum(np.abs(df['region_land_mm_anom']) - df['zero_line'])
outregion_land_Aus_mm_anom_area = sum(np.abs(df['outregion_land_Aus_mm_anom']) - df['zero_line'])
outregion_ocean_mm_anom_area = sum(np.abs(df['outregion_ocean_mm_anom']) - df['zero_line'])


prop_from_land = region_land_mm_anom_area/P_mm_anom_area*100
prop_from_outregion_land = outregion_land_Aus_mm_anom_area/P_mm_anom_area*100
prop_from_ocean = outregion_ocean_mm_anom_area/P_mm_anom_area*100
check_sum = (P_mm_anom_area - \
    np.nansum([region_land_mm_anom_area,outregion_ocean_mm_anom_area,outregion_land_Aus_mm_anom_area]))/P_mm_anom_area*100

# Put this plot in its own script.
plt.plot(df['Year'],df['P_mm_anom'],'k');plt.plot(df['Year'],df['region_land_mm_anom'],'g');plt.plot(df['Year'],df['outregion_ocean_mm_anom'],'b');plt.plot(df['Year'],dataf['outregion_land_Aus_mm'],'c') 
plt.savefig(dir_out+region+'_1979-2013_'+timeblock+'_anomalies.png',bbox_inches='tight') 

#==============================================================================
# CDF and PDF of ocean and land % contributions (as anomalies) on an annual basis
#==============================================================================
# Histogram
df.plot(y='RR_anom',kind='hist'); df.plot(y='outregion_ocean_pct_anom',kind='hist'); df.plot(y='P_anom',kind='hist')
df.plot(y='P_total',kind='hist'); df.plot(y='region_land_mm',kind='hist'); df.plot(y='outregion_ocean_mm',kind='hist')
# Normalised to 1
df.plot(y='RR_anom',kind='hist',normed=True,bins=10)
# CDF
df.plot(y='RR_anom',kind='hist',normed=True,bins=10,cumulative=True)


#==============================================================================
# Variance and standard deviation of ocean and land contribution samples
#==============================================================================
df.describe()
#RR_var = df['RR'].std()**2
#outregion_ocean_pct_var = df['outregion_ocean_pct'].std()**2
#P_anom_var = df['P_anom'].std()**2
P_mm_var = df['P_total'].std()**2
land_mm_var = df['region_land_mm'].std()**2
ocean_mm_var = df['outregion_ocean_mm'].std()**2

prop_var_from_land = land_mm_var/P_mm_var*100
prop_var_from_ocean = ocean_mm_var/P_mm_var*100
