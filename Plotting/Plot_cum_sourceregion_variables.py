# -*- coding: utf-8 -*-
"""
Created on Tue Apr  2 10:46:14 2019

This script plots cumulative monthly totals for NARCliM variables following 
"Roy, T., Martinez, J. A., Herrera-Estrada, J. E., Zhang, Y., Dominguez, F., Berg, A., et al. (2018). 
Role of Moisture Transport and Recycling in Characterizing Droughts: Perspectives from Two Recent 
U.S. Droughts and the CFSv2 System. Journal of Hydrometeorology, 20(1), 139–154. 
https://doi.org/10.1175/JHM-D-18-0159.1"

Variables include
(1) Evaporation converted from latent heat 
(2) Total precipitable water 
 
 

@author: z3131380
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from netCDF4 import Dataset 
from mpl_toolkits.basemap import Basemap

#==============================================================================
# Definitions
#==============================================================================
region = 'MurrayDarlingBasin'

dir_in = r'/srv/ccrc/data19/z3131380/PartB/Output/Australia/100parcels/TS10min/exp01/Processed/Sourceregion_variable_timeseries/'

focus_years = [1999,2002,2006,2010]

fsize = 14

#==============================================================================
# Load data
#==============================================================================
# Evaporation & evapotranspiration 
E_regions = ['Arafura','Pacific','Southern','Indian']
df_E_land = pd.read_csv(dir_in+'Monthly_E_landregions_1979-2013.csv')
E_land_groups = df_E_land.groupby('Year')
df_E_oceans = pd.read_csv(dir_in+'Monthly_E_oceanregions_1979-2013.csv')
E_ocean_groups = df_E_oceans.groupby('Year')

# Total precipitable water
df_TPW_land = pd.read_csv(dir_in+'Monthly_TPW_landregions_1979-2013.csv')
TPW_land_groups = df_TPW_land.groupby('Year')
df_TPW_oceans = pd.read_csv(dir_in+'Monthly_TPW_oceanregions_1979-2013.csv')
TPW_ocean_groups = df_TPW_oceans.groupby('Year')

# Rainfall
dir_rain = r'/srv/ccrc/data19/z3131380/PartB/Output/Australia/100parcels/TS10min/exp01/Processed/Yearly/'
fh = Dataset(dir_rain+region+'_1980_wvcont.nc', mode='r') 
rain = fh.variables['pre'][:]
fh.close()

#E_monthly_1979_CC = df_E_oceans.ix[:10,'CarpentariaCoast']

#==============================================================================
# Configure data
#==============================================================================
# Get sum of E over each ocean area for each focus year
#E_Arafura = []; E_Pacific = []; E_Indian = []; E_Southern = []
#for y in focus_years:
#    E_Arafura.append(E_ocean_groups.get_group(y)['Arafura'].sum())
#    E_Pacific.append(E_ocean_groups.get_group(y)['Pacific'].sum())
#    E_Indian.append(E_ocean_groups.get_group(y)['Indian'].sum())
#    E_Southern.append(E_ocean_groups.get_group(y)['Southern'].sum())

E_1999 = []; E_2002 = []; E_2006 = []; E_2010 = []
for r in E_regions:
    E_1999.append(E_ocean_groups.get_group(1999)[r].sum()/1000)
    E_2002.append(E_ocean_groups.get_group(2002)[r].sum()/1000)
    E_2006.append(E_ocean_groups.get_group(2006)[r].sum()/1000)
    E_2010.append(E_ocean_groups.get_group(2010)[r].sum()/1000)
E_1999.append(E_land_groups.get_group(1999)[region].sum()/1000*10) # to scale on same axis
E_2002.append(E_land_groups.get_group(2002)[region].sum()/1000*10)
E_2006.append(E_land_groups.get_group(2006)[region].sum()/1000*10)
E_2010.append(E_land_groups.get_group(2010)[region].sum()/1000*10)

#TPW_1999 = []; TPW_2002 = []; TPW_2006 = []; TPW_2010 = []
#for r in E_regions:
#    TPW_1999.append(TPW_ocean_groups.get_group(1999)[r].sum()/1000/1e3)
#    TPW_2002.append(TPW_ocean_groups.get_group(2002)[r].sum()/1000/1e3)
#    TPW_2006.append(TPW_ocean_groups.get_group(2006)[r].sum()/1000/1e3)
#    TPW_2010.append(TPW_ocean_groups.get_group(2010)[r].sum()/1000/1e3)
#TPW_1999.append(TPW_land_groups.get_group(1999)[region].sum()/1000/1e3*10) # to scale on same axis
#TPW_2002.append(TPW_land_groups.get_group(2002)[region].sum()/1000/1e3*10)
#TPW_2006.append(TPW_land_groups.get_group(2006)[region].sum()/1000/1e3*10)
#TPW_2010.append(TPW_land_groups.get_group(2010)[region].sum()/1000/1e3*10)
TPW_region = np.zeros([12,4])*np.nan
for y in focus_years:
    TPW_region[:,focus_years.index(y)] = TPW_land_groups.get_group(y)[region]/1000

#==============================================================================
# Figure
#==============================================================================
fig, ax1 = plt.subplots(figsize=(16,11.7))        
ax1 = plt.subplot2grid((3, 1), (0, 0), rowspan=1)
plt.annotate('(a)',xy=(-0.2,4500),xycoords='data',fontsize=fsize)
#plt.plot(E_groups.get_group(2002)['Month'],np.cumsum(E_groups.get_group(2002)[region]),label='2002')
#plt.legend(bbox_to_anchor=(1.3,1))
index = np.arange(5)
bar_width = 0.2
plt.bar(index, E_1999 , bar_width, color='none', hatch='.', align='center', label='1999')
plt.bar(index+bar_width*1, E_2002 , bar_width, color='none', hatch='-', align='center', label='2002')
plt.bar(index+bar_width*2, E_2006 , bar_width, color='none', hatch='x', align='center', label='2006')
plt.bar(index+bar_width*3, E_2010 , bar_width, color='none', hatch='*', align='center', label='2010')
plt.xlim(-0.5,5.5)
plt.xticks(index + bar_width, ('Arafura','Pacific','Southern','Indian','MDB x 10'))
plt.ylabel('Evaporation [m]')
plt.legend()

ax2 = plt.subplot2grid((3, 1), (1, 0), rowspan=1)
plt.annotate('(b)',xy=(0.75,4100),xycoords='data',fontsize=fsize)
#plt.bar(index, TPW_1999 , bar_width, color='none', hatch='.', align='center', label='1999')
#plt.bar(index+bar_width*1, TPW_2002 , bar_width, color='none', hatch='-', align='center', label='2002')
#plt.bar(index+bar_width*2, TPW_2006 , bar_width, color='none', hatch='x', align='center', label='2006')
#plt.bar(index+bar_width*3, TPW_2010 , bar_width, color='none', hatch='*', align='center', label='2010')
#plt.xlim(-0.5,5.5)
#plt.xticks(index + bar_width, ('Arafura','Pacific','Southern','Indian','MDB x 10'))
#plt.ylabel('TPW [m] x 10$^3$')

#plt.plot(range(12),TPW_region)
#plt.legend(focus_years)
#plt.ylabel('TPW [m]')

plt.plot(TPW_land_groups.get_group(1999)['Month'],np.cumsum(TPW_land_groups.get_group(1999)[region]),label='1999')
plt.plot(TPW_land_groups.get_group(2002)['Month'],np.cumsum(TPW_land_groups.get_group(2002)[region]),label='2002')
plt.plot(TPW_land_groups.get_group(2006)['Month'],np.cumsum(TPW_land_groups.get_group(2006)[region]),label='2006')
plt.plot(TPW_land_groups.get_group(2010)['Month'],np.cumsum(TPW_land_groups.get_group(2010)[region]),label='2010')
plt.legend(loc='best')


file = dir_in+region+'_1979-2013_seasonal_climatology.nc'
fh = Dataset(file, mode='r') 
# Get rid of boundaries for contour sake
row_start = 5; row_end = -5
col_start = 5; col_end = -50
#wvcont_pct_seasonal_climatology = fh.variables['wvcont_pct_seasonal_climatology']\
#    [:,row_start:row_end,col_start:col_end]
wvcont_pct_annual_climatology = fh.variables['wvcont_pct_annual_climatology'][5:-5,5:-50]
latitcrs = fh.variables['latitcrs'][row_start:row_end,col_start:col_end]
longicrs = fh.variables['longicrs'][row_start:row_end,col_start:col_end]
latitcrs_original = fh.variables['latitcrs'][:]
longicrs_original = fh.variables['longicrs'][:]
fh.close()      

file = '/srv/ccrc/data19/z3131380/PartB/Output/Australia/100parcels/TS10min/exp01/Processed/Climatology/Ocean_polygons/Ocean_source_regions.nc'
fh = Dataset(file, mode='r') 
wsmask_Pacific = fh.variables['wsmask_Pacific']
wsmask_Indian = fh.variables['wsmask_Indian']
wsmask_Southern = fh.variables['wsmask_Southern']
wsmask_Arafura = fh.variables['wsmask_Arafura']

plt.clf()    
fig = plt.figure(figsize=(16,11.7)) 
ax = fig.add_subplot(111)
m = Basemap(projection='stere',lon_0=135,lat_0=-25.,\
                llcrnrlat=-45,urcrnrlat=0,\
                llcrnrlon=85,urcrnrlon=175,\
                resolution='l',area_thresh=10000)
x1, y1 = m(longicrs,latitcrs)
x,y = m(longicrs_original,latitcrs_original) 
m.drawcoastlines()
h = m.pcolormesh(x1,y1,wvcont_pct_annual_climatology)
a=m.contourf(x,y,wsmask_Pacific)
b=m.pcolormesh(x,y,wsmask_Indian)
c=m.pcolormesh(x,y,wsmask_Southern)
d=m.pcolormesh(x,y,wsmask_Arafura)