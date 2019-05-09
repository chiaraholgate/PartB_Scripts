# -*- coding: utf-8 -*-
"""
Created on Tue Apr  9 13:46:06 2019

Figure 4. Interannual terrestrial and oceanic contributions to Australian precipitation in (a) relative terms 
[%] and (b) absolute terms [mm].

@author: z3131380
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

#==============================================================================
# Definitions
#==============================================================================
dir_in = r'/srv/ccrc/data19/z3131380/PartB/Output/Australia/100parcels/TS10min/exp01/Processed/Rainfall_recycling/'
dir_out = r'/home/z3131380/hdrive/PhD/PartB/Reporting/Paper_1/Figures/'

region = 'Australia'
df = pd.read_csv(dir_in+'Yearly/'+region+'_annual_rainfall_recycling_1979-2013.csv')

fsize = 14

#==============================================================================
# Functions
#==============================================================================
def make_patch_spines_invisible(ax):
    ax.set_frame_on(True)
    ax.patch.set_visible(False)
    for sp in ax.spines.values():
        sp.set_visible(False)
        

#==============================================================================
# Figure 4. Interannual terrestrial and oceanic contributions to Australian precipitation in (a) relative terms 
# [%] and (b) absolute terms [mm].
#==============================================================================

fig, ax2 = plt.subplots(figsize=(16,11.7))        
ax2 = plt.subplot2grid((2, 1), (0, 0), rowspan=1)
plt.annotate('(a)',xy=(1980,0.4),xycoords='data',fontsize=fsize)
#p3, = ax2.plot(df['Year'],df['P_anom'],'-k^', label="Rainfall anomaly [%]",markersize=4,linewidth=1,alpha=0.7)
ax2.fill_between(df['Year'],df['P_anom'],0,alpha=0.5,edgecolor='None',facecolor='gray')
ax2.set_ylim(-0.5,0.5)
ax2.grid(axis='x')
ax2.tick_params(labelbottom=False,labeltop=True) 
ax2.axhline(0,color='gray',linestyle=':')
ax2.set_yticklabels(ax2.get_yticks(),fontsize=fsize)
ax2.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))   
plt.annotate('Precipitation anomaly [%]',xy=(-0.05,0.7),xycoords='axes fraction',fontsize=fsize,rotation=90)
ax2a = ax2.twinx()
ax2b = ax2.twinx()
p2a, = ax2a.plot(df['Year'],df['outregion_ocean_pct'],'-bo', label="Ocean water vapour contribution [%]",markersize=4,markeredgecolor='b',markerfacecolor='b')#,alpha=0.5)
p2b, = ax2b.plot(df['Year'],df['RR'],'-gs', label="Land water vapour contribution [%]",markersize=4,markeredgecolor='g',markerfacecolor='g')#,alpha=0.5)
ax2a.spines["right"].set_position(("axes", 1.01))
ax2b.spines["right"].set_position(("axes", 1.07))
make_patch_spines_invisible(ax2a)
make_patch_spines_invisible(ax2b)
ax2a.spines["right"].set_visible(True)
ax2b.spines["right"].set_visible(True)
tkw = dict(size=4, width=1.5)
ax2a.tick_params(axis='y', colors=p2a.get_color(), **tkw)
ax2b.tick_params(axis='y', colors=p2b.get_color(), **tkw)
ax2a.set_yticklabels(ax2a.get_yticks(),fontsize=fsize)
ax2a.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))   
ax2b.set_yticklabels(ax2b.get_yticks(),fontsize=fsize)
ax2b.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))  
ax2a.set_ylim(70,95)
ax2b.set_ylim(5,25)
plt.annotate("Oceanic contribution [%]",xy=(1.04,0.7),xycoords = 'axes fraction',rotation=90,color='blue',fontsize=fsize)
plt.annotate("Terrestrialcontribution [%]",xy=(1.1,0.74),xycoords = 'axes fraction',rotation=90,color='green',fontsize=fsize)
ax2.set_xlim(1979,2014)    
ax2.set_xticks(np.arange(min(df['Year']), max(df['Year'])+1, 1))
ax2.set_xticklabels(ax2.get_xticks(),rotation=90,fontsize=fsize)
ax2.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))    

 
ax1 = plt.subplot2grid((2, 1), (1, 0), sharex=ax2)
plt.annotate('(b)',xy=(1980,2700),xycoords='data',fontsize=fsize)
#p1, = ax1.plot(df['Year'],df['P_total']/1000,'-k^', label="Rainfall total [m]",markersize=4,linewidth=1,alpha=0.7)
ax1.fill_between(df['Year'],df['P_total']/1000,0,alpha=0.5,edgecolor='None',facecolor='gray')
#ax1.set_ylim(0,2000)
ax1.grid(axis='x')
#ax1.tick_params(labelbottom=False,labeltop=True) 
ax1.set_xlim(1979,2014) 
ax1.set_ylim(500,3000)   
ax1.set_xticks(np.arange(min(df['Year']), max(df['Year'])+1, 1))
ax1.set_xticklabels(ax1.get_xticks(),rotation=90,fontsize=fsize)
ax1.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))   
ax1.set_yticklabels(ax1.get_yticks(),fontsize=fsize)
ax1.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))   
plt.annotate('Precipitation total [m]',xy=(-0.06,0.65),xycoords='axes fraction',fontsize=fsize,rotation=90)
ax1a = ax1.twinx()
ax1b = ax1.twinx()
p1a, = ax1a.plot(df['Year'],df['outregion_ocean_mm']/1000,'-bo', label="Ocean water vapour contribution [m]",markersize=4,markeredgecolor='b',markerfacecolor='b')
p1b, = ax1b.plot(df['Year'],df['region_land_mm']/1000,'-gs', label="Land water vapour contribution [m]",markersize=4,markeredgecolor='g',markerfacecolor='g')
ax1a.spines["right"].set_position(("axes", 1.01))
ax1b.spines["right"].set_position(("axes", 1.07))
make_patch_spines_invisible(ax1a)
make_patch_spines_invisible(ax1b)
ax1a.spines["right"].set_visible(True)
ax1b.spines["right"].set_visible(True)
tkw = dict(size=4, width=1.5)
ax1a.tick_params(axis='y', colors=p1a.get_color(), **tkw)
ax1b.tick_params(axis='y', colors=p1b.get_color(), **tkw)
ax1a.set_yticklabels(ax1a.get_yticks(),fontsize=fsize)
ax1a.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))   
ax1b.set_yticklabels(ax1b.get_yticks(),fontsize=fsize)
ax1b.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))   
ax1a.set_ylim(600,1800)
ax1b.set_ylim(100,400)
plt.annotate("Oceanic contribution [m]",xy=(1.05,0.7),xycoords = 'axes fraction',rotation=90,color='blue',fontsize=fsize)
plt.annotate("Terrestrial contribution [m]",xy=(1.105,0.75),xycoords = 'axes fraction',rotation=90,color='green',fontsize=fsize)
   
   
plt.tight_layout()
fig.savefig(dir_out+'Fig4.png',bbox_inches='tight') 
plt.close()
    
