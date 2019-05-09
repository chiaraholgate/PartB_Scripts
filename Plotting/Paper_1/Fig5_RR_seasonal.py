# -*- coding: utf-8 -*-
"""
Created on Tue Apr  9 13:46:06 2019

Figure 5. Seasonal terrestrial and oceanic contributions to Australian precipitation in (a) relative terms [%] 
and (b) absolute terms [mm].

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
df = pd.read_csv(dir_in+'Seasonal/'+region+'_seasonal_rainfall_recycling_1979-2013.csv')
groups = df.groupby('Season')
#yr_groups = df.groupby(df.Year.astype(int))
#RR_seas_mean = [groups.get_group('DJF')['RR'].mean(),\
#                            groups.get_group('MAM')['RR'].mean(),\
#                            groups.get_group('JJA')['RR'].mean(),\
#                            groups.get_group('SON')['RR'].mean()]
#RR_seas_min = [groups.get_group('DJF')['RR'].min(),\
#                            groups.get_group('MAM')['RR'].min(),\
#                            groups.get_group('JJA')['RR'].min(),\
#                            groups.get_group('SON')['RR'].min()]
#RR_seas_max = [groups.get_group('DJF')['RR'].max(),\
#                            groups.get_group('MAM')['RR'].max(),\
#                            groups.get_group('JJA')['RR'].max(),\
#                            groups.get_group('SON')['RR'].max()]     
#RR_seas_25th = [groups.get_group('DJF')['RR'].quantile(q=0.25),\
#                            groups.get_group('MAM')['RR'].quantile(q=0.25),\
#                            groups.get_group('JJA')['RR'].quantile(q=0.25),\
#                            groups.get_group('SON')['RR'].quantile(q=0.25)]  
#RR_seas_75th = [groups.get_group('DJF')['RR'].quantile(q=0.75),\
#                            groups.get_group('MAM')['RR'].quantile(q=0.75),\
#                            groups.get_group('JJA')['RR'].quantile(q=0.75),\
#                            groups.get_group('SON')['RR'].quantile(q=0.75)]                              
                            
# Open annual df to plot comparison
df2 = pd.read_csv(dir_in+'Yearly/'+region+'_annual_rainfall_recycling_1979-2013.csv')

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
# Plot seasonal RR for selected years
#==============================================================================
#
#
## Pick wet and dry years
#sel_yrs = [1986,1989,1994,1999,2002,2010]
#
#fig, ax1 = plt.subplots()#figsize=(16,11.7))      
#for y in sel_yrs:#np.arange(1979,2014):
#    ax1.plot(range(4),yr_groups.get_group(y)['RR'],':o', label=y,markersize=4,linewidth=2)
#    plt.annotate(y,xy=(3.0,yr_groups.get_group(y)['RR'].reset_index(drop=True)[3]),fontsize=12)#xycoords='data',

        
#==============================================================================
# Figure 5. Seasonal terrestrial and oceanic contributions to Australian precipitation in (a) relative terms 
# [%] and (b) absolute terms [mm], (c) mean seasonal precipitation recycling.
#==============================================================================

fig, ax2 = plt.subplots(figsize=(16,11.7))        
ax2 = plt.subplot2grid((3, 1), (0, 0), rowspan=1)
plt.annotate('(a)',xy=(1980,1.3),xycoords='data',fontsize=fsize)
#p3, = ax2.plot(df['Year'],df['P_anom'],'-k', label="Rainfall anomaly [%]",linewidth=1.5)#,alpha=0.7)
ax2.fill_between(df['Year'],df['P_anom'],0,alpha=0.5,edgecolor='None',facecolor='gray')
# NOTE: the annual time series has been moved to show its values in the MIDDLE of each year. This differs
# from the purely annual plots, with show values at the START of each year. Careful when analysing.
#p2, = ax1.plot(df2['Year']+0.5,df2['P_anom'],':', label="Rainfall anomaly [%]",markersize=4,linewidth=2,color='gray')
ax2.set_ylim(-1.5,1.5)
ax2.grid(axis='x')
ax2.tick_params(labelbottom=False,labeltop=True) 
ax2.axhline(0,color='gray',linestyle=':')
ax2.set_yticklabels(ax2.get_yticks(),fontsize=fsize)
ax2.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))   
plt.annotate('Precipitation anomaly [%]',xy=(-0.06,0.85),xycoords='axes fraction',fontsize=fsize,rotation=90)
ax2a = ax2.twinx()
ax2b = ax2.twinx()
p2a, = ax2a.plot(df['Year'],df['outregion_ocean_pct'],'-bo', label="Ocean water vapour contribution [%]",markersize=3,markeredgecolor='b',markerfacecolor='b',alpha=0.4)
p2b, = ax2b.plot(df['Year'],df['RR'],'-gs', label="Land water vapour contribution [%]",markersize=3,markeredgecolor='g',markerfacecolor='g')#,alpha=0.5)
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
plt.annotate("Oceanic contribution [%]",xy=(1.04,0.85),xycoords = 'axes fraction',rotation=90,color='blue',fontsize=fsize)
plt.annotate("Terrestrialcontribution [%]",xy=(1.1,0.9),xycoords = 'axes fraction',rotation=90,color='green',fontsize=fsize)
ax2.set_xlim(1979,2014)    
ax2.set_xticks(np.arange(min(df['Year']), max(df['Year'])+1, 1))
ax2.set_xticklabels(ax2.get_xticks(),rotation=90,fontsize=fsize)
ax2.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))    

ax1 = plt.subplot2grid((3, 1), (1, 0), sharex=ax2)
plt.annotate('(b)',xy=(1980,1800),xycoords='data',fontsize=fsize)
#p1, = ax1.plot(df['Year'],df['P_total']/1000,'-k', label="Rainfall total [m]",linewidth=1.5)#,alpha=0.7)
ax1.fill_between(df['Year'],df['P_total']/1000,0,alpha=0.5,edgecolor='None',facecolor='gray')
# NOTE: the annual time series has been moved to show its values in the MIDDLE of each year. This differs
# from the purely annual plots, with show values at the START of each year. Careful when analysing.
#p2, = ax1.plot(df2['Year']+0.5,df2['P_total']/1000,':', label="Rainfall total [m]",markersize=4,linewidth=2,color='gray')
ax1.set_ylim(0,2000)
ax1.grid(axis='x')
ax1.set_xlim(1979,2014)    
ax1.set_xticks(np.arange(min(df['Year']), max(df['Year'])+1, 1))
ax1.set_xticklabels(ax1.get_xticks(),rotation=90,fontsize=fsize)
ax1.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))   
ax1.set_yticklabels(ax1.get_yticks(),fontsize=fsize)
ax1.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))   
plt.annotate('Precipitation total [m]',xy=(-0.06,0.75),xycoords='axes fraction',fontsize=fsize,rotation=90)
ax1a = ax1.twinx()
ax1b = ax1.twinx()
p1a, = ax1a.plot(df['Year'],df['outregion_ocean_mm']/1000,'-bo', label="Ocean water vapour contribution [m]",markersize=3,markeredgecolor='b',markerfacecolor='b',alpha=0.4)
p1b, = ax1b.plot(df['Year'],df['region_land_mm']/1000,'-gs', label="Land water vapour contribution [m]",markersize=3,markeredgecolor='g',markerfacecolor='g')#,alpha=0.5)
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
plt.annotate("Oceanic contribution [m]",xy=(1.05,0.85),xycoords = 'axes fraction',rotation=90,color='blue',fontsize=fsize)
plt.annotate("Terrestrial contribution [m]",xy=(1.105,0.9),xycoords = 'axes fraction',rotation=90,color='green',fontsize=fsize)
    
#ax3  = fig.add_axes([0.5,0.42,0.25,0.15])
ax3 = plt.subplot2grid((3, 1), (2, 0))#, sharex=ax2)
plt.annotate('(c)',xy=(0.6,22),xycoords='data',fontsize=fsize)
#ax3.bar(range(4),RR_seas_mean,width=0.6,align='center',facecolor='g',alpha=0.8)
data = np.column_stack([groups.get_group('DJF')['RR'][1:],\
                        groups.get_group('MAM')['RR'][1:],\
                        groups.get_group('JJA')['RR'][1:],\
                        groups.get_group('SON')['RR'][1:]])
bp = ax3.boxplot(data,whis='range')
plt.setp(bp['boxes'], color='black')
plt.setp(bp['whiskers'], color='black')
plt.setp(bp['medians'], color='black')
for line in bp['medians']:
    # get position data for median line
    x = bp['medians'].index(line)+1#line.get_xydata()[1][0]-line.get_xydata()[0][0]
    y = line.get_xydata()[0,1] + 0.1
    mn = line.get_xydata()[0,1] 
    # overlay median value
    plt.text(x, y, '%.1f' % mn,
         horizontalalignment='center',fontsize=12) # draw above, centered
ax3.set_ylim(5,25)
#plt.tick_params(
#    axis='both',          # changes apply to the x-axis
#    which='both',      # both major and minor ticks are affected
#    bottom=False,      # ticks along the bottom edge are off
#    top=False,         # ticks along the top edge are off
#    left=False,
#    labelleft=False,
#    labelbottom=False) # labels along the bottom edge are off
ax3.set_xticklabels(['Summer','Fall','Winter','Spring'],fontsize=fsize)  
ax3.set_yticklabels(ax3.get_yticks(),fontsize=fsize)
ax3.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))  
plt.annotate('Precipitation recycling [%]',xy=(-0.06,0.9),xycoords='axes fraction',fontsize=fsize,rotation=90)  
#plt.annotate('Summer',xy=(0,11),xycoords='data',fontsize=12,rotation=90)    
#plt.annotate('Fall',xy=(1,4),xycoords='data',fontsize=12,rotation=90)    
#plt.annotate('Winter',xy=(2,8),xycoords='data',fontsize=12,rotation=90)    
#plt.annotate('Spring',xy=(3,9),xycoords='data',fontsize=12,rotation=90)    
#plt.annotate(str(round(RR_seas_mean[0],1))+'%',xy=(-0.2,20),xycoords='data',fontsize=12)  
#plt.annotate(str(round(RR_seas_mean[1],1))+'%',xy=(0.8,17),xycoords='data',fontsize=12)  
#plt.annotate(str(round(RR_seas_mean[2],1))+'%',xy=(1.8,13),xycoords='data',fontsize=12) 
#plt.annotate(str(round(RR_seas_mean[3],1))+'%',xy=(2.8,17.9),xycoords='data',fontsize=12) 

    
plt.tight_layout()
fig.savefig(dir_out+'Fig5.png',bbox_inches='tight') 
plt.close()
    
