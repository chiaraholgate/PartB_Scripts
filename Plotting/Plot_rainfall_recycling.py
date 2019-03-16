# -*- coding: utf-8 -*-
"""
Created on Fri Nov  9 14:35:42 2018

This script creates plots of rainfall recyling (equiv. to land vapour contribution) and ocean vapour contribution 
per time based on .csv's made in ../Analysis/Rainfall_recycling_region.py.

Plots are made for a defined region:
- on an annual basis and seasonal basis
- as percentage vapour contributions and rainfall anomaly
- as mm depth contributions and total rainfall

*** NOTE: annual values are plotted at the START of each year. Seasonal values are plotted at an increment
within each year. The annual line shown on the top seasonal panel has had the annual values moved to the 
MIDDLE of the year.
>> So be careful when comparing seasonal plots to annual ones.

TO DO: CHECK SEASONAL.

@author: z3131380
"""
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

dir_in = '/srv/ccrc/data19/z3131380/PartB/Output/Australia/100parcels/TS10min/exp01/Processed/'
dir_out = '/srv/ccrc/data19/z3131380/PartB/Output/Australia/100parcels/TS10min/exp01/Processed/Rainfall_recycling/'
domain = 'Australia'
region = 'MDB'
timeblock = 'seasonal'
n_i,n_j = 134,205 # QIBT model dimensions
fsize = 14

def make_patch_spines_invisible(ax):
    ax.set_frame_on(True)
    ax.patch.set_visible(False)
    for sp in ax.spines.values():
        sp.set_visible(False)

# Set plot limits        
if region == 'Australia':
    y1 = 75; y2 = 85; y3 = 14; y4 = 22; y5 = 0; y6 = 2000; y7 = 600; y8 = 1800; y9 = 100; y10 = 400
if region == 'MDB':
    y1 = 75; y2 = 95; y3 = 0; y4 = 14; y5 = 0; y6 = 500; y7 = 50; y8 = 300; y9 = 0; y10 = 50      
if region == 'SWWA':
    y1 = 75; y2 = 100; y3 = 0; y4 = 10; y5 = 0; y6 = 100; y7 = 0; y8 = 50; y9 = 0; y10 = 10
    
dir_out = '/srv/ccrc/data19/z3131380/PartB/Output/Australia/100parcels/TS10min/exp01/Processed/Rainfall_recycling/'     


#==============================================================================
# Make plots
#==============================================================================

if timeblock == 'annual':
    df = pd.read_csv(dir_out+'Yearly/'+region+'_'+timeblock+'_rainfall_recycling_1979-2013.csv')
    
    #==============================================================================
    # Plot time series of  land/ocean contributions and rainfall as %
    #==============================================================================
    fig, host = plt.subplots(figsize=(16,11.7))
    fig.subplots_adjust(right=0.75)
    
    par1 = host.twinx()
    par2 = host.twinx()
    par3 = host.twinx()
    
    # Offset the right spine of par2.  The ticks and label have already been
    # placed on the right by twinx above.
    par2.spines["right"].set_position(("axes", 1.1))
    par3.spines["right"].set_position(("axes", 1.2))
    # Having been created by twinx, par2 has its frame off, so the line of its
    # detached spine is invisible.  First, activate the frame but make the patch
    # and spines invisible.
    make_patch_spines_invisible(par2)
    make_patch_spines_invisible(par3)
    # Second, show the right spine.
    par2.spines["right"].set_visible(True)
    par3.spines["right"].set_visible(True)
    
    p1, = host.plot(df['Year'],df['P_anom'],':ko', label="Rainfall anomaly [%]",markersize=2,linewidth=2)
    p2, = par1.plot(df['Year'],df['outregion_ocean_pct'],'-bo', label="Ocean water vapour contribution [%]",markersize=2)
    p3, = par2.plot(df['Year'],df['RR'],'-go', label="Land water vapour contribution [%]",markersize=2)
    p4, = par3.plot(df['Year'],df['outregion_land_Aus_pct'],'-co', label="Out of region land water vapour contribution [%]",markersize=2)
     
    host.set_xticks(np.arange(min(df['Year']), max(df['Year'])+1, 1))
    host.set_xticklabels(host.get_xticks(),rotation=90)
    
    host.set_xlim(1979,2013)
    host.set_ylim(-1, 1)
    par1.set_ylim(y1,y2)
    par2.set_ylim(y3,y4)
    par3.set_ylim(y3,y4)
    
    host.set_xlabel("Year",fontsize=fsize)
    host.set_ylabel("Rainfall anomaly [%]",fontsize=fsize)
    par1.set_ylabel("Ocean water vapour contribution [%]",fontsize=fsize)
    par2.set_ylabel("Land water vapour contribution [%]",fontsize=fsize)
    par3.set_ylabel("Out of region (Aus) land water vapour contribution [%]",fontsize=fsize)
    
    host.yaxis.label.set_color(p1.get_color())
    par1.yaxis.label.set_color(p2.get_color())
    par2.yaxis.label.set_color(p3.get_color())
    par3.yaxis.label.set_color(p4.get_color())
    
    tkw = dict(size=4, width=1.5)
    host.tick_params(axis='y', colors=p1.get_color(), **tkw)
    par1.tick_params(axis='y', colors=p2.get_color(), **tkw)
    par2.tick_params(axis='y', colors=p3.get_color(), **tkw)
    par3.tick_params(axis='y', colors=p4.get_color(), **tkw)
    host.tick_params(axis='x', **tkw)
    
    host.axhline(0,color='k',linestyle=':')
    host.grid(axis='x')
    
    plt.annotate('Rainfall vs land contribution = '+str(round(stats.pearsonr(df['P_anom'],df['RR'])[0],2)),xy=(0.15,0.05),xycoords='axes fraction')
    plt.annotate('Rainfall vs out of region (Aus) land contribution = '+str(round(stats.pearsonr(df['P_anom'],df['outregion_land_Aus_pct'])[0],2)),xy=(0.15,0.03),xycoords='axes fraction')
    plt.annotate('Rainfall vs ocean contribution = '+str(round(stats.pearsonr(df['P_anom'],df['outregion_ocean_pct'])[0],2)),xy=(0.15,0.01),xycoords='axes fraction')
    
    fig.savefig(dir_out+region+'_1979-2013_'+timeblock+'_wvcont_pct.png',bbox_inches='tight') 
    plt.close()
    
    
    #==============================================================================
    # Plot time series of annual land/ocean contributions and rainfall as depth [m]
    # ...with ocean and out of region land contributions seperate
    #==============================================================================
    fig, host = plt.subplots(figsize=(16,11.7))
    fig.subplots_adjust(right=0.75)
    par1 = host.twinx()
    par2 = host.twinx()
    par3 = host.twinx()
    
    # Offset the right spine of par2.  The ticks and label have already been
    # placed on the right by twinx above.
    par2.spines["right"].set_position(("axes", 1.1))
    par3.spines["right"].set_position(("axes", 1.2))
    # Having been created by twinx, par2 has its frame off, so the line of its
    # detached spine is invisible.  First, activate the frame but make the patch
    # and spines invisible.
    make_patch_spines_invisible(par2)
    make_patch_spines_invisible(par3)
    # Second, show the right spine.
    par2.spines["right"].set_visible(True)
    par3.spines["right"].set_visible(True)
    
    p1, = host.plot(df['Year'],df['P_total']/1000,':ko', label="Annual rainfall [m]",markersize=2,linewidth=2)
    p2, = par1.plot(df['Year'],df['outregion_ocean_mm']/1000,'-bo', label="Ocean water vapour contribution [m]",markersize=2)
    p3, = par2.plot(df['Year'],df['region_land_mm']/1000,'-go', label="Land water vapour contribution [m]",markersize=2)
    p4, = par3.plot(df['Year'],df['outregion_land_Aus_mm']/1000,'-co', label="Out of region (Aus) land water vapour contribution [m]",markersize=2)
    
    host.set_xticks(np.arange(min(df['Year']), max(df['Year'])+1, 1))
    host.set_xticklabels(host.get_xticks(),rotation=90)
    
    host.set_xlim(1979,2013)
    host.set_ylim(y5,y6)
    par1.set_ylim(y7,y8)
    par2.set_ylim(y9,y10)
    par3.set_ylim(y9,y10)
    
    host.set_xlabel("Year",fontsize=fsize)
    host.set_ylabel("Annual rainfall [m]",fontsize=fsize)
    par1.set_ylabel("Ocean water vapour contribution [m]",fontsize=fsize)
    par2.set_ylabel("Land water vapour contribution [m]",fontsize=fsize)
    par3.set_ylabel("Out of region (Aus) land water vapour contribution [m]",fontsize=fsize)
    
    host.yaxis.label.set_color(p1.get_color())
    par1.yaxis.label.set_color(p2.get_color())
    par2.yaxis.label.set_color(p3.get_color())
    par3.yaxis.label.set_color(p4.get_color())
    
    tkw = dict(size=4, width=1.5)
    host.tick_params(axis='y', colors=p1.get_color(), **tkw)
    par1.tick_params(axis='y', colors=p2.get_color(), **tkw)
    par2.tick_params(axis='y', colors=p3.get_color(), **tkw)
    par3.tick_params(axis='y', colors=p4.get_color(), **tkw)
    host.tick_params(axis='x', **tkw)
    
    host.axhline(0,color='k',linestyle=':')
    host.grid(axis='x')
    
    plt.annotate('Rainfall vs land contribution = '+str(round(stats.pearsonr(df['P_total'],df['region_land_mm'])[0],2)),xy=(0.15,0.05),xycoords='axes fraction')
    plt.annotate('Rainfall vs out of region (Aus) land contribution = '+str(round(stats.pearsonr(df['P_total'],df['outregion_land_Aus_mm'])[0],2)),xy=(0.15,0.03),xycoords='axes fraction')
    plt.annotate('Rainfall vs ocean contribution = '+str(round(stats.pearsonr(df['P_total'],df['outregion_ocean_mm'])[0],2)),xy=(0.15,0.01),xycoords='axes fraction')
    
    fig.savefig(dir_out+region+'_1979-2013_annual_wvcont_m_2.png',bbox_inches='tight') 
    plt.close()
    
    #==============================================================================
    # Plot time series of annual land/ocean contributions and rainfall as depth [m]
    # ...with ocean and out of region land contributions combined
    #==============================================================================
    fig, host = plt.subplots(figsize=(16,11.7))
    fig.subplots_adjust(right=0.75)
    par1 = host.twinx()
    par2 = host.twinx()
    
    # Offset the right spine of par2.  The ticks and label have already been
    # placed on the right by twinx above.
    par2.spines["right"].set_position(("axes", 1.1))
    #par3.spines["right"].set_position(("axes", 1.2))
    # Having been created by twinx, par2 has its frame off, so the line of its
    # detached spine is invisible.  First, activate the frame but make the patch
    # and spines invisible.
    make_patch_spines_invisible(par2)
    make_patch_spines_invisible(par3)
    # Second, show the right spine.
    par2.spines["right"].set_visible(True)
    
    p1, = host.plot(df['Year'],df['P_total']/1000,':ko', label="Annual rainfall [m]",markersize=3,linewidth=2)
    p2, = par1.plot(df['Year'],df['outregion_ocean_mm']/1000+df['outregion_land_outside_Aus_mm']/1000,'-bo', label="Out of region water vapour contribution [m]",markersize=3)
    p3, = par2.plot(df['Year'],df['region_land_mm']/1000,'-go', label="Land water vapour contribution [m]",markersize=3)
    
    host.set_xticks(np.arange(min(df['Year']), max(df['Year'])+1, 1))
    host.set_xticklabels(host.get_xticks(),rotation=90)
    
    host.set_xlim(1979,2013)
    host.set_ylim(y5,y6)
    par1.set_ylim(y7,y8)
    par2.set_ylim(y9,y10)
    #par3.set_ylim(100,400)
    
    host.set_xlabel("Year",fontsize=fsize)
    host.set_ylabel("Annual rainfall [m]",fontsize=fsize)
    par1.set_ylabel("Out of region water vapour contribution [m]",fontsize=fsize)
    par2.set_ylabel("Land water vapour contribution [m]",fontsize=fsize)
    
    host.yaxis.label.set_color(p1.get_color())
    par1.yaxis.label.set_color(p2.get_color())
    par2.yaxis.label.set_color(p3.get_color())
    
    tkw = dict(size=4, width=1.5)
    host.tick_params(axis='y', colors=p1.get_color(), **tkw)
    par1.tick_params(axis='y', colors=p2.get_color(), **tkw)
    par2.tick_params(axis='y', colors=p3.get_color(), **tkw)
    host.tick_params(axis='x', **tkw)
    
    host.axhline(0,color='k',linestyle=':')
    host.grid(axis='x')
    
#    plt.annotate('Rainfall vs land contribution = '+str(round(stats.pearsonr(df['P_total'],df['region_land_mm'])[0],2)),xy=(0.15,0.05),xycoords='axes fraction')
#    plt.annotate('Rainfall vs out of region Australian land contribution = '+str(round(stats.pearsonr(df['P_total'],df['outregion_land_Aus_mm'])[0],2)),xy=(0.15,0.03),xycoords='axes fraction')
#    plt.annotate('Rainfall vs ocean contribution = '+str(round(stats.pearsonr(df['P_total'],df['outregion_ocean_mm'])[0],2)),xy=(0.15,0.01),xycoords='axes fraction')
    
    fig.savefig(dir_out+region+'_1979-2013_annual_wvcont_m.png',bbox_inches='tight') 
    plt.close()
        
    
        
    #==============================================================================
    # Plot time series of annual land/ocean contributions and rainfall [%]
    # ...with dL/dt
    #==============================================================================
    fig, ax1 = plt.subplots(figsize=(16,11.7))        
    ax1 = plt.subplot2grid((4, 1), (0, 0), colspan=1)
    #p1, = ax1.plot(df['Year'],df['dL/dt_2'],':ko', label="dL/dt_5",markersize=3)
    #p2, = ax1.plot(df['Year'],df['dL/dt_5'],':ko', label="dL/dt_5",markersize=3)
    ax1.bar(df['Year'],df['dL/dt_2'],color='gray',align='center',label='2-year',edgecolor='gray')
    ax1.bar(df['Year'],df['dL/dt_5'],color='k',align='center',label='5-year',fill=False)
    ax1.tick_params(labelbottom=False) 
    ax1.set_xticks(np.arange(min(df['Year']), max(df['Year'])+1, 1))
    ax1.set_xlim(1979,2013)
    ax1.set_ylabel("dL/dt",fontsize=fsize)
    plt.legend(loc='center left',bbox_to_anchor=(1.01,0.5),frameon=False)
    
    ax2 = plt.subplot2grid((4, 1), (1, 0), rowspan=3)
    plt.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0.1, hspace=0.2)
    p3, = ax2.plot(df['Year'],df['P_anom'],'--o', label="Precipitation anomaly [%]",markersize=3,color='gray')
    ax2.set_xticks(np.arange(min(df['Year']), max(df['Year'])+1, 1))
    ax2.set_xticklabels(ax2.get_xticks(),rotation=90)
    ax2.grid(axis='x')
    ax3 = ax2.twinx()
    ax4 = ax2.twinx()
    #ax5 = ax2.twinx()
    p4, = ax3.plot(df['Year'],df['outregion_ocean_pct'],'-bo', label="Ocean water vapour contribution [%]",markersize=3)
    p5, = ax4.plot(df['Year'],df['RR'],'-go', label="Land water vapour contribution [%]",markersize=3)
    #p6, = ax5.plot(df['Year'],df['outregion_land_pct'],'-co', label="Out of region land water vapour contribution [%]",markersize=3)
    ax2.set_xlim(1979,2013)
    #host.set_xlim(1979,2013)
    host.set_ylim(-1, 1)
    par1.set_ylim(y1,y2)
    par2.set_ylim(y3,y4)
    par3.set_ylim(y3,y4)

    ax4.spines["right"].set_position(("axes", 1.05))
    #ax5.spines["right"].set_position(("axes", 1.1))
    make_patch_spines_invisible(ax4)
    #make_patch_spines_invisible(ax5)
    ax4.spines["right"].set_visible(True)
    #ax5.spines["right"].set_visible(True)
    ax2.set_xlabel("Year",fontsize=fsize)
    ax2.set_ylabel("Rainfall anomaly [%]",fontsize=fsize)
    ax3.set_ylabel("Ocean water vapour contribution [%]",fontsize=fsize)
    ax4.set_ylabel("Land water vapour contribution [%]",fontsize=fsize)
    #ax5.set_ylabel("Out of region land water vapour contribution [%]",fontsize=fsize)
    #ax1.yaxis.label.set_color(p1.get_color())
    #ax2.yaxis.label.set_color(p3.get_color())
    ax3.yaxis.label.set_color(p4.get_color())
    ax4.yaxis.label.set_color(p5.get_color())
    #ax5.yaxis.label.set_color(p6.get_color())
    tkw = dict(size=4, width=1.5)
    #ax1.tick_params(axis='y', colors=p1.get_color(), **tkw)
    #ax2.tick_params(axis='y', colors=p3.get_color(), **tkw)
    ax3.tick_params(axis='y', colors=p4.get_color(), **tkw)
    ax4.tick_params(axis='y', colors=p5.get_color(), **tkw)
    #ax5.tick_params(axis='y', colors=p6.get_color(), **tkw)
    ax1.tick_params(axis='x', **tkw)
    ax1.axhline(0,color='k',linestyle=':')
    ax2.axhline(0,color='gray',linestyle=':')
    ax1.grid(axis='x')
    plt.annotate('Rainfall vs land contribution = '+str(round(stats.pearsonr(df['P_anom'],df['RR'])[0],2)),xy=(0.15,0.05),xycoords='axes fraction')
    plt.annotate('Rainfall vs ocean contribution = '+str(round(stats.pearsonr(df['P_anom'],df['outregion_ocean_pct'])[0],2)),xy=(0.15,0.03),xycoords='axes fraction')
    plt.annotate('Rainfall vs out of region land contribution = '+str(round(stats.pearsonr(df['P_anom'],df['outregion_land_outside_Aus_pct'])[0],2)),xy=(0.15,0.01),xycoords='axes fraction')
    
    fig.savefig(dir_out+region+'_1979-2013_annual_wvcont_pct_dLdt.png',bbox_inches='tight') 
    plt.close()
    
 
 
 
elif timeblock == 'seasonal':
    df = pd.read_csv(dir_out+'Seasonal/'+region+'_'+timeblock+'_rainfall_recycling_1979-2013.csv')
    
    # Open annual df to plot comparison
    df2 = pd.read_csv(dir_out+'Yearly/'+region+'_annual_rainfall_recycling_1979-2013.csv')
    
    # Group data by season
    groups = df.groupby('Season')
    
    #==============================================================================
    # Plot time series of  land/ocean contributions and rainfall as %
    #==============================================================================
    fig, ax1 = plt.subplots(figsize=(16,11.7))        
    ax1 = plt.subplot2grid((6, 1), (0, 0), rowspan=2)
    p1, = ax1.plot(df['Year'],df['P_anom'],'-ko', label="Rainfall anomaly [%]",markersize=4,linewidth=1)
    # NOTE: the annual time series has been moved to show its values in the MIDDLE of each year. This differs
    # from the purely annual plots, with show values at the START of each year. Careful when analysing.
    p2, = ax1.plot(df2['Year']+0.5,df2['P_anom'],':o', label="Rainfall anomaly [%]",markersize=4,linewidth=2,color='gray')
    ax1.grid(axis='x')
    ax1.axhline(0,color='gray',linestyle=':')
    ax1.tick_params(labelbottom=False) 
    ax1.tick_params(labelbottom=False) 
    plt.legend(['Seasonal P anom','Annual P anom'],loc=0,frameon=False) #loc='center left', bbox_to_anchor=(1.01,0.5)
    plt.annotate('All seasons',xy=(1980,0.8),xycoords='data',fontsize=fsize)
    ax1a = ax1.twinx()
    ax1b = ax1.twinx()
    p1a, = ax1a.plot(df['Year'],df['outregion_ocean_pct'],'-bo', label="Ocean water vapour contribution [%]",markersize=4,alpha=0.5)
    p1b, = ax1b.plot(df['Year'],df['RR'],'-go', label="Land water vapour contribution [%]",markersize=4,alpha=0.5)
    ax1a.spines["right"].set_position(("axes", 1.01))
    ax1b.spines["right"].set_position(("axes", 1.07))
    make_patch_spines_invisible(ax1a)
    make_patch_spines_invisible(ax1b)
    ax1a.spines["right"].set_visible(True)
    ax1b.spines["right"].set_visible(True)
    tkw = dict(size=4, width=1.5)
    ax1a.tick_params(axis='y', colors=p1a.get_color(), **tkw)
    ax1b.tick_params(axis='y', colors=p1b.get_color(), **tkw)
    ax1a.set_ylim(70,95)
    ax1b.set_ylim(5,25)

    
    ax2 = plt.subplot2grid((6, 1), (2, 0), sharex=ax1)
    p3, = ax2.plot(groups.get_group('DJF')['Year'],groups.get_group('DJF')['P_anom'],':o', label="Rainfall anomaly [%]",markersize=4,linewidth=2,color='gray')
    ax2.grid(axis='x')
    ax2.axhline(0,color='gray',linestyle=':')
    ax2.tick_params(labelbottom=False) 
    plt.annotate('Summer',xy=(1980,0.7),xycoords='data',fontsize=fsize)
    plt.annotate('Rainfall anomaly [%]',xy=(-0.05,0.45),xycoords='axes fraction',fontsize=fsize,rotation=90)
    ax2a = ax2.twinx()
    ax2b = ax2.twinx()
    p3a, = ax2a.plot(groups.get_group('DJF')['Year'],groups.get_group('DJF')['outregion_ocean_pct'],'-bo', label="Ocean water vapour contribution [%]",markersize=4)
    p3b, = ax2b.plot(groups.get_group('DJF')['Year'],groups.get_group('DJF')['RR'],'-go', label="Land water vapour contribution [%]",markersize=4)
    ax2a.spines["right"].set_position(("axes", 1.01))
    ax2b.spines["right"].set_position(("axes", 1.07))
    make_patch_spines_invisible(ax2a)
    make_patch_spines_invisible(ax2b)
    ax2a.spines["right"].set_visible(True)
    ax2b.spines["right"].set_visible(True)
    tkw = dict(size=4, width=1.5)
    ax2a.tick_params(axis='y', colors=p3a.get_color(), **tkw)
    ax2b.tick_params(axis='y', colors=p3b.get_color(), **tkw)
    ax2a.set_ylim(70,95)
    ax2b.set_ylim(5,25)
   
    ax3 = plt.subplot2grid((6, 1), (3, 0), sharex=ax1)
    p4, = ax3.plot(groups.get_group('MAM')['Year'],groups.get_group('MAM')['P_anom'],':o', label="Rainfall anomaly [%]",markersize=4,linewidth=2,color='gray')
    ax3.grid(axis='x')
    ax3.axhline(0,color='gray',linestyle=':')
    ax3.tick_params(labelbottom=False) 
    plt.annotate('Autumn',xy=(1980,0.7),xycoords='data',fontsize=fsize)
    ax3a = ax3.twinx()
    ax3b = ax3.twinx()
    p4a, = ax3a.plot(groups.get_group('MAM')['Year'],groups.get_group('MAM')['outregion_ocean_pct'],'-bo', label="Ocean water vapour contribution [%]",markersize=4)
    p4b, = ax3b.plot(groups.get_group('MAM')['Year'],groups.get_group('MAM')['RR'],'-go', label="Land water vapour contribution [%]",markersize=4)
    ax3a.spines["right"].set_position(("axes", 1.01))
    ax3b.spines["right"].set_position(("axes", 1.07))
    make_patch_spines_invisible(ax3a)
    make_patch_spines_invisible(ax3b)
    ax3a.spines["right"].set_visible(True)
    ax3b.spines["right"].set_visible(True)
    tkw = dict(size=4, width=1.5)
    ax3a.tick_params(axis='y', colors=p4a.get_color(), **tkw)
    ax3b.tick_params(axis='y', colors=p4b.get_color(), **tkw)
    ax3a.set_ylim(70,95)
    ax3b.set_ylim(5,25)
   
    ax4 = plt.subplot2grid((6, 1), (4, 0), sharex=ax1)
    p5, = ax4.plot(groups.get_group('JJA')['Year'],groups.get_group('JJA')['P_anom'],':o', label="Rainfall anomaly [%]",markersize=4,linewidth=2,color='gray')
    ax4.grid(axis='x')
    ax4.axhline(0,color='gray',linestyle=':')
    ax4.tick_params(labelbottom=False) 
    plt.annotate('Winter',xy=(1980,0.7),xycoords='data',fontsize=fsize)
    ax4a = ax4.twinx()
    ax4b = ax4.twinx()
    p5a, = ax4a.plot(groups.get_group('JJA')['Year'],groups.get_group('JJA')['outregion_ocean_pct'],'-bo', label="Ocean water vapour contribution [%]",markersize=4)
    p5b, = ax4b.plot(groups.get_group('JJA')['Year'],groups.get_group('JJA')['RR'],'-go', label="Land water vapour contribution [%]",markersize=4)
    ax4a.spines["right"].set_position(("axes", 1.01))
    ax4b.spines["right"].set_position(("axes", 1.07))
    make_patch_spines_invisible(ax4a)
    make_patch_spines_invisible(ax4b)
    ax4a.spines["right"].set_visible(True)
    ax4b.spines["right"].set_visible(True)
    plt.annotate("Ocean water vapour contribution [%]",xy=(1.04,2.5),xycoords = 'axes fraction',rotation=90,color='blue',fontsize=fsize)
    plt.annotate("Land water vapour contribution [%]",xy=(1.1,2.5),xycoords = 'axes fraction',rotation=90,color='green',fontsize=fsize)
    ax4a.yaxis.label.set_color(p5a.get_color())
    ax4b.yaxis.label.set_color(p5b.get_color())
    tkw = dict(size=4, width=1.5)
    ax4a.tick_params(axis='y', colors=p5a.get_color(), **tkw)
    ax4b.tick_params(axis='y', colors=p5b.get_color(), **tkw)
    ax4a.set_ylim(70,95)
    ax4b.set_ylim(5,25)
    
    ax5 = plt.subplot2grid((6, 1), (5, 0), sharex=ax1)
    p6, = ax5.plot(groups.get_group('SON')['Year'],groups.get_group('SON')['P_anom'],':o', label="Rainfall anomaly [%]",markersize=4,linewidth=2,color='gray')
    ax5.grid(axis='x')
    ax5.axhline(0,color='gray',linestyle=':')
    ax5.set_xlim(1979,2014)    
    ax5.set_xticks(np.arange(min(df['Year']), max(df['Year'])+1, 1))
    ax5.set_xticklabels(ax5.get_xticks(),rotation=90)
    ax5.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
    plt.annotate('Spring',xy=(1980,0.7),xycoords='data',fontsize=fsize)
    ax5a = ax5.twinx()
    ax5b = ax5.twinx()
    p6a, = ax5a.plot(groups.get_group('SON')['Year'],groups.get_group('SON')['outregion_ocean_pct'],'-bo', label="Ocean water vapour contribution [%]",markersize=4)
    p6b, = ax5b.plot(groups.get_group('SON')['Year'],groups.get_group('SON')['RR'],'-go', label="Land water vapour contribution [%]",markersize=4)
    ax5a.spines["right"].set_position(("axes", 1.01))
    ax5b.spines["right"].set_position(("axes", 1.07))
    make_patch_spines_invisible(ax5a)
    make_patch_spines_invisible(ax5b)
    ax5a.spines["right"].set_visible(True)
    ax5b.spines["right"].set_visible(True)
    tkw = dict(size=4, width=1.5)
    ax5a.tick_params(axis='y', colors=p6a.get_color(), **tkw)
    ax5b.tick_params(axis='y', colors=p6b.get_color(), **tkw)
    ax5a.set_ylim(70,95)
    ax5b.set_ylim(5,25)
    
    ax1.set_ylim(-1,1); ax2.set_ylim(-1,1); ax3.set_ylim(-1,1); ax4.set_ylim(-1,1); ax5.set_ylim(-1,1)
    
    fig.savefig(dir_out+'Seasonal/'+region+'_1979-2013_'+timeblock+'_wvcont_pct.png',bbox_inches='tight') 
    plt.close()
    
    
    #==============================================================================
    # Plot time series of  land/ocean contributions and rainfall as depth
    #==============================================================================
    fig, ax1 = plt.subplots(figsize=(16,11.7))        
    ax1 = plt.subplot2grid((6, 1), (0, 0), rowspan=2)
    p1, = ax1.plot(df['Year'],df['P_total']/1000,'-ko', label="Rainfall total [m]",markersize=4,linewidth=1)
    # NOTE: the annual time series has been moved to show its values in the MIDDLE of each year. This differs
    # from the purely annual plots, with show values at the START of each year. Careful when analysing.
    p2, = ax1.plot(df2['Year']+0.5,df2['P_total']/1000,':o', label="Rainfall total [m]",markersize=4,linewidth=2,color='gray')
    ax1.grid(axis='x')
    ax1.axhline(0,color='gray',linestyle=':')
    ax1.tick_params(labelbottom=False) 
    ax1.tick_params(labelbottom=False) 
    plt.legend(['Seasonal','Annual'],loc=4,frameon=False) #loc='center left', bbox_to_anchor=(1.01,0.5)
    plt.annotate('All seasons',xy=(1981,0.8),xycoords='data',fontsize=fsize)
    ax1a = ax1.twinx()
    ax1b = ax1.twinx()
    p1a, = ax1a.plot(df['Year'],df['outregion_ocean_mm']/1000,'-bo', label="Ocean water vapour contribution [m]",markersize=4,alpha=0.5)
    p1b, = ax1b.plot(df['Year'],df['outregion_land_Aus_mm']/1000,'-go', label="Land water vapour contribution [m]",markersize=4,alpha=0.5)
    ax1a.spines["right"].set_position(("axes", 1.01))
    ax1b.spines["right"].set_position(("axes", 1.07))
    make_patch_spines_invisible(ax1a)
    make_patch_spines_invisible(ax1b)
    ax1a.spines["right"].set_visible(True)
    ax1b.spines["right"].set_visible(True)
    tkw = dict(size=4, width=1.5)
    ax1a.tick_params(axis='y', colors=p1a.get_color(), **tkw)
    ax1b.tick_params(axis='y', colors=p1b.get_color(), **tkw)
#    ax1a.set_ylim(700,1700)
#    ax1b.set_ylim(100,4000)

    
    ax2 = plt.subplot2grid((6, 1), (2, 0), sharex=ax1)
    p3, = ax2.plot(groups.get_group('DJF')['Year'],groups.get_group('DJF')['P_total']/1000,':o', label="Rainfall total [m]",markersize=4,linewidth=2,color='gray')
    ax2.grid(axis='x')
    ax2.axhline(0,color='gray',linestyle=':')
    ax2.tick_params(labelbottom=False) 
    plt.annotate('Summer',xy=(1981,0.7),xycoords='data',fontsize=fsize)
    plt.annotate('Rainfall total [m]',xy=(-0.05,0.45),xycoords='axes fraction',fontsize=fsize,rotation=90)
    ax2a = ax2.twinx()
    ax2b = ax2.twinx()
    p3a, = ax2a.plot(groups.get_group('DJF')['Year'],groups.get_group('DJF')['outregion_ocean_mm']/1000,'-bo', label="Ocean water vapour contribution [m]",markersize=4)
    p3b, = ax2b.plot(groups.get_group('DJF')['Year'],groups.get_group('DJF')['outregion_land_Aus_mm']/1000,'-go', label="Land water vapour contribution [m]",markersize=4)
    ax2a.spines["right"].set_position(("axes", 1.01))
    ax2b.spines["right"].set_position(("axes", 1.07))
    make_patch_spines_invisible(ax2a)
    make_patch_spines_invisible(ax2b)
    ax2a.spines["right"].set_visible(True)
    ax2b.spines["right"].set_visible(True)
    tkw = dict(size=4, width=1.5)
    ax2a.tick_params(axis='y', colors=p3a.get_color(), **tkw)
    ax2b.tick_params(axis='y', colors=p3b.get_color(), **tkw)
#    ax2a.set_ylim(70,95)
#    ax2b.set_ylim(5,25)
   
    ax3 = plt.subplot2grid((6, 1), (3, 0), sharex=ax1)
    p4, = ax3.plot(groups.get_group('MAM')['Year'],groups.get_group('MAM')['P_total']/1000,':o', label="Rainfall total [m]",markersize=4,linewidth=2,color='gray')
    ax3.grid(axis='x')
    ax3.axhline(0,color='gray',linestyle=':')
    ax3.tick_params(labelbottom=False) 
    plt.annotate('Autumn',xy=(1981,0.7),xycoords='data',fontsize=fsize)
    ax3a = ax3.twinx()
    ax3b = ax3.twinx()
    p4a, = ax3a.plot(groups.get_group('MAM')['Year'],groups.get_group('MAM')['outregion_ocean_mm']/1000,'-bo', label="Ocean water vapour contribution [m]",markersize=4)
    p4b, = ax3b.plot(groups.get_group('MAM')['Year'],groups.get_group('MAM')['outregion_land_Aus_mm']/1000,'-go', label="Land water vapour contribution [m]",markersize=4)
    ax3a.spines["right"].set_position(("axes", 1.01))
    ax3b.spines["right"].set_position(("axes", 1.07))
    make_patch_spines_invisible(ax3a)
    make_patch_spines_invisible(ax3b)
    ax3a.spines["right"].set_visible(True)
    ax3b.spines["right"].set_visible(True)
    tkw = dict(size=4, width=1.5)
    ax3a.tick_params(axis='y', colors=p4a.get_color(), **tkw)
    ax3b.tick_params(axis='y', colors=p4b.get_color(), **tkw)
#    ax3a.set_ylim(70,95)
#    ax3b.set_ylim(5,25)
   
    ax4 = plt.subplot2grid((6, 1), (4, 0), sharex=ax1)
    p5, = ax4.plot(groups.get_group('JJA')['Year'],groups.get_group('JJA')['P_total']/1000,':o', label="Rainfall total [m]",markersize=4,linewidth=2,color='gray')
    ax4.grid(axis='x')
    ax4.axhline(0,color='gray',linestyle=':')
    ax4.tick_params(labelbottom=False) 
    plt.annotate('Winter',xy=(1981,0.7),xycoords='data',fontsize=fsize)
    ax4a = ax4.twinx()
    ax4b = ax4.twinx()
    p5a, = ax4a.plot(groups.get_group('JJA')['Year'],groups.get_group('JJA')['outregion_ocean_mm']/1000,'-bo', label="Ocean water vapour contribution [m]",markersize=4)
    p5b, = ax4b.plot(groups.get_group('JJA')['Year'],groups.get_group('JJA')['outregion_land_Aus_mm']/1000,'-go', label="Land water vapour contribution [m]",markersize=4)
    ax4a.spines["right"].set_position(("axes", 1.01))
    ax4b.spines["right"].set_position(("axes", 1.07))
    make_patch_spines_invisible(ax4a)
    make_patch_spines_invisible(ax4b)
    ax4a.spines["right"].set_visible(True)
    ax4b.spines["right"].set_visible(True)
    plt.annotate("Ocean water vapour contribution [m]",xy=(1.04,2.5),xycoords = 'axes fraction',rotation=90,color='blue',fontsize=fsize)
    plt.annotate("Land water vapour contribution [m]",xy=(1.1,2.5),xycoords = 'axes fraction',rotation=90,color='green',fontsize=fsize)
    ax4a.yaxis.label.set_color(p5a.get_color())
    ax4b.yaxis.label.set_color(p5b.get_color())
    tkw = dict(size=4, width=1.5)
    ax4a.tick_params(axis='y', colors=p5a.get_color(), **tkw)
    ax4b.tick_params(axis='y', colors=p5b.get_color(), **tkw)
#    ax4a.set_ylim(70,95)
#    ax4b.set_ylim(5,25)
    
    ax5 = plt.subplot2grid((6, 1), (5, 0), sharex=ax1)
    p6, = ax5.plot(groups.get_group('SON')['Year'],groups.get_group('SON')['P_total']/1000,':o', label="Rainfall total [m]",markersize=4,linewidth=2,color='gray')
    ax5.grid(axis='x')
    ax5.axhline(0,color='gray',linestyle=':')
    ax5.set_xlim(1979,2014)    
    ax5.set_xticks(np.arange(min(df['Year']), max(df['Year'])+1, 1))
    ax5.set_xticklabels(ax5.get_xticks(),rotation=90)
    ax5.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
    plt.annotate('Spring',xy=(1981,0.7),xycoords='data',fontsize=fsize)
    ax5a = ax5.twinx()
    ax5b = ax5.twinx()
    p6a, = ax5a.plot(groups.get_group('SON')['Year'],groups.get_group('SON')['outregion_ocean_mm']/1000,'-bo', label="Ocean water vapour contribution [m]",markersize=4)
    p6b, = ax5b.plot(groups.get_group('SON')['Year'],groups.get_group('SON')['outregion_land_Aus_mm']/1000,'-go', label="Land water vapour contribution [m]",markersize=4)
    ax5a.spines["right"].set_position(("axes", 1.01))
    ax5b.spines["right"].set_position(("axes", 1.07))
    make_patch_spines_invisible(ax5a)
    make_patch_spines_invisible(ax5b)
    ax5a.spines["right"].set_visible(True)
    ax5b.spines["right"].set_visible(True)
    tkw = dict(size=4, width=1.5)
    ax5a.tick_params(axis='y', colors=p6a.get_color(), **tkw)
    ax5b.tick_params(axis='y', colors=p6b.get_color(), **tkw)
#    ax5a.set_ylim(70,95)
#    ax5b.set_ylim(5,25)
    
    #ax1.set_ylim(-1,1); ax2.set_ylim(-1,1); ax3.set_ylim(-1,1); ax4.set_ylim(-1,1); ax5.set_ylim(-1,1)
    
    fig.savefig(dir_out+'Seasonal/'+region+'_1979-2013_'+timeblock+'_wvcont_m.png',bbox_inches='tight') 
    plt.close()