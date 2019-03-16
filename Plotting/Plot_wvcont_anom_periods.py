# -*- coding: utf-8 -*-
"""
Created on Tue Oct  9 14:38:37 2018

This script creates maps of the evaporative source regions for rain falling in anomalous periods.
It also creates difference maps: the difference between the annual anomaly map and the annual climatology.
Differences should be in RELATIVE terms...otherwise you'll just see the obvious reduction/increase in moisture
depth during dry/wet periods...

## TO DO:
- plot seasonal deviations
-plot deviations, one subplot for each year
- SCRIPT CREATING SEAS ANOM NETCDF STILL NEEDS TO BE CHECKED!

@author: z3131380
"""
import numpy as np
from netCDF4 import Dataset 
import matplotlib.pyplot as plt
import os
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import BoundaryNorm
import math


dir_in = '/srv/ccrc/data19/z3131380/PartB/Output/Australia/100parcels/TS10min/exp01/Processed/Anom_periods/'
domain = 'Australia'
region = 'Australia'
anom_period_name = 'MilDrought'

if region =='MDB':
    lon_0=145
    lat_0=-30
    llcrnrlat=-40
    urcrnrlat=-20
    llcrnrlon=135
    urcrnrlon=155

fsize=12
fsize_coords = 12

##==============================================================================
## Annual deviation from climatology
##==============================================================================
#file = dir_in+region+'_annual_anom_'+anom_period_name+'.nc'
#if os.path.isfile(file) == True:
#    fh = Dataset(file, mode='r') 
##    wvcont_mm_anom = fh.variables['wvcont_mm_anom'][:]    
##    wvcont_pct_anom = fh.variables['wvcont_pct_anom'][:] 
##    pre_anom = fh.variables['pre_anom'][:]  
#    wvcont_mm_departure = fh.variables['wvcont_mm_departure'][:] 
#    wvcont_mm_departure = np.ma.array(wvcont_mm_departure,mask=(np.isnan(wvcont_mm_departure)))
#    wvcont_pct_departure = fh.variables['wvcont_pct_departure'][:] 
#    wvcont_pct_departure = np.ma.array(wvcont_pct_departure,mask=(np.isnan(wvcont_pct_departure)))
#    pre_mm_departure = fh.variables['pre_mm_departure'][:] 
#    pre_mm = fh.variables['pre_mm'][:]
#    pre_mm = np.ma.array(pre_mm,mask=(np.isnan(pre_mm)))
#    yrlist = list(fh.variables['time'][:])
#    latitcrs = fh.variables['latitcrs'][:]  
#    longicrs = fh.variables['longicrs'][:]  
#    fh.close() 
#
## Plot wvcont departure from climatology as a percentage
#levs = [-0.05,-0.04,-0.03,-0.02,-0.01,0,0.01]
#cMap = plt.get_cmap('RdYlBu_r',len(levs))
#cmaplist = [cMap(i) for i in range(cMap.N)]
#cMap = cMap.from_list('Custom cmap', cmaplist, cMap.N)
#norm = BoundaryNorm(levs, ncolors=cMap.N, clip=True)
#
## Subplots
#plt.clf()    
#fig = plt.figure(figsize=(16,11.7)) # A4, portrait=8.27,11.7, else 16,11.7     
#for yr in yrlist:
#    ax = fig.add_subplot(3,3,yrlist.index(yr)+1)
#    ax.set_title(str(int(yr)),fontsize=fsize)
##    ax.spines['top'].set_visible(False)
##    ax.spines['right'].set_visible(False)
##    ax.spines['bottom'].set_visible(False)
##    ax.spines['left'].set_visible(False)        
#    m = Basemap(projection='stere',lon_0=135,lat_0=-25.,\
#                    llcrnrlat=-45,urcrnrlat=0,\
#                    llcrnrlon=85,urcrnrlon=175,\
#                    resolution='i',area_thresh=10000)
#    x, y = m(longicrs,latitcrs) #m(XLONG[5:-5,5:-5],XLAT[5:-5,5:-5]) 
#    m.drawcoastlines()
#    h=m.pcolormesh(x,y,wvcont_pct_departure[yrlist.index(yr),:,:],norm=norm,cmap=cMap)
#    cax = plt.axes([0.125,0.125,0.78,0.02])
#    cbar = plt.colorbar(h,cax=cax,ticks=levs,orientation='horizontal',boundaries=levs,extend='max',spacing='proportional')
#    cbar.set_label(label='Annual departure from climatological water vapour contribution to rainfall in '+region+' [%]',size=fsize)
#    cbar.ax.tick_params(labelsize=16)
#fig.savefig(dir_in+region+'_'+anom_period_name+'_wvcont_departure_pct.png',bbox_inches='tight') 
#plt.close()    
#
## Plot wvcont departure from climatology as a depth
#levs = [-20,-15,-10,-5,0,5,10,15,20]
#cMap = plt.get_cmap('RdYlBu_r',len(levs))
#cmaplist = [cMap(i) for i in range(cMap.N)]
#cmaplist[3] = (0.9,0.9,0.9,1)
#cmaplist[4] = (0.9,0.9,0.9,1)
#cMap = cMap.from_list('Custom cmap', cmaplist, cMap.N)
#norm = BoundaryNorm(levs, ncolors=cMap.N, clip=True)
#
## Subplots
#plt.clf()    
#fig = plt.figure(figsize=(16,11.7)) # A4, portrait=8.27,11.7, else 16,11.7     
#for yr in yrlist:
#    ax = fig.add_subplot(3,3,yrlist.index(yr)+1)
#    ax.set_title(str(int(yr)),fontsize=fsize)
##    ax.spines['top'].set_visible(False)
##    ax.spines['right'].set_visible(False)
##    ax.spines['bottom'].set_visible(False)
##    ax.spines['left'].set_visible(False)        
#    m = Basemap(projection='stere',lon_0=135,lat_0=-25.,\
#                    llcrnrlat=-45,urcrnrlat=0,\
#                    llcrnrlon=85,urcrnrlon=175,\
#                    resolution='i',area_thresh=10000)
#    x, y = m(longicrs,latitcrs) #m(XLONG[5:-5,5:-5],XLAT[5:-5,5:-5]) 
#    m.drawcoastlines()
#    h=m.pcolormesh(x,y,wvcont_mm_departure[yrlist.index(yr),:,:],norm=norm,cmap=cMap)#,levels=levs)
#    cax = plt.axes([0.125,0.125,0.78,0.02])
#    cbar = plt.colorbar(h,cax=cax,ticks=levs,orientation='horizontal',boundaries=levs,extend='both',spacing='proportional')
#    cbar.set_label(label='Annual departure from climatological water vapour contribution to rainfall in '+region+' [mm]',size=fsize)
#    cbar.ax.tick_params(labelsize=16)    
#fig.savefig(dir_in+region+'_'+anom_period_name+'_wvcont_departure_mm_greyed.png',bbox_inches='tight') 
#plt.close()
#
## Plot rainfall as a depth
#levs = [0,100,250,350,500,750,1000]
#cMap = plt.get_cmap('RdYlBu',len(levs))
##cmaplist = [cMap(i) for i in range(cMap.N)]
##cMap = cMap.from_list('Custom cmap', cmaplist, cMap.N)
#norm = BoundaryNorm(levs, ncolors=cMap.N, clip=True)
#
## Subplots
#plt.clf()    
#fig = plt.figure(figsize=(16,11.7)) # A4, portrait=8.27,11.7, else 16,11.7     
#for yr in yrlist:
#    ax = fig.add_subplot(3,3,yrlist.index(yr)+1)
#    ax.set_title(str(int(yr)),fontsize=fsize)
##    ax.spines['top'].set_visible(False)
##    ax.spines['right'].set_visible(False)
##    ax.spines['bottom'].set_visible(False)
##    ax.spines['left'].set_visible(False)        
#    m = Basemap(projection='stere',lon_0=135,lat_0=-25.,\
#                    llcrnrlat=-45,urcrnrlat=0,\
#                    llcrnrlon=85,urcrnrlon=175,\
#                    resolution='i',area_thresh=10000)
#    x, y = m(longicrs,latitcrs) 
#    m.drawcoastlines()
#    h=m.contourf(x,y,pre_mm[yrlist.index(yr),:,:],norm=norm,cmap=cMap,levels=levs)
#    cax = plt.axes([0.125,0.125,0.78,0.02])
#    cbar = plt.colorbar(h,cax=cax,orientation='horizontal',extend='max',spacing='proportional',ticks=levs,boundaries=levs)
#    cbar.set_label(label='Annual rainfall [mm]',size=fsize)
#    cbar.ax.tick_params(labelsize=16)    
#fig.savefig(dir_in+region+'_'+anom_period_name+'_rainfall_mm.png',bbox_inches='tight') 
#plt.close()



#==============================================================================
# Seasonal deviation from climatology
#==============================================================================
file = dir_in+region+'_seasonal_anom_'+anom_period_name+'.nc'
if os.path.isfile(file) == True:
    fh = Dataset(file, mode='r') 
    wvcont_mm_seasonal_departure = fh.variables['wvcont_mm_seasonal_departure'][:] 
    wvcont_mm_seasonal_departure = np.ma.array(wvcont_mm_seasonal_departure,mask=(np.isnan(wvcont_mm_seasonal_departure)))
    wvcont_pct_seasonal_departure = fh.variables['wvcont_pct_seasonal_departure'][:] 
    wvcont_pct_seasonal_departure = np.ma.array(wvcont_pct_seasonal_departure,mask=(np.isnan(wvcont_pct_seasonal_departure)))
    pre_seasonal_departure = fh.variables['pre_seasonal_departure'][:] 
    pre_seasonal_departure = np.ma.array(pre_seasonal_departure,mask=(np.isnan(pre_seasonal_departure)))
    yrlist = list(fh.variables['time'][:])
    latitcrs = fh.variables['latitcrs'][:]  
    longicrs = fh.variables['longicrs'][:]  
    fh.close() 

# mm    
if region=='Australia':
    levs = [-40,-30,-20,-10,0,10,20,30,40]
    cMap = plt.get_cmap('RdYlBu',len(levs))
    cmaplist = [cMap(i) for i in range(cMap.N)]
    cmaplist[3] = (0.9,0.9,0.9,1)
    cmaplist[4] = (0.9,0.9,0.9,1)
    cMap = cMap.from_list('Custom cmap', cmaplist, cMap.N)
    norm = BoundaryNorm(levs, ncolors=cMap.N, clip=True)
if region=='MDB':
    levs = [-20,-15,-10,-5,0,5,10,15,20]
    cMap = plt.get_cmap('RdYlBu',len(levs))
    cmaplist = [cMap(i) for i in range(cMap.N)]
    cmaplist[3] = (0.9,0.9,0.9,1)
    cmaplist[4] = (0.9,0.9,0.9,1)
    cMap = cMap.from_list('Custom cmap', cmaplist, cMap.N)
    norm = BoundaryNorm(levs, ncolors=cMap.N, clip=True)
    
plt.clf()    
fig = plt.figure(figsize=(8.27,11.7)) # A4, portrait=8.27,11.7, else 16,11.7     
for yr in yrlist:
    ax = fig.add_subplot(11,4,yrlist.index(yr)+1)
    if yrlist.index(yr)==0:
        ax.set_title('DJF',fontsize=fsize)
        plt.ylabel(str(int(yr)),rotation=90)
    elif yrlist.index(yr)==1:
        ax.set_title('MAM',fontsize=fsize)
    elif yrlist.index(yr)==2:
        ax.set_title('JJA',fontsize=fsize)
    elif yrlist.index(yr)==3:
        ax.set_title('SON',fontsize=fsize)
    elif math.modf(yr)[0]==0:
        plt.ylabel(str(int(yr)),rotation=90)

    m = Basemap(projection='stere',lon_0=135,lat_0=-25.,\
                    llcrnrlat=-45,urcrnrlat=0,\
                    llcrnrlon=85,urcrnrlon=175,\
                    resolution='l',area_thresh=10000)
    x, y = m(longicrs,latitcrs) #m(XLONG[5:-5,5:-5],XLAT[5:-5,5:-5]) 
    m.drawcoastlines()
    h=m.pcolormesh(x,y,wvcont_mm_seasonal_departure[yrlist.index(yr),:,:],norm=norm,cmap=cMap)
    cax = plt.axes([0.125,0.08,0.78,0.015])
    cbar = plt.colorbar(h,cax=cax,ticks=levs,orientation='horizontal',boundaries=levs,extend='both',spacing='proportional')
    cbar.set_label(label='Seasonal departure from climatological water vapour contribution to rainfall in '+region+' [mm]')#,size=fsize)
    cbar.ax.tick_params(labelsize=16)
fig.savefig(dir_in+region+'_'+anom_period_name+'_wvcont_seasonal_departure_mm.png',bbox_inches='tight') 
plt.close()        


# pct    
levs = [-0.1,-0.05,-0.025,-0.01,-0.005,0,0.005,0.01,0.025,0.05,0.1]
cMap = plt.get_cmap('RdYlBu',len(levs))
cmaplist = [cMap(i) for i in range(cMap.N)]
cmaplist[4] = (0.9,0.9,0.9,1)
cmaplist[5] = (0.9,0.9,0.9,1)
cMap = cMap.from_list('Custom cmap', cmaplist, cMap.N)
norm = BoundaryNorm(levs, ncolors=cMap.N, clip=True)
    
plt.clf()    
fig = plt.figure(figsize=(8.27,11.7)) # A4, portrait=8.27,11.7, else 16,11.7     
for yr in yrlist:
    ax = fig.add_subplot(11,4,yrlist.index(yr)+1)
    if yrlist.index(yr)==0:
        ax.set_title('DJF',fontsize=fsize)
        plt.ylabel(str(int(yr)),rotation=90)
    elif yrlist.index(yr)==1:
        ax.set_title('MAM',fontsize=fsize)
    elif yrlist.index(yr)==2:
        ax.set_title('JJA',fontsize=fsize)
    elif yrlist.index(yr)==3:
        ax.set_title('SON',fontsize=fsize)
    elif math.modf(yr)[0]==0:
        plt.ylabel(str(int(yr)),rotation=90)

    m = Basemap(projection='stere',lon_0=135,lat_0=-25.,\
                    llcrnrlat=-45,urcrnrlat=0,\
                    llcrnrlon=85,urcrnrlon=175,\
                    resolution='l',area_thresh=10000)
    x, y = m(longicrs,latitcrs) #m(XLONG[5:-5,5:-5],XLAT[5:-5,5:-5]) 
    m.drawcoastlines()
    h=m.pcolormesh(x,y,wvcont_pct_seasonal_departure[yrlist.index(yr),:,:],norm=norm,cmap=cMap)
    
    # Add wind vectors - NOTE THAT THESE MAY NEED TO BE ROTATED. CHECK PY DOCS    
    
    
    cax = plt.axes([0.125,0.08,0.78,0.015])
    cbar = plt.colorbar(h,cax=cax,ticks=levs,orientation='horizontal',boundaries=levs,extend='both')#,spacing='proportional')
    cbar.set_label(label='Seasonal departure from climatological water vapour contribution to '+region+' rainfall [%]')#,size=fsize)
    #cbar.ax.tick_params(labelsize=16)
#plt.tight_layout()    
fig.savefig(dir_in+region+'_'+anom_period_name+'_wvcont_seasonal_departure_pct.png',bbox_inches='tight') 
plt.close()        