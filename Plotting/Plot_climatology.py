# -*- coding: utf-8 -*-
"""
Created on Tue Oct  9 14:38:37 2018

*** TO DO:
- Plot the mean wind vectors over the climatology. 

@author: z3131380
"""
import numpy as np
from netCDF4 import Dataset 
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import BoundaryNorm

fsize=18
fsize_coords = 12

dir_in= r'/srv/ccrc/data19/z3131380/PartB/Output/Australia/100parcels/TS10min/exp01/Processed/Climatology/'
domain = 'Australia'
region = 'MDB'
seasons = ['DJF','MAM','JJA','SON']

#==============================================================================
# Plot contributions as depth total annual rainfall over the domain
#==============================================================================
#levs = [0,10,20,30,40,50,100,250,500]
#levs = [0,10,20,30,40,50,60,70,80,90,100]
if region == 'Australia':
    levs = [0,5,10,20,40,60,80,100,200]    
else:
    levs = [0,2,5,10,20,30,40,50] 
cMap = plt.get_cmap('RdYlBu_r',len(levs))
cmaplist = [cMap(i) for i in range(cMap.N)]
cMap = cMap.from_list('Custom cmap', cmaplist, cMap.N)
norm = BoundaryNorm(levs, ncolors=cMap.N, clip=True)

file = dir_in+region+'_1979-2013_climatology.nc'
#file = dir_in+domain+'/100parcels/TS10min/Processed/Yearly/'+fname
fh = Dataset(file, mode='r') 
wvcont_mm_seasonal_climatology = fh.variables['wvcont_mm_seasonal_climatology'][:]
wvcont_pct_seasonal_climatology = fh.variables['wvcont_pct_seasonal_climatology'][:]
pre_seasonal_climatology = fh.variables['pre_seasonal_climatology'][:]
latitcrs = fh.variables['latitcrs'][:]
longicrs = fh.variables['longicrs'][:]
fh.close()      
    
for s in seasons:
    plt.clf()    
    fig = plt.figure(figsize=(16,11.7)) # A4, portrait=8.27,11.7, else 16,11.7            
    m = Basemap(projection='stere',lon_0=135,lat_0=-25.,\
                    llcrnrlat=-45,urcrnrlat=0,\
                    llcrnrlon=85,urcrnrlon=175,\
                    resolution='l',area_thresh=10000)
    x, y = m(longicrs,latitcrs) #m(XLONG[5:-5,5:-5],XLAT[5:-5,5:-5]) 
    m.drawcoastlines()
    #m.readshapefile('/home/z3131380/hdrive/PhD/PartB/Masks/mdb_sw_bdys/MDB-WRPA-SurfaceWater/Murray-Darling Basin Water Resource Plan Areas - Surface Water','wsmask',color='g',linewidth=2)
    m.drawmeridians(np.arange(-180,180,10),labels=[0,0,1,0],fontsize=fsize)
    m.drawparallels(np.arange(-90,90,5),labels=[1,0,0,0],fontsize=fsize)
    h=m.pcolormesh(x,y,wvcont_mm_seasonal_climatology[seasons.index(s),:,:],norm=norm,cmap=cMap)
    #i=m.contour(x,y,pre_seasonal_climatology[seasons.index(s),:,:],linewidths=1.5,cmap=plt.cm.jet_r)#,colors='gray'
    ##plt.clabel(i,fontsize=9,inline=1)
    cax = plt.axes([0.125,0.125,0.78,0.02])
    cbar = plt.colorbar(h,cax=cax,ticks=levs,orientation='horizontal',boundaries=levs,extend='max',spacing='proportional')
    cbar.set_label(label='Mean seasonal water vapour contribution [mm]',size=fsize)
    cbar.ax.tick_params(labelsize=16)
#    cax2 = plt.axes([0.91,0.2,0.02,0.6])
#    cbar2 = plt.colorbar(i,cax=cax2)#,ticks=[250,500,750,1000,1250,1500,2000])
#    #cbar2.set_ticklabels([250,500,750,1000,1250,1500,2000])
#    cbar2.outline.set_visible(False)
#    cbar2.set_label(label='Precipitation [mm]',size=fsize)
#    cbar2.ax.tick_params(labelsize=16)
    fig.savefig(dir_in+region+'_1979-2013_wvcont_mm_'+s+'_climatology.png',bbox_inches='tight') 
    plt.close()

# Subplots
plt.clf()    
fig = plt.figure(figsize=(16,11.7)) # A4, portrait=8.27,11.7, else 16,11.7     
for s in seasons:
    ax = fig.add_subplot(2,2,seasons.index(s)+1)
    ax.set_title(s,fontsize=fsize)
#    ax.spines['top'].set_visible(False)
#    ax.spines['right'].set_visible(False)
#    ax.spines['bottom'].set_visible(False)
#    ax.spines['left'].set_visible(False)        
    m = Basemap(projection='stere',lon_0=135,lat_0=-25.,\
                    llcrnrlat=-45,urcrnrlat=0,\
                    llcrnrlon=85,urcrnrlon=175,\
                    resolution='i',area_thresh=10000)
    x, y = m(longicrs,latitcrs) #m(XLONG[5:-5,5:-5],XLAT[5:-5,5:-5]) 
    m.drawcoastlines()
    #m.readshapefile('/srv/ccrc/data03/z3131380/PartB/Masks/mdb_boundary/mdb_boundary','wsmask',color='g',linewidth=1)
    h=m.pcolormesh(x,y,wvcont_mm_seasonal_climatology[seasons.index(s),:,:],norm=norm,cmap=cMap)
    i=m.contour(x,y,pre_seasonal_climatology[seasons.index(s),:,:],[50,100,200,300,400,500,600],linewidths=1.5,cmap=plt.cm.jet_r)#,colors='gray'
    if s=='DJF': 
        m.drawparallels(np.arange(-90,90,5),labels=[1,0,0,0],fontsize=fsize_coords)
        m.drawmeridians(np.arange(-180,180,10),labels=[0,0,0,1],fontsize=fsize_coords)
    if s == 'MAM':
        m.drawmeridians(np.arange(-180,180,10),labels=[0,0,0,1],fontsize=fsize_coords)
        m.drawparallels(np.arange(-90,90,5),labels=[0,1,0,0],fontsize=fsize_coords)
    if s=='JJA': 
        m.drawparallels(np.arange(-90,90,5),labels=[1,0,0,0],fontsize=fsize_coords)
        m.drawmeridians(np.arange(-180,180,10),labels=[0,0,0,1],fontsize=fsize_coords)
    if s == 'SON':
        m.drawmeridians(np.arange(-180,180,10),labels=[0,0,0,1],fontsize=fsize_coords)
        m.drawparallels(np.arange(-90,90,5),labels=[0,1,0,0],fontsize=fsize_coords)
    
#plt.clabel(i,fontsize=9,inline=1)
cax = plt.axes([0.01,-0.05,0.97,0.03])
cbar = plt.colorbar(h,cax=cax,ticks=levs,orientation='horizontal',boundaries=levs,extend='max',spacing='proportional')
cbar.set_label(label='Mean seasonal water vapour contribution [mm]',size=fsize)
cbar.ax.tick_params(labelsize=16)
#cax2 = plt.axes([0.91,0.2,0.02,0.6])
cax2 = plt.axes([0.01,-0.15,0.97,0.03])
cbar2 = plt.colorbar(i,cax=cax2,orientation='horizontal')#,ticks=[250,500,750,1000,1250,1500,2000])
#cbar2.set_ticklabels([10,250,500,750,1000,1250,1500,2000])
cbar2.outline.set_visible(False)
cbar2.set_label(label='Precipitation [mm]',size=fsize)
cbar2.ax.tick_params(labelsize=16)
plt.tight_layout()
#fig.subplots_adjust(wspace=0)#, hspace=0)
fig.savefig(dir_in+region+'_1979-2013_wvcont_mm_seasonal_climatology.png',bbox_inches='tight') 
plt.close()

#==============================================================================
# Plot annual and seasonal on same figure [mm]
#==============================================================================
if region == 'Australia':
    levs = [0,5,10,20,40,60,80,100,200]    
else:
    levs = [0,2,5,10,20,30,40,50] 
cMap = plt.get_cmap('RdYlBu_r',len(levs))
cmaplist = [cMap(i) for i in range(cMap.N)]
cMap = cMap.from_list('Custom cmap', cmaplist, cMap.N)
norm = BoundaryNorm(levs, ncolors=cMap.N, clip=True)

file = dir_in+region+'_1979-2013_climatology.nc'
fh = Dataset(file, mode='r') 
wvcont_mm_seasonal_climatology = fh.variables['wvcont_mm_seasonal_climatology'][:]
pre_seasonal_climatology = fh.variables['pre_seasonal_climatology'][:]
wvcont_mm_AprNov_climatology = fh.variables['wvcont_mm_AprNov_climatology'][:]
pre_AprNov_climatology = fh.variables['pre_AprNov_climatology'][:]
wvcont_mm_annual_climatology = fh.variables['wvcont_mm_annual_climatology'][:]
pre_annual_climatology = fh.variables['pre_annual_climatology'][:]
latitcrs = fh.variables['latitcrs'][:]
longicrs = fh.variables['longicrs'][:]
fh.close()     

plt.clf()    
fig = plt.figure(figsize=(12,11.7)) # A4, portrait=8.27,11.7, else 16,11.7     
ax = fig.add_subplot(3,2,1)
ax.set_title('Annual',fontsize=fsize)
m = Basemap(projection='stere',lon_0=135,lat_0=-25.,\
                    llcrnrlat=-45,urcrnrlat=0,\
                    llcrnrlon=85,urcrnrlon=175,\
                    resolution='i',area_thresh=10000)
x, y = m(longicrs,latitcrs) 
m.drawcoastlines()
h=m.pcolormesh(x,y,wvcont_mm_annual_climatology,norm=norm,cmap=cMap)
i=m.contour(x,y,pre_annual_climatology,[50,100,200,300,400,500,600],linewidths=1.5,cmap=plt.cm.jet_r)#,colors='gray'
    
ax = fig.add_subplot(3,2,2)
ax.set_title('April - November',fontsize=fsize)
m = Basemap(projection='stere',lon_0=135,lat_0=-25.,\
                    llcrnrlat=-45,urcrnrlat=0,\
                    llcrnrlon=85,urcrnrlon=175,\
                    resolution='i',area_thresh=10000)
x, y = m(longicrs,latitcrs) 
m.drawcoastlines()
h=m.pcolormesh(x,y,wvcont_mm_AprNov_climatology,norm=norm,cmap=cMap)
i=m.contour(x,y,pre_AprNov_climatology,[50,100,200,300,400,500,600],linewidths=1.5,cmap=plt.cm.jet_r)#,colors='gray'
    
ax = fig.add_subplot(3,2,3)
ax.set_title('Summer',fontsize=fsize)
m = Basemap(projection='stere',lon_0=135,lat_0=-25.,\
                    llcrnrlat=-45,urcrnrlat=0,\
                    llcrnrlon=85,urcrnrlon=175,\
                    resolution='i',area_thresh=10000)
x, y = m(longicrs,latitcrs) 
m.drawcoastlines()
h=m.pcolormesh(x,y,wvcont_mm_seasonal_climatology[0,:,:],norm=norm,cmap=cMap)
i=m.contour(x,y,pre_seasonal_climatology[0,:,:],[50,100,200,300,400,500,600],linewidths=1.5,cmap=plt.cm.jet_r)#,colors='gray'

ax = fig.add_subplot(3,2,4)
ax.set_title('Autumn',fontsize=fsize)
m = Basemap(projection='stere',lon_0=135,lat_0=-25.,\
                    llcrnrlat=-45,urcrnrlat=0,\
                    llcrnrlon=85,urcrnrlon=175,\
                    resolution='i',area_thresh=10000)
x, y = m(longicrs,latitcrs)
m.drawcoastlines()
h=m.pcolormesh(x,y,wvcont_mm_seasonal_climatology[1,:,:],norm=norm,cmap=cMap)
i=m.contour(x,y,pre_seasonal_climatology[1,:,:],[50,100,200,300,400,500,600],linewidths=1.5,cmap=plt.cm.jet_r)#,colors='gray'


ax = fig.add_subplot(3,2,5)
ax.set_title('Winter',fontsize=fsize)
m = Basemap(projection='stere',lon_0=135,lat_0=-25.,\
                    llcrnrlat=-45,urcrnrlat=0,\
                    llcrnrlon=85,urcrnrlon=175,\
                    resolution='i',area_thresh=10000)
x, y = m(longicrs,latitcrs) 
m.drawcoastlines()
h=m.pcolormesh(x,y,wvcont_mm_seasonal_climatology[2,:,:],norm=norm,cmap=cMap)
i=m.contour(x,y,pre_seasonal_climatology[2,:,:],[50,100,200,300,400,500,600],linewidths=1.5,cmap=plt.cm.jet_r)#,colors='gray'


ax = fig.add_subplot(3,2,6)
ax.set_title('Spring',fontsize=fsize)
m = Basemap(projection='stere',lon_0=135,lat_0=-25.,\
                    llcrnrlat=-45,urcrnrlat=0,\
                    llcrnrlon=85,urcrnrlon=175,\
                    resolution='i',area_thresh=10000)
x, y = m(longicrs,latitcrs)
m.drawcoastlines()
h=m.pcolormesh(x,y,wvcont_mm_seasonal_climatology[3,:,:],norm=norm,cmap=cMap)
i=m.contour(x,y,pre_seasonal_climatology[3,:,:],[50,100,200,300,400,500,600],linewidths=1.5,cmap=plt.cm.jet_r)#,colors='gray'

cax = plt.axes([0.05,-0.03,0.9,0.03])
cbar = plt.colorbar(h,cax=cax,ticks=levs,orientation='horizontal',boundaries=levs,extend='max',spacing='proportional')
cbar.set_label(label='Mean water vapour contribution [mm]',size=fsize)
cbar.ax.tick_params(labelsize=16)
#cax2 = plt.axes([0.91,0.2,0.02,0.6])
cax2 = plt.axes([0.05,-0.13,0.9,0.03])
cbar2 = plt.colorbar(i,cax=cax2,orientation='horizontal')#,ticks=[250,500,750,1000,1250,1500,2000])
##cbar2.set_ticklabels([10,250,500,750,1000,1250,1500,2000])
cbar2.outline.set_visible(False)
cbar2.set_label(label='Precipitation [mm]',size=fsize)
cbar2.ax.tick_params(labelsize=16)
plt.tight_layout()
fig.subplots_adjust(wspace=0)#, hspace=0)
fig.savefig(dir_in+region+'_1979-2013_wvcont_mm_seasonal_climatology.png',bbox_inches='tight') 
plt.close()

#==============================================================================
# Plot annual and seasonal on same figure [%]
#==============================================================================
if region == 'Australia':
    levs = [0,0.0025,0.005,0.01,0.015,0.02,0.03,0.04,0.05,0.06] 
else:
    levs = [0,0.0025,0.005,0.01,0.015,0.02,0.03,0.04,0.05,0.06] 
cMap = plt.get_cmap('RdYlBu_r',len(levs))
cmaplist = [cMap(i) for i in range(cMap.N)]
cMap = cMap.from_list('Custom cmap', cmaplist, cMap.N)
norm = BoundaryNorm(levs, ncolors=cMap.N, clip=True)

file = dir_in+region+'_1979-2013_climatology.nc'
fh = Dataset(file, mode='r') 
wvcont_pct_seasonal_climatology = fh.variables['wvcont_pct_seasonal_climatology'][:]
pre_seasonal_climatology = fh.variables['pre_seasonal_climatology'][:]
wvcont_pct_AprNov_climatology = fh.variables['wvcont_pct_AprNov_climatology'][:]
pre_AprNov_climatology = fh.variables['pre_AprNov_climatology'][:]
wvcont_pct_annual_climatology = fh.variables['wvcont_pct_annual_climatology'][:]
pre_annual_climatology = fh.variables['pre_annual_climatology'][:]
latitcrs = fh.variables['latitcrs'][:]
longicrs = fh.variables['longicrs'][:]
fh.close()  

plt.clf()    
fig = plt.figure(figsize=(12,11.7)) # A4, portrait=8.27,11.7, else 16,11.7     
ax = fig.add_subplot(3,2,1)
ax.set_title('Annual',fontsize=fsize)
m = Basemap(projection='stere',lon_0=135,lat_0=-25.,\
                    llcrnrlat=-45,urcrnrlat=0,\
                    llcrnrlon=85,urcrnrlon=175,\
                    resolution='i',area_thresh=10000)
x, y = m(longicrs,latitcrs) 
m.drawcoastlines()
h=m.pcolormesh(x,y,wvcont_pct_annual_climatology,norm=norm,cmap=cMap)
i=m.contour(x,y,pre_annual_climatology,[50,100,200,300,400,500,600],linewidths=1.5,cmap=plt.cm.jet_r)#,colors='gray'
    
ax = fig.add_subplot(3,2,2)
ax.set_title('April - November',fontsize=fsize)
m = Basemap(projection='stere',lon_0=135,lat_0=-25.,\
                    llcrnrlat=-45,urcrnrlat=0,\
                    llcrnrlon=85,urcrnrlon=175,\
                    resolution='i',area_thresh=10000)
x, y = m(longicrs,latitcrs) 
m.drawcoastlines()
h=m.pcolormesh(x,y,wvcont_pct_AprNov_climatology,norm=norm,cmap=cMap)
i=m.contour(x,y,pre_AprNov_climatology,[50,100,200,300,400,500,600],linewidths=1.5,cmap=plt.cm.jet_r)#,colors='gray'
    
ax = fig.add_subplot(3,2,3)
ax.set_title('Summer',fontsize=fsize)
m = Basemap(projection='stere',lon_0=135,lat_0=-25.,\
                    llcrnrlat=-45,urcrnrlat=0,\
                    llcrnrlon=85,urcrnrlon=175,\
                    resolution='i',area_thresh=10000)
x, y = m(longicrs,latitcrs) 
m.drawcoastlines()
h=m.pcolormesh(x,y,wvcont_pct_seasonal_climatology[0,:,:],norm=norm,cmap=cMap)
i=m.contour(x,y,pre_seasonal_climatology[0,:,:],[50,100,200,300,400,500,600],linewidths=1.5,cmap=plt.cm.jet_r)#,colors='gray'

ax = fig.add_subplot(3,2,4)
ax.set_title('Autumn',fontsize=fsize)
m = Basemap(projection='stere',lon_0=135,lat_0=-25.,\
                    llcrnrlat=-45,urcrnrlat=0,\
                    llcrnrlon=85,urcrnrlon=175,\
                    resolution='i',area_thresh=10000)
x, y = m(longicrs,latitcrs)
m.drawcoastlines()
h=m.pcolormesh(x,y,wvcont_pct_seasonal_climatology[1,:,:],norm=norm,cmap=cMap)
i=m.contour(x,y,pre_seasonal_climatology[1,:,:],[50,100,200,300,400,500,600],linewidths=1.5,cmap=plt.cm.jet_r)#,colors='gray'


ax = fig.add_subplot(3,2,5)
ax.set_title('Winter',fontsize=fsize)
m = Basemap(projection='stere',lon_0=135,lat_0=-25.,\
                    llcrnrlat=-45,urcrnrlat=0,\
                    llcrnrlon=85,urcrnrlon=175,\
                    resolution='i',area_thresh=10000)
x, y = m(longicrs,latitcrs) 
m.drawcoastlines()
h=m.pcolormesh(x,y,wvcont_pct_seasonal_climatology[2,:,:],norm=norm,cmap=cMap)
i=m.contour(x,y,pre_seasonal_climatology[2,:,:],[50,100,200,300,400,500,600],linewidths=1.5,cmap=plt.cm.jet_r)#,colors='gray'


ax = fig.add_subplot(3,2,6)
ax.set_title('Spring',fontsize=fsize)
m = Basemap(projection='stere',lon_0=135,lat_0=-25.,\
                    llcrnrlat=-45,urcrnrlat=0,\
                    llcrnrlon=85,urcrnrlon=175,\
                    resolution='i',area_thresh=10000)
x, y = m(longicrs,latitcrs)
m.drawcoastlines()
h=m.pcolormesh(x,y,wvcont_pct_seasonal_climatology[3,:,:],norm=norm,cmap=cMap)
i=m.contour(x,y,pre_seasonal_climatology[3,:,:],[50,100,200,300,400,500,600],linewidths=1.5,cmap=plt.cm.jet_r)#,colors='gray'

cax = plt.axes([0.05,-0.03,0.9,0.03])
cbar = plt.colorbar(h,cax=cax,orientation='horizontal',boundaries=levs,extend='max',\
spacing='proportional',format='%.3f')
cbar.set_label(label='Mean water vapour contribution [%]',size=fsize)
cbar.ax.tick_params(labelsize=16)
cbar.ax.set_xticklabels(levs,rotation=45)
#cax2 = plt.axes([0.91,0.2,0.02,0.6])
cax2 = plt.axes([0.05,-0.18,0.9,0.03])
cbar2 = plt.colorbar(i,cax=cax2,orientation='horizontal')#,ticks=[250,500,750,1000,1250,1500,2000])
##cbar2.set_ticklabels([10,250,500,750,1000,1250,1500,2000])
cbar2.outline.set_visible(False)
cbar2.set_label(label='Precipitation [mm]',size=fsize)
cbar2.ax.tick_params(labelsize=16)
plt.tight_layout()
fig.subplots_adjust(wspace=0)#, hspace=0)
fig.savefig(dir_in+region+'_1979-2013_wvcont_pct_seasonal_climatology.png',bbox_inches='tight') 
plt.close()