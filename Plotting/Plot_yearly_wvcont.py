# -*- coding: utf-8 -*-
"""
Created on Tue Oct  9 14:38:37 2018

@author: z3131380
"""
import numpy as np
from netCDF4 import Dataset 
import matplotlib.pyplot as plt
import os
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import BoundaryNorm

fsize=18

dir_in = '/srv/ccrc/data19/z3131380/PartB/Output/'
domain = 'Australia'
region = 'SWWA'

#==============================================================================
# Plot contributions as percent total annual rainfall over the domain
#==============================================================================
levs = [10**-7,10**-6,
        10**-5,2.5*10**-5,5*10**-5,7.5*10**-5,
        10**-4,2.5*10**-4,5*10**-4,7.5*10**-4]
cMap = plt.get_cmap('RdYlBu_r',len(levs))
cmaplist = [cMap(i) for i in range(cMap.N)]
cMap = cMap.from_list('Custom cmap', cmaplist, cMap.N)
norm = BoundaryNorm(levs, ncolors=cMap.N, clip=True)

for year in np.arange(1979,2014):
    fname = region+'_'+str(year)+'_wvcont.nc'
    file = dir_in+domain+'/100parcels/TS10min/Processed/Yearly/'+fname
    if os.path.isfile(file) == True:
        fh = Dataset(file, mode='r') 
        latitcrs = fh.variables['latitcrs'][:]
        longicrs = fh.variables['longicrs'][:]
        pre = fh.variables['pre'][:]
        wv_cont_sum_yearly_mm = fh.variables['wv_cont_sum_yearly_mm'][:]
        wv_cont_sum_yearly_pct = fh.variables['wv_cont_sum_yearly_pct'][:]
        fh.close()      
        pre_cells = np.ma.array(pre,mask=np.isnan(pre))
        # If it rained that day, plot the contribution.
        if pre_cells.count() > 0:
            plt.clf()    
            fig = plt.figure(figsize=(16,11.7)) # A4, portrait=8.27,11.7, else 16,11.7            
            m = Basemap(projection='stere',lon_0=135,lat_0=-25.,\
                            llcrnrlat=-45,urcrnrlat=0,\
                            llcrnrlon=85,urcrnrlon=175,\
                            resolution='l',area_thresh=10000)
            x, y = m(longicrs,latitcrs) #m(XLONG[5:-5,5:-5],XLAT[5:-5,5:-5]) 
            m.drawcoastlines()
            m.drawmeridians(np.arange(-180,180,10),labels=[0,0,1,0],fontsize=fsize)
            m.drawparallels(np.arange(-90,90,5),labels=[1,0,0,0],fontsize=fsize)
            h=m.pcolormesh(x,y,wv_cont_sum_yearly_pct,norm=norm,cmap=cMap)
            i=m.contour(x,y,np.nansum(pre,axis=0),linewidths=1.9,cmap=plt.cm.jet_r) #colors='gray',
            #plt.clabel(i,fontsize=9,inline=1)
            cax = plt.axes([0.125,0.125,0.78,0.02])
            cbar = plt.colorbar(h,cax=cax,ticks=levs,orientation='horizontal',boundaries=levs,format='%.1e')#,spacing='proportional')
            cbar.set_label(label='Water Vapour Contribution [%] \n Year = '+str(year),size=fsize)
            cbar.ax.tick_params(labelsize=16)
            cax2 = plt.axes([0.91,0.2,0.02,0.6])
            cbar2 = plt.colorbar(i,cax=cax2)
            cbar2.set_ticklabels([250,500,750,1000,1250,1500,2000])
            cbar2.outline.set_visible(False)
            cbar2.set_label(label='Precipitation [mm]',size=fsize)
            cbar2.ax.tick_params(labelsize=16)
            fig.savefig(dir_in+domain+'/100parcels/TS10min/Processed/Yearly/Plots/'+region+'_exp01_wv_cont_pct_'+str(year)+'.png',bbox_inches='tight') #automate this line
            plt.close()
#==============================================================================
# Plot contributions as depth total annual rainfall over the domain
#==============================================================================
if region == 'Australia':
    levs = [0,10,20,30,40,50,100,250,500,1000]
else:
    levs = [0,20,40,60,80,100,200]    
cMap = plt.get_cmap('RdYlBu_r',len(levs))
cmaplist = [cMap(i) for i in range(cMap.N)]
cMap = cMap.from_list('Custom cmap', cmaplist, cMap.N)
norm = BoundaryNorm(levs, ncolors=cMap.N, clip=True)

for year in np.arange(1979,2014):
    fname = region+'_'+str(year)+'_wvcont.nc'
    file = dir_in+domain+'/100parcels/TS10min/Processed/Yearly/'+fname
    if os.path.isfile(file) == True:
        fh = Dataset(file, mode='r') 
        latitcrs = fh.variables['latitcrs'][:]
        longicrs = fh.variables['longicrs'][:]
        pre_annual_total = fh.variables['pre_annual_total'][:]
        wv_cont_sum_yearly_mm = fh.variables['wv_cont_sum_yearly_mm'][:]
        wv_cont_sum_yearly_pct = fh.variables['wv_cont_sum_yearly_pct'][:]
        fh.close()      
        pre_cells = np.ma.array(pre,mask=np.isnan(pre))
        # If it rained that day, plot the contribution.
        if pre_cells.count() > 0:    
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
            h=m.pcolormesh(x,y,wv_cont_sum_yearly_mm,norm=norm,cmap=cMap)
            i=m.contour(x,y,pre_annual_total,linewidths=1.9,cmap=plt.cm.jet_r)#,colors='gray'
            #plt.clabel(i,fontsize=9,inline=1)
            cax = plt.axes([0.125,0.125,0.78,0.02])
            cbar = plt.colorbar(h,cax=cax,ticks=levs,orientation='horizontal',boundaries=levs,extend='max',spacing='proportional')
            cbar.set_label(label='Water Vapour Contribution [mm] \n Year = '+str(year),size=fsize)
            cbar.ax.tick_params(labelsize=16)
            cax2 = plt.axes([0.91,0.2,0.02,0.6])
            cbar2 = plt.colorbar(i,cax=cax2)#,ticks=[250,500,750,1000,1250,1500,2000])
            cbar2.set_ticklabels([250,500,750,1000,1250,1500,2000])
            cbar2.outline.set_visible(False)
            cbar2.set_label(label='Precipitation [mm]',size=fsize)
            cbar2.ax.tick_params(labelsize=16)
            fig.savefig(dir_in+domain+'/100parcels/TS10min/Processed/Yearly/Plots/'+region+'_exp01_wv_cont_mm_'+str(year)+'.png',bbox_inches='tight') 
            plt.close()

#==============================================================================
# Create gif of images
#==============================================================================
#import imageio
#import glob
#
#out_dir = r'/srv/ccrc/data19/z3131380/PartB/Output/Australia/100parcels/TS10min/Processed/Yearly/gif/'
#
#file_names = glob.glob('/srv/ccrc/data19/z3131380/PartB/Output/Australia/100parcels/TS10min/Processed/Yearly/Plots/*_mm_*')
#
#with imageio.get_writer(out_dir+'wv_cont_mm_1979-2013.gif', mode='I',duration=2) as writer:
#    for filename in file_names:
#        image = imageio.imread(filename)
#        writer.append_data(image)

#==============================================================================
# Place in subplots
#==============================================================================
import pandas as pd

Yearlist = pd.date_range('1979','2013',freq='AS').year

if region == 'Australia':
    levs = [0,10,20,30,40,50,100,150,200,250,300,500]
else:
    levs = [0,20,40,60,80,100,200]  
cMap = plt.get_cmap('RdYlBu_r',len(levs))
cmaplist = [cMap(i) for i in range(cMap.N)]
cMap = cMap.from_list('Custom cmap', cmaplist, cMap.N)
norm = BoundaryNorm(levs, ncolors=cMap.N, clip=True)

fig = plt.figure(figsize=(18,16.5)) # A4, portrait=8.27,11.7, else 16,11.7
for year in range(len(Yearlist)):
    fname = region+'_'+str(Yearlist[year])+'_wvcont.nc'
    file = dir_in+domain+'/100parcels/TS10min/Processed/Yearly/'+fname
    if os.path.isfile(file) == True:
        fh = Dataset(file, mode='r') 
        latitcrs = fh.variables['latitcrs'][:]
        longicrs = fh.variables['longicrs'][:]
        pre = fh.variables['pre'][:]
        wv_cont_sum_yearly_mm = fh.variables['wv_cont_sum_yearly_mm'][:]
        #wv_cont_sum_yearly_pct = fh.variables['wv_cont_sum_yearly_pct'][:]
        fh.close()      
        #else: print 'File for '+Yearlist[year]+' does not exist'
        pre_cells = np.ma.array(pre,mask=np.isnan(pre))
        
        ax = fig.add_subplot(7,5,year+1)
        ax.set_title(Yearlist[year],fontsize=14)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(False)        
        m = Basemap(projection='stere',lon_0=135,lat_0=-25.,\
                        llcrnrlat=-45,urcrnrlat=0,\
                        llcrnrlon=85,urcrnrlon=175,\
                        resolution='i',area_thresh=10000)
        x, y = m(longicrs,latitcrs) #m(XLONG[5:-5,5:-5],XLAT[5:-5,5:-5]) 
        m.drawcoastlines()
    #    m.drawmeridians(np.arange(-180,180,10),labels=[0,0,1,0],fontsize=fsize)
    #    m.drawparallels(np.arange(-90,90,5),labels=[1,0,0,0],fontsize=fsize)
        h=m.pcolormesh(x,y,wv_cont_sum_yearly_mm,norm=norm,cmap=cMap)
    #    i=m.contour(x,y,np.nansum(pre,axis=0),linewidths=1.9,cmap=plt.cm.jet_r)#,colors='gray'
    #    #plt.clabel(i,fontsize=9,inline=1)
cax = plt.axes([0.05,-0.015,0.9,0.02])
cbar = plt.colorbar(h,cax=cax,ticks=levs,orientation='horizontal',boundaries=levs,extend='max')#,spacing='proportional')
cbar.set_label(label='Water Vapour Contribution [mm]',size=fsize)
cbar.ax.tick_params(labelsize=16)
##cax2 = plt.axes([0.99,0.2,0.02,0.6])
#cax2 = plt.axes([0.05,-0.09,0.9,0.02])
#cbar2 = plt.colorbar(i,cax=cax2,orientation='horizontal')#,ticks=[250,500,750,1000,1250,1500,2000])
##cbar2.set_ticklabels([250,500,750,1000,1250,1500,2000])
#cbar2.outline.set_visible(False)
#cbar2.set_label(label='Precipitation [mm]',size=fsize)
#cbar2.ax.tick_params(labelsize=16)
plt.tight_layout()
fig.subplots_adjust(wspace=0)#, hspace=0)
fig.savefig(dir_in+domain+'/100parcels/TS10min/Processed/Yearly/Plots/'+region+'_exp01_wvcont_mm_1979-2013.png',bbox_inches='tight') #automate this line
plt.close()