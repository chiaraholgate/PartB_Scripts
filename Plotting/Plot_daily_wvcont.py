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
import pandas as pd


dir_in = '/srv/ccrc/data19/z3131380/PartB/Output/Australia/100parcels/TS10min/Processed/Daily/'
dir_out = r'/srv/ccrc/data19/z3131380/PartB/Output/Australia/100parcels/TS10min/Processed/Daily/Plots/'
region = 'Australia'
fsize=18

Start_date = '19790101' ; Start_year = Start_date[0:4] ; Start_month = Start_date[4:6]; Start_day = Start_date[6:8]
End_date = '20131231' ; End_year = End_date[0:4] ; End_month = End_date[4:6]; End_day = End_date[6:8]
daylist = pd.date_range(Start_date,End_date,freq='d') 

if region == 'Australia':
    levs = [0.1,0.25,0.5,1,1.5,2,4,6,8,10,12,14,16,18,20]
if region == 'MDB':
    levs = [0,0.001,0.01,0.1,0.25,0.5,1,1.5,2,4,6,8,10]
cMap = plt.get_cmap('RdYlBu_r',len(levs))
cmaplist = [cMap(i) for i in range(cMap.N)]
cMap = cMap.from_list('Custom cmap', cmaplist, cMap.N)
norm = BoundaryNorm(levs, ncolors=cMap.N, clip=True)

for i in range(len(daylist)):
    yyyy = daylist[i].year; mm = '%02u' % daylist[i].month; dd = daylist[i].day
    fname = region+'_'+str(yyyy)+mm+'_'+str(dd)+'_wvcont.nc'
    file = dir_in+fname
    if os.path.isfile(file) == True:
        fh = Dataset(file, mode='r') 
        latitcrs = fh.variables['latitcrs'][:]
        longicrs = fh.variables['longicrs'][:]
        pre = fh.variables['pre'][:]
        wv_cont_sum_daily_mm = fh.variables['wv_cont_sum_daily_mm'][:]
        #wv_cont_sum_daily_pct = fh.variables['wv_cont_sum_daily_pct'][:]
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
            m.drawmeridians(np.arange(-180,180,10),labels=[0,0,1,0])
            m.drawparallels(np.arange(-90,90,5),labels=[1,0,0,0])
            h=m.pcolormesh(x,y,wv_cont_sum_daily_mm,norm=norm,cmap=cMap)
            j=m.contour(x,y,pre,cmap=plt.cm.jet_r,linewidths=1.5)
            #plt.clabel(j,fontsize=9,inline=1)
            cax = plt.axes([0.125,0.125,0.78,0.02])
            cbar = plt.colorbar(h,cax=cax,orientation='horizontal',boundaries=levs,extend='max',ticks=levs)#,spacing='proportional')
            labels = [item.get_text() for item in cbar.ax.get_xticklabels()]
            if region == 'Australia':
                labels = levs
            if region == 'MDB':
                labels = [0,0.001,0.01,0.1,0.25,0.5,1,1.5,2,4,6,8]
            cbar.ax.set_xticklabels(labels)
            cbar.set_label(label='Water Vapour Contribution [mm] \n Day = '+str(yyyy)+mm+str(dd),size=fsize)
            cbar.ax.tick_params(labelsize=16)
            cax2 = plt.axes([0.91,0.2,0.02,0.6])
            cbar2 = plt.colorbar(j,cax=cax2)#,ticks=[250,500,750,1000,1250,1500,2000])
            #cbar2.set_ticklabels([250,500,750,1000,1250,1500,2000])
            cbar2.outline.set_visible(False)
            cbar2.set_label(label='Precipitation [mm]',size=fsize)
            cbar2.ax.tick_params(labelsize=16)
            fig.savefig(dir_out+region+'_exp01_wvcont_mm_'+str(yyyy)+mm+'_'+str(dd)+'.png',bbox_inches='tight')  
            plt.close()
        
#==============================================================================
# Create gif of images
#==============================================================================
#import imageio
#import glob
#
#file_names = glob.glob('/srv/ccrc/data19/z3131380/PartB/Output/Australia/100parcels/TS10min/Processed/Daily/Plots/'+region+'*_mm_*.png')
#
#with imageio.get_writer(dir_out+region+'_exp01_wvcont_mm_1979-2013.gif', mode='I',duration=2) as writer:
#    for filename in file_names:
#        image = imageio.imread(filename)
#        writer.append_data(image)    
