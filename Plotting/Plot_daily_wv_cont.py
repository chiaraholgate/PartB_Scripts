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

levs = [0.1,0.25,0.5,1,1.5,2,4,6,8,10,12,14,16,18,20]
cMap = plt.get_cmap('RdYlBu_r',len(levs))
cmaplist = [cMap(i) for i in range(cMap.N)]
cMap = cMap.from_list('Custom cmap', cmaplist, cMap.N)
norm = BoundaryNorm(levs, ncolors=cMap.N, clip=True)

year=1979;month='01';day=31
dir_in = '/srv/ccrc/data19/z3131380/PartB/Output/'
region = 'Australia'
fname = region+'_'+str(year)+month+'_'+str(day)+'_wvcont.nc'
file = dir_in+region+'/100parcels/TS10min/Processed/Daily/'+fname
if os.path.isfile(file) == True:
    print fname
    fh = Dataset(file, mode='r') 
    latitcrs = fh.variables['latitcrs'][:]
    longicrs = fh.variables['longicrs'][:]
    pre = fh.variables['pre'][:]
    wv_cont_sum_daily_mm = fh.variables['wv_cont_sum_daily_mm'][:]
    wv_cont_sum_daily_pct = fh.variables['wv_cont_sum_daily_pct'][:]
    fh.close()      
pre_cells = np.ma.array(pre,mask=np.isnan(pre))
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
i=m.contour(x,y,pre,colors='gray',linewidths=1.5)
#plt.clabel(i,fontsize=9,inline=1)
cax = plt.axes([0.03,0.03,0.95,0.02])
cbar = plt.colorbar(h,cax=cax,ticks=levs,orientation='horizontal',boundaries=levs,label='Water Vapour Contribution [%]')
cbar.ax.set_xticklabels(levs)
#fig.savefig(dir_out+str(day)+'_June1981_30day_backtrack_'+str(p)+'_parcels.png',bbox_inches='tight') #automate this line
