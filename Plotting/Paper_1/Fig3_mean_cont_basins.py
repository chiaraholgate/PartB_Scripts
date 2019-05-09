# -*- coding: utf-8 -*-
"""
Created on Tue Apr  2 13:18:32 2019

Figure 3. Mean moisture contribution [%] for precipitation in each season and sub-basin. 
[Have one basin to a row, one season to a column. 13 rows x 4 columns. Show [mm] equivalent in 
supplementary?]

**TO DO:
- Need to have seasonal climatology processed for each sub-basin.
- Add coordinate labels to map
- Draw basin polygon on each subplot?
- put distance scale

@author: z3131380
"""

from netCDF4 import Dataset 
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import BoundaryNorm
import os.path

#==============================================================================
# 

#==============================================================================
dir_in= r'/srv/ccrc/data19/z3131380/PartB/Output/Australia/100parcels/TS10min/exp01/Processed/Climatology/'
dir_out = r'/home/z3131380/hdrive/PhD/PartB/Reporting/Paper_1/Figures/'

domain = 'Australia'
seasons = ['DJF','MAM','JJA','SON']

# MAKE SURE THE NAMES MATCH!
regions = ['CarpentariaCoast','LakeEyreBasin','NorthEastCoast','NorthWesternPlateau','MurrayDarlingBasin',\
                'SouthWestCoast','PilbaraGascoyne','SouthAustralianGulf','SouthEastCoastNSW','SouthEastCoastVictoria',\
                'SouthWesternPlateau','TanamiTimorSeaCoast','Tasmania']

regions_names = ['Carpentaria\nCoast','Lake Eyre\nBasin','North East\nCoast',\
                'North Western\nPlateau','Murray-Darling\nBasin','South West\nCoast','Pilbara\nGascoyne',\
                'South Australian\nGulf','South East\nCoast NSW','South East\nCoast Victoria',\
                'South Western\nPlateau','Tanami-Timor\nSea Coast','Tasmania']          

fsize=10
fsize_coords = 12

levs_mm = [0,2,5,10,15,20,30,40,50] 
levs_pct = [0,0.0025,0.005,0.01,0.02,0.03,0.04,0.05,0.06] 

# build a rectangle in axes coords
left, width = -0.07, .5
bottom, height = 0.05, 0.9
right = left + width
top = bottom + height*1.5
centre = left + width/2

#==============================================================================
# Plot [pct]
#==============================================================================
cMap = plt.get_cmap('RdYlBu_r',len(levs_pct))
cmaplist = [cMap(i) for i in range(cMap.N)]
cMap = cMap.from_list('Custom cmap', cmaplist, cMap.N)
norm = BoundaryNorm(levs_pct, ncolors=cMap.N, clip=True)

fig, axs = plt.subplots(len(regions), 4, sharex='col', sharey='row', figsize=(7,10))  #, #8.27,11.7
axs[0,0].set_title('Summer',size=fsize)
axs[0,1].set_title('Fall',size=fsize)
axs[0,2].set_title('Winter',size=fsize)
axs[0,3].set_title('Spring',size=fsize)      
# Delete once all results finished
for region in regions:
    file = dir_in+region+'_1979-2013_seasonal_climatology.nc'
    if os.path.isfile(file) == True:
        fh = Dataset(file, mode='r') 
        wvcont_pct_seasonal_climatology = fh.variables['wvcont_pct_seasonal_climatology'][:]
        latitcrs = fh.variables['latitcrs'][:]
        longicrs = fh.variables['longicrs'][:]
        fh.close()  
        for xs in range(4):
            ax = axs[regions.index(region),xs]
            if xs == 0:
                ax.text(left,0,regions_names[regions.index(region)], size=fsize,horizontalalignment='right',\
                    verticalalignment='bottom',transform=ax.transAxes)#,rotation=90)
            # Plot the seasonal climatology for each region
            m = Basemap(projection='stere',lon_0=135,lat_0=-25.,\
                                llcrnrlat=-45,urcrnrlat=0,\
                                llcrnrlon=85,urcrnrlon=175,\
                                resolution='i',area_thresh=10000,\
                                ax = ax)
            m.drawcoastlines()
            x, y = m(longicrs,latitcrs) 
            h = ax.pcolormesh(x,y,wvcont_pct_seasonal_climatology[xs],norm=norm,cmap=cMap)

cax = plt.axes([0.01,-0.02,0.98,0.02] )
cbar = plt.colorbar(h,cax=cax,orientation='horizontal',boundaries=levs_pct,extend='max',\
spacing='proportional')
cbar.ax.tick_params(labelsize=fsize)
cbar.ax.set_xticklabels(levs_pct)#,rotation=45)      
cbar.set_label(label='Moisture contribution [%]',size=fsize)#,rotation=90)

#plt.subplot_tool()
#plt.tight_layout()
#fig.subplots_adjust(hspace=0) #wspace=0
fig.savefig(dir_out+'Fig3_pct.png',bbox_inches='tight') 
plt.close()

#==============================================================================
# Plot [mm]
#==============================================================================
cMap = plt.get_cmap('RdYlBu_r',len(levs_mm))
cmaplist = [cMap(i) for i in range(cMap.N)]
cMap = cMap.from_list('Custom cmap', cmaplist, cMap.N)
norm = BoundaryNorm(levs_mm, ncolors=cMap.N, clip=True)

fig, axs = plt.subplots(len(regions), 4, sharex='col', sharey='row', figsize=(8.27,11.7))  #, #8.27,11.7
axs[0,0].set_title('Summer',size=fsize)
axs[0,1].set_title('Fall',size=fsize)
axs[0,2].set_title('Winter',size=fsize)
axs[0,3].set_title('Spring',size=fsize)      
for region in regions:
    file = dir_in+region+'_1979-2013_seasonal_climatology.nc'
    if os.path.isfile(file) == True:
        fh = Dataset(file, mode='r') 
        wvcont_mm_seasonal_climatology = fh.variables['wvcont_mm_seasonal_climatology'][:]
        latitcrs = fh.variables['latitcrs'][:]
        longicrs = fh.variables['longicrs'][:]
        fh.close()  
        for xs in range(4):
            ax = axs[regions.index(region),xs]
            if xs == 0:
                ax.text(left,0,regions_names[regions.index(region)], size=fsize,horizontalalignment='right',\
                    verticalalignment='bottom',transform=ax.transAxes)#,rotation=90)
            # Plot the seasonal climatology for each region
            m = Basemap(projection='stere',lon_0=135,lat_0=-25.,\
                                llcrnrlat=-45,urcrnrlat=0,\
                                llcrnrlon=85,urcrnrlon=175,\
                                resolution='i',area_thresh=10000,\
                                ax = ax)
            m.drawcoastlines()
            x, y = m(longicrs,latitcrs) 
            h = ax.pcolormesh(x,y,wvcont_mm_seasonal_climatology[xs],norm=norm,cmap=cMap)

cax = plt.axes([0.01,-0.02,0.98,0.02] )
cbar = plt.colorbar(h,cax=cax,orientation='horizontal',boundaries=levs_mm,extend='max',\
spacing='proportional')
cbar.ax.tick_params(labelsize=fsize)
cbar.ax.set_xticklabels(levs_mm)#,rotation=45)      
cbar.set_label(label='Moisture contribution [mm]',size=fsize)#,rotation=90)

plt.tight_layout()
#fig.subplots_adjust(hspace=0) #wspace=0
fig.savefig(dir_out+'Fig3_mm.png',bbox_inches='tight') 
plt.close()