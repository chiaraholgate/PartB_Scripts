# -*- coding: utf-8 -*-
"""
Created on Thu Apr  4 09:24:15 2019

@author: z3131380
"""
import numpy as np
from netCDF4 import Dataset 
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import BoundaryNorm
import shapefile
from matplotlib.collections import PatchCollection
from matplotlib.patches import PathPatch
from matplotlib.patches import Polygon
import matplotlib.patheffects as PathEffects

#==============================================================================
# Definitions
#==============================================================================
dir_clim_zones = r'/srv/ccrc/data03/z3131380/PartA/Data/BOM/Climate_Zones/'
dir_shp = r'/srv/ccrc/data03/z3131380/PartB/Masks/Aus_Drainage_Divisions/geofabric/Divisions/Single/'
dir_out = r'/home/z3131380/hdrive/PhD/PartB/Reporting/Paper_1/Figures/'

# MAKE SURE THE NAMES MATCH!
regions = ['CarpentariaCoast','LakeEyreBasin','NorthEastCoast','NorthWesternPlateau','MurrayDarlingBasin',\
                'SouthWestCoast','PilbaraGascoyne','SouthAustralianGulf','SouthEastCoastNSW','SouthEastCoastVictoria',\
                'SouthWesternPlateau','TanamiTimorSeaCoast','Tasmania']

regions_names = ['Carpentaria\nCoast','Lake Eyre\nBasin','North East\nCoast',\
                'North Western\nPlateau','Murray\nDarling\nBasin','South\nWest\nCoast','Pilbara\nGascoyne',\
                'South\nAustralian\nGulf','South East\nCoast (NSW)','South East\nCoast (VIC)',\
                'South Western\nPlateau','Tanami-Timor\nSea Coast','Tasmania']   
       
region_centres = [[136.03,-17.33],\
                            [136.46,-26.26],\
                            [146.88,-22.47],\
                            [120.75,-22.5],\
                            [145.24,-31.61],\
                            [116.11,-33.5],\
                            [113.68,-25.62],\
                            [133.03,-35.61],\
                            [150.88,-33.26],\
                            [140.67,-40.33],\
                            [122.82,-30.54],\
                            [125.03,-17.62],\
                            [144.20,-42.5]]
ocean_names = ['Indian\nOcean','Arafura Sea','Coral Sea','South\nPacific\nOcean','Tasman Sea','Southern Ocean']
ocean_centres = [[107.18,-23.81],\
                            [134.38,-9.92],\
                            [149.70,-16.21],\
                            [157.58,-28.08],\
                            [152.81,-36.32],\
                            [125.22,-42.04]]              
              
fsize=12
fsize_coords = 12
#plt.rcParams["font.family"] = "Times New Roman"

#==============================================================================
# Plot                
#==============================================================================

plt.clf()    
fig = plt.figure(figsize=(16,11.7)) # A4, portrait=8.27,11.7, else 16,11.7     
ax = fig.add_subplot(1,1,1)
ax.axis('off')
m = Basemap(projection='stere',lon_0=135,lat_0=-25.,\
                    llcrnrlat=-45,urcrnrlat=0,\
                    llcrnrlon=85,urcrnrlon=175,\
                    resolution='i',area_thresh=10000)
m.readshapefile('/srv/ccrc/data03/z3131380/PartA/Data/BOM/Climate_Zones/kpngrp/kprngrp_vector_clip','kpngrp_vector_clip')         
temperate   = []; grassland = []; desert = []; subtropical = []; tropical = []; equatorial = []
for info, shape in zip(m.kpngrp_vector_clip_info, m.kpngrp_vector_clip):
    if info['DN'] < 12:
        temperate.append(Polygon(np.array(shape),True))
    elif np.logical_and(info['DN']>=21,info['DN']<30):
        desert.append(Polygon(np.array(shape),True))
    elif np.logical_and(info['DN']>30,info['DN']<35):
        subtropical.append(Polygon(np.array(shape),True))
    elif np.logical_and(info['DN']>=12,info['DN']<21):
        grassland.append(Polygon(np.array(shape),True))
    elif np.logical_and(info['DN']>34,info['DN']<41):
        tropical.append(Polygon(np.array(shape),True))
    elif info['DN']>=41:
        equatorial.append(Polygon(np.array(shape),True))
ax.add_collection(PatchCollection(temperate, facecolor= 'green', linewidths=0,edgecolor='none',zorder=2))
ax.add_collection(PatchCollection(desert, facecolor= 'navajowhite', edgecolor='none',linewidths=0, zorder=2))
ax.add_collection(PatchCollection(subtropical, facecolor= 'skyblue', edgecolor='none',linewidths=0, zorder=2))
ax.add_collection(PatchCollection(tropical, facecolor= 'steelblue', edgecolor='none',linewidths=0, zorder=2))
ax.add_collection(PatchCollection(equatorial, facecolor= 'b', edgecolor='none',linewidths=0, zorder=2))
ax.add_collection(PatchCollection(grassland, facecolor= 'greenyellow', edgecolor='none',linewidths=0, zorder=2))          
for region in regions:
    poly = m.readshapefile(dir_shp+region,'wsmask')
    text_x, text_y = m(region_centres[regions.index(region)][0], region_centres[regions.index(region)][1])
    t = plt.text(text_x,text_y,regions_names[regions.index(region)],fontsize=fsize,ha='left',\
        va='bottom',color='k',fontweight='bold') #m.wsmask_info[0]['Division'] 
    t.set_path_effects([PathEffects.withStroke(linewidth=5, foreground='w')])
for ocean in ocean_names:
    tx, ty = m(ocean_centres[ocean_names.index(ocean)][0], ocean_centres[ocean_names.index(ocean)][1])
    tt = plt.text(tx,ty,ocean,fontsize=fsize,ha='left',va='bottom',color='k',fontweight='bold',style='italic')

#plt.legend(handles=temperate[0])
from matplotlib.lines import Line2D 
customlines = [Line2D([0], [0], color='navajowhite',lw=7),\
                        Line2D([0], [0], color='greenyellow',lw=7), \
                        Line2D([0], [0], color='green',lw=7),\
                        Line2D([0], [0], color='skyblue',lw=7),\
                       Line2D([0], [0], color='steelblue',lw=7),\
                       Line2D([0], [0], color='b',lw=7)]
leg = ax.legend(customlines,['Desert','Grassland','Temperate','Subtropical','Tropical','Equatorial'],\
    loc='upper center', bbox_to_anchor=(0.44, 0.2),  ncol=6)
leg.get_frame().set_linewidth(0.0)
plt.tight_layout()
fig.savefig(dir_out+'Fig1.png',bbox_inches='tight') 
plt.close()    