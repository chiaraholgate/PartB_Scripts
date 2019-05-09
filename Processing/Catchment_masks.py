# -*- coding: utf-8 -*-
"""
Created on Wed Nov 28 17:25:34 2018

Steps to create catchment masks:
(1) Get Australia's drainage divisions
        * Map: http://www.bom.gov.au/water/geofabric/documents/BOM002_Map_Poster_A3_Web.pdf
        * Shp:/g/data/xc0/original/GIS/Australia/geofabric/HR_Regions_river_region_Divisions.shp
(2) Combine individual catchments into their corresponding basin
        * In QGIS: Vector>Geoprocessing Tools>Dissolve>Dissolve field = 'Division'
(3) Match projection to shp file previously provided by MDBA (doesn't seem to open in python correctly otw)        
        * In QGIS: vector>Data management tools>Define current projection
(4) Save basins as individual shp files. These basins still contain multiple catchments
        * In QGIS: vector>Data management tools>Split vector layer>ID field = 'Division'          
(5) Combine catchments of a basin into a single polygon
        * In QGIS: Open files from (4), then vector>Geometry tools>Multiparts to singlepart
        * Read in this shapefile below!

** Issue:
QGIS outputs the shapefile with multiple shapes (looks like 1, but len(all_shapes)>1!)        .
I'm currently using the all_shapes[i].points of the i with the most coord points. Combining
them gives an incorrect polygon.

**Note:
This script was originally written to produce masks that align to QIBT model output.
However I need to apply the catchment masks to SM data for analysis alongside the QIBT output.
So I'm making both types of mask!

@author: z3131380
"""
import os
import numpy as np
from netCDF4 import Dataset 
import shapefile #https://pypi.org/project/pyshp/
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon

# Mask type
#grid = 'WATERDYN_050deg_grid'
grid = 'QIBT_grid'


dir_in = r'/srv/ccrc/data03/z3131380/PartB/Masks/Aus_Drainage_Divisions/geofabric/Divisions/Single/'
dir_out = r'/srv/ccrc/data03/z3131380/PartB/Masks/Aus_Drainage_Divisions/geofabric/Divisions/netcdf/'+grid+'/test/'

# Get lons,lats from first model output file
if grid == 'QIBT_grid':    
    file = '/srv/ccrc/data19/z3131380/PartB/Output/Australia/100parcels/TS10min/exp01/bt.197901_31.nc'
    fh = Dataset(file, mode='r', format='NETCDF3_CLASSIC') 
    lats = fh.variables['latitcrs'][:] # latitudes of curvilinear grid
    lons = fh.variables['longicrs'][:] # longitude of curvilinear grid
    fh.close() 
if grid == 'WATERDYN_050deg_grid':
    file = '/srv/ccrc/data03/z3131380/PartA/Data/WATERDYN/050deg/netcdf/WRel1_daily_050deg_1979.nc'
    fh = Dataset(file, mode='r') 
    lats = fh.variables['lat'][:] 
    lons = fh.variables['lon'][:] 
    fh.close()    

shp_list = os.listdir(dir_in)
divisions = [i for i in shp_list if i.endswith('.shp')]

for d in divisions:
    basin = d[:-4]
        
    # Read in division shapefile 
    shp = shapefile.Reader(dir_in+d)    
    all_shapes = shp.shapes() # get all the polygons
    all_points = []
    if len(all_shapes)>1:
        a = []
        for i in range(len(all_shapes)):
            a.append(len(all_shapes[i].points))
        loc = np.where(a==np.max(a))[0][0]
        all_points = all_shapes[loc].points  
    lats_vect = []; lons_vect = []
    for i in range(len(all_points)):
        lats_vect.append(all_points[i][1])
        lons_vect.append(all_points[i][0])
    lons_lats_vect = np.column_stack((lons_vect,lats_vect))    # Reshape coordinates
    polygon = Polygon(lons_lats_vect) # create polygon


    # Create basin mask
    if grid=='QIBT_grid':
        nrows,ncols = np.shape(lats)
        wsmask = np.zeros_like(lats)
        for r in range(nrows):
            for c in range(ncols):
                point = Point(lons[r,c],lats[r,c]) # create point
                if polygon.contains(point): # OR point.within(polygon)
                    wsmask[r,c] = 1
    
        # Save to netcdf 
        ofile = dir_out+basin+'.nc'    
        with Dataset(ofile, 'w', format='NETCDF4_CLASSIC') as of: 
            of.createDimension('i_cross', nrows)
            of.createDimension('j_cross', ncols)   
        
            of.createVariable('latitcrs', 'f4', ('i_cross', 'j_cross'))
            of['latitcrs'].long_name = 'LATITUDE (SOUTH NEGATIVE)'
            of['latitcrs'].units = 'degrees'
            of['latitcrs'][:] = lats
            
            of.createVariable('longicrs', 'f4', ('i_cross', 'j_cross'))
            of['longicrs'].long_name = 'LONGITUDE (WEST NEGATIVE)'
            of['longicrs'].units = 'degrees'
            of['longicrs'][:] = lons
            
            of.createVariable('wsmask', 'f4', ('i_cross', 'j_cross'))
            of['wsmask'].long_name = 'Australian drainage division mask, QIBT_grid'
            of['wsmask'][:] = wsmask           
    
    if grid == 'WATERDYN_050deg_grid':
        nrows,ncols = len(lats),len(lons)
        wsmask = np.zeros([nrows,ncols])
        for r in range(nrows):
            for c in range(ncols):
                point = Point(lons[c],lats[r]) # create point
                if polygon.contains(point): # OR point.within(polygon)
                    wsmask[r,c] = 1
                    
        # Save to netcdf 
        ofile = dir_out+basin+'.nc'    
        with Dataset(ofile, 'w', format='NETCDF4_CLASSIC') as of: 
            of.createDimension('i', nrows)
            of.createDimension('j', ncols)   
        
            of.createVariable('lats', 'f4', ('i',))
            of['lats'].long_name = 'LATITUDE'
            of['lats'].units = 'degrees'
            of['lats'][:] = lats
            
            of.createVariable('lons', 'f4', ('j',))
            of['lons'].long_name = 'LONGITUDE'
            of['lons'].units = 'degrees'
            of['lons'][:] = lons
            
            of.createVariable('wsmask', 'f4', ('i', 'j'))
            of['wsmask'].long_name = 'Australian drainage division mask, WATERDYN_050deg_grid'
            of['wsmask'][:] = wsmask    

    print basin
# Attempted alternatives...
#shp = shapefile.Reader(dir_in+d)

#    all_shapes = shp.shapes() # get all the polygons
#    all_points = all_shapes[0].points  
#    lats_vect = []; lons_vect = []
#    for i in range(len(all_points)):
#        lats_vect.append(all_points[i][1])
#        lons_vect.append(all_points[i][0])
#    lons_lats_vect = np.column_stack((lons_vect,lats_vect))    # Reshape coordinates
#    polygon = Polygon(lons_lats_vect) # create polygon

#polygon = shp.shapes()
#shpfilePoints = [ shape.points for shape in polygon ]
#polygons = shpfilePoints
#for polygon in polygons:
#    poly = Polygon(polygon)

#feature = shp.shapeRecords()[0]
#first = feature.shape.__geo_interface__  
##from shapely.geometry import shape
##shp_geom = shape(first['geometry']) 
#poly_coords = first['coordinates']
#coords = []
#for i in range(len(poly_coords)):
#    coords.extend(poly_coords[i][0])
#coords = np.array(coords)    
#polygon = Polygon(coords) # create polygon
        
#shpfilePoints = [ shape.points for shape in polygon ]
#polygons = shpfilePoints
#for polygon in polygons:
#poly = Polygon(polygon)

