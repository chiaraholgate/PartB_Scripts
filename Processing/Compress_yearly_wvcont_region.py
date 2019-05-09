# -*- coding: utf-8 -*-
"""
Created on Thu Oct 18 11:16:20 2018

This script creates yearly grids of water vapour contribution from individual cells to rainfall
wihtin a region of interest, defined by a netcdf mask. 
Contribution is given in mm and percent.

Checked and OK 28/10/18.

@author: z3131380
"""
import numpy as np
from netCDF4 import Dataset 
import pandas
import os.path

#==============================================================================
# Definitions
#==============================================================================
dir_in = '/srv/ccrc/data19/z3131380/PartB/Output/Australia/100parcels/TS10min/exp01/'
dir_mask = r'/srv/ccrc/data03/z3131380/PartB/Masks/Aus_Drainage_Divisions/geofabric/Divisions/netcdf/QIBT_grid/'
dir_out = '/srv/ccrc/data19/z3131380/PartB/Output/Australia/100parcels/TS10min/exp01/Processed/Yearly/'
domain = 'Australia'
regions = ['SouthWesternPlateau','TanamiTimorSeaCoast','Tasmania']
                # already done: 'MurrayDarlingBasin', 'SouthWestCoast', SouthEastCoastVIC, 'Australia'
#               'CarpentariaCoast','LakeEyreBasin','NorthEastCoast','NorthWesternPlateau',\
#                'PilbaraGascoyne','SouthAustralianGulf','SouthEastCoastNSW',\

n_i,n_j = 134,205 # QIBT model dimensions               
#==============================================================================
# Dates               
#==============================================================================
Start_date = '19790101' ; Start_year = Start_date[0:4] ; Start_month = Start_date[4:6]; Start_day = Start_date[6:8]
End_date = '20131231' ; End_year = End_date[0:4] ; End_month = End_date[4:6]; End_day = End_date[6:8]
Yearlist = pandas.date_range(Start_year,End_year,freq='AS').year # Note: 'A' gives end of year, 'AS' start of year. Must use full years in this script.

#==============================================================================
# Compress raw results by region
#==============================================================================
for region in regions:
    fh = Dataset(dir_mask+region+'.nc',mode='r')
    wsmask = fh.variables['wsmask'][:]
    fh.close
    

#if region == 'Australia':
#    fh = Dataset(dir_mask+'NARCliM_AUS_land_sea_mask.nc', mode='r') 
#    wsmask = fh.variables['wsmask'][4:138,4:209] # match shape to results output
#    fh.close()

# Load shapefile of region to use as a mask
#if region=='Australia':
#    shp = shapefile.Reader('/srv/ccrc/data03/z3131380/PartB/Masks/Aus_Coastline/australia/cstauscd_r.shp') 
#elif region=='MDB':
#    shp = shapefile.Reader('/srv/ccrc/data03/z3131380/PartB/Masks/mdb_boundary/mdb_boundary.shp') 
#all_shapes = shp.shapes() # get all the polygons
#all_points = all_shapes[0].points  
#
#lats_vect = []; lons_vect = []
#for i in range(len(all_points)):
#    lats_vect.append(all_points[i][1])
#    lons_vect.append(all_points[i][0])
#lons_lats_vect = np.column_stack((lons_vect,lats_vect))    # Reshape coordinates
#polygon = Polygon(lons_lats_vect) # create polygon


    for year in Yearlist:
        daylist = pandas.date_range(str(year)+'0101',str(year)+'1231',freq='d') 
        wv_cont_sum_daily_mm = np.zeros([len(daylist),n_i,n_j])*np.nan
        wv_cont_apbl_sum_daily_mm = np.zeros([len(daylist),n_i,n_j])*np.nan
        pre_grid_daily = np.zeros([len(daylist),n_i,n_j])*np.nan
        pre_total_daily = np.zeros([len(daylist)])
        wv_cont_sum_daily_pct = np.zeros([len(daylist),n_i,n_j])*np.nan
        wv_cont_apbl_sum_daily_pct = np.zeros([len(daylist),n_i,n_j])*np.nan
        for d in daylist:
            idx = np.where(daylist==d)[0][0]
            mm = '%02u' % d.month
            dd = d.day
            fname = 'bt.'+str(year)+mm+'_'+str(dd)+'.nc'
            file = dir_in+fname
            if os.path.isfile(file) == True:
                fh = Dataset(file, mode='r', format='NETCDF3_CLASSIC') 
                if fh.dimensions['gridcell_wvc'].size > 0: 
                    x_loc = fh.variables['x_loc'][:] #x index of precipitation in no. cells from [0,0]
                    y_loc = fh.variables['y_loc'][:] #y index of precipitation in no. cells from [0,0]
                    latitcrs = fh.variables['latitcrs'][:] # latitudes of curvilinear grid
                    longicrs = fh.variables['longicrs'][:] # longitude of curvilinear grid
                    wv_cont = fh.variables['wv_cont'][:] # water vapour contribution of domain cells to each cell where it rained
                    wv_cont_apbl = fh.variables['wv_cont_apbl'][:] # as above, but above the PBL
                    pre = fh.variables['pre'][:] # rainfall in each cell on that day
                    #day = fh.variables['day'][:]    
                fh.close()    
                latitcrs_wsmask = np.ma.array(latitcrs,mask=(wsmask==0))
                longicrs_wsmask = np.ma.array(longicrs,mask=(wsmask==0))
                
                # Find where it rained that day
                rain_lat = []; rain_lon = []; 
                rain_coords = []; rain_loc = []
                for i in range(len(pre)):
                    rain_yloc,rain_xloc =  y_loc[i],x_loc[i]
                    rain_lat = latitcrs[rain_yloc,rain_xloc]
                    rain_lon = longicrs[rain_yloc,rain_xloc] 
                    rain_coords.append([rain_lon,rain_lat])#rain_coords.append(['%.2f' % rain_lon,'%.2f' % rain_lat])     
                    rain_loc.append([rain_yloc,rain_xloc])
                
                # identify those rain locations within the region of interest
                rain_coords_region = []; rain_loc_region = []; rain_idx = []
                for r  in range(len(rain_coords)):
                    if np.logical_and(rain_coords[r][0] in longicrs_wsmask,rain_coords[r][1] in latitcrs_wsmask):
                        rain_coords_region.append(rain_coords[r])
                        rain_loc_region.append(rain_loc[r])
                        rain_idx.append(r)
    #            if region=='Australia':
    #                rain_coords_region = rain_coords; rain_loc_region = rain_loc; rain_idx = range(len(pre))
    #            else: 
    #                rain_coords_region = []; rain_loc_region = []; rain_idx = []
    #                for r  in range(len(rain_coords)):
    #                    point = Point(rain_coords[r]) # create point
    #                    if polygon.contains(point): # OR point.within(polygon)
    #                        rain_coords_region.append(rain_coords[r])
    #                        rain_loc_region.append(rain_loc[r])
    #                        rain_idx.append(r)
    
                # Find contribution (in mm) of each cell's vapour to each rain cell, and the total contribution 
                # of each cell to all rain in region each day
                wv_mm = np.zeros([len(rain_idx),n_i,n_j])*np.nan
                wv_apbl_mm = np.zeros([len(rain_idx),n_i,n_j])*np.nan            
                for k in range(len(rain_idx)):
                    wv_mm[k,:,:] = wv_cont[rain_idx[k],:,:]*pre[rain_idx[k]] # Only consider contributions to region rainfall
                    wv_apbl_mm[k,:,:] = wv_cont_apbl[rain_idx[k],:,:]*pre[rain_idx[k]]
                wv_cont_sum_daily_mm[idx,:,:] = np.nansum(wv_mm,axis=0) #Total contribution of each cell to rainfall in region each day
                wv_cont_apbl_sum_daily_mm[idx,:,:] = np.nansum(wv_apbl_mm,axis=0)
    
                # Place rain cells on the x,y grid for plotting purposes/checking
                for m in range(len(rain_loc_region)):
                    pre_grid_daily[idx,rain_loc_region[m][0],rain_loc_region[m][1]] = pre[rain_idx[m]]
                    
                # Total rain in the region each day
                pre_total_daily[idx] = np.nansum(pre_grid_daily[idx,:,:])
                
                # Find the proportion of vapour each cell contributed to total rainfall in the region each day
                wv_cont_sum_daily_pct[idx,:,:] = wv_cont_sum_daily_mm[idx,:,:] / pre_total_daily[idx]
                wv_cont_apbl_sum_daily_pct[idx,:,:] = wv_cont_apbl_sum_daily_mm[idx,:,:] / pre_total_daily[idx]
                
            else: continue # if file does not exist, go to next day.
        
        pre_total_yearly = np.nansum(pre_total_daily)   
        pre_yearly = np.nansum(pre_grid_daily,axis=0)
        wv_cont_sum_yearly_mm = np.nansum(wv_cont_sum_daily_mm,axis=0)
        wv_cont_sum_yearly_pct = wv_cont_sum_yearly_mm / pre_total_yearly
        wv_cont_apbl_sum_yearly_mm = np.nansum(wv_cont_apbl_sum_daily_mm,axis=0)
        wv_cont_apbl_sum_yearly_pct = wv_cont_apbl_sum_yearly_mm / pre_total_yearly
        
                        
        # Save to netcdf (save in same NETCDF3 as model output??)
        ofile = dir_out+region+'_'+str(year)+'_wvcont.nc'    
        with Dataset(ofile, 'w', format='NETCDF4_CLASSIC') as of: 
            of.createDimension('i_cross', n_i)
            of.createDimension('j_cross', n_j)    
            of.createDimension('time', len(daylist))
            
            of.createVariable('time', 'f4', ('time'))
            of['time'].long_name = 'days since '+str(daylist[0])[0:10]
            of['time'].units = 'days'
            of['time'][:] = np.arange(1,len(daylist)+1)
    
            of.createVariable('latitcrs', 'f4', ('i_cross', 'j_cross'))
            of['latitcrs'].long_name = 'LATITUDE (SOUTH NEGATIVE)'
            of['latitcrs'].units = 'degrees'
            of['latitcrs'][:] = latitcrs
            
            of.createVariable('longicrs', 'f4', ('i_cross', 'j_cross'))
            of['longicrs'].long_name = 'LONGITUDE (WEST NEGATIVE)'
            of['longicrs'].units = 'degrees'
            of['longicrs'][:] = longicrs
            
            of.createVariable('pre', 'f4', ('time','i_cross', 'j_cross'))
            of['pre'].long_name = 'precipitation'
            of['pre'].units = 'mm'
            of['pre'][:] = pre_grid_daily
            of['pre'].fill_value = -9999.0
            
            of.createVariable('pre_annual_total', 'f4', ('i_cross', 'j_cross'))
            of['pre_annual_total'].long_name = 'annual precipitation'
            of['pre_annual_total'].units = 'mm'
            of['pre_annual_total'][:] = pre_yearly
            of['pre_annual_total'].fill_value = -9999.0
            
            of.createVariable('wv_cont_sum_daily_mm', 'f4', ('time','i_cross', 'j_cross'))
            of['wv_cont_sum_daily_mm'].long_name = 'mm contribution of each cell to all rain in the region'
            of['wv_cont_sum_daily_mm'].units = 'mm'
            of['wv_cont_sum_daily_mm'][:] = wv_cont_sum_daily_mm
            of['wv_cont_sum_daily_mm'].fill_value = -9999.0
            
            of.createVariable('wv_cont_sum_yearly_mm', 'f4', ('i_cross', 'j_cross'))
            of['wv_cont_sum_yearly_mm'].long_name = 'mm contribution of each cell to all rain in the region in year'
            of['wv_cont_sum_yearly_mm'].units = 'mm'
            of['wv_cont_sum_yearly_mm'][:] = wv_cont_sum_yearly_mm
            of['wv_cont_sum_yearly_mm'].fill_value = -9999.0
            
            of.createVariable('wv_cont_sum_daily_pct', 'f4', ('time','i_cross', 'j_cross'))
            of['wv_cont_sum_daily_pct'].long_name = '% contribution of each cell to all rain in the region'
            of['wv_cont_sum_daily_pct'].units = '%'
            of['wv_cont_sum_daily_pct'][:] = wv_cont_sum_daily_pct
            of['wv_cont_sum_daily_pct'].fill_value = -9999.0
            
            of.createVariable('wv_cont_sum_yearly_pct', 'f4', ('i_cross', 'j_cross'))
            of['wv_cont_sum_yearly_pct'].long_name = '% contribution of each cell to all rain in the region in year'
            of['wv_cont_sum_yearly_pct'].units = '%'
            of['wv_cont_sum_yearly_pct'][:] = wv_cont_sum_yearly_pct
            of['wv_cont_sum_yearly_pct'].fill_value = -9999.0
            
    #        of.createVariable('wv_cont_apbl_sum_yearly_mm', 'f4', ('i_cross', 'j_cross'))
    #        of['wv_cont_apbl_sum_yearly_mm'].long_name = 'mm contribution above PBL of each cell to all rain in the region'
    #        of['wv_cont_apbl_sum_yearly_mm'].units = 'mm'
    #        of['wv_cont_apbl_sum_yearly_mm'][:] = wv_cont_apbl_sum_yearly_mm
    #        of['wv_cont_apbl_sum_yearly_mm'].fill_value = -9999.0
    #        
    #        of.createVariable('wv_cont_apbl_sum_yearly_pct', 'f4', ('i_cross', 'j_cross'))
    #        of['wv_cont_apbl_sum_yearly_pct'].long_name = '% contribution above PBL of each cell to all rain in the region'
    #        of['wv_cont_apbl_sum_yearly_pct'].units = '%'
    #        of['wv_cont_apbl_sum_yearly_pct'][:] = wv_cont_apbl_sum_yearly_pct
    #        of['wv_cont_apbl_sum_yearly_pct'].fill_value = -9999.0
            
            print region,year
                
        
    
