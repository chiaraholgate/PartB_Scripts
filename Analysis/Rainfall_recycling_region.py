# -*- coding: utf-8 -*-
"""
Created on Fri Nov 23 10:03:54 2018

This script computes the rainfall recycling ratio for a specified region based on back-trajectory model output.
It outputs the rainfall recycling and ocean vapour contribution per time as .csv files.

Calulations are made on an annual, seasonal and daily basis. Daily currently unused.

Checked and OK 17/1/19, 10/4/19.


@author: z3131380
"""
import sys
import numpy as np
from netCDF4 import Dataset 
import pandas
import os.path

#==============================================================================
# Definitions
#==============================================================================
dir_in = '/srv/ccrc/data19/z3131380/PartB/Output/Australia/100parcels/TS10min/exp04/Processed/Yearly/'
dir_out = '/srv/ccrc/data19/z3131380/PartB/Output/Australia/100parcels/TS10min/exp04/Processed/Rainfall_recycling/'

timeblock = sys.argv[1]

#regions = ['Australia','MurrayDarlingBasin','SouthWesternPlateau','TanamiTimorSeaCoast','Tasmania', \
#                 'SouthWestCoast', 'SouthEastCoastVictoria','CarpentariaCoast','LakeEyreBasin',\
#                 'NorthEastCoast','NorthWesternPlateau','PilbaraGascoyne','SouthAustralianGulf',\
#                 'SouthEastCoastNSW']
regions = ['Australia']
                
n_i,n_j = 134,205 # QIBT model dimensions
wvcont_mm_threshold = 5 #mm
cell_area = 50*50 # km2

#==============================================================================
# Dates
#==============================================================================
Start_date = '19790101' ; Start_year = Start_date[0:4] ; Start_month = Start_date[4:6]; Start_day = Start_date[6:8]
End_date = '20131231' ; End_year = End_date[0:4] ; End_month = End_date[4:6]; End_day = End_date[6:8]
Yearlist = pandas.date_range(Start_year,End_year,freq='AS').year # Note: 'A' gives end of year, 'AS' start of year. Must use full years in this script.
Seasonlist = []
for i in range(len(Yearlist)):
    Seasonlist.append(Yearlist[i]+0)
    Seasonlist.append(Yearlist[i]+0.25)
    Seasonlist.append(Yearlist[i]+0.5)
    Seasonlist.append(Yearlist[i]+0.75)

#==============================================================================
# Load region mask
#==============================================================================
#dir_mask = r'/srv/ccrc/data03/z3131380/PartB/Masks/Aus_Drainage_Divisions/geofabric/Divisions/netcdf/QIBT_grid/'
dir_mask = '/srv/ccrc/data03/z3131380/PartB/Masks/n_s_basin/'

# Get land and ocean outside Australia, as this is applicable regardless of region choice
fh = Dataset('/srv/ccrc/data03/z3131380/PartB/Masks/NARCliM_AUS_land_sea_mask.nc', mode='r') 
wsmask_Aus = fh.variables['wsmask'][4:138,4:209] # match shape to results output
fh.close()
fh = Dataset('/srv/ccrc/data03/z3131380/PartB/Masks/NARCLiM_land_sea_mask.nc', mode='r') 
wsmask_1 = fh.variables['topo'][4:138,4:209] 
fh.close()
wsmask_2 = np.ma.array(wsmask_1,mask=(wsmask_Aus==1))
land_outside_Aus  = np.ma.array(wsmask_2,mask=(wsmask_1==0))
ocean_outside_Aus = np.ma.array(wsmask_1,mask=(wsmask_1==1))
ocean_outside_Aus[~ocean_outside_Aus.mask] = 1
del  wsmask_1, wsmask_2

# Get region-specific land
for region in regions:    
    if region == 'Australia':
        fh = Dataset('/srv/ccrc/data03/z3131380/PartB/Masks/NARCliM_AUS_land_sea_mask.nc', mode='r') 
        wsmask_in = fh.variables['wsmask'][4:138,4:209] # match shape to results output
        fh.close()
        wsmask = np.ma.array(wsmask_in,mask=(wsmask_in==0))
        outregion_land_Aus = np.ma.zeros([n_i,n_j])
    else: 
        #fh = Dataset(dir_mask+region+'.nc', mode='r') 
        #wsmask = fh.variables['wsmask'][:]
        fh = Dataset(dir_mask+'nbasin_rotpole.nc',mode='r')
        wsmask = fh.variables['wsmask'][:-10,:-10]
        fh.close()       
        wsmask_1=np.ma.array(wsmask_Aus,mask=wsmask)
        outregion_land_Aus=np.ma.array(wsmask_1,mask=wsmask_Aus==0) 
        del wsmask_1           

   
    if timeblock == 'annual':
        #==============================================================================
        # Vapour contributions and rainfall recycling by year
        #==============================================================================
        df = pandas.DataFrame()
        RR = np.ma.zeros([len(Yearlist),n_i,n_j])*np.nan
        
        for year in Yearlist:
            # Load QIBT model results of water vapour contribution, processed by region and year
            fname = region+'_'+str(year)+'_wvcont.nc'
            file = dir_in+fname
            if os.path.isfile(file) == True:
                fh = Dataset(file, mode='r') 
                latitcrs = fh.variables['latitcrs'][:] # latitudes of curvilinear grid
                longicrs = fh.variables['longicrs'][:] # longitude of curvilinear grid
                wv_cont_sum_yearly_mm = fh.variables['wv_cont_sum_yearly_mm'][:] # water vapour contribution of domain cells to each cell where it rained
                pre_annual_total = fh.variables['pre_annual_total'][:] 
                fh.close()    
                
                # Cut up processed results by area
                wv_cont_region = np.ma.array(wv_cont_sum_yearly_mm,mask=(wsmask==0))
                wv_cont_outregion_land_Aus = np.ma.array(wv_cont_sum_yearly_mm,mask=(outregion_land_Aus==0))
                wv_cont_outregion_land_outside_Aus = np.ma.array(wv_cont_sum_yearly_mm,mask=(land_outside_Aus==0))
                wv_cont_outregion_ocean = np.ma.array(wv_cont_sum_yearly_mm,mask=(ocean_outside_Aus==0))
                wv_cont_outregion_ocean_area = np.count_nonzero(np.ma.array(wv_cont_outregion_ocean,\
                    mask=(wv_cont_outregion_ocean<wvcont_mm_threshold)))*cell_area #[km2]
                pre_region = np.ma.array(pre_annual_total,mask=(wsmask==0))
                P_total = np.nansum(pre_region) 
                
                # Calculate total water vapour contributions [%] by area
                outregion_ocean_cont = round(np.nansum(wv_cont_outregion_ocean)/np.nansum(pre_region),4)*100  
                outregion_land_Aus_cont = round(np.nansum(wv_cont_outregion_land_Aus)/np.nansum(pre_region),4)*100  
                outregion_land_outside_Aus_cont = round(np.nansum(wv_cont_outregion_land_outside_Aus)/np.nansum(pre_region),4)*100
                RR_region = round(np.nansum(wv_cont_region)/np.nansum(pre_region),4)*100  
                check_sum = np.nansum([RR_region,outregion_ocean_cont,outregion_land_Aus_cont, outregion_land_outside_Aus_cont])
                model_error = (np.nansum(P_total)-np.nansum(wv_cont_sum_yearly_mm))/np.nansum(P_total)*100
            
                df = df.append(pandas.DataFrame({'Year':year,'RR':RR_region,\
                'region_land_mm':np.nansum(wv_cont_region),\
                'outregion_land_Aus_pct':outregion_land_Aus_cont,\
                'outregion_land_Aus_mm':np.nansum(wv_cont_outregion_land_Aus),\
                'outregion_land_outside_Aus_pct':outregion_land_outside_Aus_cont,\
                'outregion_land_outside_Aus_mm':np.nansum(wv_cont_outregion_land_outside_Aus),\
                'outregion_ocean_pct':outregion_ocean_cont,\
                'outregion_ocean_mm':np.nansum(wv_cont_outregion_ocean),\
                'outregion_ocean_km':wv_cont_outregion_ocean_area,\
                'P_total':P_total,\
                'check_sum_pct':check_sum,\
                'model_error_pct':model_error},\
                index=[np.where(Yearlist==year)[0][0]]))
                
                # Replace any Nan's with zero, as pandas fills cells ignoring Nans (so ordering messed up unless you make them zero)
                df.fillna(0, inplace=True)
                
                # Calculate RR in the region
                rr = (wv_cont_region/P_total)*100
                RR[np.where(Yearlist==year)[0][0],:,:] =  np.ma.array(rr,mask=rr.mask)
        
        # Add extra results to data frame
        P_mean = df['P_total'].mean()
        df['P_anom'] = pandas.Series(df['P_total']*np.nan) # make a blank column
        df['P_anom'] = pandas.Series((df['P_total']-P_mean)/P_mean)
        # Calc the gradient of RR over different lag times
#        df['dL/dt_2'] = pandas.Series(df['P_total']*np.nan) # make a blank column
#        df['dL/dt_5'] = pandas.Series(df['P_total']*np.nan) # make a blank column
#        for row in np.arange(2,len(df)):
#            df.iloc[row,14] = (df.iloc[row,1] - df.iloc[row-2,1])/(df.iloc[row,2] - df.iloc[row-2,2])
#        for row in np.arange(5,len(df)):    
#            df.iloc[row,15] = (df.iloc[row,1] - df.iloc[row-5,1])/(df.iloc[row,2] - df.iloc[row-5,2]) 
        outname = dir_out+region+'_annual_rainfall_recycling_1979-2013.csv'
        df.to_csv(outname)      
        
        # Save RR grids to netcdf
        ofile = dir_out+region+'_1979-2013_annual_RR.nc'    
        with Dataset(ofile, 'w') as of: 
            of.createDimension('i_cross', n_i)
            of.createDimension('j_cross', n_j)    
            of.createDimension('time', len(Yearlist))
            
            of.createVariable('time', 'f4', ('time'))
            of['time'].long_name = 'Year'
            of['time'].units = 'years'
            of['time'][:] = Yearlist
        
            of.createVariable('latitcrs', 'f4', ('i_cross', 'j_cross'))
            of['latitcrs'].long_name = 'LATITUDE (SOUTH NEGATIVE)'
            of['latitcrs'].units = 'degrees'
            of['latitcrs'][:] = latitcrs
            
            of.createVariable('longicrs', 'f4', ('i_cross', 'j_cross'))
            of['longicrs'].long_name = 'LONGITUDE (WEST NEGATIVE)'
            of['longicrs'].units = 'degrees'
            of['longicrs'][:] = longicrs
            
            of.createVariable('RR', 'f4', ('time','i_cross', 'j_cross'))
            of['RR'].long_name = 'Ratio of grid cell water vapour contributed to rainfall anywhere in the region to total rainfall in the region'
            of['RR'].units = '%'
            of['RR'][:] = RR.filled(fill_value=np.nan) # convert to a non-masked array so netcdf doesn't screw up
            
    elif timeblock == 'seasonal':
        #==============================================================================
        # Rainfall recycling by season
        #==============================================================================
        df = pandas.DataFrame()
        RR = np.ma.zeros([len(Yearlist)*4,n_i,n_j])*np.nan
        
        for year in Yearlist:
            start_day = str(year)+'0101' ; end_day = str(year)+'1231'
            daylist = pandas.date_range(start_day,end_day,freq='d')  
            daylist_MAM = []; daylist_JJA = []; daylist_SON = []; daylist_DJF = []; daylist_AprNov = []
            for d in daylist:
                if d.month==3 or d.month==4 or d.month==5:
                    daylist_MAM.append(d)
                if d.month==6 or d.month==7 or d.month==8:
                    daylist_JJA.append(d)
                if d.month==9 or d.month==10 or d.month==11:
                    daylist_SON.append(d)
                if np.logical_and(d.month >=4,d.month<=11):
                    daylist_AprNov.append(d)
                    
                # DJF must include the Dec dates from the previous year        
                if year > 1979: 
                    if d.month==1 or d.month==2:
                        daylist_DJF.append(d)
                    if d.month==12:
                        daylist_DJF.append(d-pandas.DateOffset(years=1))
        
            daylist_DJF = pandas.to_datetime(daylist_DJF); daylist_DJF = daylist_DJF.sort_values()    
            daylist_MAM = pandas.to_datetime(daylist_MAM) 
            daylist_JJA = pandas.to_datetime(daylist_JJA)
            daylist_SON = pandas.to_datetime(daylist_SON)
            daylist_AprNov = pandas.to_datetime(daylist_AprNov)
            
            wv_cont_sum_daily_mm_DJF = np.zeros([len(daylist_DJF),n_i,n_j])*np.nan; pre_daily_DJF = np.zeros([len(daylist_DJF),n_i,n_j])*np.nan
            wv_cont_sum_daily_mm_MAM = np.zeros([len(daylist_MAM),n_i,n_j])*np.nan; pre_daily_MAM = np.zeros([len(daylist_MAM),n_i,n_j])*np.nan
            wv_cont_sum_daily_mm_JJA = np.zeros([len(daylist_JJA),n_i,n_j])*np.nan; pre_daily_JJA = np.zeros([len(daylist_JJA),n_i,n_j])*np.nan
            wv_cont_sum_daily_mm_SON = np.zeros([len(daylist_SON),n_i,n_j])*np.nan; pre_daily_SON = np.zeros([len(daylist_SON),n_i,n_j])*np.nan
            
        
            fname = region+'_'+str(year)+'_wvcont.nc'
            file = dir_in+fname
            if os.path.isfile(file) == True:
                fh = Dataset(file, mode='r') 
                latitcrs = fh.variables['latitcrs'][:] 
                longicrs = fh.variables['longicrs'][:] 
                wv_cont_sum_daily_mm = fh.variables['wv_cont_sum_daily_mm'][:]  
                pre = fh.variables['pre'][:] 
                fh.close()    
            #else: continue
            
            if year > 1979:
                # DJF must include the Dec data from the previous year   
                fname = region+'_'+str(year-1)+'_wvcont.nc'
                file = dir_in+fname
                if os.path.isfile(file) == True:
                    fh = Dataset(file, mode='r') 
                    wv_cont_sum_daily_mm_prevyr = fh.variables['wv_cont_sum_daily_mm'][:] 
                    pre_prevyr = fh.variables['pre'][:]  
                    fh.close()
                
            ## DJF ##
            daylist_prevyr = pandas.date_range(str(year-1)+'0101', str(year-1)+'1231',freq='d')
            for d in daylist_DJF:
                if d.month < 3:
                    idx_daylist = np.where(daylist==d)[0][0]
                    idx_daylist_DJF = np.where(daylist_DJF==d)[0][0]
                    wv_cont_sum_daily_mm_DJF[idx_daylist_DJF,:,:] = wv_cont_sum_daily_mm[idx_daylist,:,:]
                    pre_daily_DJF[idx_daylist_DJF,:,:] = pre[idx_daylist,:,:]
                elif d.month==12:
                    idx_daylist = np.where(daylist_prevyr==d)[0][0]
                    idx_daylist_DJF = np.where(daylist_DJF==d)[0][0]
                    wv_cont_sum_daily_mm_DJF[idx_daylist_DJF,:,:] = wv_cont_sum_daily_mm_prevyr[idx_daylist,:,:]
                    pre_daily_DJF[idx_daylist_DJF,:,:] = pre_prevyr[idx_daylist,:,:]
            wvcont_mm_DJF = np.nansum(wv_cont_sum_daily_mm_DJF,axis=0)
            pre_DJF = np.nansum(pre_daily_DJF, axis=0)
            
            # Cut up processed results by area
            wv_cont_region = np.ma.array(wvcont_mm_DJF,mask=(wsmask==0))
            wv_cont_outregion_land_Aus = np.ma.array(wvcont_mm_DJF,mask=(outregion_land_Aus==0))
            wv_cont_outregion_land_outside_Aus = np.ma.array(wvcont_mm_DJF,mask=(land_outside_Aus==0))
            wv_cont_outregion_ocean = np.ma.array(wvcont_mm_DJF,mask=(ocean_outside_Aus==0))
            wv_cont_outregion_ocean_area = np.count_nonzero(np.ma.array(wv_cont_outregion_ocean,\
                mask=(wv_cont_outregion_ocean<wvcont_mm_threshold)))*cell_area #[km2]
            pre_region = np.ma.array(pre_DJF,mask=(wsmask==0))
            P_total = np.nansum(pre_region)     
            
            # Calculate total water vapour contributions [%] by area
            outregion_ocean_cont = round(np.nansum(wv_cont_outregion_ocean)/np.nansum(pre_region),4)*100  
            outregion_land_Aus_cont = round(np.nansum(wv_cont_outregion_land_Aus)/np.nansum(pre_region),4)*100  
            outregion_land_outside_Aus_cont = round(np.nansum(wv_cont_outregion_land_outside_Aus)/np.nansum(pre_region),4)*100  
            RR_region = round(np.nansum(wv_cont_region)/np.nansum(pre_region),4)*100     
            check_sum = np.nansum([RR_region,outregion_ocean_cont,outregion_land_Aus_cont, outregion_land_outside_Aus_cont])
            model_error = (np.nansum(P_total)-np.nansum(wvcont_mm_DJF))/np.nansum(P_total)*100
            
            df = df.append(pandas.DataFrame({'Year':year,'Season':'DJF',\
            'RR':RR_region,\
            'region_land_mm':np.nansum(wv_cont_region),\
            'outregion_land_Aus_pct':outregion_land_Aus_cont,\
            'outregion_land_Aus_mm':np.nansum(wv_cont_outregion_land_Aus),\
            'outregion_land_outside_Aus_pct':outregion_land_outside_Aus_cont,\
            'outregion_land_outside_Aus_mm':np.nansum(wv_cont_outregion_land_outside_Aus),\
            'outregion_ocean_pct':outregion_ocean_cont,\
            'outregion_ocean_mm':np.nansum(wv_cont_outregion_ocean),\
            'outregion_ocean_km':wv_cont_outregion_ocean_area,\
            'P_total':P_total,\
            'check_sum_pct':check_sum,\
            'model_error_pct':model_error},\
            index=[4*np.where(Yearlist==year)[0][0]+0]))
            
            # Calculate RR in the region
            rr = (wv_cont_region/P_total)*100
            RR[4*np.where(Yearlist==year)[0][0]+0,:,:] =  np.ma.array(rr,mask=rr.mask)
                        
            ## MAM ##                
            for d in daylist_MAM:
                idx_daylist = np.where(daylist==d)[0][0]
                idx_daylist_MAM = np.where(daylist_MAM==d)[0][0]
                wv_cont_sum_daily_mm_MAM[idx_daylist_MAM,:,:] = wv_cont_sum_daily_mm[idx_daylist,:,:]
                pre_daily_MAM[idx_daylist_MAM,:,:] = pre[idx_daylist,:,:]
            wvcont_mm_MAM = np.nansum(wv_cont_sum_daily_mm_MAM,axis=0)
            pre_MAM = np.nansum(pre_daily_MAM, axis=0)
            
            # Cut up processed results by area
            wv_cont_region = np.ma.array(wvcont_mm_MAM,mask=(wsmask==0))
            wv_cont_outregion_land_Aus = np.ma.array(wvcont_mm_MAM,mask=(outregion_land_Aus==0))
            wv_cont_outregion_land_outside_Aus = np.ma.array(wvcont_mm_MAM,mask=(land_outside_Aus==0))
            wv_cont_outregion_ocean = np.ma.array(wvcont_mm_MAM,mask=(ocean_outside_Aus==0))
            wv_cont_outregion_ocean_area = np.count_nonzero(np.ma.array(wv_cont_outregion_ocean,\
                mask=(wv_cont_outregion_ocean<wvcont_mm_threshold)))*cell_area #[km2]
            pre_region = np.ma.array(pre_MAM,mask=(wsmask==0))
            P_total = np.nansum(pre_region)     
            
            # Calculate total water vapour contributions [%] by area
            outregion_ocean_cont = round(np.nansum(wv_cont_outregion_ocean)/np.nansum(pre_region),4)*100  
            outregion_land_Aus_cont = round(np.nansum(wv_cont_outregion_land_Aus)/np.nansum(pre_region),4)*100  
            outregion_land_outside_Aus_cont = round(np.nansum(wv_cont_outregion_land_outside_Aus)/np.nansum(pre_region),4)*100  
            RR_region = round(np.nansum(wv_cont_region)/np.nansum(pre_region),4)*100     
            check_sum = np.nansum([RR_region,outregion_ocean_cont,outregion_land_Aus_cont, outregion_land_outside_Aus_cont])
            model_error = (np.nansum(P_total)-np.nansum(wvcont_mm_MAM))/np.nansum(P_total)*100
            
            df = df.append(pandas.DataFrame({'Year':year+0.25,'Season':'MAM',\
            'RR':RR_region,\
            'region_land_mm':np.nansum(wv_cont_region),\
            'outregion_land_Aus_pct':outregion_land_Aus_cont,\
            'outregion_land_Aus_mm':np.nansum(wv_cont_outregion_land_Aus),\
            'outregion_land_outside_Aus_pct':outregion_land_outside_Aus_cont,\
            'outregion_land_outside_Aus_mm':np.nansum(wv_cont_outregion_land_outside_Aus),\
            'outregion_ocean_pct':outregion_ocean_cont,\
            'outregion_ocean_mm':np.nansum(wv_cont_outregion_ocean),\
            'outregion_ocean_km':wv_cont_outregion_ocean_area,\
            'P_total':P_total,\
            'check_sum_pct':check_sum,\
            'model_error_pct':model_error},\
            index=[4*np.where(Yearlist==year)[0][0]+1]))
            
            # Calculate RR in the region
            rr = (wv_cont_region/P_total)*100
            RR[4*np.where(Yearlist==year)[0][0]+1,:,:] =  np.ma.array(rr,mask=rr.mask)
        
            ## JJA ##    
            for d in daylist_JJA:
                idx_daylist = np.where(daylist==d)[0][0]
                idx_daylist_JJA = np.where(daylist_JJA==d)[0][0]
                wv_cont_sum_daily_mm_JJA[idx_daylist_JJA,:,:] = wv_cont_sum_daily_mm[idx_daylist,:,:]
                pre_daily_JJA[idx_daylist_JJA,:,:] = pre[idx_daylist,:,:]
            wvcont_mm_JJA = np.nansum(wv_cont_sum_daily_mm_JJA,axis=0)
            pre_JJA = np.nansum(pre_daily_JJA, axis=0)
            
            # Cut up processed results by area
            wv_cont_region = np.ma.array(wvcont_mm_JJA,mask=(wsmask==0))
            wv_cont_outregion_land_Aus = np.ma.array(wvcont_mm_JJA,mask=(outregion_land_Aus==0))
            wv_cont_outregion_land_outside_Aus = np.ma.array(wvcont_mm_JJA,mask=(land_outside_Aus==0))
            wv_cont_outregion_ocean = np.ma.array(wvcont_mm_JJA,mask=(ocean_outside_Aus==0))
            wv_cont_outregion_ocean_area = np.count_nonzero(np.ma.array(wv_cont_outregion_ocean,\
                mask=(wv_cont_outregion_ocean<wvcont_mm_threshold)))*cell_area #[km2]
            pre_region = np.ma.array(pre_JJA,mask=(wsmask==0))
            P_total = np.nansum(pre_region)     
            
            # Calculate total water vapour contributions [%] by area
            outregion_ocean_cont = round(np.nansum(wv_cont_outregion_ocean)/np.nansum(pre_region),4)*100  
            outregion_land_Aus_cont = round(np.nansum(wv_cont_outregion_land_Aus)/np.nansum(pre_region),4)*100  
            outregion_land_outside_Aus_cont = round(np.nansum(wv_cont_outregion_land_outside_Aus)/np.nansum(pre_region),4)*100  
            RR_region = round(np.nansum(wv_cont_region)/np.nansum(pre_region),4)*100     
            check_sum = np.nansum([RR_region,outregion_ocean_cont,outregion_land_Aus_cont, outregion_land_outside_Aus_cont])
            model_error = (np.nansum(P_total)-np.nansum(wvcont_mm_JJA))/np.nansum(P_total)*100
            
            df = df.append(pandas.DataFrame({'Year':year+0.5,'Season':'JJA',\
            'RR':RR_region,\
            'region_land_mm':np.nansum(wv_cont_region),\
            'outregion_land_Aus_pct':outregion_land_Aus_cont,\
            'outregion_land_Aus_mm':np.nansum(wv_cont_outregion_land_Aus),\
            'outregion_land_outside_Aus_pct':outregion_land_outside_Aus_cont,\
            'outregion_land_outside_Aus_mm':np.nansum(wv_cont_outregion_land_outside_Aus),\
            'outregion_ocean_pct':outregion_ocean_cont,\
            'outregion_ocean_mm':np.nansum(wv_cont_outregion_ocean),\
            'outregion_ocean_km':wv_cont_outregion_ocean_area,\
            'P_total':P_total,\
            'check_sum_pct':check_sum,\
            'model_error_pct':model_error},\
            index=[4*np.where(Yearlist==year)[0][0]+2]))
            
            # Calculate RR in the region
            rr = (wv_cont_region/P_total)*100
            RR[4*np.where(Yearlist==year)[0][0]+2,:,:] =  np.ma.array(rr,mask=rr.mask)
                
            ## SON ##    
            for d in daylist_SON:
                idx_daylist = np.where(daylist==d)[0][0]
                idx_daylist_SON = np.where(daylist_SON==d)[0][0]
                wv_cont_sum_daily_mm_SON[idx_daylist_SON,:,:] = wv_cont_sum_daily_mm[idx_daylist,:,:]
                pre_daily_SON[idx_daylist_SON,:,:] = pre[idx_daylist,:,:]
            wvcont_mm_SON = np.nansum(wv_cont_sum_daily_mm_SON,axis=0)
            pre_SON = np.nansum(pre_daily_SON, axis=0)
            
            # Cut up processed results by area
            wv_cont_region = np.ma.array(wvcont_mm_SON,mask=(wsmask==0))
            wv_cont_outregion_land_Aus = np.ma.array(wvcont_mm_SON,mask=(outregion_land_Aus==0))
            wv_cont_outregion_land_outside_Aus = np.ma.array(wvcont_mm_SON,mask=(land_outside_Aus==0))
            wv_cont_outregion_ocean = np.ma.array(wvcont_mm_SON,mask=(ocean_outside_Aus==0))
            wv_cont_outregion_ocean_area = np.count_nonzero(np.ma.array(wv_cont_outregion_ocean,\
                mask=(wv_cont_outregion_ocean<wvcont_mm_threshold)))*cell_area #[km2]
            pre_region = np.ma.array(pre_SON,mask=(wsmask==0))
            P_total = np.nansum(pre_region)     
            
            # Calculate total water vapour contributions [%] by area
            outregion_ocean_cont = round(np.nansum(wv_cont_outregion_ocean)/np.nansum(pre_region),4)*100  
            outregion_land_Aus_cont = round(np.nansum(wv_cont_outregion_land_Aus)/np.nansum(pre_region),4)*100  
            outregion_land_outside_Aus_cont = round(np.nansum(wv_cont_outregion_land_outside_Aus)/np.nansum(pre_region),4)*100  
            RR_region = round(np.nansum(wv_cont_region)/np.nansum(pre_region),4)*100     
            check_sum = np.nansum([RR_region,outregion_ocean_cont,outregion_land_Aus_cont, outregion_land_outside_Aus_cont])
            model_error = (np.nansum(P_total)-np.nansum(wvcont_mm_SON))/np.nansum(P_total)*100
            
            df = df.append(pandas.DataFrame({'Year':year+0.75,'Season':'SON',\
            'RR':RR_region,\
            'region_land_mm':np.nansum(wv_cont_region),\
            'outregion_land_Aus_pct':outregion_land_Aus_cont,\
            'outregion_land_Aus_mm':np.nansum(wv_cont_outregion_land_Aus),\
            'outregion_land_outside_Aus_pct':outregion_land_outside_Aus_cont,\
            'outregion_land_outside_Aus_mm':np.nansum(wv_cont_outregion_land_outside_Aus),\
            'outregion_ocean_pct':outregion_ocean_cont,\
            'outregion_ocean_mm':np.nansum(wv_cont_outregion_ocean),\
            'outregion_ocean_km':wv_cont_outregion_ocean_area,\
            'P_total':P_total,\
            'check_sum_pct':check_sum,\
            'model_error_pct':model_error},\
            index=[4*np.where(Yearlist==year)[0][0]+3]))
            
            # Calculate RR in the region
            rr = (wv_cont_region/P_total)*100
            RR[4*np.where(Yearlist==year)[0][0]+3,:,:] =  np.ma.array(rr,mask=rr.mask)
            
               
        P_mean_DJF = df.groupby(by=(df['Season']=='DJF'))['P_total'].mean()[1]
        P_mean_MAM = df.groupby(by=(df['Season']=='MAM'))['P_total'].mean()[1]
        P_mean_JJA = df.groupby(by=(df['Season']=='JJA'))['P_total'].mean()[1]
        P_mean_SON = df.groupby(by=(df['Season']=='SON'))['P_total'].mean()[1]
        df['P_anom'] = pandas.Series(df['P_total']*np.nan) # make a blank column
        for i in range(len(df)):
            if df.ix[i,'Season'] == 'DJF':
                df.ix[i,'P_anom'] = (df.iloc[i,0]-P_mean_DJF)/P_mean_DJF
            if df.ix[i,'Season'] == 'MAM':
                df.ix[i,'P_anom'] = (df.iloc[i,0]-P_mean_MAM)/P_mean_MAM
            if df.ix[i,'Season'] == 'JJA':
                df.ix[i,'P_anom'] = (df.iloc[i,0]-P_mean_JJA)/P_mean_JJA
            if df.ix[i,'Season'] == 'SON':
                df.ix[i,'P_anom'] = (df.iloc[i,0]-P_mean_SON)/P_mean_SON     
        # Change first row values from zero to nan so they don't plot            
        df.iloc[0,0]=np.nan; df.iloc[0,6]=np.nan; df.iloc[0,8]=np.nan; df.iloc[0,10]=np.nan; df.iloc[0,-1]=np.nan
        
        outname = dir_out+region+'_seasonal_rainfall_recycling_1979-2013.csv'
        df.to_csv(outname)     
        
        # Save RR grids to netcdf
        ofile = dir_out+region+'_1979-2013_seasonal_RR.nc'    
        with Dataset(ofile, 'w') as of: 
            of.createDimension('i_cross', n_i)
            of.createDimension('j_cross', n_j)    
            of.createDimension('time', len(Seasonlist))
            
            of.createVariable('time', 'f4', ('time'))
            of['time'].long_name = 'Year+0=DJF,Year+0.25=MAM,Year+0.5=JJA,Year+0.75=SON'
            of['time'].units = 'years'
            of['time'][:] = Seasonlist
        
            of.createVariable('latitcrs', 'f4', ('i_cross', 'j_cross'))
            of['latitcrs'].long_name = 'LATITUDE (SOUTH NEGATIVE)'
            of['latitcrs'].units = 'degrees'
            of['latitcrs'][:] = latitcrs
            
            of.createVariable('longicrs', 'f4', ('i_cross', 'j_cross'))
            of['longicrs'].long_name = 'LONGITUDE (WEST NEGATIVE)'
            of['longicrs'].units = 'degrees'
            of['longicrs'][:] = longicrs
            
            of.createVariable('RR', 'f4', ('time','i_cross', 'j_cross'))
            of['RR'].long_name = 'Ratio of grid cell water vapour contributed to rainfall anywhere in the region to total rainfall in the region'
            of['RR'].units = '%'
            of['RR'][:] = RR.filled(fill_value=np.nan) # convert to a non-masked array so netcdf doesn't screw up

    elif timeblock == 'daily':
        #==============================================================================
        # Rainfall recycling by day
        #==============================================================================
       
        for year in Yearlist:
            start_day = str(year)+'0101' ; end_day = str(year)+'1231'
            Daylist = pandas.date_range(start_day,end_day,freq='d')  
            fname = region+'_'+str(year)+'_wvcont.nc'
            file = dir_in+fname
            if os.path.isfile(file) == True:
                fh = Dataset(file, mode='r') 
                latitcrs = fh.variables['latitcrs'][:] # latitudes of curvilinear grid
                longicrs = fh.variables['longicrs'][:] # longitude of curvilinear grid
                wv_cont_sum_daily_mm = fh.variables['wv_cont_sum_daily_mm'][:] # water vapour contribution of domain cells to each cell where it rained
                pre = fh.variables['pre'][:] 
                fh.close()    
                
                #Daylist = pandas.date_range(start='1/1/'+str(year),end='31/12/'+str(year),freq='d')
                
                df = pandas.DataFrame() # Make a new data frame for every year
                RR = np.ma.zeros([len(Daylist),n_i,n_j])*np.nan
                
                for day in Daylist:#range(len(wv_cont_sum_daily_mm)):
                    date = day.strftime('%Y%m%d')
                    dd = np.where(Daylist==day)[0][0]
        
                    # Cut up processed results by area
                    wv_cont_region = np.ma.array(wv_cont_sum_daily_mm[dd,:,:],mask=(wsmask==0))
                    wv_cont_outregion_land_Aus = np.ma.array(wv_cont_sum_daily_mm[dd,:,:],mask=(outregion_land_Aus==0))
                    wv_cont_outregion_land_outside_Aus = np.ma.array(wv_cont_sum_daily_mm[dd,:,:],mask=(land_outside_Aus==0))
                    wv_cont_outregion_ocean = np.ma.array(wv_cont_sum_daily_mm[dd,:,:],mask=(ocean_outside_Aus==0))
                    wv_cont_outregion_ocean_area = np.count_nonzero(np.ma.array(wv_cont_outregion_ocean,\
                        mask=(wv_cont_outregion_ocean<wvcont_mm_threshold)))*cell_area #[km2]
                    pre_region = np.ma.array(pre[dd,:,:],mask=(wsmask==0))
                    P_total = np.nansum(pre_region) 
                    
                    # Calculate total water vapour contributions [%] by area
                    outregion_ocean_cont = round(np.nansum(wv_cont_outregion_ocean)/np.nansum(pre_region),4)*100  
                    outregion_land_Aus_cont = round(np.nansum(wv_cont_outregion_land_Aus)/np.nansum(pre_region),4)*100  
                    outregion_land_outside_Aus_cont = round(np.nansum(wv_cont_outregion_land_outside_Aus)/np.nansum(pre_region),4)*100
                    RR_region = round(np.nansum(wv_cont_region)/np.nansum(pre_region),4)*100  
                    check_sum = np.nansum([RR_region,outregion_ocean_cont,outregion_land_Aus_cont, outregion_land_outside_Aus_cont])
                    if P_total > 0:
                        model_error = (P_total-np.nansum(wv_cont_sum_daily_mm[dd,:,:]))/P_total*100
                    else: model_error = np.nan
                
                    df = df.append(pandas.DataFrame({'Day':[day],'RR':[RR_region],\
                    'region_land_mm':[np.nansum(wv_cont_region)],\
                    'outregion_land_Aus_pct':[outregion_land_Aus_cont],\
                    'outregion_land_Aus_mm':[np.nansum(wv_cont_outregion_land_Aus)],\
                    'outregion_land_outside_Aus_pct':[outregion_land_outside_Aus_cont],\
                    'outregion_land_outside_Aus_mm':[np.nansum(wv_cont_outregion_land_outside_Aus)],\
                    'outregion_ocean_pct':[outregion_ocean_cont],\
                    'outregion_ocean_mm':[np.nansum(wv_cont_outregion_ocean)],\
                    'outregion_ocean_km':[wv_cont_outregion_ocean_area],\
                    'P_total':[P_total],\
                    'check_sum_pct':[check_sum],\
                    'model_error_pct':[model_error]}))#,\
                    #index=date))
                    
                    # Calculate RR in the region
                    rr = (wv_cont_region/P_total)*100
                    RR[dd,:,:] =  np.ma.array(rr,mask=rr.mask)

            outname = dir_out+region+'_daily_rainfall_recycling_'+str(year)+'.csv'
            df.to_csv(outname)      
            
            # Save RR grids to netcdf
            ofile = dir_out+region+'_daily_RR_'+str(year)+'.nc'    
            with Dataset(ofile, 'w') as of: 
                of.createDimension('i_cross', n_i)
                of.createDimension('j_cross', n_j)    
                of.createDimension('time', len(Daylist))
                
                of.createVariable('time', 'f4', ('time'))
                of['time'].long_name = 'Day'
                of['time'].units = 'days'
                of['time'][:] = Daylist.strftime('%Y%m%d')
            
                of.createVariable('latitcrs', 'f4', ('i_cross', 'j_cross'))
                of['latitcrs'].long_name = 'LATITUDE (SOUTH NEGATIVE)'
                of['latitcrs'].units = 'degrees'
                of['latitcrs'][:] = latitcrs
                
                of.createVariable('longicrs', 'f4', ('i_cross', 'j_cross'))
                of['longicrs'].long_name = 'LONGITUDE (WEST NEGATIVE)'
                of['longicrs'].units = 'degrees'
                of['longicrs'][:] = longicrs
                
                of.createVariable('RR', 'f4', ('time','i_cross', 'j_cross'))
                of['RR'].long_name = 'Ratio of grid cell water vapour contributed to rainfall anywhere in the region to total rainfall in the region'
                of['RR'].units = '%'
                of['RR'][:] = RR.filled(fill_value=np.nan) # convert to a non-masked array so netcdf doesn't screw up
                
            print year