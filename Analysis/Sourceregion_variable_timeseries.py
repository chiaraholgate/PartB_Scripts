# -*- coding: utf-8 -*-
"""
Created on Thu Mar 28 15:25:58 2019

This script calculates daily and monthly timeseries of selected NARCliM variables over each source area. 
The script does calculations for ocean polygons separate from land regions. Could combine if desired.
Output variables include:
    (1) Evap [mm]
    (2) TPW [mm]
    (3) CAPE [J/kg]
    (4) CIN [J/kg]

** Daily evaporation is calculated by first taking the daily mean latent heat flux over the region, then 
converting it to daily evaporation in mm. If desired, could replace daily mean of LH with a weighted
average, weighted by cell area.

** Daily TPW is calculated by subtracting the first TPW value for the day from the last. 

** CAPE and CIN are taken as means.


@author: z3131380
"""
import numpy as np
from netCDF4 import Dataset 
import pandas as pd

#==============================================================================
# Definitions
#==============================================================================
dir_E = r'/srv/ccrc/data03/z3131380/PartB/NARCliM_postprocess/'
dir_TPW = r'/srv/ccrc/data03/z3131380/PartB/NARCliM_postprocess/TPW/'
dir_CAPE = r'/srv/ccrc/data03/z3131380/PartB/NARCliM_postprocess/CAPE/'
dir_oceanpolys = r'/srv/ccrc/data19/z3131380/PartB/Output/Australia/100parcels/TS10min/exp01/Processed/Climatology/Ocean_polygons/'
dir_out = r'/srv/ccrc/data19/z3131380/PartB/Output/Australia/100parcels/TS10min/exp01/Processed/Sourceregion_variable_timeseries/'

Lv = 2.25E6   #latent heat of vaporization of water (Jkg-1)

#==============================================================================
# Dates
#==============================================================================
Start_date = '19790201' ; Start_year = Start_date[0:4] ; Start_month = Start_date[4:6]; Start_day = Start_date[6:8]
End_date = '20131230' ; End_year = End_date[0:4] ; End_month = End_date[4:6]; End_day = End_date[6:8]
Daylist = pd.date_range(Start_date,End_date,freq='d') # Note: 'A' gives end of year, 'AS' start of year. Must use full years in this script.
      
#==============================================================================
# Integrate over ocean source areas
#==============================================================================
fh = Dataset(dir_oceanpolys+'Ocean_source_regions.nc', mode='r') 
wsmask_Pacific = fh.variables['wsmask_Pacific'][:]
wsmask_Indian = fh.variables['wsmask_Indian'][:]
wsmask_Southern = fh.variables['wsmask_Southern'][:]
wsmask_Arafura = fh.variables['wsmask_Arafura'][:]
fh.close()      

# Make masks same time length shape as input data (which has 8 timesteps)
wsmask_Pacific = np.repeat(wsmask_Pacific[np.newaxis,:,:],8,axis=0)
wsmask_Indian = np.repeat(wsmask_Indian[np.newaxis,:,:],8,axis=0)
wsmask_Southern = np.repeat(wsmask_Southern[np.newaxis,:,:],8,axis=0)
wsmask_Arafura = np.repeat(wsmask_Arafura[np.newaxis,:,:],8,axis=0)

df_E_ocean = pd.DataFrame()
df_TPW_ocean = pd.DataFrame()
df_CAPE_ocean = pd.DataFrame()
df_CIN_ocean = pd.DataFrame()

for d in Daylist:
    yyyy = str(d.year)
    mm = '%02u' % d.month
    dd = '%02u' % d.day
    
    fh = Dataset(dir_E+'wrfhrly_d01_'+yyyy+'-'+mm+'-'+dd+'_00:00:00_LH.nc', mode='r') 
    LH = fh.variables['LH'][:,4:138,4:209] # match shape to QIBT results output
    fh.close()
    
    fh = Dataset(dir_TPW+'TPW_'+yyyy+'_'+mm+'_'+dd+'.nc', mode='r') 
    TPW = fh.variables['TPW'][:,4:138,4:209] # match shape to QIBT results output
    fh.close()
    
    fh = Dataset(dir_CAPE+'CAPE_'+yyyy+'_'+mm+'_'+dd+'.nc', mode='r') 
    CAPE = fh.variables['CAPE'][:,4:138,4:209] # match shape to QIBT results output
    CIN = fh.variables['CIN'][:,4:138,4:209] # match shape to QIBT results output
    fh.close()
    
    # Mask array
    LH_Pacific = np.ma.array(LH,mask=(wsmask_Pacific.mask))
    LH_Indian = np.ma.array(LH,mask=(wsmask_Indian.mask))
    LH_Southern = np.ma.array(LH,mask=(wsmask_Southern.mask))
    LH_Arafura = np.ma.array(LH,mask=(wsmask_Arafura.mask))
    
    TPW_Pacific = np.ma.array(TPW,mask=(wsmask_Pacific.mask))
    TPW_Indian = np.ma.array(TPW,mask=(wsmask_Indian.mask))
    TPW_Southern = np.ma.array(TPW,mask=(wsmask_Southern.mask))
    TPW_Arafura = np.ma.array(TPW,mask=(wsmask_Arafura.mask))
    
    CAPE_Pacific = np.ma.array(CAPE,mask=(wsmask_Pacific.mask))
    CAPE_Indian = np.ma.array(CAPE,mask=(wsmask_Indian.mask))
    CAPE_Southern = np.ma.array(CAPE,mask=(wsmask_Southern.mask))
    CAPE_Arafura = np.ma.array(CAPE,mask=(wsmask_Arafura.mask))
    
    CIN_Pacific = np.ma.array(CIN,mask=(wsmask_Pacific.mask))
    CIN_Indian = np.ma.array(CIN,mask=(wsmask_Indian.mask))
    CIN_Southern = np.ma.array(CIN,mask=(wsmask_Southern.mask))
    CIN_Arafura = np.ma.array(CIN,mask=(wsmask_Arafura.mask))
    
    # Convert average daily from latent heat [W/m2] to total daily evap [mm]
    evap_Pacific_daytotal = LH_Pacific.mean()*(60*60*24)/Lv
    evap_Indian_daytotal = LH_Indian.mean()*(60*60*24)/Lv
    evap_Southern_daytotal = LH_Southern.mean()*(60*60*24)/Lv
    evap_Arafura_daytotal = LH_Arafura.mean()*(60*60*24)/Lv
    
    # Find change in TPW storage over the day [mm]
    TPW_Pacific_daytotal = np.sum(TPW_Pacific[7,:,:] - TPW_Pacific[0,:,:])
    TPW_Indian_daytotal = np.sum(TPW_Indian[7,:,:] - TPW_Indian[0,:,:])
    TPW_Southern_daytotal = np.sum(TPW_Southern[7,:,:] - TPW_Southern[0,:,:])
    TPW_Arafura_daytotal = np.sum(TPW_Arafura[7,:,:] - TPW_Arafura[0,:,:])
    
    # Find mean daily CAPE 
    CAPE_Pacific_daymean = CAPE_Pacific.mean()
    CAPE_Indian_daymean = CAPE_Indian.mean()
    CAPE_Southern_daymean = CAPE_Southern.mean()
    CAPE_Arafura_daymean = CAPE_Arafura.mean()
    
    # Find mean daily CIN 
    CIN_Pacific_daymean = CIN_Pacific.mean()
    CIN_Indian_daymean = CIN_Indian.mean()
    CIN_Southern_daymean = CIN_Southern.mean()
    CIN_Arafura_daymean = CIN_Arafura.mean()
         
    df_E_ocean = df_E_ocean.append(pd.DataFrame({'Date':d,'Month':d.month,'Year':d.year,\
        'Pacific':[evap_Pacific_daytotal],\
        'Indian':[evap_Indian_daytotal],\
        'Southern':[evap_Southern_daytotal],\
        'Arafura':[evap_Arafura_daytotal]}))
    
    df_TPW_ocean = df_TPW_ocean.append(pd.DataFrame({'Date':d,'Month':d.month,'Year':d.year,\
        'Pacific':[TPW_Pacific_daytotal],\
        'Indian':[TPW_Indian_daytotal],\
        'Southern':[TPW_Southern_daytotal],\
        'Arafura':[TPW_Arafura_daytotal]}))
    
    df_CAPE_ocean = df_CAPE_ocean.append(pd.DataFrame({'Date':d,'Month':d.month,'Year':d.year,\
        'Pacific':[CAPE_Pacific_daymean],\
        'Indian':[CAPE_Indian_daymean],\
        'Southern':[CAPE_Southern_daymean],\
        'Arafura':[CAPE_Arafura_daymean]}))
    
    df_CIN_ocean = df_CIN_ocean.append(pd.DataFrame({'Date':d,'Month':d.month,'Year':d.year,\
        'Pacific':[CIN_Pacific_daymean],\
        'Indian':[CIN_Indian_daymean],\
        'Southern':[CIN_Southern_daymean],\
        'Arafura':[CIN_Arafura_daymean]}))
        
df_E_ocean.index = df_E_ocean['Date']
df_E_ocean.to_csv(dir_out+'Daily_E_mm_oceanregions_1979-2013.csv')  
df_E_ocean_grouped = df_E_ocean.groupby(['Year','Month']).sum()
df_E_ocean_grouped.to_csv(dir_out+'Monthly_E_mm_oceanregions_1979-2013.csv')  
#df.plot()
#df_grouped.plot()

df_TPW_ocean.index = df_TPW_ocean['Date']
df_TPW_ocean.to_csv(dir_out+'Daily_TPW_mm_oceanregions_1979-2013.csv')  
df_TPW_ocean_grouped = df_TPW_ocean.groupby(['Year','Month']).sum()
df_TPW_ocean_grouped.to_csv(dir_out+'Monthly_TPW_mm_oceanregions_1979-2013.csv')  

df_CAPE_ocean.index = df_CAPE_ocean['Date']
df_CAPE_ocean.to_csv(dir_out+'Daily_CAPE_Jkg-1_oceanregions_1979-2013.csv')  
df_CAPE_ocean_grouped = df_CAPE_ocean.groupby(['Year','Month']).sum()
df_CAPE_ocean_grouped.to_csv(dir_out+'Monthly_CAPE_Jkg-1_oceanregions_1979-2013.csv')  

df_CIN_ocean.index = df_CIN_ocean['Date']
df_CIN_ocean.to_csv(dir_out+'Daily_CIN_Jkg-1_oceanregions_1979-2013.csv')  
df_CIN_ocean_grouped = df_CIN_ocean.groupby(['Year','Month']).sum()
df_CIN_ocean_grouped.to_csv(dir_out+'Monthly_CIN_Jkg-1_oceanregions_1979-2013.csv')  

#==============================================================================
# Integrate evap and TPW over land sub-basins
#==============================================================================
regions = ['CarpentariaCoast','LakeEyreBasin','MurrayDarlingBasin','NorthEastCoast','NorthWesternPlateau',\
                'PilbaraGascoyne','SouthAustralianGulf','SouthEastCoastNSW','SouthEastCoastVictoria',\
                'SouthWestCoast','SouthWesternPlateau','TanamiTimorSeaCoast']

dir_mask = r'/srv/ccrc/data03/z3131380/PartB/Masks/Aus_Drainage_Divisions/geofabric/Divisions/netcdf/QIBT_grid/'
df_E_regions = pd.DataFrame(); df_E_regions['Date'] = Daylist
df_T_regions = pd.DataFrame(); df_T_regions['Date'] = Daylist
df_CAPE_regions = pd.DataFrame(); df_CAPE_regions['Date'] = Daylist
df_CIN_regions = pd.DataFrame(); df_CIN_regions['Date'] = Daylist

for r in regions:
    E = []; T = []; cp = []; cn = []
    
    fh = Dataset(dir_mask+r+'.nc',mode='r')
    wsmask = fh.variables['wsmask'][:]
    wsmask_region = np.repeat(wsmask[np.newaxis,:,:],8,axis=0)
    fh.close; del wsmask
    
    for d in Daylist:
        yyyy = str(d.year)
        mm = '%02u' % d.month
        dd = '%02u' % d.day
        
        fh = Dataset(dir_E+'wrfhrly_d01_'+yyyy+'-'+mm+'-'+dd+'_00:00:00_LH.nc', mode='r') 
        LH = fh.variables['LH'][:,4:138,4:209] # match shape to QIBT results output
        fh.close()
        
        fh = Dataset(dir_TPW+'TPW_'+yyyy+'_'+mm+'_'+dd+'.nc', mode='r') 
        TPW = fh.variables['TPW'][:,4:138,4:209] # match shape to QIBT results output
        fh.close()
        
        fh = Dataset(dir_CAPE+'CAPE_'+yyyy+'_'+mm+'_'+dd+'.nc', mode='r') 
        CAPE = fh.variables['CAPE'][:,4:138,4:209] # match shape to QIBT results output
        CIN = fh.variables['CIN'][:,4:138,4:209] # match shape to QIBT results output
        fh.close()
              
        # Mask array
        LH_region = np.ma.array(LH,mask=(wsmask_region==0))        
        TPW_region = np.ma.array(TPW,mask=(wsmask_region==0))
        CAPE_region = np.ma.array(CAPE,mask=(wsmask_region==0))
        CIN_region = np.ma.array(CIN,mask=(wsmask_region==0))
        
        # Convert average daily from latent heat [W/m2] to total daily evap [mm]
        evap_region_daytotal = LH_region.mean()*(60*60*24)/Lv
        
        # Find change in TPW storage over the day [mm]
        TPW_region_daytotal = np.sum(TPW_region[7,:,:] - TPW_region[0,:,:])
        
        # Find mean daily CAPE 
        CAPE_region_daymean = CAPE_region.mean()
        
        # Find mean daily CIN 
        CIN_region_daymean = CIN_region.mean()
        
        # Append to variable
        E.append(evap_region_daytotal) #[mm]
        T.append(TPW_region_daytotal) #[mm]
        cp.append(CAPE_region_daymean) #[J/kg]
        cn.append(CIN_region_daymean) #[J/kg]

        
    df_E_regions[r] = E
    df_T_regions[r] = T
    df_CAPE_regions[r] = cp
    df_CIN_regions[r] = cn
    
df_E_regions.index = df_E_regions['Date']
df_E_regions['Year'] = Daylist.year
df_E_regions['Month'] = Daylist.month
df_E_regions.to_csv(dir_out+'Daily_E_mm_landregions_1979-2013.csv')    
df_E_regions_grouped = df_E_regions.groupby(['Year','Month']).sum()
df_E_regions_grouped.to_csv(dir_out+'Monthly_E_mm_landregions_1979-2013.csv')    
#df_E_regions.plot()
#df_E_regions_grouped.plot()

df_T_regions.index = df_T_regions['Date']
df_T_regions['Year'] = Daylist.year
df_T_regions['Month'] = Daylist.month  
df_T_regions.to_csv(dir_out+'Daily_TPW_mm_landregions_1979-2013.csv')    
df_T_regions_grouped = df_T_regions.groupby(['Year','Month']).sum()
df_T_regions_grouped.to_csv(dir_out+'Monthly_TPW_mm_landregions_1979-2013.csv')    

df_CAPE_regions.index = df_CAPE_regions['Date']
df_CAPE_regions['Year'] = Daylist.year
df_CAPE_regions['Month'] = Daylist.month  
df_CAPE_regions.to_csv(dir_out+'Daily_CAPE_Jkg-1_landregions_1979-2013.csv')    
df_CAPE_regions_grouped = df_CAPE_regions.groupby(['Year','Month']).sum()
df_CAPE_regions_grouped.to_csv(dir_out+'Monthly_CAPE_Jkg-1_landregions_1979-2013.csv')   

df_CIN_regions.index = df_CIN_regions['Date']
df_CIN_regions['Year'] = Daylist.year
df_CIN_regions['Month'] = Daylist.month  
df_CIN_regions.to_csv(dir_out+'Daily_CIN_Jkg-1_landregions_1979-2013.csv')    
df_CIN_regions_grouped = df_CIN_regions.groupby(['Year','Month']).sum()
df_CIN_regions_grouped.to_csv(dir_out+'Monthly_CIN_Jkg-1_landregions_1979-2013.csv')   
