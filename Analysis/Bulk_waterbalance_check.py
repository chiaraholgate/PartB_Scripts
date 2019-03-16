# -*- coding: utf-8 -*-
"""
Created on Thu Feb 28 09:56:22 2019

This script calculates P-E=Q, where P and E are taken from the QIBT input (ie NARCliM output) and Q
is taken from (a) NARCliM output and (b) AWRA-L. The idea is to get an approximation of the water 
balance error of QIBT, although it's obviously going to be affected when Q comes from a different model.

P and E have been remapped to Q, when Q is from AWRA. See:
/home/z3131380/hdrive/PhD/PartB/Scripts/Processing/NARCliM_yearly_evap.sh
/home/z3131380/hdrive/PhD/PartB/Scripts/Processing/NARCliM_yearly_rain.sh
/home/z3131380/hdrive/PhD/PartB/Scripts/Processing/NARCliM_yearly_remap_to_AWRAqtot.sh

## NOTE!
- Surface runoff (SFROFF) from NARCliM seems to be an accumulative variable, even though the metadata
doesn't indicate this. I have processed is (NARCliM_yearly_runoff.sh) as accumulative.
- The yearly rainfall files loaded in this script have been made using NARCliM_yearly_rain.sh. There are some
cells where rainfall is very negative (e.g. -1000mm). This is because the convective rainfall (RAINC) variable 
from NARCliM, which is accumulative, sometimes drops between timesteps. So when you subtract one time
from a later time, the result is negative. Need to ask Jason about this! In this script I just force P<0 to be zero.
Fix!
- In addition to SFROFF, NARCliM/wrf also outputs several variables that may also need to be included in the
water balance to be correct. E.g. SH2O the soil liquid water, SMOIS soil moisture



@author: z3131380
"""
from netCDF4 import Dataset 
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import pandas

dir_P = '/srv/ccrc/data03/z3131380/PartB/NARCliM_postprocess/Yearly_rain/'
dir_E = '/srv/ccrc/data03/z3131380/PartB/NARCliM_postprocess/Yearly_evap/'
dir_out = '/srv/ccrc/data03/z3131380/PartB/NARCliM_postprocess/Bulk_waterbalance_check/'

#==============================================================================
# P and E from NARCliM, Q from NARCliM
#==============================================================================
dir_Q = '/srv/ccrc/data03/z3131380/PartB/NARCliM_postprocess/Yearly_runoff/'
dir_SM = '/srv/ccrc/data03/z3131380/PartB/NARCliM_postprocess/Yearly_SM/'
dir_mask = '/srv/ccrc/data03/z3131380/PartB/Masks/'
file = dir_mask+'NARCliM_AUS_land_sea_mask.nc'
fh = Dataset(file, mode='r') 
wsmask = fh.variables['wsmask'][:] 
fh.close() 

df = pandas.DataFrame()

for year in range(1979,2014):
    # Load Q
    file = dir_Q+'wrfout_d01_'+str(year)+'_runoff_mm.nc'
    fh = Dataset(file, mode='r') 
    Q = fh.variables['SFROFF'][:] 
    Q = Q[0,:,:] 
    fh.close() 
    
    # Load SM
    file = dir_SM+'wrfout_d01_'+str(year)+'_delta_SMOIS_mm.nc'
    fh = Dataset(file, mode='r') 
    SM = fh.variables['SMOIS_mm'][:] 
    SM = SM[0,0,:,:] 
    fh.close() 
          
    # Load P
    file = dir_P+'wrfhrly_d01_'+str(year)+'_rain_mm.nc'
    fh = Dataset(file, mode='r') 
    P = fh.variables['RAIN'][:] 
    P = P[0,:,:] 
    P[P<0]=0 ## NOTE!! ERROR IN RAIN FILES- SEE NOTES 11/7/17 & 1/3/19!!!
    fh.close() 
    
    # Load E
    file = dir_E+'wrfhrly_d01_'+str(year)+'_evap_mm.nc'
    fh = Dataset(file, mode='r') 
    E = fh.variables['evap'][:]  
    E = E[0,:,:]
    fh.close() 
    
    # Mask to Aus land
    Qm = np.ma.array(Q,mask=(wsmask==0))  
    SMm = np.ma.array(SM,mask=(wsmask==0)) 
    Pm = np.ma.array(P,mask=(wsmask==0))  
    Em = np.ma.array(E,mask=(wsmask==0))  
    
    # Calc sum
    Q_sum = np.nansum(Qm)
    SM_sum = np.nansum(SMm)
    P_sum = np.nansum(Pm)
    E_sum = np.nansum(Em)
    
    df = df.append(pandas.DataFrame({'Year':[year],'Q_sum':[Q_sum],\
        'SM_sum':[SM_sum],\
        'P_sum':[P_sum],\
        'E_sum':[E_sum],\
        'P-E':[P_sum-E_sum],\
        'Q+deltaSM':[Q_sum+SM_sum],\
        'P-E-Q-deltaSM':[P_sum-E_sum-Q_sum-SM_sum]},\
        index=[range(1979,2014).index(year)]))
outname = dir_out+'NARCliM_P_E_Q_SM_annual_waterbalance_1979-2013.csv'
df.to_csv(outname)    

#Plot
plt.plot(df['Year'],df['P_sum'],'b',label='P')
plt.plot(df['Year'],df['E_sum'],'g',label='E')
plt.plot(df['Year'],df['Q_sum'],'r',label='Q')
plt.plot(df['Year'],df['SM_sum'],'c',label='SM')
plt.legend()

>> THE RUNOFF MUST BE WRONG. IT IS TOO SMALL!




#==============================================================================
# P and E from NARCliM, remapped to AWRA; Q from AWRA
#==============================================================================
dir_Q = '/srv/ccrc/data03/z3131380/PartA/Data/AWRA/Runoff/'
# Load runoff estimates
file = dir_Q+'qtot.nc'
fh = Dataset(file, mode='r') 
Q = fh.variables['qtot'][:]  
lat = fh.variables['latitude'][:]  
lon = fh.variables['longitude'][:]  
fh.close() 

df = pandas.DataFrame()

for year in range(1979,2014):
    # Load P
    file = dir_P+'wrfhrly_d01_'+str(year)+'_rain_mm_remap_AWRAqtot.nc'
    fh = Dataset(file, mode='r') 
    P = fh.variables['RAIN'][:] 
    P = P[0,:,:] 
    P[P<0]=0 ## NOTE!! ERROR IN RAIN FILES- SEE NOTES 11/7/17 & 1/3/19!!!
    fh.close() 
    # Mask according to Q
    P_masked = np.ma.array(P,mask=Q.mask[0,:,:])
    P_sum = np.ma.sum(P_masked)
    
    # Load E
    file = dir_E+'wrfhrly_d01_'+str(year)+'_evap_mm_remap_AWRAqtot.nc'
    fh = Dataset(file, mode='r') 
    E = fh.variables['evap'][:]  
    E = E[0,:,:]
    fh.close() 
    # Mask according to Q
    E_masked = np.ma.array(E,mask=Q.mask[0,:,:])
    E_sum = np.ma.sum(E_masked)
    
    Q_sum = np.ma.sum(Q[68+range(1979,2014).index(year),:,:]) # Since 1979 is the 68th index
    
    df = df.append(pandas.DataFrame({'Year':[year],'Q_sum':[Q_sum],\
        'P_sum':[P_sum],\
        'E_sum':[E_sum],\
        'P-E':[P_sum-E_sum],\
        'P-E-Q':[P_sum-E_sum-Q_sum]},\
        index=[range(1979,2014).index(year)]))
outname = dir_out+'NARCliM_P_E_AWRA_Q_annual_waterbalance_1979-2013.csv'
df.to_csv(outname)             

#Plot
plt.plot(df['Year'],df['P_sum'],'b',df['Year'],df['E_sum'],'g',df['Year'],df['Q_sum'],'r')




## Calc P-E
#P_E = P_masked - E_masked
#residual = Q[68,:,:] - P_E
#P_E_tot_annual = np.ma.sum(P_masked) - np.ma.sum(E_masked)
#Q_tot_annual = np.ma.sum(Q[68,:,:])

# Count number of cells that have a low residual
total_cells = float(np.count_nonzero(residual))
low_residual = float(np.ma.count(residual[(residual<500) & (residual > 0)]))
low_pct = low_residual/total_cells*100

# Plot
plt.clf()    
fig = plt.figure()
m = Basemap(projection='cyl',llcrnrlat=lat.min(),urcrnrlat=lat.max(),\
             llcrnrlon=lon.min(),urcrnrlon=lon.max(),resolution='c')
x, y = m(lon,lat) 
m.drawcoastlines()
m.drawmeridians(np.arange(-180,180,10),labels=[0,0,1,0])
m.drawparallels(np.arange(-90,90,5),labels=[1,0,0,0])
h=m.pcolor(x,y,residual)
cbar = plt.colorbar(h)

