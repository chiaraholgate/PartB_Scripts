############################################################################

# This script creates input data files from NARCliM output for QIBT analysis over Australia. 
# Necessary atmospheric variables from daily files of 3 hourly output sourced directly from server.
# Precipitation and evaporation data extracted from monthly files of 1 hourly output and accumulated to daily files of 3 hourly to match atmospheric variables.

# Data sourced from NARCliM output, version R2_nudging. This version uses wrf output to constrain upper tropospheric variables on a scale of ~500km to ensure simulated large-scale processes are similar to observed/reanalysis.

###########################################################################

dir_in='/srv/ccrc/data33/z3481416/CCRC-WRF3.6.0.5-SEB/ERA-Interim/R2_nudging/out/'
dir_out='/srv/ccrc/data03/z3131380/PartB/NARCliM_postprocess/'

##### Extract precipitation from monthly files of hourly timestep, and save as daily files of 3-hourly timestep


# start with non-convective rainfall
# get month
# get number of days in month
# split file total timesteps into daily increments of 24
# convert cumulative rainfall to rainfall in the 3-hour period
# repeat for convective rainfall
# add non-convective and convective rainfall
# save in daily files of 3-hourly timestep.



#### Extract precipitation from 1-hourly monthly files 

# for each year in period of interest...
for year in {1979..2013};do
# for each month in the period of interest....


#### Process "RAINC": cumulative convective rainfall
# Extract relevant rainfall variables from hourly data to work with
cdo selname,RAINC,Times,XLAT,XLONG ${dir_in}wrfhrly_d01_2007-08-01_00:00:00 ${dir_out}RAINC.nc 
ncrename -v RAINC,RAIN ${dir_out}RAINC.nc
# Find number of days in this month by reading number of time steps
cdo ntime ${dir_out}RAINC.nc > ${dir_out}timesteps_in_month << EOF
EOF
days_in_month=$((`cat timesteps_in_month` / 24))
timesteps_in_month=$(($days_in_month*24))
# Convert cumulative rainfall to amount that fell in the timestep; then merge into 3-hourly daily file
# For each day in the month
for day in `seq 1 $days_in_month`;do
period=0
# For each 3-hourly period in the day
for t in `seq 1 3 24`;do
time1=$(($t + 24*$(($day - 1))))
cdo seltimestep,$time1 ${dir_out}RAINC.nc ${dir_out}temp_${t}_T1.nc 
time2=$(($time1 + 2))
cdo seltimestep,$time2 ${dir_out}RAINC.nc ${dir_out}temp_${t}_T2.nc 
period=$(($period + 1))
cdo sub ${dir_out}temp_${t}_T2.nc ${dir_out}temp_${t}_T1.nc ${dir_out}temp_${period}_P.nc 
rm ${dir_out}temp_*_T*.nc
done
if [ $day -lt 10 ] 
then 
today=0$day
else
today=$day
fi 
cdo mergetime ${dir_out}temp_*_P.nc ${dir_out}wrfhrly_d01_2007-08-${today}_00:00:00_RAINC.nc
rm ${dir_out}temp_*_P.nc
done

#### Process "RAINNC": cumulative non-convective rainfall
# Extract relevant rainfall variables from hourly data to work with
cdo selname,RAINC,Times,XLAT,XLONG ${dir_in}wrfhrly_d01_2007-08-01_00:00:00 ${dir_out}RAINNC.nc 
ncrename -v RAINNC,RAIN ${dir_out}RAINNC.nc
# Convert non-cumulative rainfall to amount that fell in the timestep; then merge into 3-hourly daily file
# For each day in the month
for day in `seq 1 $days_in_month`;do
period=0
# For each 3-hourly period in the day
for t in `seq 1 3 24`;do
time1=$(($t + 24*$(($day - 1))))
cdo seltimestep,$time1 ${dir_out}RAINNC.nc ${dir_out}temp_${t}_T1.nc
time2=$(($time1 + 2))
cdo seltimestep,$time2 ${dir_out}RAINNC.nc ${dir_out}temp_${t}_T2.nc 
period=$(($period + 1))
cdo sub ${dir_out}temp_${t}_T2.nc ${dir_out}temp_${t}_T1.nc ${dir_out}temp_${period}_P.nc 
rm ${dir_out}temp_*_T*.nc
done
if [ $day -lt 10 ] 
then 
today=0$day
else
today=$day
fi 
cdo mergetime ${dir_out}temp_*_P.nc ${dir_out}wrfhrly_d01_2007-08-${today}_00:00:00_RAINNC.nc
rm ${dir_out}temp_*_P.nc
done


#### Sum convective and non-convective rainfall into one file
cdo sum ${dir_out}wrfhrly_d01_2007-08-${today}_00:00:00_RAINC.nc ${dir_out}wrfhrly_d01_2007-08-${today}_00:00:00_RAINNC.nc




















