# Try to save .nc as classic- see cdo 1.2.1 to match wrfout
#check evap mean/sum
# change variable description in RAIN file
# Try saving netcdf (RAIN) to a local place and rerunning.



############################################################################

# This script creates input data files from NARCliM output for QIBT analysis over Australia. 
# Necessary atmospheric variables from daily files of 3 hourly output sourced directly from server.
# Precipitation and evaporation data extracted from monthly files of 1 hourly output and accumulated to daily files of 3 hourly to match atmospheric variables.

# Data sourced from NARCliM output, version R2_nudging. This version uses wrf output to constrain upper tropospheric variables on a scale of ~500km to ensure simulated large-scale processes are similar to observed/reanalysis.

dir_in='/srv/ccrc/data33/z3481416/CCRC-WRF3.6.0.5-SEB/ERA-Interim/R2_nudging/out/'
dir_out='/srv/ccrc/data03/z3131380/PartB/NARCliM_postprocess/'

###########################################################################
# RAINFALL

# Take difference between first and third timesteps to convert hourly cumulative rainfall to a total over 3-hours

for year in {1979..2013};do
for month in {01..12};do
# "RAINC": cumulative convective rainfall
cdo selname,RAINC,Times,XLAT,XLONG ${dir_in}wrfhrly_d01_${year}-${month}-01_00:00:00 ${dir_out}RAINC.nc 
ncrename -v RAINC,RAIN ${dir_out}RAINC.nc
# "RAINNC": cumulative non-convective rainfall
cdo selname,RAINNC,Times,XLAT,XLONG ${dir_in}wrfhrly_d01_${year}-${month}-01_00:00:00 ${dir_out}RAINNC.nc 
ncrename -v RAINNC,RAIN ${dir_out}RAINNC.nc
# Find number of days in this month by reading number of time steps
cdo ntime ${dir_out}RAINC.nc > ${dir_out}timesteps_in_month << EOF
EOF
days_in_month=$((`cat ${dir_out}timesteps_in_month` / 24))
timesteps_in_month=`cat ${dir_out}timesteps_in_month` #$(($days_in_month*24))
# Convert cumulative rainfall to amount that fell in the timestep; then merge into 3-hourly daily file
# For each day in the month
for day in `seq 1 $days_in_month`;do
period=0
# For each 3-hourly period in the day
for t in `seq 1 3 24`;do
time1=$(($t + 24*$(($day - 1))))
cdo seltimestep,$time1 ${dir_out}RAINC.nc ${dir_out}temp_${t}_T1_RAINC.nc 
cdo seltimestep,$time1 ${dir_out}RAINNC.nc ${dir_out}temp_${t}_T1_RAINNC.nc 
time2=$(($time1 + 2))
cdo seltimestep,$time2 ${dir_out}RAINC.nc ${dir_out}temp_${t}_T2_RAINC.nc 
cdo seltimestep,$time2 ${dir_out}RAINNC.nc ${dir_out}temp_${t}_T2_RAINNC.nc 
period=$(($period + 1))
cdo sub ${dir_out}temp_${t}_T2_RAINC.nc ${dir_out}temp_${t}_T1_RAINC.nc ${dir_out}temp_${period}_P_RAINC.nc 
cdo sub ${dir_out}temp_${t}_T2_RAINNC.nc ${dir_out}temp_${t}_T1_RAINNC.nc ${dir_out}temp_${period}_P_RAINNC.nc 
rm ${dir_out}temp_*_T*.nc
done
if [ $day -lt 10 ] 
then 
today=0$day
else
today=$day
fi 
cdo mergetime ${dir_out}temp_*_P_RAINC.nc ${dir_out}wrfhrly_d01_${year}-${month}-${today}_00:00:00_RAINC.nc
cdo mergetime ${dir_out}temp_*_P_RAINNC.nc ${dir_out}wrfhrly_d01_${year}-${month}-${today}_00:00:00_RAINNC.nc
rm ${dir_out}temp_*_P*.nc
# Sum convective and non-convective rainfall into one file
cdo -f nc4c add ${dir_out}wrfhrly_d01_${year}-${month}-${today}_00:00:00_RAINC.nc ${dir_out}wrfhrly_d01_${year}-${month}-${today}_00:00:00_RAINNC.nc ${dir_out}wrfhrly_d01_${year}-${month}-${today}_00:00:00_RAIN.nc
ncatted -O -a description,,o,c,"Sum of total cumulus and non-cumulus precipitation" ${dir_out}wrfhrly_d01_${year}-${month}-${today}_00:00:00_RAIN.nc
rm ${dir_out}wrfhrly_d01_${year}-${month}-${today}_00:00:00_RAINC.nc
rm ${dir_out}wrfhrly_d01_${year}-${month}-${today}_00:00:00_RAINNC.nc
done
rm ${dir_out}RAINC.nc 
rm ${dir_out}RAINNC.nc 
done
done

###########################################################################


###########################################################################
# EVAPORATION

# Take average over the 3 one-hourly timesteps
# [LH given in W/m2 or J/s/m2...so you are averaging the J/s over the three hours.]

for year in {1979..2013};do
for month in {01..12};do
# "LH": latent heat flux at the surface
cdo selname,LH,Times,XLAT,XLONG ${dir_in}wrfhrly_d01_${year}-${month}-01_00:00:00 ${dir_out}LH.nc 
# Find number of days in this month by reading number of time steps
cdo ntime ${dir_out}LH.nc > ${dir_out}timesteps_in_month << EOF
EOF
days_in_month=$((`cat ${dir_out}timesteps_in_month` / 24))
timesteps_in_month=`cat ${dir_out}timesteps_in_month` #$(($days_in_month*24))

# For each day in the month
for day in `seq 1 $days_in_month`;do
period=0
# For each 3-hourly period in the day
for t in `seq 1 3 24`;do
time1=$(($t + 24*$(($day - 1))))
time2=$(($time1 + 1))
time3=$(($time1 + 2))
cdo seltimestep,$time1 ${dir_out}LH.nc ${dir_out}temp_${t}_T1_LH.nc 
cdo seltimestep,$time2 ${dir_out}LH.nc ${dir_out}temp_${t}_T2_LH.nc 
cdo seltimestep,$time3 ${dir_out}LH.nc ${dir_out}temp_${t}_T3_LH.nc 
period=$(($period + 1))
#cdo ensmean ${dir_out}temp_${t}_T1_LH.nc ${dir_out}temp_${t}_T2_LH.nc ${dir_out}temp_${t}_T3_LH.nc ${dir_out}temp_${period}_P_LH.nc 
cdo enssum ${dir_out}temp_${t}_T1_LH.nc ${dir_out}temp_${t}_T2_LH.nc ${dir_out}temp_${t}_T3_LH.nc ${dir_out}temp_${period}_P_LH.nc
rm ${dir_out}temp_*_T*.nc
done
if [ $day -lt 10 ] 
then 
today=0$day
else
today=$day
fi 
cdo mergetime ${dir_out}temp_*_P_LH.nc ${dir_out}wrfhrly_d01_${year}-${month}-${today}_00:00:00_LH.nc
rm ${dir_out}temp_*_P*.nc
done
rm ${dir_out}LH.nc 
done
done

###########################################################################


###########################################################################
# SURFACE PRESSURE REFERENCE

# Take average over the 3 one-hourly timesteps

# for year in {1979..2013};do
# for month in {01..12};do
# # "PSFC": pressure at the surface
# cdo selname,PSFC,Times,XLAT,XLONG ${dir_in}wrfhrly_d01_${year}-${month}-01_00:00:00 ${dir_out}PSFC.nc 
# # Find number of days in this month by reading number of time steps
# cdo ntime ${dir_out}PSFC.nc > ${dir_out}timesteps_in_month << EOF
# EOF
# days_in_month=$((`cat ${dir_out}timesteps_in_month` / 24))
# timesteps_in_month=`cat ${dir_out}timesteps_in_month` #$(($days_in_month*24))
# 
# # For each day in the month
# for day in `seq 1 $days_in_month`;do
# period=0
# # For each 3-hourly period in the day
# for t in `seq 1 3 24`;do
# time1=$(($t + 24*$(($day - 1))))
# time2=$(($time1 + 1))
# time3=$(($time1 + 2))
# cdo seltimestep,$time1 ${dir_out}PSFC.nc ${dir_out}temp_${t}_T1_PSFC.nc 
# cdo seltimestep,$time2 ${dir_out}PSFC.nc ${dir_out}temp_${t}_T2_PSFC.nc 
# cdo seltimestep,$time3 ${dir_out}PSFC.nc ${dir_out}temp_${t}_T3_PSFC.nc 
# period=$(($period + 1))
# cdo ensmean ${dir_out}temp_${t}_T1_PSFC.nc ${dir_out}temp_${t}_T2_PSFC.nc ${dir_out}temp_${t}_T3_PSFC.nc ${dir_out}temp_${period}_P_PSFC.nc 
# rm ${dir_out}temp_*_T*.nc
# done
# if [ $day -lt 10 ] 
# then 
# today=0$day
# else
# today=$day
# fi 
# cdo mergetime ${dir_out}temp_*_P_PSFC.nc ${dir_out}wrfhrly_d01_${year}-${month}-${today}_00:00:00_PSFC.nc
# rm ${dir_out}temp_*_P*.nc
# done
# rm ${dir_out}PSFC.nc 
# done
# done

###########################################################################