
############################################################################

# This script creates input data files from NARCliM output for QIBT analysis over Australia. 
# Necessary atmospheric variables from daily files of 3 hourly output sourced directly from server.
# Precipitation and evaporation data extracted from monthly files of 1 hourly output and accumulated to daily files of 3 hourly to match atmospheric variables.

# Data sourced from NARCliM output, version R2_nudging. This version uses wrf output to constrain upper tropospheric variables on a scale of ~500km to ensure simulated large-scale processes are similar to observed/reanalysis.

# March 2019: Found error in RAIN processing. Last day of month was taking first timestep of next month for all "time2" of that day, when it should have only been taking "time2" as the next day for the last period of the last day of the month. Fixed and re-processed 14/3/19. Previous RAIN processed files moved to ${dir_out}/Erroneous_RAIN.

dir_in='/srv/ccrc/data33/z3481416/CCRC-WRF3.6.0.5-SEB/ERA-Interim/R2_nudging/out/'
dir_out='/srv/ccrc/data03/z3131380/PartB/NARCliM_postprocess/'

###########################################################################
# RAINFALL

# Take difference between first and third timesteps to convert hourly cumulative rainfall to a total over 3-hours

for year in {1979..2013};do
for month in {1..12};do
if [ $month -le 8 ]; then
this_month=0$month
next_month=0$(($month+1))
next_year=$year
elif [ $month -eq 9 ];then
this_month=0$month
next_month=$(($month+1))
next_year=$year
elif [ $month -eq 10 ];then
this_month=$month
next_month=$(($month+1))
next_year=$year
elif [ $month -eq 11 ];then
this_month=$month
next_month=$(($month+1))
next_year=$year
elif [ $month -eq 12 ];then
this_month=$month
next_month=01
next_year=$(($year+1))
fi

# "RAINC": cumulative convective rainfall ; "RAINNC": cumulative non-convective rainfall
# > this month
cdo -s selname,RAINC,Times,XLAT,XLONG ${dir_in}wrfhrly_d01_${year}-${this_month}-01_00:00:00 ${dir_out}RAINC1.nc 
ncrename -v RAINC,RAIN ${dir_out}RAINC1.nc
cdo -s selname,RAINNC,Times,XLAT,XLONG ${dir_in}wrfhrly_d01_${year}-${this_month}-01_00:00:00 ${dir_out}RAINNC1.nc 
ncrename -v RAINNC,RAIN ${dir_out}RAINNC1.nc
# > next month
cdo -s selname,RAINC,Times,XLAT,XLONG ${dir_in}wrfhrly_d01_${next_year}-${next_month}-01_00:00:00 ${dir_out}RAINC2.nc 
ncrename -v RAINC,RAIN ${dir_out}RAINC2.nc
cdo -s selname,RAINNC,Times,XLAT,XLONG ${dir_in}wrfhrly_d01_${next_year}-${next_month}-01_00:00:00 ${dir_out}RAINNC2.nc 
ncrename -v RAINNC,RAIN ${dir_out}RAINNC2.nc

# Find number of days in this month by reading number of time steps
cdo -s ntime ${dir_out}RAINC1.nc > ${dir_out}timesteps_in_month << EOF
EOF
days_in_month=$((`cat ${dir_out}timesteps_in_month` / 24))
timesteps_in_month=`cat ${dir_out}timesteps_in_month` #$(($days_in_month*24))

# Convert cumulative rainfall to amount that fell in the timestep
# For each day in the month
for day in `seq 1 $days_in_month`;do
period=0
# For each 3-hourly period in the day
for t in `seq 1 3 24`;do
time1=$(($t + 24*$(($day - 1))))
cdo -s seltimestep,$time1 ${dir_out}RAINC1.nc ${dir_out}temp_${t}_T1_RAINC.nc 
cdo -s seltimestep,$time1 ${dir_out}RAINNC1.nc ${dir_out}temp_${t}_T1_RAINNC.nc 
time2=$(($time1 + 3))
if [ $day -lt $days_in_month ]; then
cdo -s seltimestep,$time2 ${dir_out}RAINC1.nc ${dir_out}temp_${t}_T2_RAINC.nc 
cdo -s seltimestep,$time2 ${dir_out}RAINNC1.nc ${dir_out}temp_${t}_T2_RAINNC.nc 
else
	if [ $t -lt 22 ]; then # if it's the last day of month, but not last period, take T2 as usual
	cdo -s seltimestep,$time2 ${dir_out}RAINC1.nc ${dir_out}temp_${t}_T2_RAINC.nc 
	cdo -s seltimestep,$time2 ${dir_out}RAINNC1.nc ${dir_out}temp_${t}_T2_RAINNC.nc 
	elif [ $t -eq 22 ]; then # if this is the last period in the day, of last day in month, use 1st TS of next month-day as T2
	cdo -s seltimestep,1 ${dir_out}RAINC2.nc ${dir_out}temp_${t}_T2_RAINC.nc 
	cdo -s seltimestep,1 ${dir_out}RAINNC2.nc ${dir_out}temp_${t}_T2_RAINNC.nc
	fi
fi
period=$(($period + 1))
cdo -s sub ${dir_out}temp_${t}_T2_RAINC.nc ${dir_out}temp_${t}_T1_RAINC.nc ${dir_out}temp_${period}_P_RAINC.nc 
cdo -s sub ${dir_out}temp_${t}_T2_RAINNC.nc ${dir_out}temp_${t}_T1_RAINNC.nc ${dir_out}temp_${period}_P_RAINNC.nc 
rm ${dir_out}temp_*_T*.nc
done
if [ $day -lt 10 ]; then 
today=0$day
else
today=$day
fi 
# Merge into 3-hourly daily file
cdo -s mergetime ${dir_out}temp_*_P_RAINC.nc ${dir_out}wrfhrly_d01_${year}-${this_month}-${today}_00:00:00_RAINC.nc
cdo -s mergetime ${dir_out}temp_*_P_RAINNC.nc ${dir_out}wrfhrly_d01_${year}-${this_month}-${today}_00:00:00_RAINNC.nc
rm ${dir_out}temp_*_P*.nc

# Sum convective and non-convective rainfall into one file
cdo -s -z szip -f nc4c add ${dir_out}wrfhrly_d01_${year}-${this_month}-${today}_00:00:00_RAINC.nc ${dir_out}wrfhrly_d01_${year}-${this_month}-${today}_00:00:00_RAINNC.nc ${dir_out}wrfhrly_d01_${year}-${this_month}-${today}_00:00:00_RAIN_temp.nc
ncatted -O -a description,RAIN,o,c,"Sum of total cumulus and non-cumulus precipitation" ${dir_out}wrfhrly_d01_${year}-${this_month}-${today}_00:00:00_RAIN_temp.nc
rm ${dir_out}wrfhrly_d01_${year}-${this_month}-${today}_00:00:00_RAINC.nc
rm ${dir_out}wrfhrly_d01_${year}-${this_month}-${today}_00:00:00_RAINNC.nc

# Find negative values from use of rainfall_bucket in climate model, and add bucket depth (1000mm) to those values
cdo -mulc,1000 -ltc,0 ${dir_out}wrfhrly_d01_${year}-${this_month}-${today}_00:00:00_RAIN_temp.nc ${dir_out}wrfhrly_d01_${year}-${this_month}-${today}_00:00:00_RAIN_temp_lt0.nc
cdo add ${dir_out}wrfhrly_d01_${year}-${this_month}-${today}_00:00:00_RAIN_temp_lt0.nc ${dir_out}wrfhrly_d01_${year}-${this_month}-${today}_00:00:00_RAIN_temp.nc ${dir_out}wrfhrly_d01_${year}-${this_month}-${today}_00:00:00_RAIN.nc
rm ${dir_out}wrfhrly_d01_${year}-${this_month}-${today}_00:00:00_RAIN_temp*


done
rm ${dir_out}RAINC1.nc 
rm ${dir_out}RAINNC1.nc 
rm ${dir_out}RAINC2.nc 
rm ${dir_out}RAINNC2.nc
done
done

###########################################################################


###########################################################################
# EVAPORATION

# Take average over the 3 one-hourly timesteps
# [LH given in W/m2 or J/s/m2...so you are averaging the J/s rate over the three hours.]

# for year in {1979..2013};do
# for month in {01..12};do
# # "LH": latent heat flux at the surface
# cdo selname,LH,Times,XLAT,XLONG ${dir_in}wrfhrly_d01_${year}-${month}-01_00:00:00 ${dir_out}LH.nc 
# # Find number of days in this month by reading number of time steps
# cdo ntime ${dir_out}LH.nc > ${dir_out}timesteps_in_month << EOF
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
# cdo seltimestep,$time1 ${dir_out}LH.nc ${dir_out}temp_${t}_T1_LH.nc 
# cdo seltimestep,$time2 ${dir_out}LH.nc ${dir_out}temp_${t}_T2_LH.nc 
# cdo seltimestep,$time3 ${dir_out}LH.nc ${dir_out}temp_${t}_T3_LH.nc 
# period=$(($period + 1))
# cdo ensmean ${dir_out}temp_${t}_T1_LH.nc ${dir_out}temp_${t}_T2_LH.nc ${dir_out}temp_${t}_T3_LH.nc ${dir_out}temp_${period}_P_LH.nc 
# rm ${dir_out}temp_*_T*.nc
# done
# if [ $day -lt 10 ] 
# then 
# today=0$day
# else
# today=$day
# fi 
# cdo mergetime ${dir_out}temp_*_P_LH.nc ${dir_out}wrfhrly_d01_${year}-${month}-${today}_00:00:00_LH.nc
# rm ${dir_out}temp_*_P*.nc
# done
# rm ${dir_out}LH.nc 
# done
# done

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