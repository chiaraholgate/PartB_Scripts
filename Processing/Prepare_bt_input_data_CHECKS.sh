test_data=/srv/ccrc/data03/z3131380/PartB/test_data/
dir_out=/srv/ccrc/data03/z3131380/PartB/NARCliM_postprocess/Processing_checks/

#cdo seltimestep,1/49 ${test_data}wrfhrly_d01_2007-08-01_00:00:00.nc ${test_data}wrfhrly_d01_2007-08-01_00:00:00_test.nc

# this month
cdo -s selname,RAINC,Times,XLAT,XLONG ${test_data}wrfhrly_d01_2007-08-01_00:00:00.nc ${dir_out}RAINC1.nc 
ncrename -v RAINC,RAIN ${dir_out}RAINC1.nc
cdo -s selname,RAINNC,Times,XLAT,XLONG ${test_data}wrfhrly_d01_2007-08-01_00:00:00.nc ${dir_out}RAINNC1.nc 
ncrename -v RAINNC,RAIN ${dir_out}RAINNC1.nc
# next month
cdo -s selname,RAINC,Times,XLAT,XLONG ${test_data}wrfhrly_d01_2007-09-01_00:00:00.nc ${dir_out}RAINC2.nc 
ncrename -v RAINC,RAIN ${dir_out}RAINC2.nc
cdo -s selname,RAINNC,Times,XLAT,XLONG ${test_data}wrfhrly_d01_2007-09-01_00:00:00.nc ${dir_out}RAINNC2.nc 
ncrename -v RAINNC,RAIN ${dir_out}RAINNC2.nc

# Convert cumulative rainfall to amount that fell in the timestep
# For each day in the month
days_in_month=$((`cat /srv/ccrc/data03/z3131380/PartB/NARCliM_postprocess/timesteps_in_month` / 24))
day=30
period=0
# For each 3-hourly period in the day
# (Note that the min timestep CDO allows is 1, not zero.)
for t in `seq 1 3 24`;do
time1=$(($t + 24*$(($day - 1))))
cdo -s seltimestep,$time1 ${dir_out}RAINC1.nc ${dir_out}temp_${t}_T1_RAINC.nc 
cdo -s seltimestep,$time1 ${dir_out}RAINNC1.nc ${dir_out}temp_${t}_T1_RAINNC.nc 
echo $t,$time1,'time1 as usual'

time2=$(($time1 + 3))
if [ $day -lt $days_in_month ]; then
cdo -s seltimestep,$time2 ${dir_out}RAINC1.nc ${dir_out}temp_${t}_T2_RAINC.nc 
cdo -s seltimestep,$time2 ${dir_out}RAINNC1.nc ${dir_out}temp_${t}_T2_RAINNC.nc 
echo $time2,'time 2 as usual, not the last day of month'
else
	if [ $t -lt 22 ]; then #otw if it's the last day of month, but not last period, take T2 as usual
	cdo -s seltimestep,$time2 ${dir_out}RAINC1.nc ${dir_out}temp_${t}_T2_RAINC.nc 
	cdo -s seltimestep,$time2 ${dir_out}RAINNC1.nc ${dir_out}temp_${t}_T2_RAINNC.nc 
	echo  $time2,'last day of month but not last period, so time 2 is T2 as usual' 
	elif [ $t -eq 22 ]; then # if this is the last period in the day, of last day in month, use 1st TS of next month-day as T2
	cdo -s seltimestep,1 ${dir_out}RAINC2.nc ${dir_out}temp_${t}_T2_RAINC.nc 
	cdo -s seltimestep,1 ${dir_out}RAINNC2.nc ${dir_out}temp_${t}_T2_RAINNC.nc
	echo  $time2,'last day and last period of month, so time 2 is first TS of next month'
	fi
fi
period=$(($period + 1))
cdo -s sub ${dir_out}temp_${t}_T2_RAINC.nc ${dir_out}temp_${t}_T1_RAINC.nc ${dir_out}temp_${period}_P_RAINC.nc 
cdo -s sub ${dir_out}temp_${t}_T2_RAINNC.nc ${dir_out}temp_${t}_T1_RAINNC.nc ${dir_out}temp_${period}_P_RAINNC.nc 
done

if [ $day -lt 10 ]; then 
today=0$day
else
today=$day
fi 
# Merge into 3-hourly daily file
cdo -s mergetime ${dir_out}temp_*_P_RAINC.nc ${dir_out}wrfhrly_d01_${year}-${this_month}-${today}_00:00:00_RAINC.nc
cdo -s mergetime ${dir_out}temp_*_P_RAINNC.nc ${dir_out}wrfhrly_d01_${year}-${this_month}-${today}_00:00:00_RAINNC.nc


# Sum convective and non-convective rainfall into one file
cdo -s -f nc4c add ${dir_out}wrfhrly_d01_${year}-${this_month}-${today}_00:00:00_RAINC.nc ${dir_out}wrfhrly_d01_${year}-${this_month}-${today}_00:00:00_RAINNC.nc ${dir_out}wrfhrly_d01_${year}-${this_month}-${today}_00:00:00_RAIN_temp.nc
ncatted -O -a description,RAIN,o,c,"Sum of total cumulus and non-cumulus precipitation" ${dir_out}wrfhrly_d01_${year}-${this_month}-${today}_00:00:00_RAIN_temp.nc

# Find negative values from use of rainfall_bucket in climate model, and add bucket depth (1000mm) to those values
cdo -mulc,1000 -ltc,0 ${dir_out}wrfhrly_d01_${year}-${this_month}-${today}_00:00:00_RAIN_temp.nc ${dir_out}wrfhrly_d01_${year}-${this_month}-${today}_00:00:00_RAIN_temp_lt0.nc
cdo add ${dir_out}wrfhrly_d01_${year}-${this_month}-${today}_00:00:00_RAIN_temp_lt0.nc ${dir_out}wrfhrly_d01_${year}-${this_month}-${today}_00:00:00_RAIN_temp.nc ${dir_out}wrfhrly_d01_${year}-${this_month}-${today}_00:00:00_RAIN.nc
rm ${dir_out}wrfhrly_d01_${year}-${this_month}-${today}_00:00:00_RAIN_temp*


# # Calculate histogram of wind speeds based on pre-defined bins. Outputs histogram per grid cell.
bins=-1100,-1000,-900,-5,0,5,100,1000
cdo histcount,$bins ${dir_out}wrfhrly_d01_${year}-${this_month}-${today}_00:00:00_RAIN.nc ${dir_out}wrfhrly_d01_${year}-${this_month}-${today}_00:00:00_RAIN_histcount_temp.nc
cdo setattribute,description="Histogram count of days when rainfall depth is in each bin. Bin size -1100mm to +10000mm." ${dir_out}wrfhrly_d01_${year}-${this_month}-${today}_00:00:00_RAIN_histcount_temp.nc ${dir_out}wrfhrly_d01_${year}-${this_month}-${today}_00:00:00_RAIN_histcount.nc
rm ${dir_out}wrfhrly_d01_${year}-${this_month}-${today}_00:00:00_RAIN_histcount_temp.nc
