## This script creates grids of mean daily and mean monthly U and V wind components at NARCliM model level 10.
## Checked and OK 27/2/19.

dirdata_atm=/srv/ccrc/data33/z3481416/CCRC-WRF3.6.0.5-SEB/ERA-Interim/R2_nudging/out/
dir_out=/srv/ccrc/data03/z3131380/PartB/NARCliM_postprocess/Prevailing_winds/

for year in {1990..2013};do
for filename in /${dirdata_atm}wrfout_d01_${year}*;do 
# Extract variables
ncks -v XLAT_U,XLAT_V,XLONG_U,XLONG_V,XLAT,XLONG,U,V,W ${filename} ${dir_out}$(basename ${filename})_temp1.nc
for level in {2,10};do
# Select model level  
cdo -s sellevel,${level} ${dir_out}$(basename ${filename})_temp1.nc ${dir_out}$(basename ${filename})_temp2.nc
# Translate into mean daily values
cdo -s timmean ${dir_out}$(basename ${filename})_temp2.nc ${dir_out}$(basename ${filename})_dailymean_L${level}.nc
done
rm ${dir_out}$(basename ${filename})_temp*.nc
done

for level in {2,10};do
# Translate into mean monthly values
for month in {01..12};do
if [ $month -eq 01 ]
then ed=31
elif [ $month -eq 02 ]
then
if [ $(( year % 400 == 0 || (year % 4 == 0 && year % 100 != 0) )) -eq 0 ]
then ed=28
else ed=29
fi
elif [ $month -eq 03 ]
then ed=31
elif [ $month -eq 04 ]
then ed=30
elif [ $month -eq 05 ]
then ed=31
elif [ $month -eq 06 ]
then ed=30
elif [ $month -eq 07 ]
then ed=31
elif [ $month -eq 08 ]
then ed=31
elif [ $month -eq 09 ]
then ed=30
elif [ $month -eq 10 ]
then ed=31
elif [ $month -eq 11 ]
then ed=30
elif [ $month -eq 12 ]
then ed=31
fi
eval cdo -s -z szip ensmean ${dir_out}wrfout_d01_${year}-${month}-{01..$ed}_00:00:00_dailymean_L${level}.nc ${dir_out}wrfout_d01_${year}-${month}_00:00:00_monthlymean_L${level}.nc
done
done
done
