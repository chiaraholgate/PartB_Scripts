dirdata_atm=/srv/ccrc/data33/z3481416/CCRC-WRF3.6.0.5-SEB/ERA-Interim/R2_nudging/out/
dir_out=/srv/ccrc/data03/z3131380/PartB/NARCliM_postprocess/wrf_var_extraction/

for year in {1979..2013};do
for filename in /${dirdata_atm}wrfout_d01_${year}*;do 
# Extract variables
ncks -v XLAT_U,XLAT_V,XLONG_U,XLONG_V,XLAT,XLONG,U,V,W ${filename} ${dir_out}$(basename ${filename})_temp1.nc
# Select first 10 model levels 
cdo sellevel,1/10 ${dir_out}$(basename ${filename})_temp1.nc ${dir_out}$(basename ${filename})_temp2.nc
# Find mean over model levels
cdo vertmean ${dir_out}$(basename ${filename})_temp2.nc ${dir_out}$(basename ${filename})_temp3.nc
# Translate into mean daily values
cdo timmean ${dir_out}$(basename ${filename})_temp3.nc ${dir_out}$(basename ${filename})_temp4.nc
rm ${dir_out}$(basename ${filename})_temp1.nc 
rm ${dir_out}$(basename ${filename})_temp2.nc 
rm ${dir_out}$(basename ${filename})_temp3.nc
done
cdo -z szip mergetime ${dir_out}wrfout_d01_${year}*_temp4.nc ${dir_out}wrfout_d01_${year}_daily_mean_WS_L1-10.nc
rm ${dir_out}wrfout_d01_${year}*_temp4.nc
done

To get prevailing wind direction in each model layer:


Python wind speed and directino arrows:
https://scitools.org.uk/iris/docs/latest/examples/Meteorology/wind_speed.html


## Calculate pressure and geopotential height from wrf NARCliM simulation output.

## See here: http://forum.wrfforum.com/viewtopic.php?f=8&t=1492
# and http://mailman.ucar.edu/pipermail/wrf-users/2016/004360.html

pressure = P + PB
GpH = (PH + PHB)/9.8