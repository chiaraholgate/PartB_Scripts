dirin=/srv/ccrc/data03/z3131380/PartB/NARCliM_postprocess/Yearly_rain/
dirout=/srv/ccrc/data19/z3131380/PartB/Output/Australia/100parcels/TS10min/Processed/Trends/

cdo -z szip mergetime ${dirin}wrfhrly_d01_*_rain_mm.nc ${dirin}wrfhrly_d01_1979-2013_rain_mm.nc

cdo -z szip trend ${dirin}wrfhrly_d01_1979-2013_rain_mm.nc ${dirout}wrfhrly_d01_1979-2013_rain_trend_a.nc ${dirout}wrfhrly_d01_1979-2013_rain_trend_b.nc


#TREND NOT YET WORKING