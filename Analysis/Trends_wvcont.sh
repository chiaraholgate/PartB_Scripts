dirin=/srv/ccrc/data19/z3131380/PartB/Output/Australia/100parcels/TS10min/Processed/Yearly/
dirout=/srv/ccrc/data19/z3131380/PartB/Output/Australia/100parcels/TS10min/Processed/Trends/

cdo -z szip mergetime ${dirin}Australia*wvcont.nc ${dirin}Australia_1979_2013_wvcont_temp.nc
ncks -x -v pre_annual_total,wv_cont_sum_yearly_mm,wv_cont_sum_yearly_pct ${dirin}Australia_1979_2013_wvcont_temp.nc ${dirin}Australia_1979_2013_wvcont.nc
rm Australia_1979_2013_wvcont_temp.nc

cdo -z szip trend ${dirin}Australia_1979_2013_wvcont.nc ${dirout}Australia_1979_2013_wvcont_trend_a.nc ${dirout}Australia_1979_2013_wvcont_trend_b.nc


TREND NOT YET WORKING