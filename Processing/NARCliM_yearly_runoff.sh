## This script converts daily runoff (mm) over 8 timesteps into annual runoff (mm).
## NOTE! I THINK THE WRF DATA PROVIDES SFROFF AS AN ACCUMULATED VALUE, EVEN THOUGH METADATA DOESN'T INICATE THIS! 
# When I initially summed the time steps over a day, and then summed the daily totals over a year 
# (see "Bulk_waterbalance_check.py"), the graph looked like Q was accumulative. Then looking the the wrf output, 
# it seems that each time step is an accumulated value.
# THEREFORE I am changing this script to deal with Q as an accumulative value.
# Checked and ok 1/3/19.

dir_Q=/srv/ccrc/data33/z3481416/CCRC-WRF3.6.0.5-SEB/ERA-Interim/R2_nudging/out/
dir_out=/srv/ccrc/data03/z3131380/PartB/NARCliM_postprocess/Yearly_runoff/

# Taking Q as a NON-accumulative value
# for year in {1979..2013};do
# for filename in ${dir_Q}wrfout_d01_${year}*; do
# cdo -s select,name=SFROFF ${filename} ${dir_out}$(basename ${filename})_SFROFF.nc
# cdo -s timsum ${dir_out}$(basename ${filename})_SFROFF.nc ${dir_out}$(basename ${filename})_runoff_daily.nc
# done
# cdo -s -z szip enssum ${dir_out}wrfout_d01_${year}*_runoff_daily.nc ${dir_out}wrfout_d01_${year}_runoff_mm.nc
# rm ${dir_out}*_SFROFF.nc
# rm ${dir_out}*_runoff_daily.nc
# done

# Taking Q as an accumulative value
for year in {1979..2013};do
for filename in ${dir_Q}wrfout_d01_${year}*; do
cdo -s select,name=SFROFF ${filename} ${dir_out}$(basename ${filename})_SFROFF_acc.nc
cdo -s seltimestep,1 ${dir_out}$(basename ${filename})_SFROFF_acc.nc ${dir_out}$(basename ${filename})_SFROFF_acc_TS1.nc
cdo -s seltimestep,8 ${dir_out}$(basename ${filename})_SFROFF_acc.nc ${dir_out}$(basename ${filename})_SFROFF_acc_TS8.nc
cdo -s sub ${dir_out}$(basename ${filename})_SFROFF_acc_TS8.nc ${dir_out}$(basename ${filename})_SFROFF_acc_TS1.nc ${dir_out}$(basename ${filename})_runoff_daily.nc
done
cdo -s -z szip enssum ${dir_out}wrfout_d01_${year}*_runoff_daily.nc ${dir_out}wrfout_d01_${year}_runoff_mm.nc
rm ${dir_out}$(basename ${filename})_SFROFF_acc.nc
rm ${dir_out}$(basename ${filename})_SFROFF_acc_TS1.nc
rm ${dir_out}$(basename ${filename})_SFROFF_acc_TS8.nc
rm ${dir_out}$(basename ${filename})_runoff_daily.nc
done
