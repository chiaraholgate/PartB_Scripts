# Find maximum U,V winds in NARCLiM data set to help determine QIBT model time step.
# This script finds the maximum u and v components as output from NARCliM.
# The max windspeed (assuming that these max components occur at the same time and place) is then calculated in /home/z3131380/hdrive/PhD/PartB/Scripts/Model_timestep.py

dir_data='/srv/ccrc/data33/z3481416/CCRC-WRF3.6.0.5-SEB/ERA-Interim/R2_nudging/out/'
dir_out='/srv/ccrc/data03/z3131380/PartB/Model_testing/Max_windspeed/'

for year in {1979..2013};do
for month in {01..12};do
# max grid over all files
cdo ensmax ${dir_data}wrfout_d01_${year}-${month}* ${dir_out}wrfout_d01_${year}-${month}_max.nc
ncks -v U,V ${dir_out}wrfout_d01_${year}-${month}_max.nc ${dir_out}wrfout_d01_${year}-${month}_UV_max_1.nc

# max of levels
cdo vertmax ${dir_out}wrfout_d01_${year}-${month}_UV_max_1.nc ${dir_out}wrfout_d01_${year}-${month}_UV_max_2.nc

# max of timesteps
cdo timmax ${dir_out}wrfout_d01_${year}-${month}_UV_max_2.nc ${dir_out}wrfout_d01_${year}-${month}_UV_max_3.nc

# max of grid
cdo fldmax ${dir_out}wrfout_d01_${year}-${month}_UV_max_3.nc ${dir_out}wrfout_d01_${year}-${month}_UV_max.nc

rm ${dir_out}wrfout_d01_${year}-${month}_UV_max_*.nc
rm ${dir_out}wrfout_d01_${year}-${month}_max.nc
done
done
