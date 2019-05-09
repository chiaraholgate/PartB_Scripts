# Find maximum U,V winds in NARCLiM data set to help determine QIBT model time step.
# This script finds the maximum u and v components as output from NARCliM. It only considers wind speeds below 
# approx 300hPa, which is approx model level = (300-50)/(1000-50)=0.263~~model level 22 from wrfout variable "z".
# The max windspeed (assuming that these max components occur at the same time and place) is then calculated in /home/z3131380/hdrive/PhD/PartB/Scripts/Model_timestep.py

dir_data='/srv/ccrc/data33/z3481416/CCRC-WRF3.6.0.5-SEB/ERA-Interim/R2_nudging/out/'
dir_out='/srv/ccrc/data03/z3131380/PartB/Model_testing/Max_windspeed/'

for year in {1979..2013};do
for month in {01..12};do

# max grid over all files
cdo -s ensmax ${dir_data}wrfout_d01_${year}-${month}* ${dir_out}wrfout_d01_${year}-${month}_max_1.nc

# Extract U,V variables
cdo -s -z szip selname,U,V ${dir_out}wrfout_d01_${year}-${month}_max_1.nc ${dir_out}wrfout_d01_${year}-${month}_UV_max_2.nc

# Max of time and max of levels from level 22 downward, since we're only interested in the max windspeed of the troposphere 
# (assuming there's little moisture in the stratosphere so QIBT won't be moving parcels that high).
cdo -s sellevel,1/22 ${dir_out}wrfout_d01_${year}-${month}_UV_max_2.nc ${dir_out}wrfout_d01_${year}-${month}_UV_max_3.nc
cdo -s -z szip vertmax -timmax ${dir_out}wrfout_d01_${year}-${month}_UV_max_3.nc ${dir_out}wrfout_d01_${year}-${month}_UV_max_L0-22.nc

rm ${dir_out}wrfout_d01_${year}-${month}_max_1.nc
rm ${dir_out}wrfout_d01_${year}-${month}_UV_max_2.nc
rm ${dir_out}wrfout_d01_${year}-${month}_UV_max_3.nc

done
done

# for year in {1979..2013};do
# for month in {01..12};do
# # max grid over all files
# cdo ensmax ${dir_data}wrfout_d01_${year}-${month}* ${dir_out}wrfout_d01_${year}-${month}_max.nc
# ncks -v U,V ${dir_out}wrfout_d01_${year}-${month}_max.nc ${dir_out}wrfout_d01_${year}-${month}_UV_max_1.nc
# 
# # max of levels from level 22 downward, since we're only interested in the max windspeed of the troposphere 
# # (assuming there's little moisture in the stratosphere so QIBT won't be moving parcels that high).
# cdo sellev,0/22 ${dir_out}wrfout_d01_${year}-${month}_UV_max_1.nc ${dir_out}wrfout_d01_${year}-${month}_UV_max_1_levs.nc
# cdo vertmax ${dir_out}wrfout_d01_${year}-${month}_UV_max_1_levs.nc ${dir_out}wrfout_d01_${year}-${month}_UV_max_2.nc
# 
# # max of timesteps
# cdo timmax ${dir_out}wrfout_d01_${year}-${month}_UV_max_2.nc ${dir_out}wrfout_d01_${year}-${month}_UV_max_3.nc
# 
# # max of grid
# cdo fldmax ${dir_out}wrfout_d01_${year}-${month}_UV_max_3.nc ${dir_out}wrfout_d01_${year}-${month}_UV_max_L0_22.nc
# 
# rm ${dir_out}wrfout_d01_${year}-${month}_UV_max_*.nc
# rm ${dir_out}wrfout_d01_${year}-${month}_max.nc
# done
# done



