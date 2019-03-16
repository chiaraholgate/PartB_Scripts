## This script converts daily soil moisture (m3/m3) over 8 timesteps into annual soil moisture (mm).
## Note that NARCliM also outputs SH2O (soil liquid water) but for Aus they should be the same.
# To convert soil moisture in volumetric units to a depth, the wrf soil_layers_stag values are used.
# i.e. soil depth [m] x volumetric SM [m3/m3] x 1000 = soil moisture [mm]
# NARCliM output gives 4 soil depths: 0.1, 0.3, 0.7 and 1.5.
# Since we want the total storage change over the course of one year, we only determine the SM on 
# the first and last days of the year, and take the difference. 
# Checked and OK 4/3/19.

dir_Q=/srv/ccrc/data33/z3481416/CCRC-WRF3.6.0.5-SEB/ERA-Interim/R2_nudging/out/
dir_out=/srv/ccrc/data03/z3131380/PartB/NARCliM_postprocess/Yearly_SM/

for year in {1979..2013};do

# First day of year
# Select soil layer
cdo -s select,name=SMOIS,level=1 ${dir_Q}wrfout_d01_${year}-01-01* ${dir_out}wrfout_d01_${year}-01-01_SMOIS_layer1.nc
cdo -s select,name=SMOIS,level=2 ${dir_Q}wrfout_d01_${year}-01-01* ${dir_out}wrfout_d01_${year}-01-01_SMOIS_layer2.nc
cdo -s select,name=SMOIS,level=3 ${dir_Q}wrfout_d01_${year}-01-01* ${dir_out}wrfout_d01_${year}-01-01_SMOIS_layer3.nc
cdo -s select,name=SMOIS,level=4 ${dir_Q}wrfout_d01_${year}-01-01* ${dir_out}wrfout_d01_${year}-01-01_SMOIS_layer4.nc
# Convert from m3/m3 to mm
cdo -s expr,'SMOIS_mm=SMOIS*0.1*1000;' ${dir_out}wrfout_d01_${year}-01-01_SMOIS_layer1.nc ${dir_out}wrfout_d01_${year}-01-01_SMOIS_layer1_mm.nc
cdo -s expr,'SMOIS_mm=SMOIS*0.3*1000;' ${dir_out}wrfout_d01_${year}-01-01_SMOIS_layer2.nc ${dir_out}wrfout_d01_${year}-01-01_SMOIS_layer2_mm.nc
cdo -s expr,'SMOIS_mm=SMOIS*0.7*1000;' ${dir_out}wrfout_d01_${year}-01-01_SMOIS_layer3.nc ${dir_out}wrfout_d01_${year}-01-01_SMOIS_layer3_mm.nc
cdo -s expr,'SMOIS_mm=SMOIS*1.5*1000;' ${dir_out}wrfout_d01_${year}-01-01_SMOIS_layer4.nc ${dir_out}wrfout_d01_${year}-01-01_SMOIS_layer4_mm.nc
# Average over the day for that soil layer
cdo -s timmean ${dir_out}wrfout_d01_${year}-01-01_SMOIS_layer1_mm.nc ${dir_out}wrfout_d01_${year}-01-01_SMOIS_layer1_mm_daily.nc
cdo -s timmean ${dir_out}wrfout_d01_${year}-01-01_SMOIS_layer2_mm.nc ${dir_out}wrfout_d01_${year}-01-01_SMOIS_layer2_mm_daily.nc
cdo -s timmean ${dir_out}wrfout_d01_${year}-01-01_SMOIS_layer3_mm.nc ${dir_out}wrfout_d01_${year}-01-01_SMOIS_layer3_mm_daily.nc
cdo -s timmean ${dir_out}wrfout_d01_${year}-01-01_SMOIS_layer4_mm.nc ${dir_out}wrfout_d01_${year}-01-01_SMOIS_layer4_mm_daily.nc
# Sum the depth of water over all soil layers
cdo -s enssum ${dir_out}wrfout_d01_${year}-01-01_SMOIS_layer*_mm_daily.nc ${dir_out}wrfout_d01_${year}-01-01_SMOIS_mm_daily_alllayers.nc

# Last day of year
# Select soil layer
cdo -s select,name=SMOIS,level=1 ${dir_Q}wrfout_d01_${year}-12-31* ${dir_out}wrfout_d01_${year}-12-31_SMOIS_layer1.nc
cdo -s select,name=SMOIS,level=2 ${dir_Q}wrfout_d01_${year}-12-31* ${dir_out}wrfout_d01_${year}-12-31_SMOIS_layer2.nc
cdo -s select,name=SMOIS,level=3 ${dir_Q}wrfout_d01_${year}-12-31* ${dir_out}wrfout_d01_${year}-12-31_SMOIS_layer3.nc
cdo -s select,name=SMOIS,level=4 ${dir_Q}wrfout_d01_${year}-12-31* ${dir_out}wrfout_d01_${year}-12-31_SMOIS_layer4.nc
# Convert from m3/m3 to mm
cdo -s expr,'SMOIS_mm=SMOIS*0.1*1000;' ${dir_out}wrfout_d01_${year}-12-31_SMOIS_layer1.nc ${dir_out}wrfout_d01_${year}-12-31_SMOIS_layer1_mm.nc
cdo -s expr,'SMOIS_mm=SMOIS*0.3*1000;' ${dir_out}wrfout_d01_${year}-12-31_SMOIS_layer2.nc ${dir_out}wrfout_d01_${year}-12-31_SMOIS_layer2_mm.nc
cdo -s expr,'SMOIS_mm=SMOIS*0.7*1000;' ${dir_out}wrfout_d01_${year}-12-31_SMOIS_layer3.nc ${dir_out}wrfout_d01_${year}-12-31_SMOIS_layer3_mm.nc
cdo -s expr,'SMOIS_mm=SMOIS*1.5*1000;' ${dir_out}wrfout_d01_${year}-12-31_SMOIS_layer4.nc ${dir_out}wrfout_d01_${year}-12-31_SMOIS_layer4_mm.nc
# Average over the day for that soil layer
cdo -s timmean ${dir_out}wrfout_d01_${year}-12-31_SMOIS_layer1_mm.nc ${dir_out}wrfout_d01_${year}-12-31_SMOIS_layer1_mm_daily.nc
cdo -s timmean ${dir_out}wrfout_d01_${year}-12-31_SMOIS_layer2_mm.nc ${dir_out}wrfout_d01_${year}-12-31_SMOIS_layer2_mm_daily.nc
cdo -s timmean ${dir_out}wrfout_d01_${year}-12-31_SMOIS_layer3_mm.nc ${dir_out}wrfout_d01_${year}-12-31_SMOIS_layer3_mm_daily.nc
cdo -s timmean ${dir_out}wrfout_d01_${year}-12-31_SMOIS_layer4_mm.nc ${dir_out}wrfout_d01_${year}-12-31_SMOIS_layer4_mm_daily.nc
# Sum the depth of water over all soil layers
cdo -s enssum ${dir_out}wrfout_d01_${year}-12-31_SMOIS_layer*_mm_daily.nc ${dir_out}wrfout_d01_${year}-12-31_SMOIS_mm_daily_alllayers.nc

# Find difference in soil moisture [mm] between start and end of year
cdo -s -z szip sub ${dir_out}wrfout_d01_${year}-12-31_SMOIS_mm_daily_alllayers.nc ${dir_out}wrfout_d01_${year}-01-01_SMOIS_mm_daily_alllayers.nc ${dir_out}wrfout_d01_${year}_delta_SMOIS_mm.nc

# Remove temporary files
rm ${dir_out}wrfout_d01_${year}-*_SMOIS_layer*.nc
rm ${dir_out}wrfout_d01_${year}-*_SMOIS_mm_daily_alllayers.nc

done




# #for year in {1979..2013};do
# year=1979
# for filename in ${dir_Q}wrfout_d01_${year}*; do
# # Select soil layer
# cdo -s select,name=SMOIS,level=1 ${filename} ${dir_out}$(basename ${filename})_SMOIS_layer1.nc
# cdo -s select,name=SMOIS,level=2 ${filename} ${dir_out}$(basename ${filename})_SMOIS_layer2.nc
# cdo -s select,name=SMOIS,level=3 ${filename} ${dir_out}$(basename ${filename})_SMOIS_layer3.nc
# cdo -s select,name=SMOIS,level=4 ${filename} ${dir_out}$(basename ${filename})_SMOIS_layer4.nc
# # Convert from m3/m3 to mm
# cdo -s expr,'SMOIS_mm=SMOIS*0.1*1000;' ${dir_out}$(basename ${filename})_SMOIS_layer1.nc ${dir_out}$(basename ${filename})_SMOIS_layer1_mm.nc
# cdo -s expr,'SMOIS_mm=SMOIS*0.3*1000;' ${dir_out}$(basename ${filename})_SMOIS_layer2.nc ${dir_out}$(basename ${filename})_SMOIS_layer2_mm.nc
# cdo -s expr,'SMOIS_mm=SMOIS*0.7*1000;' ${dir_out}$(basename ${filename})_SMOIS_layer3.nc ${dir_out}$(basename ${filename})_SMOIS_layer3_mm.nc
# cdo -s expr,'SMOIS_mm=SMOIS*1.5*1000;' ${dir_out}$(basename ${filename})_SMOIS_layer4.nc ${dir_out}$(basename ${filename})_SMOIS_layer4_mm.nc
# # Average over the day for that soil layer
# cdo -s timmean ${dir_out}$(basename ${filename})_SMOIS_layer1_mm.nc ${dir_out}$(basename ${filename})_SMOIS_layer1_mm_daily.nc
# cdo -s timmean ${dir_out}$(basename ${filename})_SMOIS_layer2_mm.nc ${dir_out}$(basename ${filename})_SMOIS_layer2_mm_daily.nc
# cdo -s timmean ${dir_out}$(basename ${filename})_SMOIS_layer3_mm.nc ${dir_out}$(basename ${filename})_SMOIS_layer3_mm_daily.nc
# cdo -s timmean ${dir_out}$(basename ${filename})_SMOIS_layer4_mm.nc ${dir_out}$(basename ${filename})_SMOIS_layer4_mm_daily.nc
# # Sum the depth of water over all soil layers
# cdo -s enssum ${dir_out}$(basename ${filename})_SMOIS_layer*_mm_daily.nc ${dir_out}$(basename ${filename})_SMOIS_mm_daily_alllayers.nc
# done
# cdo -s -z szip enssum ${dir_out}wrfout_d01_${year}*_SMOIS_mm_daily_alllayers.nc ${dir_out}wrfout_d01_${year}_SM_mm.nc
# #rm ${dir_out}*_SMOIS_layer*.nc
# #rm ${dir_out}*_SMOIS_mm_daily_alllayers.nc
# #done
