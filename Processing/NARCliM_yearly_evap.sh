## This script converts daily latent heat (W/m2) over 8 timesteps into annual evaporation (mm).

dir_E=/srv/ccrc/data03/z3131380/PartB/NARCliM_postprocess/
dir_out=/srv/ccrc/data03/z3131380/PartB/NARCliM_postprocess/Yearly_evap/


for year in {1980..2013};do
for filename in ${dir_E}wrfhrly_d01_${year}*LH.nc; do
cdo -s expr,'evap=LH*(1440/8)*60/2.25E6;' ${filename} ${dir_out}$(basename ${filename})_evap_TS.nc
cdo -s timsum ${dir_out}$(basename ${filename})_evap_TS.nc ${dir_out}$(basename ${filename})_evap_daily.nc
cdo -s setattribute,evap@units=mm ${dir_out}$(basename ${filename})_evap_daily.nc ${dir_out}$(basename ${filename})_evap_daily2.nc
done
cdo -s -z szip enssum ${dir_out}wrfhrly_d01_${year}*_evap_daily2.nc ${dir_out}wrfhrly_d01_${year}_evap_mm.nc
rm ${dir_out}*_evap_TS.nc
rm ${dir_out}*_evap_daily.nc
rm ${dir_out}*_evap_daily2.nc
done

