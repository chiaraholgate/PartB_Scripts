## This script converts daily rainfall (mm) over 8 timesteps into annual rainfall (mm).

dir_P=/srv/ccrc/data03/z3131380/PartB/NARCliM_postprocess/
dir_out=/srv/ccrc/data03/z3131380/PartB/NARCliM_postprocess/Yearly_rain/


for year in {1979..2013};do
for filename in ${dir_P}wrfhrly_d01_${year}*RAIN.nc; do
cdo -s timsum ${filename} ${dir_out}$(basename ${filename})_rain_daily.nc
done
cdo -s -z szip enssum ${dir_out}wrfhrly_d01_${year}*_rain_daily.nc ${dir_out}wrfhrly_d01_${year}_rain_mm.nc
rm ${dir_out}*_rain_daily.nc
done

