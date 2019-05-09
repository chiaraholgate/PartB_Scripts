## This script remaps the NARCliM P and E annual totals to the grid of AWRA runoff (qtot).
## I am doing this because I want to compute P-E=Q for a bulk water balance check of QIBT.
## Note that some remap methods won't work because P and E are have a curvilinear grid, and AWRA has a regular lat/lon
## grid. So I'm using remapnn, but remapcon would technically be more correct I think.

dir_Q=/srv/ccrc/data03/z3131380/PartA/Data/AWRA/Runoff/
dir_E=/srv/ccrc/data03/z3131380/PartB/NARCliM_postprocess/Yearly_evap/
dir_P=/srv/ccrc/data03/z3131380/PartB/NARCliM_postprocess/Yearly_rain/

# cdo griddes ${dir_Q}/qtot.nc

cat > ${dir_Q}/AWRA_qtot_grid << EOF
gridtype  = lonlat
gridsize  = 572721
xsize     = 841
ysize     = 681
xname     = longitude
xlongname = "longitude"
xunits    = "degrees_east"
yname     = latitude
ylongname = "latitude"
yunits    = "degrees_north"
xfirst    = 112
xinc      = 0.05
yfirst    = -10
yinc      = -0.0500000000000007
EOF

for year in {1979..2013};do
# Regrid using first order conservative remapping for precipitation (e.g. see https://climatedataguide.ucar.edu/climate-data-tools-and-analysis/regridding-overview)

cdo -s -z zip remapnn,${dir_Q}/AWRA_qtot_grid ${dir_P}/wrfhrly_d01_${year}_rain_mm.nc ${dir_P}/wrfhrly_d01_${year}_rain_mm_remap_AWRAqtot.nc

cdo -s -z zip remapnn,${dir_Q}/AWRA_qtot_grid ${dir_E}/wrfhrly_d01_${year}_evap_mm.nc ${dir_E}/wrfhrly_d01_${year}_evap_mm_remap_AWRAqtot.nc

done