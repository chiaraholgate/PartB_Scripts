mask_dir='/srv/ccrc/data03/z3131380/PartB/Masks/'
test_file='/srv/ccrc/data03/z3131380/PartB/Output/MDB/10parcels/bt.198112_1.nc'
test_out='/srv/ccrc/data03/z3131380/PartB/Recycling_ratios/'

ncks -v latitcrs ${test_file} ${test_out}latitcrs.nc
ncks -v longicrs ${test_file} ${test_out}longicrs.nc

# latitcrs,longicrs is the same grid as XLAT,XLONG minus the boundary rows ignored by QIBT program.
# latitcrs,longicrs: 134,205
# XLAT,XLONG: 144,215

# Modify mask grid to match latitcrs,longicrs, then apply mask.
cdo selindexbox,6,210,6,139 ${mask_dir}/mdb_rotpole.nc ${test_out}/mdb_rotpole_crs_grid.nc
cdo div ${test_out}latitcrs.nc ${test_out}/mdb_rotpole_crs_grid.nc ${test_out}latitcrs_MDB.nc
cdo div ${test_out}longicrs.nc ${test_out}/mdb_rotpole_crs_grid.nc ${test_out}longicrs_MDB.nc