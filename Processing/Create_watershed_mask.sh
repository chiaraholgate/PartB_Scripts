# Create 'watershed' mask for BTA analysis.

#dir='/srv/ccrc/data03/z3131380/PartB/Masks/'
dir='/home/z3131380/hdrive/PhD/PartB/Masks/'
dir_data='/srv/ccrc/data03/z3131380/PartB/NARCliM_postprocess/'

######################################################################################
#### Create mask of Australia land area

# Create land-sea mask of whole world (land=1, sea=0)
cdo -f nc -gtc,0 -topo ${dir}land_sea_mask.nc 
#cdo sellonlatbox,110,155,-10,-45 ${dir}land_sea_mask.nc ${dir}land_sea_mask_AUS.nc 

# Remap to target grid
cdo remapnn,${dir_data}NARCliM_in_grid ${dir}land_sea_mask.nc ${dir}land_sea_mask_temp1.nc

# Set cells outside of Australia to 0
cdo setclonlatbox,0,90,180,12,-2 ${dir}land_sea_mask_temp1.nc ${dir}land_sea_mask_temp2.nc
cdo setclonlatbox,0,90,120,0,-11 ${dir}land_sea_mask_temp2.nc ${dir}land_sea_mask_temp3.nc
cdo setclonlatbox,0,160,180,0,-50 ${dir}land_sea_mask_temp3.nc ${dir}land_sea_mask_temp4.nc
cdo setclonlatbox,0,120,127,0,-10 ${dir}land_sea_mask_temp4.nc ${dir}land_sea_mask_temp5.nc
#cdo setclonlatbox,0,112,132,0,12 ${dir}land_sea_mask_temp5.nc ${dir}land_sea_mask_temp6.nc
cdo chname,topo,wsmask ${dir}land_sea_mask_temp5.nc ${dir}NARCliM_AUS_land_sea_mask.nc
rm ${dir}*temp*.nc
######################################################################################


######################################################################################
# MDB netcdf files provided by Jim Clunie (MDBA) via email 9/8/17 and 17/8/17.
#### Create mask of Murray-Darling Basin
cdo remapnn,${dir}NARCliM_in_grid ${dir}mdb.nc ${dir}mdb_rotpole_temp.nc
cdo chname,mdb,wsmask ${dir}mdb_rotpole_temp.nc ${dir}mdb_rotpole.nc

#### Create mask of southern MDB 
cdo remapnn,${dir}NARCliM_in_grid ${dir}n_s_basin/sbasin.nc ${dir}n_s_basin/sbasin_rotpole_temp1.nc
cdo chname,mdb,wsmask ${dir}n_s_basin/sbasin_rotpole_temp1.nc ${dir}n_s_basin/sbasin_rotpole_temp2.nc
cdo ifthenc,1 ${dir}n_s_basin/sbasin_rotpole_temp2.nc ${dir}n_s_basin/sbasin_rotpole_temp3.nc
cdo chname,sbasin,wsmask ${dir}n_s_basin/sbasin_rotpole_temp3.nc ${dir}n_s_basin/sbasin_rotpole.nc
rm ${dir}n_s_basin/sbasin_*temp*

#### Create mask of northern MDB 
cdo remapnn,${dir}NARCliM_in_grid ${dir}n_s_basin/nbasin.nc ${dir}n_s_basin/nbasin_rotpole_temp.nc
cdo chname,nbasin,wsmask ${dir}n_s_basin/nbasin_rotpole_temp.nc ${dir}n_s_basin/nbasin_rotpole.nc
rm ${dir}n_s_basin/nbasin_*temp*

#### Create mask of MDB surface water Water Resource Plan areas
cdo remapnn,${dir}NARCliM_in_grid ${dir}n_s_basin/swwrpa.nc ${dir}n_s_basin/swwrpa_rotpole_temp.nc
cdo chname,swwrpa,wsmask ${dir}n_s_basin/swwrpa_rotpole_temp.nc ${dir}n_s_basin/swwrpa_rotpole.nc

rm ${dir}n_s_basin/swwrpa_*temp*

#### Create mask of specific MDB surface water Water Resource Plan area

# Moonie: SW18
cdo eqc,18 ${dir}n_s_basin/swwrpa_rotpole.nc ${dir}n_s_basin/moonie_rotpole_temp1.nc
cdo setmissval,0 ${dir}n_s_basin/moonie_rotpole_temp1.nc ${dir}n_s_basin/moonie_rotpole.nc
rm ${dir}n_s_basin/moonie_*temp*

######################################################################################




######################################################################################
# Here just selecting a lat/lon box to test the BTA program over a small-ish area in SE Aust.

# f_in='/srv/ccrc/data03/z3131380/PartB/test_data/NARCliM_test_data/wrfout_d01_1950-01-01_00:00:00.nc'
# dir_out='/srv/ccrc/data03/z3131380/PartB/'
# 
# ncks -v PB $f_in ${dir_out}test_watershed_mask_temp0.nc
# ncwa -a bottom_top ${dir_out}test_watershed_mask_temp0.nc ${dir_out}test_watershed_mask_temp1.nc #Collapse level dim
# ncwa -a Times ${dir_out}test_watershed_mask_temp1.nc ${dir_out}test_watershed_mask_temp2.nc #Collape time dim
# cdo masklonlatbox,142,148,-32,-40 ${dir_out}test_watershed_mask_temp2.nc ${dir_out}test_watershed_mask_temp3.nc
# cdo setclonlatbox,1,142,148,-32,-40 ${dir_out}test_watershed_mask_temp3.nc ${dir_out}test_watershed_mask_temp4.nc
# cdo setmisstoc,0 ${dir_out}test_watershed_mask_temp4.nc ${dir_out}test_watershed_mask_temp5.nc
# cdo chname,PB,wsmask ${dir_out}test_watershed_mask_temp5.nc ${dir_out}test_watershed_mask.nc
# ncatted -O -a units,,d,,  ${dir_out}test_watershed_mask.nc 
# ncatted -O -a description,,d,,  ${dir_out}test_watershed_mask.nc 
# rm ${dir_out}*temp*.nc

######################################################################################
