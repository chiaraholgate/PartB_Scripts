#!/bin/bash
module load intel/2019.0.117 netcdf-c/4.6.1-intel netcdf-f/4.4.4-intel-netcdf-4.6.1 hdf5/1.10.1-intel
export MP_SHARED_MEMORY=yes
export OMP_NUM_THREADS=8
# need to set numthreads in .f90 as well
ulimit -s unlimited
export FC='ifort'
export NETCDFF='/share/apps/netcdf-f/intel/4.4.4'
export NETCDFFLIBDIR=$NETCDFF/lib
export NETCDFFINCDIR=$NETCDFF/include
export NETCDFC='/share/apps/netcdf-c/intel/4.4.1.1' #'/share/apps/netcdf-c/4.4.1.1'
export NETCDFCLIBDIR=$NETCDFC/lib
export HDF5='/share/apps/hdf5/intel/1.10.1'  #'/share/apps/hdf5/1.10.1-threadsafe'
export HDF5LIBDIR=$HDF5/lib
export SZIP='/share/apps/szip/2.1'
export SZIPLIBDIR=$SZIP/lib

ifort -I$NETCDFFINCDIR -c -fast -free QIBT_exp06.f90 -qopenmp -check bounds -g -traceback
ifort -o QIBT_exp06 QIBT_exp06.o -qopenmp -L$NETCDFCLIBDIR -lnetcdf -L$NETCDFFLIBDIR -lnetcdff -shared-intel -L$HDF5LIBDIR -lhdf5_hl -lhdf5 
#./QIBT_exp06 1 2 1979 2 2 1979 /srv/ccrc/data19/z3131380/PartB/Output/Australia/100parcels/TS10min/exp06/Raw/ 