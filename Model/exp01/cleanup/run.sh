#!/bin/bash
#PBS -P xc0
#PBS -q express
#PBS -l ncpus=16
#PBS -l mem=32GB
#PBS -l walltime=24:00:00
#PBS -l wd

ulimit -s unlimited
module purge
module load netcdf/4.2.1.1 openmpi/3.0.1 intel-fc/2019.0.117
mpif90 -O2 -xHost QIBT_exp02.f90 -o QIBT_exp02 -lnetcdff -fopenmp

export OMP_NUM_THREADS=16

export NMPI=$(( PBS_NCPUS / OMP_NUM_THREADS ))

mpirun --output-filename run.log --report-bindings --map-by ppr:$((16/OMP_NUM_THREADS )):node:PE=${OMP_NUM_THREADS} -n ${NMPI} ./QIBT_exp02 11 2 1979 18 2 1979
