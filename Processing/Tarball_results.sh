#!/bin/bash

####
## This script:
## (1) moves QIBT results files into subdirectories by year;
## (2) creates an archive folder for each year, and
## (3) gzip's each archive.
## To extract files: tar -xzvf <directory name>.tar.gz
## I note that it's probably better to initally move each file directly into an archive folder, instead of creating a subdirectory of year first. This is my first attempt and I'm being conservative... Next time...


results_dir=/g/data/xc0/user/Holgate/QIBT/exp03/


for year in {1979..2013};do
mkdir $year
for ff in bt.${year}*.nc;do 
mv $ff ${year}/
done
tar -cf ${year}.tar ${year}/
gzip ${year}.tar
echo ${year}
done