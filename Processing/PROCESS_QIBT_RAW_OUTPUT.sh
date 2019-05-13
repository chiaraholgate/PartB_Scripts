#!/bin/bash

# This script processes the raw QIBT results to output:
# (1) One .nc file per year of annual and daily water vapour contribution by region of interest.
# (2) One .nc file per region of seasonal climatology of water vapour contribution.
# (3) (a) One .csv file per region of annual recycling ratio, terrestrial and oceanic contributions.
#     (b) One .nc file per region of gridded annual recycling ratio.

dir_proc=/home/z3131380/hdrive/PhD/PartB/Scripts/Processing/
dir_anal=/home/z3131380/hdrive/PhD/PartB/Scripts/Analysis/

# (1) Compress raw output of 2d vapour contribution per rain grid cell, to 3d contribution per region
python ${dir_proc}Compress_yearly_wvcont_region.py

# (2) Compute seasonal climatology of source regions
python ${dir_proc}Climatology_wvcont_region.py

# (3) Compute rainfall recycling by region, annually
python ${dir_anal}Rainfall_recycling_region.py annual
# And then seasonally
python ${dir_anal}Rainfall_recycling_region.py seasonal