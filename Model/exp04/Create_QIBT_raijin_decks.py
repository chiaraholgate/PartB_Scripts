#!/usr/bin/python

"""
This script creates one deck for each model simulation week.
Each deck creates and submits a PBS jobscript to the Raijin queue.
"""

import pandas as pd
import numpy as np

#==============================================================================
# Project info
#==============================================================================
proj = 'w28'
results_dir = '/g/data/xc0/user/Holgate/QIBT/exp04'
jobfs = 3 # GB
queue = 'normalsl'

if queue == 'normal':
    ncpus = 16
    lmem = 32 # GB
elif queue == 'normalbw':
    ncpus = 28
    lmem = 128 # GB     
elif queue == 'normalsl':    
    ncpus = 32
    lmem = 190 # GB..it's good to ask for a little bit less than the available memory at any compute 
                        #node even when the job is taking the entire node so that the background `root`/system 
                        #processes can still run with their own designated memory share.

#==============================================================================
# Dates
#==============================================================================
Start_date = '19790201' 
End_date = '20131230' 

weeklist = pd.date_range(Start_date,End_date,freq='w')

#==============================================================================
# Directories
#==============================================================================
#name the input deck to use
indeck = '/home/z3131380/hdrive/PhD/PartB/Scripts/Model/exp04/QIBT_exp04.deck'

#directory you want decks to be saved in
deck_dir = '/home/z3131380/hdrive/PhD/PartB/Scripts/Model/exp04/Decks/'

#==============================================================================
# Create decks
#==============================================================================
for w in np.arange(0,len(weeklist)-1):
    startday = weeklist[w].day; endday = weeklist[w+1].day
    startmonth = weeklist[w].month; endmonth = weeklist[w+1].month
    startyear = weeklist[w].year; endyear = weeklist[w+1].year   
    #print w,startday,startmonth,startyear

    #open the sample deck 
    fin = open (indeck,"r") 

    #open the deck I am creating
    if startmonth==12 or startmonth==1 or startmonth==2:
        fout = open (deck_dir+"/Summer/QIBT_exp04_%s_%s_%s.deck"%(str(startyear),str(startmonth),str(startday)),"w")
        wallhours = 24
    if startmonth==3 or startmonth==4 or startmonth==5:    
        fout = open (deck_dir+"/Autumn/QIBT_exp04_%s_%s_%s.deck"%(str(startyear),str(startmonth),str(startday)),"w")
        wallhours = 18
    if startmonth==6 or startmonth==7 or startmonth==8:  
        fout = open (deck_dir+"/Winter/QIBT_exp04_%s_%s_%s.deck"%(str(startyear),str(startmonth),str(startday)),"w")
        wallhours = 12
    if startmonth==9 or startmonth==10 or startmonth==11:   
        fout = open (deck_dir+"/Spring/QIBT_exp04_%s_%s_%s.deck"%(str(startyear),str(startmonth),str(startday)),"w")
        wallhours = 18

    
    #Loop over the lines of the input file
    for lines in fin.readlines():
        lines = lines.replace("%proj%", proj)
        lines = lines.replace("%queue%", queue)
        lines = lines.replace("%ncpus%", str(ncpus))
        lines = lines.replace("%lmem%", str(lmem))
        lines = lines.replace("%wallhours%",str(wallhours))
        lines = lines.replace("%jobfs%",str(jobfs))
        lines = lines.replace("%nthreads%",str(ncpus))
        lines = lines.replace("%dirout%",results_dir)
        lines = lines.replace("%syear%", str(startyear))
        lines = lines.replace("%smon%", str(startmonth))
        lines = lines.replace("%sday%", str(startday))
        lines = lines.replace("%edday%", str(endday))
        lines = lines.replace("%edmon%", str(endmonth))
        lines = lines.replace("%edyear%", str(endyear))
        
        fout.write(lines)
        
    #Close input and output files
    fin.close()
    fout.close()
        