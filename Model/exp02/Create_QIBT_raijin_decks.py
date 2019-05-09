#!/usr/bin/python

"""
22/1/19: First batch of weekly sim jobs have been run on raijin. While the jobs were in the queue, some 
started, and at one point NCI noted that the efficiency was too low, so they Held the remaining jobs. 
(Efficiency was low because mpirun comman was incorrectly specifying the use of 1 CPU instead of all
16 on a node.)
413 jobs had already finished by that time, so we don't want to rerun them. This script creates the 
decks for the sim weeks that have not yet been computed, with the mpirun command fixed.

29/1/19: Winter batch of weekly sim jobs have been run on raijin, and the spring batch is still running. 
For both seasons some weeks finished and others didn't. I think the walltime needs to be increaed for 
summer and autumn, since some of the spring weeks didn't finish within 24h. 
    > Need to create a new csv of completed weeks though.


"""

import pandas as pd
import numpy as np

Start_date = '19790201' 
End_date = '20131231' 

weeklist = pd.date_range(Start_date,End_date,freq='w')

# List of completed jobs/sim weeks (taken from /home/z3131380/hdrive/PhD/PartB/Model_Testing/Parallelisation/Run_time_log.xlsx)
df = pd.read_csv('/home/z3131380/hdrive/PhD/PartB/Model_Testing/Completed_weeks_22Jan2019.csv')

#name the input deck to use
indeck = "run_QIBT_exp02.deck"

for w in np.arange(0,len(weeklist)-1):
    startday = weeklist[w].day; endday = weeklist[w+1].day
    startmonth = weeklist[w].month; endmonth = weeklist[w+1].month
    startyear = weeklist[w].year; endyear = weeklist[w+1].year   
    fname = 'QIBT_exp02_'+str(startyear)+'_'+str(startmonth)+'_'+str(startday)
    if fname not in df.values:

        #open the sample deck 
        fin = open ('/home/z3131380/hdrive/PhD/PartB/Scripts/Model/'+indeck,"r") 
    
        #open the deck I am creating
        if startmonth==12 or startmonth==1 or startmonth==2:
            fout = open ("/home/z3131380/hdrive/PhD/PartB/Scripts/Model/Decks/Summer/run_QIBT_exp02_%s_%s_%s.deck"%(str(startyear),str(startmonth),str(startday)),"w")
            wallhours = 48
        if startmonth==3 or startmonth==4 or startmonth==5:    
            fout = open ("/home/z3131380/hdrive/PhD/PartB/Scripts/Model/Decks/Autumn/run_QIBT_exp02_%s_%s_%s.deck"%(str(startyear),str(startmonth),str(startday)),"w")
            wallhours = 48
        if startmonth==6 or startmonth==7 or startmonth==8:  
            continue
            #fout = open ("/home/z3131380/hdrive/PhD/PartB/Scripts/Model/Decks/Winter/run_QIBT_exp02_%s_%s_%s.deck"%(str(startyear),str(startmonth),str(startday)),"w")
            #wallhours = 12
        if startmonth==9 or startmonth==10 or startmonth==11:    
            continue
            #fout = open ("/home/z3131380/hdrive/PhD/PartB/Scripts/Model/Decks/Spring/run_QIBT_exp02_%s_%s_%s.deck"%(str(startyear),str(startmonth),str(startday)),"w")
            #wallhours = 36
        
        #Loop over the lines of the input file
        for lines in fin.readlines():
            lines = lines.replace("%syear%", str(startyear))
            lines = lines.replace("%smon%", str(startmonth))
            lines = lines.replace("%sday%", str(startday))
            lines = lines.replace("%edday%", str(endday))
            lines = lines.replace("%edmon%", str(endmonth))
            lines = lines.replace("%edyear%", str(endyear))
            lines = lines.replace("%wallhours%",str(wallhours))
            
            fout.write(lines)
            
        #Close input and output files
        fin.close()
        fout.close()
        