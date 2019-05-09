#!/usr/bin/python

"""
22/1/19: First batch of weekly sim jobs have been run on raijin. While the jobs were in the queue, some 
started, and at one point NCI noted that the efficiency was too low, so they Held the remaining jobs. 
(Efficiency was low because mpirun comman was incorrectly specifying the use of 1 CPU instead of all
16 on a node.) >> ]13/2/19 Actually, we now think the problem was that many files were simulateneously trying to 
access the same NCI directory. Raijin's settings mean this is slow. Currently seeing if improvements made from
writing outout to local jobfs folder, and then copying to actual output directory.]
413 jobs had already finished by that time, so we don't want to rerun them. This script creates the 
decks for the sim weeks that have not yet been computed, with the mpirun command fixed.

29/1/19: Winter batch of weekly sim jobs have been run on raijin, and the spring batch is still running. 
For both seasons some weeks finished and others didn't. I think the walltime needs to be increaed for 
summer and autumn, since some of the spring weeks didn't finish within 24h. 
    > Need to create a new csv of completed weeks though.
    
This script creates .pbs scripts, not decks. This is because we need to specify raijin environment
variables, and when running a deck to create a pbs script, the environment variables are blank, as the job
is not actually running at that point. 

pbs scripts are created for those model days that are identified as incomplete, using Check_all_days_computed.py.     


"""

import pandas as pd
import numpy as np

Start_date = '19790201' 
End_date = '20131231' 

daylist = pd.date_range(Start_date,End_date,freq='d').strftime('%Y-%m-%d')

# List of incomplete days (taken from /g/data/xc0/user/Holgate/QIBT/exp02/QIBT_exp02_incomplete_days.txt)
df = pd.read_csv('/home/z3131380/hdrive/PhD/PartB/Scripts/Model/QIBT_exp02_incomplete_days.txt')
# Remove first and last days, as we know we can't run the model for those days
df = df.drop(df.index[0]); df = df.drop(df.index[-1])
# Remove column of zeros
df = df.drop('Unnamed: 0',axis=1)
# Reset the index to start from zero
df = df.reset_index()

#name the input deck to use
inscript = "QIBT_exp02_date.pbs"

for d in range(len(df)):
    startyear = df['Day'][d][3:7]
    startmonth = df['Day'][d][7:9]
    if len(df['Day'][d])==14:
        startday = '0'+df['Day'][d][-4:-3]
    elif len(df['Day'][d])==15:
        startday = df['Day'][d][-5:-3]
    
    endyear = daylist[np.where(daylist==startyear+'-'+startmonth+'-'+startday)[0][0]+1][0:4]
    endmonth = daylist[np.where(daylist==startyear+'-'+startmonth+'-'+startday)[0][0]+1][5:7]
    endday = daylist[np.where(daylist==startyear+'-'+startmonth+'-'+startday)[0][0]+1][-2:]
    
    #open the sample pbs script 
    fin = open ('/home/z3131380/hdrive/PhD/PartB/Scripts/Model/'+inscript,"r") 

    #open the script I am creating
    if int(startmonth)==12 or int(startmonth)==1 or int(startmonth)==2:
        fout = open ("/home/z3131380/hdrive/PhD/PartB/Scripts/Model/pbs_scripts/QIBT_exp02_%s_%s_%s.pbs"%(str(startyear),str(startmonth),str(startday)),"w")
        wallhours = 8
    if int(startmonth)==3 or int(startmonth)==4 or int(startmonth)==5:    
        fout = open ("/home/z3131380/hdrive/PhD/PartB/Scripts/Model/pbs_scripts/QIBT_exp02_%s_%s_%s.pbs"%(str(startyear),str(startmonth),str(startday)),"w")
        wallhours = 8
    if int(startmonth)==6 or int(startmonth)==7 or int(startmonth)==8:  
        fout = open ("/home/z3131380/hdrive/PhD/PartB/Scripts/Model/pbs_scripts/QIBT_exp02_%s_%s_%s.pbs"%(str(startyear),str(startmonth),str(startday)),"w")
        wallhours = 5
    if int(startmonth)==9 or int(startmonth)==10 or int(startmonth)==11:    
        fout = open ("/home/z3131380/hdrive/PhD/PartB/Scripts/Model/pbs_scripts/QIBT_exp02_%s_%s_%s.pbs"%(str(startyear),str(startmonth),str(startday)),"w")
        wallhours = 8
    
    #Loop over the lines of the input file
    for lines in fin.readlines():
        lines = lines.replace("%proj%", 'w28')
        lines = lines.replace("%queue%", 'normalbw')
        lines = lines.replace("%ncpus%", str(28))
        lines = lines.replace("%lmem%", str(128))
        lines = lines.replace("%wallhours%",str(wallhours))
        lines = lines.replace("%jobfs%",str(500))
        lines = lines.replace("%nthreads%",str(28))
        lines = lines.replace("%dirout%",'/g/data/xc0/user/Holgate/QIBT/exp02')
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
        