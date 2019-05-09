#!/usr/bin/python

import pandas as pd
import numpy as np

Start_date = '19790201' 
End_date = '20131231' 

weeklist = pd.date_range(Start_date,End_date,freq='w')

#name the input deck to use
indeck = "run_QIBT_exp02.deck"

for w in np.arange(0,len(weeklist)-1):
    startday = weeklist[w].day; endday = weeklist[w+1].day
    startmonth = weeklist[w].month; endmonth = weeklist[w+1].month
    startyear = weeklist[w].year; endyear = weeklist[w+1].year   

    #open the sample deck 
    fin = open ('/home/z3131380/hdrive/PhD/PartB/Scripts/Model/'+indeck,"r") 

    #open the deck I am creating
    fout = open ("/home/z3131380/hdrive/PhD/PartB/Scripts/Model/Decks/run_QIBT_exp02_%s_%s_%s.deck"%(str(startyear),str(startmonth),str(startday)),"w")
    
    #Loop over the lines of the input file
    for lines in fin.readlines():
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
        