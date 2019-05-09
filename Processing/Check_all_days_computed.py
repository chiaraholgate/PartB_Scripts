# -*- coding: utf-8 -*-
"""
Created on Fri Oct 26 10:02:23 2018

Script to check if all days in the period 31/1/1979-31/12/2013 have been computed by QIBT.

5/2/19: updated script to check for completion on NCI runs. If result day is not present, OR if filesize
differs from exp01 (storm servers), then mark as incomplete. Note that os.path.getfilesize() gives a different
file size to bash output! Therefore get file sizes from bash and then compare.
To fetch this file from hdrive: 
scp z3131380@hurricane.ccrc.unsw.edu.au:/home/z3131380/hdrive/PhD/PartB/Scripts/Processing/Check_all_days_computed.py /home/603/cxh603/PhD/PartB/Scripts/Processing/

@author: z3131380
"""
import pandas as pd
import os.path
import numpy as np

#==============================================================================
# Get file sizes from exp01 (storm servers) to compare with exp02 on NCI. 
#==============================================================================
# (1) Run on storm servers
#ls -lh /srv/ccrc/data19/z3131380/PartB/Output/Australia/100parcels/TS10min/exp01/ > exp01_results_list.txt

# (2) Fetch result list from storm servers
#scp z3131380@hurricane.ccrc.unsw.edu.au:/srv/ccrc/data19/z3131380/PartB/Output/Australia/100parcels/TS10min/exp01/exp01_results_list.txt /home/603/cxh603/PhD/PartB/Scripts/Processing/

# (3) Run on raijin
#ls -lh /g/data/xc0/user/Holgate/QIBT/exp02 > exp02_results_list.txt

#==============================================================================
# Mark as incomplete if file not present OR file size differs from exp01
#==============================================================================
# (4) Run this python script on raijin
Start_date = '19790131' ; Start_year = Start_date[0:4] ; Start_month = Start_date[4:6]; Start_day = Start_date[6:8]
End_date = '20131231' ; End_year = End_date[0:4] ; End_month = End_date[4:6]; End_day = End_date[6:8]
daylist = pd.date_range(Start_date,End_date,freq='d') 
#dirname = '/srv/ccrc/data19/z3131380/PartB/Output/Australia/100parcels/TS10min/exp01/'
dirname = '/g/data/xc0/user/Holgate/QIBT/exp03/'

#exp01_results=pd.read_fwf('/home/603/cxh603/PhD/PartB/Scripts/Processing/exp01_results_list.txt')
#exp01_results.columns = ['Permissions','No.','ID','Place','Filesize','Mon','Day','Year','Name']
#
#exp02_results=pd.read_fwf('/g/data/xc0/user/Holgate/QIBT/exp02/exp02_results_list.txt')
#exp02_results.columns = ['Permissions','No.','ID','Place','Filesize','Mon','Day','Year','Name']

df = pd.DataFrame()  
for i in range(len(daylist)):
    dd = str(daylist[i].day)
    mm = '%02u' % daylist[i].month
    yyyy = str(daylist[i].year)
    fname = 'bt.'+yyyy+mm+'_'+dd+'.nc'
    file = dirname+fname
    if os.path.isfile(file) == False:
        df = df.append(pd.DataFrame({'Day':[fname]}))
#    elif exp02_results.iloc[np.where(fname==exp02_results['Name'].values)[0][0]]['Filesize'] != exp01_results.iloc[np.where(fname==exp01_results['Name'].values)[0][0]]['Filesize']:
#        df = df.append(pd.DataFrame({'Day':[fname]}))
outname = dirname+'QIBT_exp03_incomplete_days.txt'        
df.to_csv(outname)              
# (5) Copy the outfile to /home/z3131380/hdrive/PhD/PartB/Scripts/Model
# (6) Create .pbs scripts for incomplete days list in outfile.