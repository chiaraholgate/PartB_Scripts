# -*- coding: utf-8 -*-
"""
Created on Fri Apr 12 10:36:11 2019

This script is called from a bash shell script that first uses grep to idenfify those .pbs.o*  files in 
/home/603/cxh603/PhD/PartB/Scripts/Model/exp03 that have a 271 exit code.
(Most jobs are killed due to exceeding walltime, but some are from insufficient JOBFS memory).

Each .pbs script is then opened and its walltime and jobfs is increased.

@author: z3131380
"""
#==============================================================================
# First attempt to re-run killed jobs:
#==============================================================================
dir = r'/home/603/cxh603/PhD/PartB/Scripts/Model/exp03/'
incomplete_weeks = open(dir+'Incomplete_weeks_exp03_Try1.txt','r')
files = incomplete_weeks.readlines()
incomplete_weeks.close()

pbs_to_submit = open(dir+"pbs_to_submit.txt","w")

for f in files:
    fname = f[:-10]
    if "Walltime requested: 12:00:00" in  open(dir+f.rstrip("\n")).read():
        s = open(dir+'pbs_scripts/'+fname).read()
        s = s.replace('12:00:00', '24:00:00')
        fh = open(dir+'pbs_scripts/'+fname,'w')
        fh.write(s)
        fh.close()
    if "Walltime requested: 06:00:00" in  open(dir+f.rstrip("\n")).read():
        s = open(dir+'pbs_scripts/'+fname).read()
        s = s.replace('6:00:00', '12:00:00')
        fh = open(dir+'pbs_scripts/'+fname,'w')
        fh.write(s)
        fh.close()
    if "Walltime requested: 06:00:00" in  open(dir+f.rstrip("\n")).read():
        s = open(dir+'pbs_scripts/'+fname).read()
        s = s.replace('06:00:00', '12:00:00')
        fh = open(dir+'pbs_scripts/'+fname,'w')
        fh.write(s)
        fh.close()
    if "Walltime requested: 24:00:00" in  open(dir+f.rstrip("\n")).read():
        s = open(dir+'pbs_scripts/'+fname).read()
        s = s.replace('24:00:00', '48:00:00')
        fh = open(dir+'pbs_scripts/'+fname,'w')
        fh.write(s)
        fh.close()
    if "JobFS requested:    900.0MB" in  open(dir+f.rstrip("\n")).read():
        s = open(dir+'pbs_scripts/'+fname).read()
        s = s.replace('900MB', '3GB')
        fh = open(dir+'pbs_scripts/'+fname,'w')
        fh.write(s)
        fh.close()
    if "JobFS requested:    500.0GB" in  open(dir+f.rstrip("\n")).read():
        s = open(dir+'pbs_scripts/'+fname).read()
        s = s.replace('500GB', '3GB')
        fh = open(dir+'pbs_scripts/'+fname,'w')
        fh.write(s)
        fh.close()
    pbs_to_submit.write(fname+"\n") 
pbs_to_submit.close() 
    


#==============================================================================
# Second attempt
#==============================================================================
#incomplete_weeks = open(dir+'Incomplete_weeks_exp03_Try1.txt','r')
#try1 = incomplete_weeks.readlines()
#incomplete_weeks.close()
#incomplete_weeks = open(dir+'Incomplete_weeks_exp03_Try2.txt','r')
#try2 = incomplete_weeks.readlines()
#incomplete_weeks.close()
#Remaining weeks = [x for x in try2 if x not in try1]




#==============================================================================
# ####
#==============================================================================
#fin = open(dir+'test.txt','r')
#fout = open(dir+'test2.txt','w')
#for lines in fin.readlines():
#    lines = lines.replace("gday","mate")
#    fout.write(lines)
#fin.close()
#fout.close()
####




#for f in files:
#    fname = f[:-10]
#    if "JobFS requested:    900.0MB" in  open(dir+f.rstrip("\n")).read():
#        print fname,' was 900mb'
#    elif "JobFS requested:    500.0GB" in  open(dir+f.rstrip("\n")).read():
#        print fname,' was 500gb'
#    elif "JobFS requested:    3.0GB" in  open(dir+f.rstrip("\n")).read():
#        print fname,' was 3gb'
#    else:
#        print fname,'was other'


