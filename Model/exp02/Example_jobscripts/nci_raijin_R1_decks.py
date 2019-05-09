#!/usr/bin/python

#-------------------------------------------------------------------------------
# PURPOSE
#=========
# Reads a master deck on magnus, and creates many monthly decks by substituting variables in the
# master deck. These monthly decks are then used to run WRF on magnus. These scripts
# use boundary conditions files already created and stored on the NCI system.
# This deck works for models that use leap years, and for analysis.
# This deck works for putting output to massdata storage system
#
# STATUS
#========
# Working for CCRC-WRF3.6.0.5 ERA-Interim R1 no-nudging run
# (for both Standard and projected SEB versions)
# (Do final checks on WRF output)

# METHOD
#=========
# The monthly decks are created by substituting the text in the master deck framed by % signs
# with values of variables in this script (template texts method)
# It is possible to add lines in the template csh script that use these template texts.
# The template text in these new lines will be automatically replaced. Here is a list of
# the template texts supported:
#
# Note: in the following, the date (month/year) for which we create the csh
# script is called "simulated date".
#
#%ccrcdir%:      replaced by the value of CCRC_dir
#%expdir%:       replaced by the value of exp_dir
#%ccrcuser%:     replaced by the value of CCRC_user
#%email%:        replaced by the value of email
#%project%:      replaced by the value of project
#
#%syear%:        replaced by the year of the simulated date
#%smonth%:       replaced by the month of the simulated date
#%sday%:         replaced by the first day of the simulated date
#%nyear%:        replaced by the year of the simulated date + 1 month
#%nmonth%:       replaced by the month of the simulated date + 1 month
#%eday%:         replaced by the last day of the simulation
#%rundays%:      replaced by the number of days of the simulated date.
#
#%isrestart%:    replaced by .false. if first = 1, by .true. else.
#%resint%:       replaced by 1440*rundays
#
# List of template texts with very specific meanings that should be treated carefully:
#%sexist%:       start of an optional section
#%eexist%:       end of an optional section
#%spart%:        start of an optional section
#%epart%:        end of an optional section

# NOTE: %sexist%-%eexist% lines are not printed for the first piece of the first month 
# of a new run from scratch

# Each created script (runwrf_YEAR_MONTH_DAY.deck) will create and launch 3 dependant
# PBS scripts:
# 1st to grab the boundary condition files from CCRC machine
# 2nd to create the namelist.input and run WRF
# 3rd to send the output to massdata system on raijin
#
# USAGE
#========
# 1. Update the namelist part in runwrf_nci_raijin_*.deck. It must be the same
# as what you used for creating the boundary condition files.
# 2. Create output /massdata directories explicitly USING THE SAME PROJECT
#  as in this script
# 3. Run this script
#
# 4. To run the resulting monthly scripts put them into the WRF/run directory on raijin (nci)
# and type ./runwrf_yyyy_mm_dd.deck
#
# INPUT
#========
# MANDATORY
# start_month/start_year:First month of the simulation. Will start at day 1.
# end_month/end_year:    Last month of the simulation
# first:                 Used to determine if this is a run started from scratch
#                        or a continuation run (i.e. restart files already
#                        exist for previous dates for the same WRF simulation)
# BDY_user:              Your user name at CCRC and the machine address you want to use at CCRC for
#                        boundary files
# BDY_dir/RST_dir/OUT_dir: Paths to boundary, restart, and output files, respectively. The directories
#                        should NOT end with a slash! The restart and output directories are for
#                        /massdata system.
# email:                 Your email address.
# project:               NCI project number to use for these simulations
# indeck:                Name for the master deck

# There are other optional inputs. They are not so important and are not discussed
#
# OUTPUT
#=========
# Monthly decks

# HISTORY
#=========
# Written by ?
# 17 03 14 => Commented by Roman Olson, CCRC
# 26 04 14 => Modified for decks that output to massdata system on raijin
# 28 05 2014 => Modified to work with projected CCRC-WRF3.6 no-nudging
# 19 Aug 2014 => Updated output directories to work with CCRC-WRF3.6.0.5
# 05 Dec 2014 => Updated to go to the end of 2013
#==============================================================================


import os

#======================================================================================
# INPUT
#======================================================================================

# Start month of the simulation. Will start at day 1.
start_month = 11 
start_year = 2005 

# End month of the simulation (included). 
end_month = 12
end_year = 2040

# If starting from scratch (not a continuation run), put first=1 else put first=-1
first = 1

#name the input deck to use
indeck = "runwrf_pawsey_magnus_R1.deck"

#how many pieces should the month long run be divided up into?
#Only 1 works if to have monthly runs with restarts on 00:00 the 1st day of the next month.
per_month = 1

#username on system and address of the machine containing the bdy files
#143.119.252.237 is ip address for cluster.irs.environment.nsw.gov.au
BDY_user = ""
#Path containing the boundary files.
#BDY_dir = "/g/data3/w28/jpe561/CORDEX/ACCESS1.0/historical/bdy"
BDY_dir = "/g/data3/w28/jpe561/CORDEX/ACCESS1.0/RCP8.5/bdy"
#scp flags required e.g. for port 6512 need "-P 6512"
BDY_scpflags = " "

#username for restart files
RST_user = BDY_user
#Path containing the restart files.
RST_dir = "/g/data3/hh5/tmp/w28/jpe561/CORDEX/ACCESS1.0/RCP8.5/R1/restart"
#scp flags required e.g. for port 6512 need "-P 6512"
RST_scpflags = ""

#username for moving out files
OUT_user = BDY_user
#Path containing the output files".
OUT_dir = "/g/data3/hh5/tmp/w28/jpe561/CORDEX/ACCESS1.0/RCP8.5/R1/out"
#scp flags required e.g. for port 6512 need "-P 6512"
OUT_scpflags = " "


#set the project number.
project = "n81"
#project = "w28"

#email: to receive an email at the end of each script
email = "shaoxiu.ma@unsw.edu.au"
#=========================================================================================

year = start_year
month = start_month

days_in_month = (31,28,31,30,31,30,31,31,30,31,30,31)

#set up loop to be over all months in the years mentioned
while (year < end_year or (year == end_year and month <= end_month)):
    
    #get the month as a 2 digit string
    monthstr = str(month).rjust(2,"0")
    
    
    #get number of days in each piece
    numdays = int(days_in_month[month-1]/per_month)

    # To manage optional parts
    getjobid = True
    
    #set up loop over the pieces of the month
    for dd in range(1,per_month+1):
        
        #get the start day for this run
        sday = ((dd - 1) * numdays) + 1
  
        #get the end day for this run
        if (dd < per_month):
            edaystr = str((dd * numdays) + 1).rjust(2,"0")
            nextyear = year
            nextmonthstr = monthstr
        else:
            edaystr = "01"

            numdays = (days_in_month[month-1] - sday) + 1
            
            #adjust for leap year
            if (month == 2):
                if ((year/4. == int(year/4.)) and ((year/100. != int(year/100.)) or (year/400.) == int(year/400.))):
                    numdays = numdays + 1
            
            if (month == 12):
                nextmonthstr = "01"
                nextyear = year+1
            else:
                nextmonthstr = str(month+1).rjust(2,"0")
                nextyear = year
                
        #get the day as a 2 digit string
        sdaystr = str(sday).rjust(2,"0")

        #open the sample deck 
        fin = open (indeck,"r") 

        #open the deck I am creating
        fout = open ("runwrf_%s_%s_%s.deck"%(year,monthstr,sdaystr),"w")
        
	# make sure the output file is executable by user
	#os.fchmod(fout.fileno(),0744)

        # To manage optional parts
        yesprint = True
        partprint = True
        existprint = True
        
        # Is there a restart file?
        isrestart = (first != 1) or (dd != 1)

        #Loop over the lines of the input file
        for lines in fin.readlines():

            # to manage optional parts
            if ((lines == "%spart%\n") and (dd > 1)):
                partprint = False

            if (lines == "%epart%\n"):
                partprint = True
            
	    if (lines == "%s1part%\n") and (dd == 1):
	      partprint = False
	      
	    if (lines == "%e1part%\n"):
              partprint = True
	    
            if ((lines == "%sexist%\n") and (dd == 1) and (first == 1)):
                existprint = False

            if (lines == "%eexist%\n"):
                existprint = True

            yesprint = partprint and existprint
            
            # Replace template fields by values
            lines = lines.replace("%BDYdir%", BDY_dir)
            lines = lines.replace("%BDYuser%", BDY_user)
	    lines = lines.replace("%BDYscpflags%", BDY_scpflags)
            lines = lines.replace("%RSTdir%", RST_dir)
            lines = lines.replace("%RSTuser%", RST_user)
	    lines = lines.replace("%RSTscpflags%", RST_scpflags)
            lines = lines.replace("%OUTdir%", OUT_dir)
            lines = lines.replace("%OUTuser%", OUT_user)
	    lines = lines.replace("%OUTscpflags%", OUT_scpflags)
            lines = lines.replace("%email%", email)
            lines = lines.replace("%project%", project)
            lines = lines.replace("%syear%", str(year))
            lines = lines.replace("%smonth%", monthstr)
            lines = lines.replace("%sday%", sdaystr)
            lines = lines.replace("%rundays%", str(numdays))
            lines = lines.replace("%nyear%", str(nextyear))
            lines = lines.replace("%nmonth%", nextmonthstr)
            lines = lines.replace("%eday%", edaystr)

            lines = lines.replace("%isrestart%", "."+str(isrestart).lower()+".")
            lines = lines.replace("%resint%", str(numdays*1440))
             
	    lines = lines.replace("%runhours%", str(numdays*24))
	    lines = lines.replace("%run3hours%", str(numdays*8))
           
            lines = lines.replace("%spart%", "")
            lines = lines.replace("%epart%", "")
            lines = lines.replace("%s1part%", "")
            lines = lines.replace("%e1part%", "")
	    lines = lines.replace("%sexist%","")
            lines = lines.replace("%eexist%","")

            if (yesprint):
                fout.write(lines)

        first = 0

        #Close input and output files
        fin.close()
        fout.close()
        
    #increment to next month
    if (month <= 11):
        month = month + 1
    else:
        month = 1
        year = year+1
            
