! These routines are used to perform a quasi-isentropic back trajectory
! analysis of water vapor sources based on output from an MM5 run.
!
!here I am d&b but with w not istentropes


! NOTE: INPUT DATA VERTICAL LEVELS ARE THE REVERSE OF WHAT IS EXPECTED BY THE FUNCTIONS.
! MODEL EXPECTS VERTICAL LEVELS TO RANGE FROM HIGHEST ATMOSPHERIC LEVEL TO LOWEST.
! CURRENT INPUT DATA GIVES OPPOSITE - LEVELS REVERSED WHEN LOADING INPUT DATA.
! MODEL EXPECTS INPUT DATA RANGE BETWEEN -180deg and +180deg.


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! EXPLANATION OF TIME-RELATED VARIABLES, AND TESTING VALUES
! "set"   = user defined values
! "calcd" = calculated within program

! totdays = 1                                            ! number of days to run simulation forward, based on set start/end dates(calcd)
! totbtadays = 3                                         ! number of days to back-track for (set)
! tstep = 30                                             ! number of minutes for back-track time step (set)
! daytsteps = 1440/tstep = 48                            ! number of simulation time steps in a day (calcd)
! totsteps = daytsteps*(totbtadays+1) = 192             ! total number of simulation time steps over period (calcd)
! MM5tstep = 180                                         ! input file time step in minutes (set)
! MM5daysteps = 1440/MM5tstep = 8                        ! number of input file time steps in a day (calcd)
! inMM5tsteps = MM5tstep/tstep = 6                       ! number of simulation time steps per input time step  (calcd)
! MM5totsteps = (MM5daysteps*(totbtadays+1)) + 1 = 33   ! total number of input time steps over the back-track period (calcd)
!.............we want the event day + totbtadays before it + the time step after ! (ie 0 hour time step which is the last in each file)
! ttMM5day = = ((tt-1)/inMM5tsteps) + 1 = 1              ! position of parcel time step in input file (calcd)
! ttMM5 = MM5totsteps - MM5daysteps - 1 + ttMM5day = 17  ! input file time step from the beginning of the loaded files (calcd)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


MODULE global_data

IMPLICIT NONE

SAVE

!
!*******user modified variables**********************
!
INTEGER :: sday,smon,syear    !start day for calculations
INTEGER :: edday,edmon,edyear  !end day for calculations (must be at least one day after start day)

!INTEGER :: sday = 1,smon = 2,syear = 2000    !start day for calculations
!INTEGER :: edday = 1,edmon = 3,edyear = 2000  !end day for calculations (must be at least one day after start day)
INTEGER :: totdays
INTEGER, PARAMETER :: totbtadays = 30   !number of days of data to keep for bta; i.e. how far back in time to calc.
                                       !must be less than days you have input data for

!INTEGER, PARAMETER :: totdays = 30    !total number of days to calculate for; must be <= no. days in days_of_rain.txt
INTEGER, PARAMETER :: tstep = 10   !number of minutes for back trajectory time step (simultion time step)
                      !must divide evenly into number of minutes in day 1440 and number of minutes in MM5 time step (here 180)
INTEGER, PARAMETER :: nparcels = 100   !set the number of parcels to release if it rains
REAL, PARAMETER :: minpre = 2   !min daily precip to deal with (mm)

INTEGER, PARAMETER :: bdy = 6   !boundary layers to ignore; trajectories will be tracked to this boundary

CHARACTER(LEN=50), PARAMETER :: diri = "/g/data/hh5/tmp/w28/jpe561/back_traj/"   
CHARACTER(LEN=100), PARAMETER :: diro = "/short/xc0/cxh603/QIBT_testing/Month/"  
CHARACTER(LEN=100), PARAMETER :: dirdata_atm = "/g/data/hh5/tmp/w28/jpe561/back_traj/wrfout/"
CHARACTER(LEN=100), PARAMETER :: dirdata_land = "/g/data/hh5/tmp/w28/jpe561/back_traj/wrfhrly/"  

INTEGER, PARAMETER :: numthreads = 16   !set the number of parallel openmp threads

!CHARACTER(LEN=50), PARAMETER :: fdaylist = "top300precip_days_min0.5.txt"   !file containing
!CHARACTER(LEN=50), PARAMETER :: fdaylist = "days_of_rain.txt"   !file containing list of days to do qibt on

LOGICAL, PARAMETER :: peak = .FALSE.	!does the daylist indicate storm peaks (TRUE) or whole days (FALSE)

LOGICAL, PARAMETER :: wshed = .TRUE. !only calculate trajectories for watershed

CHARACTER(LEN=50), PARAMETER :: fwshed = "NARCliM_AUS_land_sea_mask.nc"
								!"n_s_basin/sbasin_rotpole.nc"
								!"mdb_rotpole.nc"
								!"n_s_basin/moonie_rotpole.nc"
                                !set to "" if no watershed
                                !0 outside watershed, >0 inside

REAL, PARAMETER :: min_del_q = 0.0001    !the minimum change in parcel mixing ratio (kg/kg) to be considered "real"

LOGICAL, PARAMETER :: eachParcel = .TRUE.   !output the data along the trajectory of each parcel
!1
!****************************************************
!

INTEGER :: daytsteps,totsteps,inMM5tsteps,MM5daysteps,MM5totsteps
INTEGER :: dim_i,dim_j,dim_k,fdim_i,fdim_j,ssdim
!INTEGER :: sday,smon,syear
INTEGER :: mon,year,dd,totpts
!REAL :: day
INTEGER :: day


REAL, PARAMETER :: Lv = 2.25E6   !latent heat of vaporization of water (Jkg-1)
REAL, PARAMETER :: g = 9.8     !gravity (m.s-2)
REAL, PARAMETER :: P0 = 100000   !reference surface pressure (Pa)
REAL, PARAMETER :: Rd = 287.053   !ideal gas constant for dry air (J/kgK)
REAL, PARAMETER :: Cp = 1004.67   !heat capacity of air at constant pressure (J/kgK)
REAL, PARAMETER :: Rv = 461.5    !gas constant of water vapor (J/kgK)
REAL, PARAMETER :: Cl = 4400     !heat capacity of liquid water at ~-20C (J/kgK)
REAL, PARAMETER :: pi = 3.14159265
REAL, PARAMETER :: deg_dist = 111.   !average distance of 1 degree lat is assumed to be 111km

COMMON /global_vars/ daytsteps,totsteps,inMM5tsteps,MM5daysteps,MM5totsteps, &
		dim_i,dim_j,dim_k,sday,smon,syear,mon,year,day,dd,totpts, &
		fdim_i,fdim_j,ssdim 

!$OMP THREADPRIVATE(/global_vars/)

END MODULE global_data


!**************************************************************
!****************************************************************

MODULE bt_subs

IMPLICIT NONE

CONTAINS



SUBROUTINE handle_err(status)
!---------------------------------
! handle any errors from the fortran 90 netCDF interface
!---------------------------------
  
USE netcdf
  
IMPLICIT NONE
   
integer, intent (in) :: status

print *,'status = ',status

if(status /= nf90_noerr) then
  print *, trim(nf90_strerror(status))
  stop "Stopped"
end if
  
END SUBROUTINE handle_err

!***********************************************************************

INTEGER FUNCTION string_to_int(string)
!------------------------------------------
! converts a character string to an integer
!--------------------------------------------

IMPLICIT NONE

character (len=*), intent(in) :: string

! local constant
integer, parameter :: zero = iachar("0")
integer :: i,  sign, integ
character (len=50) :: str

str = trim(string)
integ = 0

select case (str(1:1))
case ("-")
  sign = -1
  str = str(2:)
case ("+")
  sign = 1
  str = str(2:)
case ("0":"9")
  sign = 1
end select

do i=len(trim(str)),1,-1
  select case (str(i:i))
    case ("0":"9")
      integ = integ + (iachar(string(i:i))-zero)*10**(len(trim(str))-i)
    case default
      print *, "cannot convert a non-integer character to an integer!!"
      return
  end select
end do

string_to_int = integ

end FUNCTION string_to_int

!******************************************************************************

CHARACTER(LEN=50) FUNCTION int_to_string(integ)
!----------------------------------------------
! converts an integer to a character string 
!----------------------------------------------

IMPLICIT NONE

integer, intent(in) :: integ

! local constant
integer, parameter :: zero = iachar("0")
character (len=50) :: str
integer :: i, inte

str="                                                  "
inte = integ

do i=1,50
  str(50-i:50-i) = achar(mod(abs(inte),10)+zero)
  inte = int(inte/10)
  if (abs(inte)<1) exit
end do

if (integ<0) str(50-i-1:50-i) = "-"

int_to_string = adjustl(str)

end FUNCTION int_to_string

!************************************************************************


CHARACTER(LEN=50) FUNCTION real_to_string(num)
!----------------------------------------------
! converts an real to a character string 
! here I allow a maximum of 10 decimal places
!----------------------------------------------

IMPLICIT NONE

!real, intent(in) :: num
INTEGER, intent(in) :: num

! local 
integer :: i, whole, frac
real :: fracnum


whole = num !AINT(num)
fracnum = num - whole

do i=1,10
  if (MOD(fracnum,1.)==0) exit
  fracnum = fracnum * 10.
end do

frac = AINT(fracnum)

real_to_string = TRIM(int_to_string(whole))//"."//TRIM(int_to_string(frac))

end FUNCTION real_to_string

!************************************************************************

INTEGER FUNCTION days_in_month(mon,year)
!-------------------------------------------
! how many days in this month?
!-----------------------------------------

IMPLICIT NONE

INTEGER, INTENT(IN) :: mon,year

INTEGER, DIMENSION(12), PARAMETER :: months = (/31,28,31,30,31,30,31,31,30,31,30,31/)

days_in_month = months(mon)

if (mon == 2 .AND. MOD(year,4)==0 ) then
  days_in_month = 29

  if (MOD(year,100)==0 .AND. MOD(year,400)/=0) then
    days_in_month = 28
  end if
  
end if

END  FUNCTION days_in_month

!**************************************************************************

INTEGER FUNCTION days_in_year(year)
!-------------------------------------------
! how many days in this year?
!-----------------------------------------

IMPLICIT NONE

INTEGER, INTENT(IN) :: year

if (MOD(year,4)==0) then
  days_in_year = 366
else
  days_in_year = 365
end if

END  FUNCTION days_in_year

!**********************************************************************

INTEGER FUNCTION julian(year,mon,day)

! From http://aa.usno.navy.mil/faq/docs/JD_Formula.php

IMPLICIT NONE
INTEGER, INTENT(IN) :: year,mon,day

julian= day-32075+1461*(year+4800+(mon-14)/12)/4+367*(mon-2-(mon-14)/12*12)/12-3*((year+4900+(mon-14)/12)/100)/4
   
END FUNCTION julian

!**********************************************************************

SUBROUTINE GREGORIAN(JD,YEAR,MONTH,DAY)!,HOUR,MINUTE,SECOND)

! From http://aa.usno.navy.mil/faq/docs/JD_Formula.php

IMPLICIT NONE
!REAL, INTENT(IN) :: JD
INTEGER :: JD
INTEGER, INTENT(OUT) :: YEAR, MONTH, DAY!, HOUR, MINUTE, SECOND
REAL :: JT
INTEGER :: I,J,K,L,N

L = INT(JD)+68569!JD= K-32075+1461*(I+4800+(J-14)/12)/4+367*(J-2-(J-14)/12*12)/12-3*((I+4900+(J-14)/12)/100)/4
N = 4*L/146097
L = L-(146097*N+3)/4
I = 4000*(L+1)/1461001
L = L-1461*I/4+31
J = 80*L/2447
K = L-2447*J/80
L = J/11
J = J+2-12*L
I = 100*(N-49)+I+L

YEAR = I
MONTH = J
DAY = K

! JT = DMOD(JD,1.D0)*24.D0
! HOUR = INT(JT)
! JT = DMOD(JT,1.D0)*60.D0
! MINUTE = INT(JT)
! JT = DMOD(JT,1.D0)*60.D0
! SECOND = NINT(JT)
! 
! IF (SECOND == 60) THEN
! 	SECOND = SECOND-60
! 	MINUTE = MINUTE+1
! END IF



END SUBROUTINE GREGORIAN

!**********************************************************************

INTEGER FUNCTION simlength(startday,startmon,startyear,endday,endmon,endyear)
!-----------------------------------------------------
! given the start and end dates of the desired simulation
! period (top of global data), calculate the number of days 
! in the period

! Functions: don't need to specify what comes out, e.g. fn_name(in,in,in)
! Subroutines: do specify in and out, e.g. sbrtn_name(in,in,in,out,out)

!-----------------------------------------------------
USE global_data

IMPLICIT NONE

INTEGER, intent(in) :: startday,startmon,startyear,endday,endmon,endyear
INTEGER :: start_jd, end_jd
! 
! INTEGER :: first_month,this_month,next_months,mm,yy,no_years,this_month_newyear,next_months_newyear
! 
! if (startyear==endyear .AND. startmon==endmon) then
! first_month = endday-startday
! else if (startday==1) then
!   first_month = days_in_month(startmon,startyear)
! else 
!   first_month = days_in_month(startmon,startyear) - startday + 1
! end if 
! 
! next_months=0
! next_months_newyear=0
! 
! if (startyear==endyear) then
!   !next_months = 0
!   do mm = startmon+1,endmon
!     this_month = days_in_month(mm,endyear)
!     next_months =  next_months + this_month 
!   end do
! else 
!   ! add rest of first year
!   do mm = startmon+1,12 
!     this_month = days_in_month(mm,startyear)
!     next_months =  next_months + this_month 
!   end do
!   ! add next years
!   no_years = endyear - startyear
!   !next_months_newyear = 0
!   if (no_years==1) then    
!     do mm=1,endmon
!       if (mm==endmon) then
!         next_months_newyear = next_months_newyear + endday
!       else
!       this_month_newyear = days_in_month(mm,endyear)
!       next_months_newyear =  next_months_newyear + this_month_newyear 
!       end if
!     end do
!   else
!     !next_months_newyear = 0
!     do yy=1,no_years
!       if (yy==no_years) then
!         do mm=1,endmon
!           if (mm==endmon) then
!             next_months_newyear = next_months_newyear + endday
!           else
!             this_month_newyear = days_in_month(mm,endyear)
!             next_months_newyear =  next_months_newyear + this_month_newyear 
!           end if
!         end do
!       else
!         next_months_newyear = next_months_newyear + days_in_year(startyear+yy)
!       end if
!     end do
!   end if
! end if
!     
! simlength = first_month + next_months + next_months_newyear
! print *,simlength

start_jd = julian(startyear,startmon,startday)
end_jd = julian(endyear,endmon,endday)
simlength = end_jd - start_jd
!print *,'simlength=',simlength

END FUNCTION simlength
!**************************************************************************


!**************************************************************************
  
SUBROUTINE days_since_start(daynum)
!-----------------------------------------
!number of days since the start of the simulation
!-------------------------------------------------

USE global_data

IMPLICIT NONE

INTEGER, INTENT(OUT) :: daynum

INTEGER :: mm,yy

!
!make sure wanted start day is legit
!
if (year<syear.OR.(year==syear.AND.mon<smon).OR.(year==syear.AND.mon==smon.AND.day<sday)) then
  print *,"Wanted day is before start of simulation!!!"
  STOP
end if

!
!now count the days
!
daynum = 0

if (year==syear) then
  if (smon==mon) then
    daynum = day-sday
  else
    do mm = smon,mon-1
      daynum = daynum + days_in_month(mm,year)
    end do
    daynum = daynum - sday + day
  end if
  
else 

  do mm = smon,12
    daynum = daynum + days_in_month(mm,syear)
  end do
  daynum = daynum - sday
  
  if (year-syear>1) then
    do yy = syear+1,year-1
      do mm = 1,12
        daynum = daynum + days_in_month(mm,yy)
      end do
    end do
  end if
  
  if (mon==1) then
    daynum = daynum + day
  else
    do mm = 1,mon-1
      daynum = daynum + days_in_month(mm,yy)
    end do
    daynum = daynum + day
  end if
  
end if

print *,'L472, daynum=', daynum
  
END SUBROUTINE days_since_start

!**********************************************************************

SUBROUTINE day_month_year(simday)
!----------------------------------------------
! given the simulation day, calculate the corresponding month and year
! simulation day one is the first day
!-----------------------------------------------------

USE global_data

IMPLICIT NONE

!REAL,INTENT(IN) :: simday	!simulation day
INTEGER,INTENT(IN) :: simday	!simulation day

INTEGER :: days_so_far 
INTEGER :: mm,yr


days_so_far = 0

!calculate the no. of days in the first year
!testing to see if wanted day is in this year
do mm = smon,12
  days_so_far = days_so_far + days_in_month(mm,syear)
  
  if (mm==smon) days_so_far = days_so_far - sday + 1
  
  if (simday.le.days_so_far) then
    year = syear
    mon = mm
    day = simday - days_so_far + days_in_month(mm,syear) 
    return
  end if
end do
  

!loop through years to find correct month and year
yr = syear + 1
do
  do mm = 1,12
    days_so_far = days_so_far + days_in_month(mm,yr)
    if (simday.le.days_so_far) then
      year = yr
      mon = mm
      day = simday - days_so_far + days_in_month(mm,yr) - sday + 1
      
      !print *,simday,day,mon,year
      return
    end if
  end do
  yr = yr + 1
end do

END SUBROUTINE day_month_year

!**************************************************************************


SUBROUTINE get_filename(d,mn,yr,filename_ext_atm,filename_ext_RAIN,filename_ext_LH,filename_ext_P)
!-----------------------------------------------
! given the month and year get the filename extension string
!---------------------------------------------------

USE global_data

IMPLICIT NONE

INTEGER, INTENT(IN) :: d
INTEGER, INTENT(IN) :: mn, yr
CHARACTER(LEN=100), INTENT(OUT) :: filename_ext_atm,filename_ext_RAIN,filename_ext_LH,filename_ext_P

if (mn<10) then
  if (d<10) then
    filename_ext_atm = "wrfout_d01_"//TRIM(int_to_string(yr))//"-0"//TRIM(int_to_string(mn))//"-0"//TRIM(int_to_string(d))//"_00:00:00"
    filename_ext_RAIN = "wrfhrly_d01_"//TRIM(int_to_string(yr))//"-0"//TRIM(int_to_string(mn))//"-0"//TRIM(int_to_string(d))//"_00:00:00_RAIN.nc"
    filename_ext_LH = "wrfhrly_d01_"//TRIM(int_to_string(yr))//"-0"//TRIM(int_to_string(mn))//"-0"//TRIM(int_to_string(d))//"_00:00:00_LH.nc"
    filename_ext_P = "wrfhrly_d01_"//TRIM(int_to_string(yr))//"-0"//TRIM(int_to_string(mn))//"-0"//TRIM(int_to_string(d))//"_00:00:00_PSFC.nc"
  else 
    filename_ext_atm = "wrfout_d01_"//TRIM(int_to_string(yr))//"-0"//TRIM(int_to_string(mn))//"-"//TRIM(int_to_string(d))//"_00:00:00"
    filename_ext_RAIN = "wrfhrly_d01_"//TRIM(int_to_string(yr))//"-0"//TRIM(int_to_string(mn))//"-"//TRIM(int_to_string(d))//"_00:00:00_RAIN.nc"
    filename_ext_LH = "wrfhrly_d01_"//TRIM(int_to_string(yr))//"-0"//TRIM(int_to_string(mn))//"-"//TRIM(int_to_string(d))//"_00:00:00_LH.nc"
    filename_ext_P = "wrfhrly_d01_"//TRIM(int_to_string(yr))//"-0"//TRIM(int_to_string(mn))//"-"//TRIM(int_to_string(d))//"_00:00:00_PSFC.nc"
  end if
else
  if (d<10) then
    filename_ext_atm = "wrfout_d01_"//TRIM(int_to_string(yr))//"-"//TRIM(int_to_string(mn))//"-0"//TRIM(int_to_string(d))//"_00:00:00"
    filename_ext_RAIN = "wrfhrly_d01_"//TRIM(int_to_string(yr))//"-"//TRIM(int_to_string(mn))//"-0"//TRIM(int_to_string(d))//"_00:00:00_RAIN.nc"
    filename_ext_LH = "wrfhrly_d01_"//TRIM(int_to_string(yr))//"-"//TRIM(int_to_string(mn))//"-0"//TRIM(int_to_string(d))//"_00:00:00_LH.nc"
    filename_ext_P = "wrfhrly_d01_"//TRIM(int_to_string(yr))//"-"//TRIM(int_to_string(mn))//"-0"//TRIM(int_to_string(d))//"_00:00:00_PSFC.nc"
  else
    filename_ext_atm = "wrfout_d01_"//TRIM(int_to_string(yr))//"-"//TRIM(int_to_string(mn))//"-"//TRIM(int_to_string(d))//"_00:00:00"
    filename_ext_RAIN = "wrfhrly_d01_"//TRIM(int_to_string(yr))//"-"//TRIM(int_to_string(mn))//"-"//TRIM(int_to_string(d))//"_00:00:00_RAIN.nc"
    filename_ext_LH = "wrfhrly_d01_"//TRIM(int_to_string(yr))//"-"//TRIM(int_to_string(mn))//"-"//TRIM(int_to_string(d))//"_00:00:00_LH.nc"
    filename_ext_P = "wrfhrly_d01_"//TRIM(int_to_string(yr))//"-"//TRIM(int_to_string(mn))//"-"//TRIM(int_to_string(d))//"_00:00:00_PSFC.nc"
  end if
end if

filename_ext_atm = ADJUSTL(filename_ext_atm)
filename_ext_RAIN = ADJUSTL(filename_ext_RAIN)
filename_ext_LH = ADJUSTL(filename_ext_LH)
filename_ext_P = ADJUSTL(filename_ext_P)

! print *,'get_filename:'
! print *,'filename_ext_atm= ',filename_ext_atm 
! print *,'filename_ext_RAIN= ',filename_ext_RAIN
! print *,'filename_ext_LH= ',filename_ext_LH


END SUBROUTINE get_filename

!**********************************************************************

SUBROUTINE new_out_file(outncid,wvcid,wvc2id,xlocid,ylocid,dayid,opreid,daynum,lat2d,lon2d)
!----------------------------------------------
!create the output netcdf file and prepare it to accept data
!--------------------------------------------------------

USE global_data
USE netcdf

IMPLICIT NONE

INTEGER, INTENT(INOUT) :: outncid,wvcid,wvc2id,xlocid,ylocid,dayid,opreid
!REAL, INTENT(IN) :: daynum
INTEGER, INTENT(IN) :: daynum
REAL,INTENT(IN),DIMENSION(:,:) :: lat2d,lon2d

INTEGER :: status,jdimid,idimid,gwvcdimid,latid,lonid


!
!create the file
!
!differentiate whether we are doing whole days or around storm peaks
if (peak) then
 print *,'we are doing peaks here!'
  if (smon<10) then
    !status = nf90_create(TRIM(diro)//"bt."//TRIM(int_to_string(syear))//"0" &
    status = nf90_create(TRIM(diro)//"bt."//TRIM(int_to_string(year))//"0" &
  		!//TRIM(int_to_string(smon))//"_"//TRIM(real_to_string(daynum-1))// &
  		!//TRIM(int_to_string(smon))//"_"//TRIM(real_to_string(daynum))// &
  		//TRIM(int_to_string(mon))//"_"//TRIM(real_to_string(daynum))// &
		".nc",nf90_clobber,outncid)
    if (status /= NF90_NOERR) call handle_err(status)
  else
    !status = nf90_create(TRIM(diro)//"bt."//TRIM(int_to_string(syear)) &
    status = nf90_create(TRIM(diro)//"bt."//TRIM(int_to_string(year)) &
  		!//TRIM(int_to_string(smon))//"_"//TRIM(real_to_string(daynum-1))// &
  		!//TRIM(int_to_string(smon))//"_"//TRIM(real_to_string(daynum))// &
  		//TRIM(int_to_string(mon))//"_"//TRIM(real_to_string(daynum))// &
		".nc",nf90_clobber,outncid)
    if (status /= NF90_NOERR) call handle_err(status)
  end if
else
  if (smon<10) then
    !status = nf90_create(TRIM(diro)//"bt."//TRIM(int_to_string(syear))//"0" &
    status = nf90_create(TRIM(diro)//"bt."//TRIM(int_to_string(year))//"0" &
  		!//TRIM(int_to_string(smon))//"_"//TRIM(int_to_string(INT(daynum-0.5)))// &
  		!//TRIM(int_to_string(smon))//"_"//TRIM(int_to_string(INT(daynum)))// &
  		//TRIM(int_to_string(mon))//"_"//TRIM(int_to_string(INT(daynum)))// &
		".nc",nf90_clobber,outncid)
    if (status /= NF90_NOERR) call handle_err(status)
  else
    !status = nf90_create(TRIM(diro)//"bt."//TRIM(int_to_string(syear)) &
    status = nf90_create(TRIM(diro)//"bt."//TRIM(int_to_string(year)) &
  		!//TRIM(int_to_string(smon))//"_"//TRIM(int_to_string(INT(daynum-0.5)))// &
  		!//TRIM(int_to_string(smon))//"_"//TRIM(int_to_string(INT(daynum)))// &
  		//TRIM(int_to_string(mon))//"_"//TRIM(int_to_string(INT(daynum)))// &
		".nc",nf90_clobber,outncid)
    if (status /= NF90_NOERR) call handle_err(status)
  end if
end if

print *,'L582, outfile=',TRIM(diro)//"bt."//TRIM(int_to_string(year))//"0" &
  		//TRIM(int_to_string(mon))//"_"//TRIM(int_to_string(INT(daynum)))// &
		".nc"

!
!define dimensions
!
status = nf90_def_dim(outncid,"j_cross",dim_j,jdimid)
if (status /= NF90_NOERR) call handle_err(status)
status = nf90_def_dim(outncid,"i_cross",dim_i,idimid)
if (status /= NF90_NOERR) call handle_err(status)
status = nf90_def_dim(outncid,"gridcell_wvc",nf90_unlimited,gwvcdimid)
if (status /= NF90_NOERR) call handle_err(status)

!
!define the variable
!
status = nf90_def_var(outncid,"wv_cont",nf90_float,(/jdimid,idimid,gwvcdimid/),wvcid)
if (status /= NF90_NOERR) call handle_err(status)
status = nf90_def_var(outncid,"wv_cont_apbl",nf90_float,(/jdimid,idimid,gwvcdimid/),wvc2id)
if (status /= NF90_NOERR) call handle_err(status)
status = nf90_def_var(outncid,"x_loc",nf90_int,(/gwvcdimid/),xlocid)
if (status /= NF90_NOERR) call handle_err(status)
status = nf90_def_var(outncid,"y_loc",nf90_int,(/gwvcdimid/),ylocid)
if (status /= NF90_NOERR) call handle_err(status)
status = nf90_def_var(outncid,"day",nf90_float,(/gwvcdimid/),dayid)
if (status /= NF90_NOERR) call handle_err(status)
status = nf90_def_var(outncid,"pre",nf90_float,(/gwvcdimid/),opreid)
if (status /= NF90_NOERR) call handle_err(status)
status = nf90_def_var(outncid,"latitcrs",nf90_float,(/jdimid,idimid/),latid)
if (status /= NF90_NOERR) call handle_err(status)
status = nf90_def_var(outncid,"longicrs",nf90_float,(/jdimid,idimid/),lonid)
if (status /= NF90_NOERR) call handle_err(status)


!
!define attributes
!
status = nf90_put_att(outncid,wvcid,"long_name","Water Vapor Contribution")
if (status /= NF90_NOERR) call handle_err(status)
status = nf90_put_att(outncid,wvcid,"units","proportion of precipitation")
if (status /= NF90_NOERR) call handle_err(status)
status = nf90_put_att(outncid,wvcid,"num_boundary_layers",bdy)
if (status /= NF90_NOERR) call handle_err(status)
status = nf90_put_att(outncid,wvcid,"parcels_per_grid_point",nparcels)
if (status /= NF90_NOERR) call handle_err(status)

!status = nf90_put_att(outncid,wvc2id,"long_name","Water Vapor Contribution above PBL")
!if (status /= NF90_NOERR) call handle_err(status)
!status = nf90_put_att(outncid,wvc2id,"units","proportion of precipitation")
!if (status /= NF90_NOERR) call handle_err(status)
!status = nf90_put_att(outncid,wvc2id,"num_boundary_layers",bdy)
!if (status /= NF90_NOERR) call handle_err(status)
!status = nf90_put_att(outncid,wvc2id,"parcels_per_grid_point",nparcels)
!if (status /= NF90_NOERR) call handle_err(status)

status = nf90_put_att(outncid,xlocid,"long_name","x index location of precipitation (from 0)")
if (status /= NF90_NOERR) call handle_err(status)
status = nf90_put_att(outncid,ylocid,"long_name","y index location of precipitation (from 0)")
if (status /= NF90_NOERR) call handle_err(status)
status = nf90_put_att(outncid,dayid,"long_name","days since "// &
	TRIM(int_to_string(sday))//"/"//TRIM(int_to_string(smon))//"/"// &
	TRIM(int_to_string(syear)))
if (status /= NF90_NOERR) call handle_err(status)
status = nf90_put_att(outncid,opreid,"long_name","precipitation")
if (status /= NF90_NOERR) call handle_err(status)
status = nf90_put_att(outncid,opreid,"units","mm")
if (status /= NF90_NOERR) call handle_err(status)
status = nf90_put_att(outncid,latid,"long_name","LATITUDE (SOUTH NEGATIVE)")
if (status /= NF90_NOERR) call handle_err(status)
status = nf90_put_att(outncid,latid,"units","degrees")
if (status /= NF90_NOERR) call handle_err(status)
status = nf90_put_att(outncid,lonid,"long_name","LONGITUDE (WEST NEGATIVE)")
if (status /= NF90_NOERR) call handle_err(status)
status = nf90_put_att(outncid,lonid,"units","degrees")
if (status /= NF90_NOERR) call handle_err(status)


!
!leave define mode
!
status = nf90_enddef(outncid)


status = nf90_put_var(outncid,latid,lat2d,start=(/1,1/),count=(/dim_j,dim_i/))
if(status /= nf90_NoErr) call handle_err(status)
status = nf90_put_var(outncid,lonid,lon2d,start=(/1,1/),count=(/dim_j,dim_i/))
if(status /= nf90_NoErr) call handle_err(status)
	


END SUBROUTINE new_out_file

!****************************************************************************

SUBROUTINE open_netcdf_files(ncid,prencid,lhncid,psfcncid,preid,lhid,uid,vid,wid,tid,qid,ppid,pbid,pblid,psfcid,filename_ext_atm,filename_ext_RAIN,filename_ext_LH,filename_ext_P)
!----------------------------------------------------------------
! open all the netcdf data files and get the variable ids
!------------------------------------------------------------------

USE netcdf
USE global_data

IMPLICIT NONE

INTEGER, INTENT(OUT) :: ncid,prencid,lhncid,psfcncid !,uncid,vncid,wncid,tncid,qncid,ppncid,pblncid
INTEGER, INTENT(OUT) :: preid,lhid,uid,vid,wid,tid,qid,ppid,pbid,pblid,psfcid 
CHARACTER(LEN=100), INTENT(IN) :: filename_ext_atm,filename_ext_RAIN,filename_ext_LH,filename_ext_P

INTEGER :: status


! open the netcdf files - ATMOSPHERIC VARIABLES
status = NF90_OPEN(TRIM(dirdata_atm)//TRIM(filename_ext_atm), NF90_NOWRITE, ncid)
if (status /= NF90_NOERR) call handle_err(status)

! open the netcdf files - PRECIP 
status = NF90_OPEN(TRIM(dirdata_land)//TRIM(filename_ext_RAIN), NF90_NOWRITE, prencid)
if (status /= NF90_NOERR) call handle_err(status)

! open the netcdf files - EVAP 
status = NF90_OPEN(TRIM(dirdata_land)//TRIM(filename_ext_LH), NF90_NOWRITE, lhncid)
if (status /= NF90_NOERR) call handle_err(status)

! open the netcdf files - SURFACE PRESSURE
status = NF90_OPEN(TRIM(dirdata_land)//TRIM(filename_ext_P), NF90_NOWRITE, psfcncid)
if (status /= NF90_NOERR) call handle_err(status)

!
!get ids for each variable
!
status = nf90_inq_varid(prencid, "RAIN", preid) 	! [mm]
if(status /= nf90_NoErr) call handle_err(status)
status = nf90_inq_varid(lhncid, "LH", lhid)			! [Wm-2] > converted to mm at end of get_data subroutine
if(status /= nf90_NoErr) call handle_err(status)
status = nf90_inq_varid(ncid, "U", uid)			! [ms-1]
if(status /= nf90_NoErr) call handle_err(status)
status = nf90_inq_varid(ncid, "V", vid)			! [ms-1]
if(status /= nf90_NoErr) call handle_err(status)
status = nf90_inq_varid(ncid, "W", wid)			! [ms-1]
if(status /= nf90_NoErr) call handle_err(status)
status = nf90_inq_varid(ncid, "T", tid) 			! [K]
if(status /= nf90_NoErr) call handle_err(status)
status = nf90_inq_varid(ncid, "QVAPOR", qid)			! [kgkg-1]
if(status /= nf90_NoErr) call handle_err(status)
status = nf90_inq_varid(ncid, "P", ppid) 		! [Pa]
if(status /= nf90_NoErr) call handle_err(status)
status = nf90_inq_varid(ncid, "PB", pbid) 		! [Pa]
if(status /= nf90_NoErr) call handle_err(status)
status = nf90_inq_varid(ncid, "PBLH", pblid)		! [m]
if(status /= nf90_NoErr) call handle_err(status)
status = nf90_inq_varid(psfcncid, "PSFC", psfcid)		! [Pa]
if(status /= nf90_NoErr) call handle_err(status)


END SUBROUTINE open_netcdf_files

!*******************************************************************


SUBROUTINE get_data(precip,evap,u,v,w,t,q,qc,qt,pp,pb,pbl_hgt,psfc)
!-----------------------------------------------
! read in the data for the first time
!-----------------------------------------------

USE global_data
USE netcdf

IMPLICIT NONE

REAL, DIMENSION(:,:,:) :: precip,evap,pbl_hgt,psfc
REAL, DIMENSION(:,:,:,:) :: u,v,w,t,q,qc,qt,pp,pb
REAL, DIMENSION(SIZE(u,1),SIZE(u,2),SIZE(u,3),MM5daysteps) :: temp

CHARACTER(LEN=100) :: filename_ext_atm,filename_ext_RAIN,filename_ext_LH,filename_ext_P

INTEGER :: ncid,prencid,lhncid,psfcncid !,uncid,vncid,wncid,tncid,qncid,ppncid,pblncid
INTEGER :: preid,lhid,uid,vid,wid,tid,qid,ppid,pbid,pblid,psfcid
INTEGER :: sind,status,sind2,i,getsteps,getsteps2
!REAL :: filedays

REAL :: dayend

INTEGER :: jd_today,jd_before,new_y,new_m,new_d


!first get the qt
call get_data_mixtot(qc,qt)

call get_filename(day,mon,year,filename_ext_atm,filename_ext_RAIN,filename_ext_LH,filename_ext_P)

!call open_netcdf_files(pncid,preid,lhid,uid,vid,wid,tid,qid,ppid,pbliqd,filename_ext)
call open_netcdf_files(ncid,prencid,lhncid,psfcncid,preid,lhid,uid,vid,wid,tid,qid,ppid,pbid,pblid,psfcid,filename_ext_atm,filename_ext_RAIN,filename_ext_LH,filename_ext_P)

!if this is a day around a storm peak we want the half day after as well
if (peak) then
  dayend = day + 0.5
else
  dayend = day
end if
  
! Since our input files only consist of one day, open all timesteps (MM5daysteps) in file (i.e. remove sind2, make it 1)

!!! SIM DAY SHOULD BE AT THE END OF THE ARRAY, DAY BEFORE JUST BEFORE THAT, ETC.
!!! LAST BACK-TRACKED DAY SHOULD BE AT THE START OF THE ARRAY.
!!! THE LAST TIME POSITION IN THE ARRAY SHOULD BE THE FIRST TIMESTEP OF SIM DAY + 1.

if (totbtadays>1) then
	! Open the first day input file
	status = nf90_get_var(prencid, preid, precip, &
	start=(/bdy,bdy,1/),count=(/dim_j,dim_i,MM5daysteps/))
	if(status /= nf90_NoErr) call handle_err(status)
	
	status = nf90_get_var(lhncid, lhid, evap(:,:,(MM5totsteps-MM5daysteps):), &
	start=(/bdy,bdy,1/),count=(/dim_j,dim_i,MM5daysteps/))
	if(status /= nf90_NoErr) call handle_err(status)
	
	status = nf90_get_var(ncid, qid, temp(:,:,:,:), &
	start=(/bdy,bdy,1,1/),count=(/dim_j,dim_i,dim_k,MM5daysteps/))
	if(status /= nf90_NoErr) call handle_err(status)
	
	q(:,:,:,(MM5totsteps-MM5daysteps):(MM5totsteps-1)) = temp(:,:,dim_k:1:-1,:)

	status = nf90_get_var(ncid, uid, temp(:,:,:,:), &
	start=(/bdy,bdy,1,1/),count=(/dim_j,dim_i,dim_k,MM5daysteps/))
	if(status /= nf90_NoErr) call handle_err(status)
	
	u(:,:,:,(MM5totsteps-MM5daysteps):(MM5totsteps-1)) = temp(:,:,dim_k:1:-1,:)

	status = nf90_get_var(ncid, vid, temp(:,:,:,:), &
	start=(/bdy,bdy,1,1/),count=(/dim_j,dim_i,dim_k,MM5daysteps/))
	if(status /= nf90_NoErr) call handle_err(status)
	
	v(:,:,:,(MM5totsteps-MM5daysteps):(MM5totsteps-1)) = temp(:,:,dim_k:1:-1,:)

	status = nf90_get_var(ncid, wid, temp(:,:,:,:), &
	start=(/bdy,bdy,1,1/),count=(/dim_j,dim_i,dim_k,MM5daysteps/))
	if(status /= nf90_NoErr) call handle_err(status)
	
	w(:,:,:,(MM5totsteps-MM5daysteps):(MM5totsteps-1)) = temp(:,:,dim_k:1:-1,:)

	status = nf90_get_var(ncid, tid, temp(:,:,:,:), &
	start=(/bdy,bdy,1,1/),count=(/dim_j,dim_i,dim_k,MM5daysteps/))
	if(status /= nf90_NoErr) call handle_err(status)
	
	t(:,:,:,(MM5totsteps-MM5daysteps):(MM5totsteps-1)) = temp(:,:,dim_k:1:-1,:)

	status = nf90_get_var(ncid, ppid, temp(:,:,:,:), &
	start=(/bdy,bdy,1,1/),count=(/dim_j,dim_i,dim_k,MM5daysteps/))
	if(status /= nf90_NoErr) call handle_err(status)
	
	pp(:,:,:,(MM5totsteps-MM5daysteps):(MM5totsteps-1)) = temp(:,:,dim_k:1:-1,:)
	
	status = nf90_get_var(ncid, pbid, temp(:,:,:,:), &
	start=(/bdy,bdy,1,1/),count=(/dim_j,dim_i,dim_k,MM5daysteps/))
	if(status /= nf90_NoErr) call handle_err(status)
	
	pb(:,:,:,(MM5totsteps-MM5daysteps):(MM5totsteps-1)) = temp(:,:,dim_k:1:-1,:)

	status = nf90_get_var(ncid, pblid, pbl_hgt(:,:,(MM5totsteps-MM5daysteps):), &
	start=(/bdy,bdy,1/),count=(/dim_j,dim_i,MM5daysteps/))
	if(status /= nf90_NoErr) call handle_err(status)
	
	status = nf90_get_var(psfcncid, psfcid, psfc(:,:,(MM5totsteps-MM5daysteps):), &
	start=(/bdy,bdy,1/),count=(/dim_j,dim_i,MM5daysteps/))
	if(status /= nf90_NoErr) call handle_err(status)

	! close the netcdf files
	status = nf90_close(ncid)
	status = nf90_close(prencid)
	status = nf90_close(lhncid)
	status = nf90_close(psfcncid)
	
	print *,'L948, Input file of first day loaded successfully:',filename_ext_atm


	! Get julian day of current day 
	jd_today = julian(year,mon,day)
	!print *,'L918, jd_today=',jd_today
	
	! Get julian day for all other totbtadays and open the corresponding input files
	do i = 1,totbtadays
		jd_before = jd_today-i 
		! Convert julidan day to gregorian
		call gregorian(jd_before,new_y,new_m,new_d)
		!print *,'L909, jd_before,new_y,new_m,new_d=',jd_before,new_y,new_m,new_d
		call get_filename(new_d,new_m,new_y,filename_ext_atm,filename_ext_RAIN,filename_ext_LH,filename_ext_P)
		call open_netcdf_files(ncid,prencid,lhncid,psfcncid,preid,lhid,uid,vid,wid,tid,qid,ppid,pbid,pblid,psfcid,filename_ext_atm,filename_ext_RAIN,filename_ext_LH,filename_ext_P)

		status = nf90_get_var(lhncid, lhid, evap(:,:,MM5totsteps-(MM5daysteps*(i+1)):(MM5totsteps-(MM5daysteps*i)-1)), &
		start=(/bdy,bdy,1/),count=(/dim_j,dim_i,MM5daysteps/))
		if(status /= nf90_NoErr) call handle_err(status)
		
		status = nf90_get_var(ncid, qid, temp(:,:,:,:), &
		start=(/bdy,bdy,1,1/),count=(/dim_j,dim_i,dim_k,MM5daysteps/))
		if(status /= nf90_NoErr) call handle_err(status)
		
		q(:,:,:,MM5totsteps-(MM5daysteps*(i+1)):(MM5totsteps-(MM5daysteps*i)-1)) = temp(:,:,dim_k:1:-1,:)

		status = nf90_get_var(ncid, uid, temp(:,:,:,:), &
		start=(/bdy,bdy,1,1/),count=(/dim_j,dim_i,dim_k,MM5daysteps/))
		if(status /= nf90_NoErr) call handle_err(status)
		
		u(:,:,:,MM5totsteps-(MM5daysteps*(i+1)):(MM5totsteps-(MM5daysteps*i)-1)) = temp(:,:,dim_k:1:-1,:)

		status = nf90_get_var(ncid, vid, temp(:,:,:,:), &
		start=(/bdy,bdy,1,1/),count=(/dim_j,dim_i,dim_k,MM5daysteps/))
		if(status /= nf90_NoErr) call handle_err(status)
		
		v(:,:,:,MM5totsteps-(MM5daysteps*(i+1)):(MM5totsteps-(MM5daysteps*i)-1)) = temp(:,:,dim_k:1:-1,:)

		status = nf90_get_var(ncid, wid, temp(:,:,:,:), &
		start=(/bdy,bdy,1,1/),count=(/dim_j,dim_i,dim_k,MM5daysteps/))
		if(status /= nf90_NoErr) call handle_err(status)
		
		w(:,:,:,MM5totsteps-(MM5daysteps*(i+1)):(MM5totsteps-(MM5daysteps*i)-1)) = temp(:,:,dim_k:1:-1,:)

		status = nf90_get_var(ncid, tid, temp(:,:,:,:), &
		start=(/bdy,bdy,1,1/),count=(/dim_j,dim_i,dim_k,MM5daysteps/))
		if(status /= nf90_NoErr) call handle_err(status)
		
		t(:,:,:,MM5totsteps-(MM5daysteps*(i+1)):(MM5totsteps-(MM5daysteps*i)-1)) = temp(:,:,dim_k:1:-1,:)

		status = nf90_get_var(ncid, ppid, temp(:,:,:,:), &
		start=(/bdy,bdy,1,1/),count=(/dim_j,dim_i,dim_k,MM5daysteps/))
		if(status /= nf90_NoErr) call handle_err(status)
		
		pp(:,:,:,MM5totsteps-(MM5daysteps*(i+1)):(MM5totsteps-(MM5daysteps*i)-1)) = temp(:,:,dim_k:1:-1,:)
		
		status = nf90_get_var(ncid, pbid, temp(:,:,:,:), &
		start=(/bdy,bdy,1,1/),count=(/dim_j,dim_i,dim_k,MM5daysteps/))
		if(status /= nf90_NoErr) call handle_err(status)
		
		pb(:,:,:,MM5totsteps-(MM5daysteps*(i+1)):(MM5totsteps-(MM5daysteps*i)-1)) = temp(:,:,dim_k:1:-1,:)

		status = nf90_get_var(ncid, pblid, pbl_hgt(:,:,MM5totsteps-(MM5daysteps*(i+1)):(MM5totsteps-(MM5daysteps*i)-1)), &
		start=(/bdy,bdy,1/),count=(/dim_j,dim_i,MM5daysteps/))
		if(status /= nf90_NoErr) call handle_err(status)
		
		status = nf90_get_var(psfcncid, psfcid, psfc(:,:,MM5totsteps-(MM5daysteps*(i+1)):(MM5totsteps-(MM5daysteps*i)-1)), &
		start=(/bdy,bdy,1/),count=(/dim_j,dim_i,MM5daysteps/))
		if(status /= nf90_NoErr) call handle_err(status)
		
		print *,'L1018, Input file of previous day loaded successfully:',i,filename_ext_atm
		!print *,'L964, :(MM5totsteps-(MM5daysteps*i)-1)=',(MM5totsteps-(MM5daysteps*i)-1)
		!print *,'shape(precip)=',shape(precip)
		!print *,'shape(evap)',shape(evap)
		!print *,'shape(q)=',shape(q)


		! close the netcdf files
		status = nf90_close(ncid)
		status = nf90_close(prencid)
		status = nf90_close(lhncid)  
		status = nf90_close(psfcncid)  

	end do
	
	! Get julian day for day after sim day (1st timestep needed) and open the corresponding input file
	jd_before = jd_today+1
	call gregorian(jd_before,new_y,new_m,new_d)
	call get_filename(new_d,new_m,new_y,filename_ext_atm,filename_ext_RAIN,filename_ext_LH,filename_ext_P)
	call open_netcdf_files(ncid,prencid,lhncid,psfcncid,preid,lhid,uid,vid,wid,tid,qid,ppid,pbid,pblid,psfcid,filename_ext_atm,filename_ext_RAIN,filename_ext_LH,filename_ext_P)
	
	status = nf90_get_var(lhncid, lhid, evap(:,:,MM5totsteps), &
	start=(/bdy,bdy,1/),count=(/dim_j,dim_i,1/))
	if(status /= nf90_NoErr) call handle_err(status)
	
	status = nf90_get_var(ncid, qid, temp(:,:,:,1), &
	start=(/bdy,bdy,1,1/),count=(/dim_j,dim_i,dim_k,1/))
	if(status /= nf90_NoErr) call handle_err(status)
	
	q(:,:,:,MM5totsteps) = temp(:,:,dim_k:1:-1,1)

	status = nf90_get_var(ncid, uid, temp(:,:,:,1), &
	start=(/bdy,bdy,1,1/),count=(/dim_j,dim_i,dim_k,1/))
	if(status /= nf90_NoErr) call handle_err(status)
	
	u(:,:,:,MM5totsteps) = temp(:,:,dim_k:1:-1,1)

	status = nf90_get_var(ncid, vid, temp(:,:,:,1), &
	start=(/bdy,bdy,1,1/),count=(/dim_j,dim_i,dim_k,1/))
	if(status /= nf90_NoErr) call handle_err(status)
	
	v(:,:,:,MM5totsteps) = temp(:,:,dim_k:1:-1,1)

	status = nf90_get_var(ncid, wid, temp(:,:,:,1), &
	start=(/bdy,bdy,1,1/),count=(/dim_j,dim_i,dim_k,1/))
	if(status /= nf90_NoErr) call handle_err(status)
	
	w(:,:,:,MM5totsteps) = temp(:,:,dim_k:1:-1,1)

	status = nf90_get_var(ncid, tid, temp(:,:,:,1), &
	start=(/bdy,bdy,1,1/),count=(/dim_j,dim_i,dim_k,1/))
	if(status /= nf90_NoErr) call handle_err(status)
	
	t(:,:,:,MM5totsteps) = temp(:,:,dim_k:1:-1,1)

	status = nf90_get_var(ncid, ppid, temp(:,:,:,1), &
	start=(/bdy,bdy,1,1/),count=(/dim_j,dim_i,dim_k,1/))
	if(status /= nf90_NoErr) call handle_err(status)
	
	pp(:,:,:,MM5totsteps) = temp(:,:,dim_k:1:-1,1)
	
	status = nf90_get_var(ncid, pbid, temp(:,:,:,1), &
	start=(/bdy,bdy,1,1/),count=(/dim_j,dim_i,dim_k,1/))
	if(status /= nf90_NoErr) call handle_err(status)
	
	pb(:,:,:,MM5totsteps) = temp(:,:,dim_k:1:-1,1)

	status = nf90_get_var(ncid, pblid, pbl_hgt(:,:,MM5totsteps), &
	start=(/bdy,bdy,1/),count=(/dim_j,dim_i,1/))
	if(status /= nf90_NoErr) call handle_err(status)
	
	status = nf90_get_var(psfcncid, psfcid, psfc(:,:,MM5totsteps), &
	start=(/bdy,bdy,1/),count=(/dim_j,dim_i,1/))
	if(status /= nf90_NoErr) call handle_err(status)
	
	! close the netcdf files
	status = nf90_close(ncid)
	status = nf90_close(prencid)
	status = nf90_close(lhncid)  
	status = nf90_close(psfcncid)
	
	print *,'L1100, Input file of next day (1st time step) loaded successfully:',filename_ext_atm
	
else
	print *,'L1036, Input file of sim day loaded successfully:',filename_ext_atm
	! Open the first day input file ONLY
	status = nf90_get_var(prencid, preid, precip, &
	start=(/bdy,bdy,1/),count=(/dim_j,dim_i,MM5daysteps/))
	if(status /= nf90_NoErr) call handle_err(status)

	status = nf90_get_var(lhncid, lhid, evap(:,:,(MM5totsteps-MM5daysteps):), &
	start=(/bdy,bdy,1/),count=(/dim_j,dim_i,MM5daysteps/))
	if(status /= nf90_NoErr) call handle_err(status)
	
	status = nf90_get_var(ncid, qid, temp(:,:,:,(MM5totsteps-MM5daysteps):), &
	start=(/bdy,bdy,1,1/),count=(/dim_j,dim_i,dim_k,MM5daysteps/))
	if(status /= nf90_NoErr) call handle_err(status)
	
	q = temp(:,:,dim_k:1:-1,:)

	status = nf90_get_var(ncid, uid, temp(:,:,:,(MM5totsteps-MM5daysteps):), &
	start=(/bdy,bdy,1,1/),count=(/dim_j,dim_i,dim_k,MM5daysteps/))
	if(status /= nf90_NoErr) call handle_err(status)
	
	u = temp(:,:,dim_k:1:-1,:)

	status = nf90_get_var(ncid, vid, temp(:,:,:,(MM5totsteps-MM5daysteps):), &
	start=(/bdy,bdy,1,1/),count=(/dim_j,dim_i,dim_k,MM5daysteps/))
	if(status /= nf90_NoErr) call handle_err(status)
	
	v = temp(:,:,dim_k:1:-1,:)

	status = nf90_get_var(ncid, wid, temp(:,:,:,(MM5totsteps-MM5daysteps):), &
	start=(/bdy,bdy,1,1/),count=(/dim_j,dim_i,dim_k,MM5daysteps/))
	if(status /= nf90_NoErr) call handle_err(status)
	
	w = temp(:,:,dim_k:1:-1,:)

	status = nf90_get_var(ncid, tid, temp(:,:,:,(MM5totsteps-MM5daysteps):), &
	start=(/bdy,bdy,1,1/),count=(/dim_j,dim_i,dim_k,MM5daysteps/))
	if(status /= nf90_NoErr) call handle_err(status)
	
	t = temp(:,:,dim_k:1:-1,:)

	status = nf90_get_var(ncid, ppid, temp(:,:,:,(MM5totsteps-MM5daysteps):), &
	start=(/bdy,bdy,1,1/),count=(/dim_j,dim_i,dim_k,MM5daysteps/))
	if(status /= nf90_NoErr) call handle_err(status)
	
	pp = temp(:,:,dim_k:1:-1,:)
	
	status = nf90_get_var(ncid, pbid, temp(:,:,:,(MM5totsteps-MM5daysteps):), &
	start=(/bdy,bdy,1,1/),count=(/dim_j,dim_i,dim_k,MM5daysteps/))
	if(status /= nf90_NoErr) call handle_err(status)
	
	pb = temp(:,:,dim_k:1:-1,:)

	status = nf90_get_var(ncid, pblid, pbl_hgt(:,:,(MM5totsteps-MM5daysteps):), &
	start=(/bdy,bdy,1/),count=(/dim_j,dim_i,MM5daysteps/))
	if(status /= nf90_NoErr) call handle_err(status)
	
	status = nf90_get_var(psfcncid, psfcid, psfc(:,:,(MM5totsteps-MM5daysteps):), &
	start=(/bdy,bdy,1/),count=(/dim_j,dim_i,MM5daysteps/))
	if(status /= nf90_NoErr) call handle_err(status)

	! close the netcdf files
	status = nf90_close(ncid)
	status = nf90_close(prencid)
	status = nf90_close(lhncid)
	status = nf90_close(psfcncid)
	

	! Get julian day of current day 
	jd_today = julian(year,mon,day)
	! Get julian day for day before and open the corresponding input files
	jd_before = jd_today-1	
	call gregorian(jd_before,new_y,new_m,new_d)
	call get_filename(new_d,new_m,new_y,filename_ext_atm,filename_ext_RAIN,filename_ext_LH,filename_ext_P)
	call open_netcdf_files(ncid,prencid,lhncid,psfcncid,preid,lhid,uid,vid,wid,tid,qid,ppid,pbid,pblid,psfcid,filename_ext_atm,filename_ext_RAIN,filename_ext_LH,filename_ext_P)
	print *,'L1110, sim day - 1:',filename_ext_atm
	
	status = nf90_get_var(lhncid, lhid, evap(:,:,(MM5totsteps-(MM5daysteps*2)):), &
	start=(/bdy,bdy,1/),count=(/dim_j,dim_i,MM5daysteps/))
	if(status /= nf90_NoErr) call handle_err(status)
	
	status = nf90_get_var(ncid, qid, temp(:,:,:,(MM5totsteps-(MM5daysteps*2)):), &
	start=(/bdy,bdy,1,1/),count=(/dim_j,dim_i,dim_k,MM5daysteps/))
	if(status /= nf90_NoErr) call handle_err(status)
	
	q(:,:,:,(MM5totsteps-(MM5daysteps*2)):) = temp(:,:,dim_k:1:-1,(MM5totsteps-(MM5daysteps*2)):)

	status = nf90_get_var(ncid, uid, temp(:,:,:,(MM5totsteps-(MM5daysteps*2)):), &
	start=(/bdy,bdy,1,1/),count=(/dim_j,dim_i,dim_k,MM5daysteps/))
	if(status /= nf90_NoErr) call handle_err(status)
	
	u(:,:,:,(MM5totsteps-(MM5daysteps*2)):) = temp(:,:,dim_k:1:-1,(MM5totsteps-(MM5daysteps*2)):)

	status = nf90_get_var(ncid, vid, temp(:,:,:,(MM5totsteps-(MM5daysteps*2)):), &
	start=(/bdy,bdy,1,1/),count=(/dim_j,dim_i,dim_k,MM5daysteps/))
	if(status /= nf90_NoErr) call handle_err(status)
	
	v(:,:,:,(MM5totsteps-(MM5daysteps*2)):) = temp(:,:,dim_k:1:-1,(MM5totsteps-(MM5daysteps*2)):)

	status = nf90_get_var(ncid, wid, temp(:,:,:,(MM5totsteps-(MM5daysteps*2)):), &
	start=(/bdy,bdy,1,1/),count=(/dim_j,dim_i,dim_k,MM5daysteps/))
	if(status /= nf90_NoErr) call handle_err(status)
	
	w(:,:,:,(MM5totsteps-(MM5daysteps*2)):) = temp(:,:,dim_k:1:-1,(MM5totsteps-(MM5daysteps*2)):)

	status = nf90_get_var(ncid, tid, temp(:,:,:,(MM5totsteps-(MM5daysteps*2)):), &
	start=(/bdy,bdy,1,1/),count=(/dim_j,dim_i,dim_k,MM5daysteps/))
	if(status /= nf90_NoErr) call handle_err(status)
	
	t(:,:,:,(MM5totsteps-(MM5daysteps*2)):) = temp(:,:,dim_k:1:-1,(MM5totsteps-(MM5daysteps*2)):)

	status = nf90_get_var(ncid, ppid, temp(:,:,:,(MM5totsteps-(MM5daysteps*2)):), &
	start=(/bdy,bdy,1,1/),count=(/dim_j,dim_i,dim_k,MM5daysteps/))
	if(status /= nf90_NoErr) call handle_err(status)
	
	pp(:,:,:,(MM5totsteps-(MM5daysteps*2)):) = temp(:,:,dim_k:1:-1,(MM5totsteps-(MM5daysteps*2)):)
	
	status = nf90_get_var(ncid, pbid, temp(:,:,:,(MM5totsteps-(MM5daysteps*2)):), &
	start=(/bdy,bdy,1,1/),count=(/dim_j,dim_i,dim_k,MM5daysteps/))
	if(status /= nf90_NoErr) call handle_err(status)
	
	pb(:,:,:,(MM5totsteps-(MM5daysteps*2)):) = temp(:,:,dim_k:1:-1,(MM5totsteps-(MM5daysteps*2)):)

	status = nf90_get_var(ncid, pblid, pbl_hgt(:,:,(MM5totsteps-(MM5daysteps*2)):), &
	start=(/bdy,bdy,1/),count=(/dim_j,dim_i,MM5daysteps/))
	if(status /= nf90_NoErr) call handle_err(status)
	
	status = nf90_get_var(psfcncid, psfcid, psfc(:,:,(MM5totsteps-(MM5daysteps*2)):), &
	start=(/bdy,bdy,1/),count=(/dim_j,dim_i,MM5daysteps/))
	if(status /= nf90_NoErr) call handle_err(status)
	
	! close the netcdf files
	status = nf90_close(ncid)
	status = nf90_close(prencid)
	status = nf90_close(lhncid)  
	status = nf90_close(psfcncid)
		
	! Get julian day for day after (1st timestep needed) and open the corresponding input file
	jd_before = jd_today+1
	call gregorian(jd_before,new_y,new_m,new_d)
	call get_filename(new_d,new_m,new_y,filename_ext_atm,filename_ext_RAIN,filename_ext_LH,filename_ext_P)
	call open_netcdf_files(ncid,prencid,lhncid,psfcncid,preid,lhid,uid,vid,wid,tid,qid,ppid,pbid,pblid,psfcid,filename_ext_atm,filename_ext_RAIN,filename_ext_LH,filename_ext_P)
	print *,'L1177, sim day + 1:',filename_ext_atm
	
	status = nf90_get_var(lhncid, lhid, evap(:,:,MM5totsteps), &
	start=(/bdy,bdy,1/),count=(/dim_j,dim_i,1/))
	if(status /= nf90_NoErr) call handle_err(status)
	
	status = nf90_get_var(ncid, qid, temp(:,:,:,1), &
	start=(/bdy,bdy,1,1/),count=(/dim_j,dim_i,dim_k,1/))
	if(status /= nf90_NoErr) call handle_err(status)
	
	q(:,:,:,MM5totsteps) = temp(:,:,dim_k:1:-1,1)

	status = nf90_get_var(ncid, uid, temp(:,:,:,1), &
	start=(/bdy,bdy,1,1/),count=(/dim_j,dim_i,dim_k,1/))
	if(status /= nf90_NoErr) call handle_err(status)
	
	u(:,:,:,MM5totsteps) = temp(:,:,dim_k:1:-1,1)

	status = nf90_get_var(ncid, vid, temp(:,:,:,1), &
	start=(/bdy,bdy,1,1/),count=(/dim_j,dim_i,dim_k,1/))
	if(status /= nf90_NoErr) call handle_err(status)
	
	v(:,:,:,MM5totsteps) = temp(:,:,dim_k:1:-1,1)

	status = nf90_get_var(ncid, wid, temp(:,:,:,1), &
	start=(/bdy,bdy,1,1/),count=(/dim_j,dim_i,dim_k,1/))
	if(status /= nf90_NoErr) call handle_err(status)
	
	w(:,:,:,MM5totsteps) = temp(:,:,dim_k:1:-1,1)

	status = nf90_get_var(ncid, tid, temp(:,:,:,1), &
	start=(/bdy,bdy,1,1/),count=(/dim_j,dim_i,dim_k,1/))
	if(status /= nf90_NoErr) call handle_err(status)
	
	t(:,:,:,MM5totsteps) = temp(:,:,dim_k:1:-1,1)

	status = nf90_get_var(ncid, ppid, temp(:,:,:,1), &
	start=(/bdy,bdy,1,1/),count=(/dim_j,dim_i,dim_k,1/))
	if(status /= nf90_NoErr) call handle_err(status)
	
	pp(:,:,:,MM5totsteps) = temp(:,:,dim_k:1:-1,1)
	
	status = nf90_get_var(ncid, pbid, temp(:,:,:,1), &
	start=(/bdy,bdy,1,1/),count=(/dim_j,dim_i,dim_k,1/))
	if(status /= nf90_NoErr) call handle_err(status)
	
	pb(:,:,:,MM5totsteps) = temp(:,:,dim_k:1:-1,1)

	status = nf90_get_var(ncid, pblid, pbl_hgt(:,:,MM5totsteps), &
	start=(/bdy,bdy,1/),count=(/dim_j,dim_i,1/))
	if(status /= nf90_NoErr) call handle_err(status)
	
	status = nf90_get_var(psfcncid, psfcid, psfc(:,:,MM5totsteps), &
	start=(/bdy,bdy,1/),count=(/dim_j,dim_i,1/))
	if(status /= nf90_NoErr) call handle_err(status)
	
	! close the netcdf files
	status = nf90_close(ncid)
	status = nf90_close(prencid)
	status = nf90_close(lhncid)  
	status = nf90_close(psfcncid)
	
! 	print *,'shape(precip)=',shape(precip)
! 	print *,'shape(evap)=',shape(evap)
! 	print *,'shape(pp)=',shape(pp)
! 	print *,'shape(psfc)=',shape(psfc)


end if
  

!evap converted to mm > Unit conversion checked and OK 18/7/17 :)
evap = evap*(1440/MM5daysteps)*60/Lv
  
qt = qt + q ! i.e. QVAPOUR + SUM(QCLD,QRAIN,QSNOW,QICE)
!qt = q
qc = qt

END SUBROUTINE get_data

!***********************************************************************


SUBROUTINE open_mixtot_netcdf_files(ncid,clwid,rnwid,snowid,iceid,filename_ext_atm)
!----------------------------------------------------------------
! open all the netcdf data files and get the variable ids
!------------------------------------------------------------------

USE netcdf
USE global_data

IMPLICIT NONE

INTEGER, INTENT(OUT) :: ncid !clwncid,rnwncid,snowncid,icencid
INTEGER, INTENT(OUT) :: clwid,rnwid,snowid,iceid
CHARACTER(LEN=100), INTENT(IN) :: filename_ext_atm

INTEGER :: status


!
! open the netcdf files - ATMOSPHERIC VARIABLES
!
status = NF90_OPEN(TRIM(dirdata_atm)//TRIM(filename_ext_atm), NF90_NOWRITE, ncid)
if (status /= NF90_NOERR) call handle_err(status)

! status = NF90_OPEN(TRIM(diri)//"clw"//TRIM(filename_ext), NF90_NOWRITE, clwncid)
! if (status /= NF90_NOERR) call handle_err(status)
! 
! status = NF90_OPEN(TRIM(diri)//"rnw"//TRIM(filename_ext), NF90_NOWRITE, rnwncid)
! if (status /= NF90_NOERR) call handle_err(status)
! 
! status = NF90_OPEN(TRIM(diri)//"snow"//TRIM(filename_ext), NF90_NOWRITE, snowncid)
! if (status /= NF90_NOERR) call handle_err(status)
! 
! status = NF90_OPEN(TRIM(diri)//"ice"//TRIM(filename_ext), NF90_NOWRITE, icencid)
! if (status /= NF90_NOERR) call handle_err(status)


!
!get ids for each variable
!
status = nf90_inq_varid(ncid, "QCLOUD", clwid)		! [kgkg-1]  ! QCLOUD
if(status /= nf90_NoErr) call handle_err(status)
status = nf90_inq_varid(ncid, "QRAIN", rnwid)		! [kgkg-1] ! QRAIN
if(status /= nf90_NoErr) call handle_err(status)
status = nf90_inq_varid(ncid, "QSNOW", snowid)		! [kgkg-1] ! QSNOW
if(status /= nf90_NoErr) call handle_err(status)
status = nf90_inq_varid(ncid, "QICE", iceid)		! [kgkg-1] ! QICE
if(status /= nf90_NoErr) call handle_err(status)
! status = nf90_inq_varid(clwncid, "clw", clwid)
! if(status /= nf90_NoErr) call handle_err(status)
! status = nf90_inq_varid(rnwncid, "rnw", rnwid)
! if(status /= nf90_NoErr) call handle_err(status)
! status = nf90_inq_varid(snowncid, "snow", snowid)
! if(status /= nf90_NoErr) call handle_err(status)
! status = nf90_inq_varid(icencid, "ice", iceid)
! if(status /= nf90_NoErr) call handle_err(status)

END SUBROUTINE open_mixtot_netcdf_files

!*******************************************************************


SUBROUTINE get_data_mixtot(qc,qt)
!-----------------------------------------------
! read in the data for the first time
!-----------------------------------------------

USE global_data
USE netcdf

IMPLICIT NONE

REAL, DIMENSION(:,:,:,:) :: qc,qt

CHARACTER(LEN=100) :: filename_ext_atm,filename_ext_RAIN,filename_ext_LH,filename_ext_P

REAL, DIMENSION(SIZE(qt,1),SIZE(qt,2),SIZE(qt,3),SIZE(qt,4)) :: clw,rnw,snow,ice,temp

INTEGER :: ncid !clwncid,rnwncid,snowncid,icencid
INTEGER :: clwid,rnwid,snowid,iceid
INTEGER :: sind,status,sind2,i,getsteps,getsteps2
!REAL :: filedays

REAL :: dayend

INTEGER :: jd_today,jd_before,new_y,new_m,new_d

call get_filename(day,mon,year,filename_ext_atm,filename_ext_RAIN,filename_ext_LH,filename_ext_P)

!call open_mixtot_netcdf_files(clwncid,rnwncid,snowncid,icencid, &
call open_mixtot_netcdf_files(ncid,clwid,rnwid,snowid,iceid,filename_ext_atm)

!if this is a day around a storm peak we want the half day after as well
if (peak) then
  dayend = day + 0.5
else
  dayend = day
end if
  
! !  
! !calculate the number of days from start of file to desired day
! !
! ! if (string_to_int(filename_ext(16:17))==mon) then
! !   filedays = dayend
! ! else
! !   filedays = days_in_month(string_to_int(filename_ext(16:17)), &
! !   		string_to_int(filename_ext(12:15))) + dayend
! ! 		
! ! end if
! ! print *,'L1080, filedays,totbtadays',filedays,totbtadays
! 
! !
! !we want the event day + totbtadays before it + the time step after
! ! (ie 0 hour time step which is the last in each file)
! !
! 
! !if filedays<totbtadays we need to load some from these files
! !and some from the previous files
! !
! ! if (filedays < (totbtadays+2)) then
! !   !
! !   !load variable from current file
! !   !  
! !   
! !   !Need to reorder dimensions below if using this section (see input data shape)
! !   
! !   getsteps = ANINT(MM5daysteps*filedays)
! !   sind = MM5totsteps - getsteps
! ! 
! !   print *,"L1099, start+count = ",bdy+dim_j,bdy+dim_i,1+dim_k,getsteps
! !   print *,"L1100, MM5daysteps,MM5totsteps,sind,getsteps",MM5daysteps,MM5totsteps,sind,getsteps
! ! 
! !   status = nf90_get_var(ncid, clwid, clw(:,:,:,sind+1:), &
! !     start=(/bdy,bdy,1,1/),count=(/dim_j,dim_i,dim_k,getsteps/)) 
! !   if(status /= nf90_NoErr) call handle_err(status)
! !   
! !   status = nf90_get_var(ncid, rnwid, rnw(:,:,:,sind+1:), &
! !     start=(/bdy,bdy,1,1/),count=(/dim_j,dim_i,dim_k,getsteps/))
! !   if(status /= nf90_NoErr) call handle_err(status)
! ! 
! !   status = nf90_get_var(ncid, snowid, snow(:,:,:,sind+1:), &
! !     start=(/bdy,bdy,1,1/),count=(/dim_j,dim_i,dim_k,getsteps/))
! !   if(status /= nf90_NoErr) call handle_err(status)
! ! 
! !   status = nf90_get_var(ncid, iceid, ice(:,:,:,sind+1:), &
! !     start=(/bdy,bdy,1,1/),count=(/dim_j,dim_i,dim_k,getsteps/))
! !   if(status /= nf90_NoErr) call handle_err(status)
! ! 
! ! !   status = nf90_get_var(clwncid, clwid, clw(:,:,:,sind+1:), &
! ! !     start=(/bdy,bdy,1,1/),count=(/dim_j,dim_i,dim_k,getsteps/))
! ! !   if(status /= nf90_NoErr) call handle_err(status)
! ! ! 
! ! !   status = nf90_get_var(rnwncid, rnwid, rnw(:,:,:,sind+1:), &
! ! !     start=(/bdy,bdy,1,1/),count=(/dim_j,dim_i,dim_k,getsteps/))
! ! !   if(status /= nf90_NoErr) call handle_err(status)
! ! ! 
! ! !   status = nf90_get_var(snowncid, snowid, snow(:,:,:,sind+1:), &
! ! !     start=(/bdy,bdy,1,1/),count=(/dim_j,dim_i,dim_k,getsteps/))
! ! !   if(status /= nf90_NoErr) call handle_err(status)
! ! ! 
! ! !   status = nf90_get_var(icencid, iceid, ice(:,:,:,sind+1:), &
! ! !     start=(/bdy,bdy,1,1/),count=(/dim_j,dim_i,dim_k,getsteps/))
! ! !   if(status /= nf90_NoErr) call handle_err(status)
! ! 
! !   
! !   !
! !   ! close the netcdf files
! !   !
! !   status = nf90_close(ncid)
! ! !   status = nf90_close(clwncid)
! ! !   status = nf90_close(rnwncid)
! ! !   status = nf90_close(snowncid)
! ! !   status = nf90_close(icencid)
! !   
! !   
! !   !
! !   !now need to get file before to get the rest of the data
! !   !
! ! !   if (string_to_int(filename_ext(16:17))==1) then
! ! !     filename_ext(16:17) = "11"
! ! !     filename_ext(12:15) = TRIM(int_to_string(year-1))
! ! !   else 
! ! !     if (string_to_int(filename_ext(16:17))==2) then
! ! !       filename_ext(16:17) = "12"
! ! !       filename_ext(12:15) = TRIM(int_to_string(year-1))    
! ! !     else 
! ! !       if (string_to_int(filename_ext(16:17))<12) then
! ! !         filename_ext(16:17) = "0"//TRIM(int_to_string(string_to_int(filename_ext(16:17))-2))
! ! !       else
! ! !         filename_ext(16:17) = TRIM(int_to_string(string_to_int(filename_ext(16:17))-2))
! ! !       end if
! ! !     end if
! ! !   end if
! !   
! !   !call open_mixtot_netcdf_files(clwncid,rnwncid,snowncid,icencid, &
! !   call open_mixtot_netcdf_files(ncid,clwid,rnwid,snowid,iceid,filename_ext)
! !   !print *,'opened mixtot_netcdf_files for get_data_mixtot, call#2 L1163'
! ! 
! ! 
! !   !
! !   !calculate the number of time steps in these new files
! !   !
! !   getsteps2 = sind 
! !   
! !   !
! !   ! calculate where to start in these new files
! !   !
! !   if (mon>2) then
! !     sind2 = (days_in_month(mon-1,year)+days_in_month(mon-2,year))*MM5daysteps - getsteps2 + 1
! !   else if (mon==1) then
! !     sind2 = (days_in_month(11,year-1)+days_in_month(12,year-1))*MM5daysteps - getsteps2 + 1
! !   else if (mon==2) then
! !     sind2 = (days_in_month(12,year-1)+days_in_month(1,year))*MM5daysteps - getsteps2 + 1
! !   end if
! !   
! !   
! !   !adjust if this is the begining of simulation (have extra 0 time step
! !   if (string_to_int(filename_ext(16:17))==smon.AND.string_to_int(filename_ext(12:15))==syear) then
! !     sind2=sind2+1
! !   end if
! ! print *, 'sind2',sind2
! !   !
! !   !load variable from current file
! !   !   
! !   !status = nf90_get_var(clwncid, clwid, clw(:,:,:,:sind), &
! !   status = nf90_get_var(ncid, clwid, clw(:,:,:,:sind), &
! !     start=(/bdy,bdy,1,sind2/),count=(/dim_j,dim_i,dim_k,getsteps2/))
! !   if(status /= nf90_NoErr) call handle_err(status) 
! ! 
! !   !status = nf90_get_var(rnwncid, rnwid, rnw(:,:,:,:sind), &
! !   status = nf90_get_var(ncid, rnwid, rnw(:,:,:,:sind), &
! !     start=(/bdy,bdy,1,sind2/),count=(/dim_j,dim_i,dim_k,getsteps2/))
! !   if(status /= nf90_NoErr) call handle_err(status) 
! ! 
! !   !status = nf90_get_var(snowncid, snowid, snow(:,:,:,:sind), &
! !   status = nf90_get_var(ncid, snowid, snow(:,:,:,:sind), &
! !     start=(/bdy,bdy,1,sind2/),count=(/dim_j,dim_i,dim_k,getsteps2/))
! !   if(status /= nf90_NoErr) call handle_err(status) 
! ! 
! !   !status = nf90_get_var(icencid, iceid, ice(:,:,:,:sind), &
! !   status = nf90_get_var(ncid, iceid, ice(:,:,:,:sind), &
! !     start=(/bdy,bdy,1,sind2/),count=(/dim_j,dim_i,dim_k,getsteps2/))
! !   if(status /= nf90_NoErr) call handle_err(status) 
! ! 
! !   
! !   !
! !   ! close the netcdf files
! !   !
! !   status = nf90_close(ncid)
! ! !   status = nf90_close(clwncid)
! ! !   status = nf90_close(rnwncid)
! ! !   status = nf90_close(snowncid)
! ! !   status = nf90_close(icencid)
! ! 
! !   
! ! else
! ! 
! !   sind2 = ANINT(filedays*MM5daysteps - MM5totsteps)
! !   
! !   print *,'L1237,sind2,filedays,MM5daysteps,MM5totsteps',sind2,filedays,MM5daysteps,MM5totsteps
! !   
! ! 
! !   !
! !   !load variable from current file
! !   !
! !   !status = nf90_get_var(ncid, clwid, clw, &
! !   !  start=(/bdy,bdy,1,sind2+1/),count=(/dim_j,dim_i,dim_k,MM5totsteps/))
! !   !if(status /= nf90_NoErr) call handle_err(status)
! !   status = nf90_get_var(ncid, clwid, clw, &
! !     start=(/bdy,bdy,1,sind2+1/),count=(/dim_j,dim_i,dim_k,MM5totsteps/))
! !   if(status /= nf90_NoErr) call handle_err(status)
! ! 
! !   !status = nf90_get_var(rnwncid, rnwid, rnw, &
! !   status = nf90_get_var(ncid, rnwid, rnw, &
! !     start=(/bdy,bdy,1,sind2+1/),count=(/dim_j,dim_i,dim_k,MM5totsteps/))
! !   if(status /= nf90_NoErr) call handle_err(status)
! ! 
! !   !status = nf90_get_var(snowncid, snowid, snow, &
! !   status = nf90_get_var(ncid, snowid, snow, &
! !     start=(/bdy,bdy,1,sind2+1/),count=(/dim_j,dim_i,dim_k,MM5totsteps/))
! !   if(status /= nf90_NoErr) call handle_err(status)
! ! 
! !   !status = nf90_get_var(icencid, iceid, ice, &
! !   status = nf90_get_var(ncid, iceid, ice, &
! !     start=(/bdy,bdy,1,sind2+1/),count=(/dim_j,dim_i,dim_k,MM5totsteps/))
! !   if(status /= nf90_NoErr) call handle_err(status)
! ! 
! !   
! !   !
! !   ! close the netcdf files
! !   !
! !   status = nf90_close(ncid)
! ! !   status = nf90_close(clwncid)
! ! !   status = nf90_close(rnwncid)
! ! !   status = nf90_close(snowncid)
! ! !   status = nf90_close(icencid)
! 
! 
! 
! 
! ! Since our input files only consist of one day, open all timesteps (MM5daysteps) in file (i.e. remove sind2, make it 1)
! 
! !load variable from current file
!   !
!   !status = nf90_get_var(ncid, clwid, clw, &
!   !  start=(/bdy,bdy,1,sind2+1/),count=(/dim_j,dim_i,dim_k,MM5totsteps/))
!   !if(status /= nf90_NoErr) call handle_err(status)
!   status = nf90_get_var(ncid, clwid, clw, &
!     start=(/bdy,bdy,1,1/),count=(/dim_j,dim_i,dim_k,MM5daysteps/))    
!   if(status /= nf90_NoErr) call handle_err(status)
! 
!   !status = nf90_get_var(rnwncid, rnwid, rnw, &
!   status = nf90_get_var(ncid, rnwid, rnw, &
!     start=(/bdy,bdy,1,1/),count=(/dim_j,dim_i,dim_k,MM5daysteps/))
!   if(status /= nf90_NoErr) call handle_err(status)
! 
!   !status = nf90_get_var(snowncid, snowid, snow, &
!   status = nf90_get_var(ncid, snowid, snow, &
!     start=(/bdy,bdy,1,1/),count=(/dim_j,dim_i,dim_k,MM5daysteps/))
!   if(status /= nf90_NoErr) call handle_err(status)
! 
!   !status = nf90_get_var(icencid, iceid, ice, &
!   status = nf90_get_var(ncid, iceid, ice, &
!     start=(/bdy,bdy,1,1/),count=(/dim_j,dim_i,dim_k,MM5daysteps/))
!   if(status /= nf90_NoErr) call handle_err(status)
! 
!   
!   !
!   ! close the netcdf files
!   !
!   status = nf90_close(ncid)
! !   status = nf90_close(clwncid)
! !   status = nf90_close(rnwncid)
! !   status = nf90_close(snowncid)
! !   status = nf90_close(icencid)
! 
! !end if



! Since our input files only consist of one day, open all timesteps (MM5daysteps) in file (i.e. remove sind2, make it 1)

! We want the event day + totbtadays before it

if (totbtadays>1) then
	! Open the first day input file
	status = nf90_get_var(ncid, clwid, temp(:,:,:,(MM5totsteps-MM5daysteps):), &
	start=(/bdy,bdy,1,1/),count=(/dim_j,dim_i,dim_k,MM5daysteps/))    
	if(status /= nf90_NoErr) call handle_err(status)
	
	clw(:,:,:,(MM5totsteps-MM5daysteps):) = temp(:,:,dim_k:1:-1,(MM5totsteps-MM5daysteps):)

	status = nf90_get_var(ncid, rnwid, temp(:,:,:,(MM5totsteps-MM5daysteps):), &
	start=(/bdy,bdy,1,1/),count=(/dim_j,dim_i,dim_k,MM5daysteps/))
	if(status /= nf90_NoErr) call handle_err(status)
	
	rnw(:,:,:,(MM5totsteps-MM5daysteps):) = temp(:,:,dim_k:1:-1,(MM5totsteps-MM5daysteps):)

	status = nf90_get_var(ncid, snowid, snow(:,:,:,(MM5totsteps-MM5daysteps):), &
	start=(/bdy,bdy,1,1/),count=(/dim_j,dim_i,dim_k,MM5daysteps/))
	if(status /= nf90_NoErr) call handle_err(status)
	
	snow(:,:,:,(MM5totsteps-MM5daysteps):) = temp(:,:,dim_k:1:-1,(MM5totsteps-MM5daysteps):)

	status = nf90_get_var(ncid, iceid, temp(:,:,:,(MM5totsteps-MM5daysteps):), &
	start=(/bdy,bdy,1,1/),count=(/dim_j,dim_i,dim_k,MM5daysteps/))
	if(status /= nf90_NoErr) call handle_err(status)
	
	ice(:,:,:,(MM5totsteps-MM5daysteps):) = temp(:,:,dim_k:1:-1,(MM5totsteps-MM5daysteps):)
	
	! close the netcdf file
	status = nf90_close(ncid)

	! Get julian day of current day 
	jd_today = julian(year,mon,day)
	
	! Get julian/gregrorian day for all other totbtadays and open the corresponding input files
	do i = 1,totbtadays
		jd_before = jd_today-i 
		! Convert julian day to gregorian
		call gregorian(jd_before,new_y,new_m,new_d)
		call get_filename(new_d,new_m,new_y,filename_ext_atm,filename_ext_RAIN,filename_ext_LH,filename_ext_P)
		call open_mixtot_netcdf_files(ncid,clwid,rnwid,snowid,iceid,filename_ext_atm)

		status = nf90_get_var(ncid, clwid, temp(:,:,:,:(MM5totsteps-(MM5daysteps*i)-1)), &
		start=(/bdy,bdy,1,1/),count=(/dim_j,dim_i,dim_k,MM5daysteps/))    
		if(status /= nf90_NoErr) call handle_err(status)
		
		clw(:,:,:,:(MM5totsteps-(MM5daysteps*i)-1)) = temp(:,:,dim_k:1:-1,:(MM5totsteps-(MM5daysteps*i)-1))

		status = nf90_get_var(ncid, rnwid, temp(:,:,:,:(MM5totsteps-(MM5daysteps*i)-1)), &
		start=(/bdy,bdy,1,1/),count=(/dim_j,dim_i,dim_k,MM5daysteps/))
		if(status /= nf90_NoErr) call handle_err(status)
		
		rnw(:,:,:,:(MM5totsteps-(MM5daysteps*i)-1)) = temp(:,:,dim_k:1:-1,:(MM5totsteps-(MM5daysteps*i)-1))

		status = nf90_get_var(ncid, snowid, temp(:,:,:,:(MM5totsteps-(MM5daysteps*i)-1)), &
		start=(/bdy,bdy,1,1/),count=(/dim_j,dim_i,dim_k,MM5daysteps/))
		if(status /= nf90_NoErr) call handle_err(status)
		
		snow(:,:,:,:(MM5totsteps-(MM5daysteps*i)-1)) = temp(:,:,dim_k:1:-1,:(MM5totsteps-(MM5daysteps*i)-1))

		status = nf90_get_var(ncid, iceid, temp(:,:,:,:(MM5totsteps-(MM5daysteps*i)-1)), &
		start=(/bdy,bdy,1,1/),count=(/dim_j,dim_i,dim_k,MM5daysteps/))
		if(status /= nf90_NoErr) call handle_err(status)
		
		ice(:,:,:,:(MM5totsteps-(MM5daysteps*i)-1)) = temp(:,:,dim_k:1:-1,:(MM5totsteps-(MM5daysteps*i)-1))
		
		! close the netcdf file
		status = nf90_close(ncid)
	end do
	
	! Get julian day for day after sim day (1st timestep needed) and open the corresponding input file
	jd_before = jd_today+1
	call gregorian(jd_before,new_y,new_m,new_d)
	call get_filename(new_d,new_m,new_y,filename_ext_atm,filename_ext_RAIN,filename_ext_LH,filename_ext_P)
	call open_mixtot_netcdf_files(ncid,clwid,rnwid,snowid,iceid,filename_ext_atm)

	status = nf90_get_var(ncid, clwid, temp(:,:,:,1), &
	start=(/bdy,bdy,1,1/),count=(/dim_j,dim_i,dim_k,1/))    
	if(status /= nf90_NoErr) call handle_err(status)
	
	clw(:,:,:,MM5totsteps) = temp(:,:,dim_k:1:-1,1)

	status = nf90_get_var(ncid, rnwid, temp(:,:,:,1), &
	start=(/bdy,bdy,1,1/),count=(/dim_j,dim_i,dim_k,1/))
	if(status /= nf90_NoErr) call handle_err(status)
	
	rnw(:,:,:,MM5totsteps) = temp(:,:,dim_k:1:-1,1)

	status = nf90_get_var(ncid, snowid, temp(:,:,:,1), &
	start=(/bdy,bdy,1,1/),count=(/dim_j,dim_i,dim_k,1/))
	if(status /= nf90_NoErr) call handle_err(status)
	
	snow(:,:,:,MM5totsteps) = temp(:,:,dim_k:1:-1,1)

	status = nf90_get_var(ncid, iceid, temp(:,:,:,1), &
	start=(/bdy,bdy,1,1/),count=(/dim_j,dim_i,dim_k,1/))
	if(status /= nf90_NoErr) call handle_err(status)
	
	ice(:,:,:,MM5totsteps) = temp(:,:,dim_k:1:-1,1)
	
	! close the netcdf file
	status = nf90_close(ncid)
	
	print *,'L1100, Input file of next day (1st time step) loaded successfully:',filename_ext_atm
	
else
	! Open the first day input file ONLY
	status = nf90_get_var(ncid, clwid, temp(:,:,:,:), &
	start=(/bdy,bdy,1,1/),count=(/dim_j,dim_i,dim_k,MM5daysteps/))    
	if(status /= nf90_NoErr) call handle_err(status)
	
	clw = temp(:,:,dim_k:1:-1,:)

	status = nf90_get_var(ncid, rnwid, temp(:,:,:,:), &
	start=(/bdy,bdy,1,1/),count=(/dim_j,dim_i,dim_k,MM5daysteps/))
	if(status /= nf90_NoErr) call handle_err(status)
	
	rnw = temp(:,:,dim_k:1:-1,:)

	status = nf90_get_var(ncid, snowid, temp(:,:,:,:), &
	start=(/bdy,bdy,1,1/),count=(/dim_j,dim_i,dim_k,MM5daysteps/))
	if(status /= nf90_NoErr) call handle_err(status)
	
	snow = temp(:,:,dim_k:1:-1,:)

	status = nf90_get_var(ncid, iceid, temp(:,:,:,:), &
	start=(/bdy,bdy,1,1/),count=(/dim_j,dim_i,dim_k,MM5daysteps/))
	if(status /= nf90_NoErr) call handle_err(status)
	
	ice = temp(:,:,dim_k:1:-1,:)
	
	! close the netcdf file
	status = nf90_close(ncid)
end if

qc = clw + rnw + snow + ice
qt = qc  

END SUBROUTINE get_data_mixtot

!***********************************************************************


REAL FUNCTION lin_interp(var,fac)
!---------------------------------------
!linearly interpolate between the values of var
!fac is the proporational distance from the first value
!------------------------------------------

IMPLICIT NONE

REAL, INTENT(IN), DIMENSION(2) :: var
REAL, INTENT(IN) :: fac

lin_interp = var(1)*(1-fac) + var(2)*fac

END FUNCTION lin_interp

!***********************************************************************

FUNCTION lin_interp2D(var,fac)
!---------------------------------------
!linearly interpolate between the values of var (last dimension must have size 2)
!fac is the proporational distance from the first value
!------------------------------------------

IMPLICIT NONE

REAL, INTENT(IN), DIMENSION(:,:,:) :: var
REAL, INTENT(IN) :: fac

REAL, DIMENSION(SIZE(var,1),SIZE(var,2)):: lin_interp2D

lin_interp2D = var(:,:,1)*(1-fac) + var(:,:,2)*fac

END FUNCTION lin_interp2D

!************************************************************************

FUNCTION lin_interp3D(var,fac)
!---------------------------------------
!linearly interpolate between the values of var (last dimension must have size 2)
!fac is the proporational distance from the first value
!------------------------------------------

IMPLICIT NONE

REAL, INTENT(IN), DIMENSION(:,:,:,:) :: var
REAL, INTENT(IN) :: fac

REAL, DIMENSION(SIZE(var,1),SIZE(var,2),SIZE(var,3)) :: lin_interp3D

lin_interp3D = var(:,:,:,1)*(1-fac) + var(:,:,:,2)*fac

END FUNCTION lin_interp3D

!************************************************************************

SUBROUTINE parcel_release_time(precip,npar,par_release)
!---------------------------------------------------
! here we calculate the times of day to release our parcels
! based on a random precipitation weighted sampling
!-----------------------------------------------------

USE global_data

IMPLICIT NONE

REAL, DIMENSION(:), INTENT(IN) :: precip
INTEGER, DIMENSION(:), INTENT(OUT) :: par_release
INTEGER, INTENT(IN) :: npar

REAL, DIMENSION(SIZE(precip)*inMM5tsteps) :: cumm_precip
REAL, DIMENSION(npar) :: rand_nums
INTEGER :: tt,rr,ss,rec


par_release = 0

call RANDOM_NUMBER(rand_nums)

cumm_precip = 0.

rec = 0

do tt = 1,SIZE(precip)
  do ss = 1,inMM5tsteps
    rec = rec + 1
    if (rec==1) then
      cumm_precip(1) = precip(1)/inMM5tsteps
    else
      cumm_precip(rec) = cumm_precip(rec-1) + precip(tt)/inMM5tsteps
    end if
  end do
end do

cumm_precip = cumm_precip/cumm_precip(SIZE(cumm_precip))

do rr = 1,npar
  do tt = 1,SIZE(cumm_precip)
    if (cumm_precip(tt)>rand_nums(rr)) then
      par_release(tt) = par_release(tt) + 1
      EXIT
    end if
  end do
end do

!print *,'L1675, par_release(tt)=',par_release(tt)

END SUBROUTINE parcel_release_time

!**********************************************************************

SUBROUTINE parcel_release_height(pw,par_lev)
!----------------------------------------------
! calculate the height to release the parcel from
! based on precipitable water weighted random sampling
!-----------------------------------------------

IMPLICIT NONE

REAL, INTENT(IN), DIMENSION(:) :: pw

INTEGER, INTENT(OUT) :: par_lev

REAL :: rand_num
INTEGER :: kk


call RANDOM_NUMBER(rand_num)
!print *,'L1609, rand_num, rand_num*size(pw)=',rand_num,rand_num*size(pw)
!rand_num = rand_num*pw(1)
rand_num = 1 + FLOOR(size(pw)*rand_num)
!print *,'L1612, rand_num=',rand_num


! do kk = SIZE(pw),1,-1 ! kk is then between 29 and 1
!   if (pw(kk)>rand_num) then
!     par_lev = kk
!     EXIT
!   end if
! end do

par_lev = rand_num
!print *,'L1623, par_lev=',par_lev

END SUBROUTINE parcel_release_height

!**************************************************************************

SUBROUTINE lin_interp_inMM5tsteps(var)
!------------------------------------------
!linearly interpolate through time inside MM5 time steps
!----------------------------------------------

USE global_data

IMPLICIT NONE

REAL, INTENT(INOUT), DIMENSION(:,:,:,:) :: var

INTEGER :: i


do i = 1,inMM5tsteps-1
  var(:,:,:,i+1::inMM5tsteps) = (1-(i*1./inMM5tsteps))*var(:,:,:,1:MM5daysteps+1-inMM5tsteps:inMM5tsteps) &
                 & + (i*1./inMM5tsteps)*var(:,:,:,1+inMM5tsteps::inMM5tsteps)
end do

END SUBROUTINE lin_interp_inMM5tsteps

!****************************************************************************

SUBROUTINE calc_pw(mix,pres,surf_pres,ptop,pw)
!------------------------------------------
! calculate the total precipitable water from MM5 fields
! only for the day of current interest
! save as accumulated field from the ground up (4D field)
!-------------------------------------------------

USE global_data

IMPLICIT NONE

REAL, INTENT(IN), DIMENSION(:,:,:,:) :: mix,pres
REAL, INTENT(IN), DIMENSION(:,:,:) :: surf_pres
REAL, INTENT(OUT), DIMENSION(:,:,:,:) :: pw

REAL, DIMENSION(dim_j,dim_i,dim_k,SIZE(mix,4)) :: dp
INTEGER :: k

REAL, INTENT(IN) :: ptop


!
! calculate the change in pressure (Pa) represented by each point
!
!for highest level
dp(:,:,1,:) = SUM(pres(:,:,:2,:),3)/2. - ptop

!for the middle levels
do k = 2,dim_k-1
  dp(:,:,k,:) = (pres(:,:,k+1,:) - pres(:,:,k-1,:)) /2. !dp(:,:,k,:) = SUM(pres(:,:,k-1:k+1:2,:),3)/2.
end do

!for the lowest level
dp(:,:,dim_k,:) = surf_pres(:,:,:) - SUM(pres(:,:,dim_k-1:,:),3)/2.

 
!mass in mm
pw(:,:,:,::inMM5tsteps) = dp*mix/g


!interpolate inside MM5 time steps
call lin_interp_inMM5tsteps(pw)

!accumulate from the bottom up   
do k = dim_k-1,1,-1             ! i.e. from level 28 to 1
  pw(:,:,k,2:) = pw(:,:,k+1,2:) + pw(:,:,k,2:)
end do


END SUBROUTINE calc_pw

!**************************************************************************

SUBROUTINE calc_tpw(mix,pres,surf_pres,ptop,tpw)
!------------------------------------------
! calculate the total precipitable water in the pbl from MM5 fields
! at every level and time (3D field)
!-------------------------------------------------

USE global_data

IMPLICIT NONE

REAL, INTENT(IN), DIMENSION(:,:,:,:) :: mix,pres
REAL, INTENT(IN), DIMENSION(:,:,:) :: surf_pres
REAL, INTENT(OUT), DIMENSION(:,:,:) :: tpw

REAL, DIMENSION(dim_j,dim_i,dim_k,MM5totsteps) :: dp
INTEGER :: j,i,k,t

REAL, INTENT(IN) :: ptop

!
! calculate the change in pressure (Pa) represented by each point
!
!for highest level
dp(:,:,1,:) = SUM(pres(:,:,:2,:),3)/2. - ptop

!for the middle levels
do k = 2,dim_k-1
  dp(:,:,k,:) = (pres(:,:,k+1,:) - pres(:,:,k-1,:)) /2.
  ! equiv to (p2+p3)/2 - (p1+p2)/2 = (p3-p1)/2
end do

!for the lowest level
dp(:,:,dim_k,:) = surf_pres(:,:,:) - SUM(pres(:,:,dim_k-1:,:),3)/2.


!mass in mmtpw
do j = 1,dim_j
  do i = 1,dim_i
    do t = 1,MM5totsteps
      tpw(j,i,t) = SUM(dp(j,i,:,t)*mix(j,i,:,t)/g)
    end do
  end do
end do

!print *,'L1836, tpw=',tpw(1,1,249), 'dp=',dp(1,1,:,249), 'pres=',pres(1,1,:,249), 'surf_pres=',surf_pres(1,1,249)
!print *,'max surf pres=',maxval(surf_pres)
! print *,'L2044, min pres=',minval(pres),'min pres loc=',minloc(pres) !(j,i,k,time)
! print *,'L2045, min dp=',minval(dp),'min dp loc=',minloc(dp)
! print *,'L2046, pres=',pres(1,1,1,:),'dp=',dp(1,1,1,:)

END SUBROUTINE calc_tpw

!**************************************************************************
SUBROUTINE calc_tpw_pbl(mix,pres,surf_pres,tpw,pbl_lev)
!------------------------------------------
! calculate the total precipitable water in the pbl from MM5 fields
! at every level and time
!-------------------------------------------------

USE global_data

IMPLICIT NONE

REAL, INTENT(IN), DIMENSION(:,:,:,:) :: mix,pres
REAL, INTENT(IN), DIMENSION(:,:,:) :: surf_pres
INTEGER, INTENT(IN), DIMENSION(:,:,:) :: pbl_lev
REAL, INTENT(OUT), DIMENSION(:,:,:) :: tpw

REAL, DIMENSION(dim_j,dim_i,dim_k,MM5totsteps) :: dp
INTEGER :: j,i,k,t


!
! calculate the change in pressure (Pa) represented by each point
!
!for highest level
dp(:,:,1,:) = SUM(pres(:,:,:2,:),3)/2.

!for the middle levels
do k = 2,dim_k-1
  !dp(:,:,k,:) = SUM(pres(:,:,k-1:k+1:2,:),3)/2.
  dp(:,:,k,:) = (pres(:,:,k+1,:) - pres(:,:,k-1,:)) /2.
end do

!for the lowest level
dp(:,:,dim_k,:) = surf_pres(:,:,:) - SUM(pres(:,:,dim_k-1:,:),3)/2.


!mass in mm
do j = 1,dim_j
  do i = 1,dim_i
    do t = 1,MM5totsteps
      tpw(j,i,t) = SUM(dp(j,i,pbl_lev(j,i,t):,t)*mix(j,i,pbl_lev(j,i,t):,t)/g)
    end do
  end do
end do


END SUBROUTINE calc_tpw_pbl

!**************************************************************************


SUBROUTINE calc_pot_temp(temp,pres,pot_temp)
!------------------------------------------------
! calculate the potential temperature at every level and time
!-----------------------------------------------------

USE global_data

IMPLICIT NONE

REAL, INTENT(IN), DIMENSION(:,:,:,:) :: temp,pres

REAL, INTENT(OUT), DIMENSION(:,:,:,:) :: pot_temp


!
! calculate the potential temperature
!
pot_temp = temp * ((P0/pres)**(Rd/Cp))

END SUBROUTINE calc_pot_temp


!**********************************************************************


SUBROUTINE calc_eq_pot_temp(mix,mixtot,temp,pres,eq_pot_temp)
!------------------------------------------------
! calculate the equivalent potential temperature at every level and time
!-----------------------------------------------------

USE global_data

IMPLICIT NONE

REAL, INTENT(IN), DIMENSION(:,:,:,:) :: mix,mixtot,temp,pres

REAL, INTENT(OUT), DIMENSION(:,:,:,:) :: eq_pot_temp


REAL, DIMENSION(dim_j,dim_i,dim_k,MM5totsteps) :: es,mix_s,pot_temp


!
! use Teten's formula to calculate the saturation vapor pressure (in stull)
!
es = 611 * exp((17.2694*(temp - 273.16))/(temp - 35.86))

!
! calculate the saturated mixing ratio everywhere
!
mix_s = (Rd/Rv)*es/(pres - es)
  
!
! calculate the equivalent potential temperature (in atmospheric convection, Emanuel)
!
eq_pot_temp = temp * ((P0/pres)**(Rd/(Cp+Cl*mix))) * &
		mix_s**(-1.*mix*Rv/(Cp+Cl*mixtot)) * &
		exp(Lv*mix/((Cp+Cl*mixtot)*temp))

!
!calculate the equivalent potential temperature (in atmospheric science..., Wallace & Hobbs)
!
!woops this isn't right - need to know mix_s at the
!lifting condensation level!!
!
!call calc_pot_temp(temp,pres,pot_temp)

!eq_pot_temp = pot_temp*exp(Lv*mix_s/(Cp*temp))


END SUBROUTINE calc_eq_pot_temp


!**********************************************************************


SUBROUTINE calc_pbl_lev(pbl_hgt,pres,surf_pres,pbl_lev)
!------------------------------------------------
! calculate the model level just above the pbl height
!-----------------------------------------------------

USE global_data

IMPLICIT NONE

REAL, INTENT(IN), DIMENSION(:,:,:,:) :: pres
REAL, INTENT(IN), DIMENSION(:,:,:) ::pbl_hgt,surf_pres

INTEGER, INTENT(OUT), DIMENSION(:,:,:) :: pbl_lev


REAL, DIMENSION(dim_j,dim_i,MM5totsteps) :: pbl_pres
INTEGER :: j,i,t
INTEGER,DIMENSION(1) :: dummy_min


!
! calculate pressure at the pbl height using the hydrostatic equation
! here I assume that the density averages 1kg m-3 in the pbl
!
pbl_pres = -1*pbl_hgt*g + surf_pres
  
!
! calculate the model level just above the pbl height
! also add the level above that as it gains moisture by detrainment from the 
! PBL and so this moisture can also be associated with the current
! location
!
do j = 1,dim_j
  do i = 1,dim_i
    do t = 1,MM5totsteps
      dummy_min = MINLOC(abs(pbl_pres(j,i,t) - pres(j,i,:,t)))
      if ((pbl_pres(j,i,t) - pres(j,i,dummy_min(1),t)) < 0.) then
        pbl_lev(j,i,t) = dummy_min(1) - 2 
      else
        pbl_lev(j,i,t) = dummy_min(1) - 1
      end if
    end do
  end do
end do

END SUBROUTINE calc_pbl_lev


!**********************************************************************

SUBROUTINE near_pt(lon2d,lat2d,lon,lat,x,y)
!---------------------------------------------------------
!calculate the grid point nearest the lat and lon location
!------------------------------------------------------

USE global_data

IMPLICIT NONE

REAL, INTENT(IN), DIMENSION(:,:) :: lon2d,lat2d
REAL, INTENT(IN) :: lon,lat
INTEGER, INTENT(OUT) :: x,y

REAL, DIMENSION(SIZE(lon2d(:,1)),SIZE(lon2d(1,:))) :: dist
INTEGER, DIMENSION(2) :: loc


!
!calculate the distance from the parcel location to every grid point
!must account for changing distance between longitude lines as latitude changes
!
dist = SQRT((lat2d-lat)**2 + (cos(lat2d*pi/180)*(lon2d-lon))**2)

loc = MINLOC(dist)

x = loc(1)
y = loc(2)

END SUBROUTINE near_pt

!*************************************************************************

SUBROUTINE bilin_interp(var2d,lon2d,lat2d,par_lon,par_lat,var)
!------------------------------------------------------------------
! find the bi-linearly interpolated value at par_lon,par_lat
! lon and lat are not regularly spaced grids
!------------------------------------------------------------------

USE global_data

IMPLICIT NONE

REAL, INTENT(IN), DIMENSION(:,:) :: var2d,lon2d,lat2d
REAL, INTENT(IN) :: par_lon,par_lat
REAL, INTENT(OUT) :: var

REAL :: fac,t,u
INTEGER :: xx,yy,xbdy,ybdy
LOGICAL :: changex,changey

changex = .FALSE.
changey = .FALSE.
xbdy = 0
ybdy = 0


call near_pt(lon2d,lat2d,par_lon,par_lat,xx,yy)

!
!check if we are currently exactly on a grid pt
!
If (lon2d(xx,yy)==par_lon.AND.lat2d(xx,yy)==par_lat) then
  var = var2d(xx,yy)
  RETURN
end if

!
!get indexes of closest grid value to south and west
!be careful of boundary
!
if (xx==1) xbdy = -1
if (yy==1) ybdy = -1

if (lon2d(xx,yy)>par_lon.AND.xbdy>-1) then
  changex = .TRUE.
end if

if (lat2d(xx,yy)>par_lat.AND.ybdy>-1) then
  changey = .TRUE.
end if

if (changex) then
  xx = xx-1
end if
if (changey) then
  yy = yy-1
end if

!
!if we are at the top or right boundary
!
if (xx==ssdim) then
  xbdy = 1
end if
if (yy==ssdim) then
  ybdy = 1
end if


!
!check to see if point in inside lower and left boundaries
!
if (xbdy<0.AND.lon2d(xx,yy)<par_lon) xbdy = 0
if (ybdy<0.AND.lat2d(xx,yy)<par_lat) ybdy = 0

!
!calculate the t and u weights
!account for possibility of being outside the boundaries
!even though this should never happen
!if outside a boundary use the value at the boundary
!
if ((xbdy/=0.AND.ybdy/=0)) then
  var = var2d(xx,yy)  
else if (xbdy/=0) then
  fac = (par_lat - lat2d(xx,yy))/(lat2d(xx,yy+1)-lat2d(xx,yy))
  var = lin_interp(var2d(xx,yy:yy+1),fac)
else if (ybdy/=0) then
  fac = (par_lon - lon2d(xx,yy))/(lon2d(xx+1,yy)-lon2d(xx,yy))
  var = lin_interp(var2d(xx:xx+1,yy),fac)
else
  t = (par_lon - lon2d(xx,yy))/(lon2d(xx+1,yy)-lon2d(xx,yy))
  u = (par_lat - lat2d(xx,yy))/(lat2d(xx,yy+1)-lat2d(xx,yy))
  var = (1-t)*(1-u)*var2d(xx,yy) + (1-t)*u*var2d(xx,yy+1) + &
         t*u*var2d(xx+1,yy+1) + t*(1-u)*var2d(xx+1,yy)
end if	 

  

END SUBROUTINE bilin_interp


!**************************************************************************************

SUBROUTINE new_parcel_level_pt(par_pot_temp,pot_temp,par_lev,lev)
!-------------------------------------------------------------------------
!calculate what level has the parcel potential temperature at this location
!--------------------------------------------------------------------------

USE global_data

IMPLICIT NONE

REAL, INTENT(IN), DIMENSION(:) :: pot_temp
INTEGER, INTENT(IN) :: par_lev
REAL, INTENT(INOUT) :: par_pot_temp
INTEGER, INTENT(OUT) :: lev

REAL, DIMENSION(dim_k) :: pot_dist
REAL :: pot_min
INTEGER :: xx,yy,kk
INTEGER, DIMENSION(1) :: lev_dummy
INTEGER, DIMENSION(2) :: lev_dis
LOGICAL :: getdownmin, getupmin


getdownmin = .TRUE.
getupmin = .TRUE.


!
!adjust parcel potential temperature if required (ie. driven into ground)
!
!print *,'1',par_pot_temp,xx,yy
pot_min = MINVAL(pot_temp)
par_pot_temp = MAX(par_pot_temp,pot_min)
!print *,'2',pot_temp(xx,yy,:)
!
!calculate nearest potential temperature in vertical column
!
pot_dist = pot_temp - par_pot_temp

!
!since equivalent pot temperature is not monotonic I need
!to find the level closest to previous level with right potential temp
!
lev_dis = dim_k*2

do kk = 1, dim_k
  if (par_lev-kk>0 .AND. getdownmin) then
    if (pot_dist(par_lev-kk)<0.AND.pot_dist(par_lev-kk+1)>0.OR. &
    	pot_dist(par_lev-kk)>0.AND.pot_dist(par_lev-kk+1)<0.) then
      lev_dummy = MINLOC(ABS(pot_dist(par_lev-kk:par_lev-kk+1)))
      lev_dis(1) = -kk-1+lev_dummy(1)
      getdownmin = .FALSE.
    end if
  end if
  
  if (par_lev+kk<dim_k+1 .AND. getupmin) then
    if (pot_dist(par_lev+kk)<0.AND.pot_dist(par_lev+kk-1)>0.OR. &
    	pot_dist(par_lev+kk)>0.AND.pot_dist(par_lev+kk-1)<0.) then
      lev_dummy = MINLOC(ABS(pot_dist(par_lev+kk-1:par_lev+kk)))
      lev_dis(2) = kk-2+lev_dummy(1)
      getupmin = .FALSE.
    end if
  end if
end do
      
!print *,lev_dis,SUM(lev_dis),par_pot_temp
      
      
if (SUM(lev_dis)==dim_k*4) then
  lev_dummy = MINLOC(ABS(pot_dist))
  lev = lev_dummy(1)
else
  lev_dummy = MINLOC(ABS(lev_dis))
  lev = lev_dis(lev_dummy(1)) + par_lev
end if

if (lev==0) then
  print *,"new_parcel_pt - pot dist",pot_dist
  print *,"new_parcel_pt - par_lev,lev",par_lev,lev
end if

END SUBROUTINE new_parcel_level_pt



!**************************************************************************************

SUBROUTINE new_parcel_level_w(par_pres,pres,w,temp,mix,lev)
!-------------------------------------------------------------------------
!calculate the new parcel level given w at this location
!--------------------------------------------------------------------------

USE global_data

IMPLICIT NONE

REAL, INTENT(IN) :: w,temp,mix
REAL, INTENT(IN), DIMENSION(:) :: pres
REAL, INTENT(INOUT) :: par_pres
INTEGER, INTENT(OUT) :: lev

INTEGER, DIMENSION(1) :: dummy_lev

!
! Here I use the hydrostatic eqn to calculate the change in pressure given w.
! deltaP = rho*g*deltaz when in hydrostatic equilibrium
!
par_pres = par_pres + -1.*(par_pres/(Rd*(1+0.61*mix)*temp))*g*w*tstep*60


! Find the model level where the difference in pressure between the parcel
! and the atmosphere is the smallest, i.e. which height in pres does the 
! smallest difference occur, where pres dims are (lat,lon,height).
dummy_lev = MINLOC(ABS(pres - par_pres))

lev = dummy_lev(1)

!if the parcel is below the lowest model level then set it to the lowest level
if (par_pres > MAXVAL(pres)) par_pres = MAXVAL(pres)

!print *,'L2188, par_pres=',par_pres

if (lev==0) then
  print *,'par_lev_w - pres_dis',(pres - par_pres),temp,w
  print *,'par_lev_w -',par_pres,pres
end if

END SUBROUTINE new_parcel_level_w


!***************************************************************************************

SUBROUTINE advect(u,v,lon,lat)
!------------------------------------
!linear advection in direction given by u and v
!u and v are m/s
!lon and lat are in degrees east and north
!time step is given by tstep from global_data
!---------------------------------------------

USE global_data

IMPLICIT NONE

REAL, INTENT(IN) :: u,v
REAL, INTENT(INOUT) :: lon,lat


!
!calculate the new lat
!
lat = lat + v*tstep*60/(deg_dist*1000)

!
!calculate new lon
!
lon = lon + u*tstep*60/(cos(lat*pi/180)*deg_dist*1000)
! This assumes the input data ranges in lon between 0deg and 360deg. Our data ranges
! between -180deg and +180deg. If calc's lon is >180deg, translate it.
if (lon>=180) then
	lon = (lon + u*tstep*60/(cos(lat*pi/180)*deg_dist*1000)) - 360
end if 

END SUBROUTINE advect

!***********************************************************************************


SUBROUTINE implicit_back_traj(u,v,w,temp,pbl_lev,pot_temp,pres,lon2d,lat2d, &
				par_lon,par_lat,par_lev, &
				par_pot_temp,par_pres,par_q,thread)
!-------------------------------------------------------------------------------
! Using Merrill's fully implicit isentropic technique
! calculate the parcels position one time step before
!
!u,v,w should only have 2 time steps in them
!pot_temp & pres should only have the end time step
!--------------------------------------------------------------------------

USE global_data

IMPLICIT NONE

REAL, INTENT(IN), DIMENSION(:,:,:,:) :: u,v,w
REAL, INTENT(IN), DIMENSION(:,:,:) :: pot_temp,pres,temp
REAL, INTENT(IN), DIMENSION(:,:) :: lon2d,lat2d
INTEGER, INTENT(IN), DIMENSION(:,:) :: pbl_lev
REAL, INTENT(INOUT) :: par_lon,par_lat,par_pot_temp,par_pres
REAL, INTENT(IN) :: par_q
INTEGER, INTENT(INOUT) :: par_lev
INTEGER, INTENT(IN) :: thread

INTEGER :: xx,yy,ll,lev
REAL :: lon,lat,u_back,v_back,w_back,temp_back,u_for,v_for,w_for,pr
REAL :: pt1,pt2,vfac,pr1,pr2
INTEGER, DIMENSION(1) :: dummy_lev


!print *,'1st',par_pres,par_pot_temp,par_lon,par_lat,par_lev,thread

!
!get u and v at parcel location
!
lon = par_lon
lat = par_lat
call bilin_interp(u(:,:,par_lev,2),lon2d,lat2d,lon,lat,u_back)
call bilin_interp(v(:,:,par_lev,2),lon2d,lat2d,lon,lat,v_back)


!
!find lat and lon after advecting back in time
!
call advect(-1.*u_back,-1.*v_back,lon,lat)

!
!calculate which vertical level has correct potential temperature 
!at new location or that w moves us to
!
call near_pt(lon2d,lat2d,lon,lat,xx,yy)

!print *,'par_lev,pbl_lev',par_lev,pbl_lev(xx,yy),xx,yy,thread

call bilin_interp(temp(:,:,par_lev),lon2d,lat2d,par_lon,par_lat,temp_back)

!if (par_lev >= pbl_lev(xx,yy)) then
if (.TRUE.) then
  call new_parcel_level_pt(par_pot_temp,pot_temp(xx,yy,:),par_lev,lev)
else
  pr = par_pres
  call bilin_interp(w(:,:,par_lev,2),lon2d,lat2d,lon,lat,w_back)
  call new_parcel_level_w(pr,pres(xx,yy,:),w_back,temp_back,par_q,lev)
end if

!print *,'2nd',par_pres,par_pot_temp,par_lon,par_lat,par_lev,thread


!
!get u and v at new location
!
call bilin_interp(u(:,:,lev,1),lon2d,lat2d,lon,lat,u_for)
call bilin_interp(v(:,:,lev,1),lon2d,lat2d,lon,lat,v_for)

!
!find new location of parcel as mean location given by back and forward trajectory
!
call advect(-1.*(u_back+u_for)/2.,-1.*(v_back+v_for)/2.,par_lon,par_lat)


!
!find final parcel potential temperature and level
!
call near_pt(lon2d,lat2d,par_lon,par_lat,xx,yy)

!print *,'par_lev,pbl_lev',par_lev,pbl_lev(xx,yy),xx,yy,thread


!if (par_lev >= pbl_lev(xx,yy)) then
if (.TRUE.) then
  call new_parcel_level_pt(par_pot_temp,pot_temp(xx,yy,:),par_lev,lev)
  
  !print *,'lev_pt1',par_pres,par_pot_temp,pot_temp(xx,yy,lev-1:lev+1)
  
  
  !
  !need to calculate the new parcel pressure
  !need to be extra careful as pot_temp may not be monotonic
  !
  call bilin_interp(pres(:,:,lev),lon2d,lat2d,par_lon,par_lat,pr1)
  if (lev == dim_k .OR. lev == 1 .OR. par_pot_temp == pot_temp(xx,yy,lev)) then
    par_pres = pr1
  else
    !
    !in case this level is a min or max in pot_temp
    !
    if (par_pot_temp>pot_temp(xx,yy,lev).AND.par_pot_temp<pot_temp(xx,yy,lev-1).AND.par_pot_temp<pot_temp(xx,yy,lev+1)) then
    
      dummy_lev = MINLOC(ABS(par_pres - pres(xx,yy,lev-1:lev+1:2)))
      if (dummy_lev(1)==1) then
        dummy_lev(1) = lev-1
      else
        dummy_lev(1) = lev+1
      end if
      call bilin_interp(pres(:,:,dummy_lev(1)),lon2d,lat2d,par_lon,par_lat,pr2)
      vfac = (pot_temp(xx,yy,dummy_lev(1))-par_pot_temp)/(pot_temp(xx,yy,dummy_lev(1))-pot_temp(xx,yy,lev))
      par_pres = exp((1-vfac)*log(pr2) + vfac*log(pr1))
      
    else if (par_pot_temp<pot_temp(xx,yy,lev).AND.par_pot_temp>pot_temp(xx,yy,lev-1).AND.par_pot_temp>pot_temp(xx,yy,lev+1)) then
    
      dummy_lev = MINLOC(ABS(par_pres - pres(xx,yy,lev-1:lev+1:2)))
      if (dummy_lev(1)==1) then
        dummy_lev(1) = lev-1
      else
        dummy_lev(1) = lev+1
      end if
      call bilin_interp(pres(:,:,dummy_lev(1)),lon2d,lat2d,par_lon,par_lat,pr2)
      vfac = (pot_temp(xx,yy,lev)-par_pot_temp)/(pot_temp(xx,yy,lev)-pot_temp(xx,yy,dummy_lev(1)))
      par_pres = exp((1-vfac)*log(pr1) + vfac*log(pr2))
      
    !
    !in other cases
    !
    else if (par_pot_temp > pot_temp(xx,yy,lev) .AND. par_pot_temp < pot_temp(xx,yy,lev-1)) then
    
      call bilin_interp(pres(:,:,lev-1),lon2d,lat2d,par_lon,par_lat,pr2)
      vfac = (pot_temp(xx,yy,lev-1)-par_pot_temp)/(pot_temp(xx,yy,lev-1)-pot_temp(xx,yy,lev))
      par_pres = exp((1-vfac)*log(pr2) + vfac*log(pr1))
      
    else if (par_pot_temp < pot_temp(xx,yy,lev) .AND. par_pot_temp > pot_temp(xx,yy,lev+1)) then
    
      call bilin_interp(pres(:,:,lev+1),lon2d,lat2d,par_lon,par_lat,pr2)
      vfac = (pot_temp(xx,yy,lev)-par_pot_temp)/(pot_temp(xx,yy,lev)-pot_temp(xx,yy,lev+1))
      par_pres = exp((1-vfac)*log(pr1) + vfac*log(pr2))
    
    else if (par_pot_temp > pot_temp(xx,yy,lev) .AND. par_pot_temp < pot_temp(xx,yy,lev+1)) then
    
      call bilin_interp(pres(:,:,lev+1),lon2d,lat2d,par_lon,par_lat,pr2)
      vfac = (pot_temp(xx,yy,lev+1)-par_pot_temp)/(pot_temp(xx,yy,lev+1)-pot_temp(xx,yy,lev))
      par_pres = exp((1-vfac)*log(pr2) + vfac*log(pr1))
      
    else if (par_pot_temp < pot_temp(xx,yy,lev) .AND. par_pot_temp > pot_temp(xx,yy,lev-1)) then
    
      call bilin_interp(pres(:,:,lev-1),lon2d,lat2d,par_lon,par_lat,pr2)
      vfac = (pot_temp(xx,yy,lev)-par_pot_temp)/(pot_temp(xx,yy,lev)-pot_temp(xx,yy,lev-1))
      par_pres = exp((1-vfac)*log(pr1) + vfac*log(pr2))
      
    end if
  end if
  
  !print *,'lev_pt2',par_pres,par_pot_temp,pot_temp(xx,yy,lev-1:lev+1)
  
else
  call bilin_interp(w(:,:,par_lev,1),lon2d,lat2d,lon,lat,w_for)
  call new_parcel_level_w(par_pres,pres(xx,yy,:),(w_back+w_for)/2.,temp_back,par_q,lev)

  !need to calculate the new parcel potential temperature
  call bilin_interp(pot_temp(:,:,lev),lon2d,lat2d,par_lon,par_lat,pt1)
  if (lev == dim_k .OR. lev == 1 .OR. par_pres == pres(xx,yy,lev)) then
    par_pot_temp = pt1
  else
    if (par_pres > pres(xx,yy,lev)) then
      call bilin_interp(pot_temp(:,:,lev+1),lon2d,lat2d,par_lon,par_lat,pt2)
      vfac = (log(pres(xx,yy,lev+1))-log(par_pres))/(log(pres(xx,yy,lev+1))-log(pres(xx,yy,lev)))
      par_pot_temp = (1-vfac)*pt2 + vfac*pt1
    else
      call bilin_interp(pot_temp(:,:,lev-1),lon2d,lat2d,par_lon,par_lat,pt2)
      vfac = (log(pres(xx,yy,lev))-log(par_pres))/(log(pres(xx,yy,lev))-log(pres(xx,yy,lev-1)))
      par_pot_temp = (1-vfac)*pt1 + vfac*pt2
    end if
  end if
end if


if (lev==0) then
  print *,'L2389, ',par_lev,lev,par_pres,par_pot_temp,temp_back,thread
  STOP
end if

par_lev = lev



END SUBROUTINE implicit_back_traj
!***********************************************************************************


SUBROUTINE implicit_back_traj_w(u,v,w,temp,pres,lon2d,lat2d, &
				par_lon,par_lat,par_lev, &
				par_pres,par_q,thread)
!-------------------------------------------------------------------------------
! Using Merrill's fully implicit  technique
! calculate the parcels position one time step before
!
! u,v,w should only have 2 time steps in them
! pres should only have the end time step

! Output parcel lat,lon, height and pressure
!--------------------------------------------------------------------------

USE global_data

IMPLICIT NONE

REAL, INTENT(IN), DIMENSION(:,:,:,:) :: u,v,w
REAL, INTENT(IN), DIMENSION(:,:,:) :: pres,temp
REAL, INTENT(IN), DIMENSION(:,:) :: lon2d,lat2d
!INTEGER, INTENT(IN), DIMENSION(:,:) :: pbl_lev
REAL, INTENT(INOUT) :: par_lon,par_lat,par_pres
REAL, INTENT(IN) :: par_q
INTEGER, INTENT(INOUT) :: par_lev
INTEGER, INTENT(IN) :: thread

INTEGER :: xx,yy,ll,lev
REAL :: lon,lat,u_back,v_back,w_par,temp_par,u_for,v_for,pr
REAL :: pt1,pt2,vfac,pr1,pr2
INTEGER, DIMENSION(1) :: dummy_lev

!
! Get u and v at parcel location, by interpolating in time (current and previous parcel timestep), and
! space (parcel lat/lon and lat/lon of the nearest grid point).
!
lon = par_lon
lat = par_lat
call bilin_interp(u(:,:,par_lev,2),lon2d,lat2d,lon,lat,u_back)
call bilin_interp(v(:,:,par_lev,2),lon2d,lat2d,lon,lat,v_back)


!
!find lat and lon after advecting back in time
!
call advect(-1.*u_back,-1.*v_back,lon,lat)

!
!calculate which vertical level that w moves us to
!
call near_pt(lon2d,lat2d,lon,lat,xx,yy)

pr = par_pres
call bilin_interp(temp(:,:,par_lev),lon2d,lat2d,lon,lat,temp_par)
call bilin_interp(w(:,:,par_lev,2),lon2d,lat2d,lon,lat,w_par)
call new_parcel_level_w(pr,pres(xx,yy,:),-1.*w_par,temp_par,par_q,lev)


!
!get u and v at new location
!
call bilin_interp(u(:,:,lev,1),lon2d,lat2d,lon,lat,u_for)
call bilin_interp(v(:,:,lev,1),lon2d,lat2d,lon,lat,v_for)

!
!find new location of parcel as mean location given by back and forward trajectory
!
call advect(-1.*(u_back+u_for)/2.,-1.*(v_back+v_for)/2.,par_lon,par_lat)

!print *,'lat,lon,pres,lev,v,u,w',par_lat,par_lon,par_pres,par_lev,(v_back+v_for)/2.,(u_back+u_for)/2.,w_par

!
!find final parcel level
!
call near_pt(lon2d,lat2d,par_lon,par_lat,xx,yy)

call bilin_interp(temp(:,:,par_lev),lon2d,lat2d,par_lon,par_lat,temp_par)
call bilin_interp(w(:,:,par_lev,1),lon2d,lat2d,par_lon,par_lat,w_par)
call new_parcel_level_w(par_pres,pres(xx,yy,:),-1.*w_par,temp_par,par_q,lev)


if (lev==0) then
  print *,'here',par_lev,lev,par_pres,temp_par,thread
  STOP
end if

par_lev = lev



END SUBROUTINE implicit_back_traj_w


  
END MODULE bt_subs


!***********************************************************************
!***********************************************************************


PROGRAM back_traj

!------------------------------------------------------------
! This program calculates the quasi-isentropic back trajectories of
! water vapor using a method based on Dirmeyer & Brubaker, 1999,
! "Contrasting evaporative moisture sources during the drought of
! 1988 and the flood of 1993", Journal of Geophysical Research, 104 D16
! pg 19,383-19,397
!
! Input data are taken from NARCliM output (precip,evap (2D), temp,wind,
! mixing ratios (3D)). Data here:
! /srv/ccrc/data33/z3481416/CCRC-WRF3.6.0.5-SEB/ERA-Interim/R2_nudging/out/
! 
! The model outputs a daily 2d grid of water vapour contribution of each 
! grid cell to the precipitation in each cell where it rained. It also
! outputs the total rainfall depth in each cell where it rained, and the 
! coordinates of each rain grid cell.
! 
! The output file is dimensioned (rec,lat,lon) where rec is the number of 
! grid cells where it rained that day, lat & lon are the sizes of the domain in 
! the lat and lon directions. Attributes include location and time and 
! amount of precip in the pixel of interest (y,x,day,month,year,precip) 
! each of which is a 1D array of length rec.
!
! Here I have amended the program to read in a list times (in days since
! start of simulation) and calculate the back trajectories for the day
! centered around the given time
!
! Here I have changed it so that parcels follow equivalent potential
! temperature surfaces within the PBL and is advected usign w above that
!-------------------------------------------------------------------

USE netcdf
USE global_data
USE bt_subs
USE omp_lib


IMPLICIT NONE


!
!netcdf id variables
!
INTEGER :: status,headncid
INTEGER :: spid,ptopid,delxid,latcrsid,loncrsid,terid,tstepid,sdayid,smonid,syearid
INTEGER :: dimjid,dimiid,sigid,dimkid
INTEGER :: outncid,wvcid,wvc2id,xlocid,ylocid,dayid,opreid,latid,lonid
INTEGER :: wsncid, wsid
CHARACTER(LEN=100):: fname
CHARACTER(LEN=10):: fname_smon, fname_sday

!
!data variables
!
INTEGER :: par_lev
REAL :: ptop,delx,par_lat,par_lon,par_pres,par_q,new_par_q,end_precip 
INTEGER :: MM5tstep
REAL,ALLOCATABLE,DIMENSION(:,:) :: lat2d,lon2d
REAL,ALLOCATABLE,DIMENSION(:,:) :: terrain,WV_cont,WV_cont_day !pstar, > make 3d
REAL,ALLOCATABLE,DIMENSION(:,:) :: WV_cont_apbl,WV_cont_day_apbl
REAL,ALLOCATABLE,DIMENSION(:,:,:) :: precip
REAL,ALLOCATABLE,DIMENSION(:,:,:) :: evap,tpw,pbl_hgt,surf_pres,pstar,psfc
REAL,ALLOCATABLE,DIMENSION(:,:,:,:) :: u,v,w,temp,mix,pp,pb,pw,mixcld,mixtot,pres
!REAL,ALLOCATABLE,DIMENSION(:,:,:,:) :: pot_temp !
REAL,ALLOCATABLE,DIMENSION(:,:,:,:) :: unow,vnow,wnow
REAL,ALLOCATABLE,DIMENSION(:,:,:) :: pres_then,tempnow
!REAL,ALLOCATABLE,DIMENSION(:,:,:) :: pot_temp_then !
INTEGER,ALLOCATABLE,DIMENSION(:,:,:) :: pbl_lev

INTEGER,ALLOCATABLE,DIMENSION(:) :: par_release
INTEGER :: xx,yy,tt,nn,mm,npar,orec,x,y,i,ttMM5,nnMM5,ttMM5day
INTEGER :: xx_omp,threadnum,torec
REAL :: ttfac,nnfac,precip_here,qfac

INTEGER,ALLOCATABLE,DIMENSION(:,:) :: wsmask

INTEGER :: ssx,ssy

INTEGER :: iii,jjj,ttt

LOGICAL :: print_test

!I want dd and totpts to persist inside and outside parallel regions and subroutines
!so they have been added to the global_data module and declared threadprivate

!for outputting the parcel stats for each parcel
REAL,ALLOCATABLE,DIMENSION(:,:) :: parcel_stats


!----------------------------------------------------------------
!
! Retrieve simulation start and end dates from the command line input
!
INTEGER :: num_args, ix
character(len=12), dimension(:), allocatable :: args
num_args = command_argument_count()
allocate(args(num_args))  ! I've omitted checking the return status of the allocation

call get_command_argument(1,args(1))
sday = string_to_int(args(1))
call get_command_argument(2,args(2))
smon = string_to_int(args(2))
call get_command_argument(3,args(3))
syear = string_to_int(args(3))
call get_command_argument(4,args(4))
edday = string_to_int(args(4))
call get_command_argument(5,args(5))
edmon = string_to_int(args(5))
call get_command_argument(6,args(6))
edyear = string_to_int(args(6))



! Find total number of days to run simulation for, given input start and end dates
totdays=simlength(sday,smon,syear,edday,edmon,edyear)

print *, 'Start and end dates: ',sday,smon,syear,edday,edmon,edyear
print *,"Total days to run analysis=",totdays
print *,"Total number of back-track days=",totbtadays
print *,"Number of parcels=",nparcels
print *,"Simulation time step (mins)=",tstep

stop

! Get header info from first input file -------------------------------------------

if (smon<10) then
fname_smon="0"//TRIM(int_to_string(smon))
else
fname_smon=TRIM(int_to_string(smon))
end if

if (sday<10) then
fname_sday="0"//TRIM(int_to_string(sday))
else
fname_sday=TRIM(int_to_string(sday))
end if

fname="wrfout_d01_"//TRIM(int_to_string(syear))//"-"//TRIM(fname_smon)//"-"//TRIM(fname_sday)//"_00:00:00"
print *,'Get header info from first input file: ',fname
status = NF90_OPEN(TRIM(dirdata_atm)//fname, NF90_NOWRITE, headncid)
if (status /= NF90_NOERR) call handle_err(status)

!
! Get ids for required variables from header

status = nf90_inq_varid(headncid, "P_TOP", ptopid)  !top pressure in model (Pa)
if(status /= nf90_NoErr) call handle_err(status)
status = nf90_inquire_attribute(headncid, nf90_global, "DX", delxid)  !grid distance (m)
if(status /= nf90_NoErr) call handle_err(status)
status = nf90_inq_varid(headncid, "XLAT", latcrsid)  !latitudes of grid points (degrees)
if(status /= nf90_NoErr) call handle_err(status)
status = nf90_inq_varid(headncid, "XLONG", loncrsid)  !longitudes of grid points (degrees)
if(status /= nf90_NoErr) call handle_err(status)
status = nf90_inq_varid(headncid, "HGT", terid)  !model terrain (m)
if(status /= nf90_NoErr) call handle_err(status)
status = nf90_inq_dimid(headncid, "Time", tstepid)  !number of time steps in file
if(status /= nf90_NoErr) call handle_err(status)

!
! Read in 1d variables
!
status = nf90_get_var(headncid, ptopid, ptop)
if(status /= nf90_NoErr) call handle_err(status)
status = nf90_get_att(headncid, NF90_GLOBAL, "DX",delx)
if(status /= nf90_NoErr) call handle_err(status)
status = nf90_inquire_dimension(headncid, tstepid,len = MM5tstep) 
MM5tstep=1440/MM5tstep ! Value must be 1440/8=180, where 8 is number of timesteps in the file (it was set up this way based on MM5 input files, where MM5 model timestep was 180mins)
if(status /= nf90_NoErr) call handle_err(status)
status = nf90_get_att(headncid, NF90_GLOBAL, "SOUTH-NORTH_GRID_DIMENSION", fdim_i) 
if(status /= nf90_NoErr) call handle_err(status)
status = nf90_get_att(headncid, NF90_GLOBAL, "WEST-EAST_GRID_DIMENSION", fdim_j) 
if(status /= nf90_NoErr) call handle_err(status)
status = nf90_get_att(headncid, NF90_GLOBAL, "BOTTOM-TOP_GRID_DIMENSION", dim_k) 
if(status /= nf90_NoErr) call handle_err(status)

dim_k=dim_k-1

!switch from dots to crosses
fdim_i = fdim_i-1
fdim_j = fdim_j-1

!print *,fdim_j,fdim_i


!get i and j dimensions when ignoring the boundaries
dim_i = fdim_i - 2*(bdy-1)
dim_j = fdim_j - 2*(bdy-1)
!print *,"dim_i=",dim_i,"dim_j=",dim_j,"dim_k=",dim_k

!Total number of grid pts inside boundaries
totpts = (dim_j-2)*(dim_i-2)

! Allocate the required arrays
ALLOCATE( terrain(dim_j,dim_i),lon2d(dim_j,dim_i),lat2d(dim_j,dim_i), STAT = status) !sigma(dim_k),pstar(dim_j,dim_i,MM5daysteps), &
 
!
! Read in more variables
!
status = nf90_get_var(headncid, latcrsid, lat2d,start=(/bdy,bdy/),count=(/dim_j,dim_i/))
if(status /= nf90_NoErr) call handle_err(status)
status = nf90_get_var(headncid, loncrsid, lon2d,start=(/bdy,bdy/),count=(/dim_j,dim_i/))
if(status /= nf90_NoErr) call handle_err(status)
status = nf90_get_var(headncid, terid, terrain,start=(/bdy,bdy/),count=(/dim_j,dim_i/))
if(status /= nf90_NoErr) call handle_err(status)

status = nf90_close(headncid)


!--------------------------------------------------------

!
! Calculate the number of trajectory time steps in a day and in input file time step
!
daytsteps = 1440/tstep                ! number of sub-daily simulation time steps
inMM5tsteps = MM5tstep/tstep          ! divide input file time step by the number of simulation time steps, as they may differ
totsteps = daytsteps*(totbtadays+1)   ! total number of simulation data time steps to remember

MM5daysteps = 1440/MM5tstep           ! number of input file time steps in day
MM5totsteps = (MM5daysteps*(totbtadays+1)) + 1 ! total number of input file time steps over the back-track period

! Allocate the variable arrays
ALLOCATE( precip(dim_j,dim_i,MM5daysteps), &
  evap(dim_j,dim_i,MM5totsteps),u(dim_j,dim_i,dim_k,MM5totsteps), &
  v(dim_j,dim_i,dim_k,MM5totsteps),temp(dim_j,dim_i,dim_k,MM5totsteps), &
  mix(dim_j,dim_i,dim_k,MM5totsteps),pp(dim_j,dim_i,dim_k,MM5totsteps),pb(dim_j,dim_i,dim_k,MM5totsteps), &
  !pot_temp(dim_j,dim_i,dim_k,MM5totsteps), &
  mixtot(dim_j,dim_i,dim_k,MM5totsteps), &
  tpw(dim_j,dim_i,MM5totsteps),pw(dim_j,dim_i,dim_k,daytsteps+1), &
  mixcld(dim_j,dim_i,dim_k,MM5totsteps),pres(dim_j,dim_i,dim_k,MM5totsteps), &
  psfc(dim_j,dim_i,MM5totsteps),surf_pres(dim_j,dim_i,MM5totsteps), w(dim_j,dim_i,dim_k,MM5totsteps), &
  pbl_hgt(dim_j,dim_i,MM5totsteps),pbl_lev(dim_j,dim_i,MM5totsteps), STAT = status )

  
!
! Read in watershed mask if required
!
if (wshed) then
  fname=TRIM(diri)//"watershed/"//TRIM(fwshed)
  print *,'using wshed from',fname 
  ALLOCATE( wsmask(dim_j,dim_i), STAT = status )
  
  status = NF90_OPEN(TRIM(diri)//"watershed/"//TRIM(fwshed), NF90_NOWRITE, wsncid)
  if (status /= NF90_NOERR) call handle_err(status)
  
  status = nf90_inq_varid(wsncid, "wsmask", wsid)  !watershed mask
  if(status /= nf90_NoErr) call handle_err(status)
  
  status = nf90_get_var(wsncid, wsid, wsmask,start=(/bdy,bdy/),count=(/dim_j,dim_i/))
  if(status /= nf90_NoErr) call handle_err(status)
  
  status = nf90_close(wsncid)
end if

! Total number of grid pts inside boundaries
totpts = (dim_j-2)*(dim_i-2)

! Set the number of threads to use in the parallel sections
call OMP_SET_NUM_THREADS(numthreads)
print *,'max threads=',omp_get_max_threads()


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! FOR EVERY DAY OF THE SIMULATION PERIOD
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

do dd = 1, totdays
  orec = 0
    
  !Get date to open correct input file   
  call day_month_year(dd) 
  print *,"L2850, day,mon,year",day,mon,year  
  
  ! Create output file (will create empty file even if it didn't rain anywhere in the domain on that day)
  call new_out_file(outncid,wvcid,wvc2id,xlocid,ylocid,dayid,opreid,day,lat2d,lon2d)
  print *,'L2860, created new out file'

  ! Get the variables required for back trajectory calculations for the current day
  call get_data(precip,evap,u,v,w,temp,mix,mixcld,mixtot,pp,pb,pbl_hgt,psfc)
  
  print *,'L2950, got data...'    

  ! Calculate the pressure field and surface pressure
  ! In previous (MM5-based) model, surface pressure had to be calculated - not so here, it is output from wrf.
  ! Pressure at heights is calculated by adding the base state pressure to the perturbation pressure at each level.
  ! pp="P"(4d)perturbation pressure. pb="PB"(4d)base state pressure. psfc="PSFC"(3d) surface pressure.
  pres = pp + pb
  surf_pres = psfc

  ! calculate the (equivalent) potential temperature field
  !call calc_pot_temp(temp,pres,pot_temp)
  !call calc_eq_pot_temp(mix,mixtot,temp,pres,pot_temp)
  !DEALLOCATE(temp)

  !calculate the model level just above the boundary layer height
  !call calc_pbl_lev(pbl_hgt,pres,surf_pres,pbl_lev)
  !print *, 'PBL level calculated'

  ! Calculate the precipitable water accumulated from the ground up on day of interest (lat,lon,height,time). 
  ! This is used to determine the parcel initial height.
  call calc_pw(mixcld(:,:,:,MM5totsteps-MM5daysteps:),pres(:,:,:,MM5totsteps-MM5daysteps:),surf_pres(:,:,MM5totsteps-MM5daysteps:),ptop,pw)

  ! Calculate the total precipitable water (lat,lon,time).
  call calc_tpw(mixtot,pres,surf_pres,ptop,tpw)  

  ! Calculate the subsection x & y dimensions, based on the max distance a parcel can travel in the sim timestep
  ssdim = (ceiling((sqrt(maxval(u)**2+maxval(v)**2)*tstep*60)/delx) *2) + 1

  !loop over x and y grid points
  !
  
  !parallelize over the grid points
  !ie each grid point will be sent to a new thread
  !program will wait until all grid points are complete before
  !moving on to next day
  
  !it remains unclear to me but it seems that I have to include
  !all calls to subroutines in critical sections even when they
  !don't do anything like changing values of shared variables etc
  !I can't actually find documentation to confirm this but it doesn't 
  !seem to work otherwise!!!????
  
  print *, 'L3011, Starting parallelisation'
  !$OMP PARALLEL DEFAULT(PRIVATE) COPYIN(daytsteps,totsteps,inMM5tsteps,MM5daysteps,MM5totsteps,dim_i,dim_j,dim_k,sday,smon,syear,mon,year,day,dd,totpts,ssdim) SHARED(pw,tpw,u,v,w,pres,temp,surf_pres,evap,precip,mix,mixtot,pbl_lev,lat2d,lon2d,orec,outncid,wvcid,wvc2id,xlocid,ylocid,dayid,opreid,wsmask)
  !,daylist)  
  
  !!!!$OMP SHARED(pw,tpw,pot_temp,u,v,w,pres,temp,surf_pres,evap,precip,mix,mixtot,pbl_lev,lat2d,lon2d,orec,outncid,wvcid,wvc2id,xlocid,ylocid,dayid,opreid,wsmask)
  !,daylist)  
  
  !print *, 'L2722, allocating arrays for each thread'
  !print *, 'array dimensions',dim_j, dim_i, ssdim, dim_k, daytsteps
  !allocate these arrays for each thread
  ALLOCATE( WV_cont(dim_j,dim_i),WV_cont_day(dim_j,dim_i), &
  		WV_cont_apbl(dim_j,dim_i),WV_cont_day_apbl(dim_j,dim_i), &
		unow(ssdim,ssdim,dim_k,2),vnow(ssdim,ssdim,dim_k,2), &
  		!pot_temp_then(ssdim,ssdim,dim_k), &
		par_release(daytsteps), &
		pres_then(ssdim,ssdim,dim_k),wnow(ssdim,ssdim,dim_k,2), &
		tempnow(ssdim,ssdim,dim_k), &
		STAT = status)
  !print *, 'L2732, arrays allocated'

  
  if (eachParcel) then
    ALLOCATE(parcel_stats(10,totsteps), STAT = status)
  end if
		
  !$OMP DO &
  !$OMP SCHEDULE (DYNAMIC)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! FOR EVERY POINT IN THE GRID
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  do xx_omp = 0, totpts-1
  
    xx = 2 + AINT(xx_omp*1./(dim_i-2))
    yy = 2 + (xx_omp - (xx-2)*(dim_i-2))
    
    threadnum = OMP_GET_THREAD_NUM()
    !print *,"start",threadnum,xx_omp,xx,yy
    
      ! Only do something if within watershed, if we care about the watershed
      if (wshed) then
        if (wsmask(xx,yy)==0) CYCLE
      end if

      ! Only do something if rain fell at this point on this day
      if (SUM(precip(xx,yy,:))>minpre) then
       
        !$OMP CRITICAL (output_index)
        orec = orec + 1
		torec = orec
	
! 	if (eachParcel) then
! 	  OPEN(unit=threadnum+10,file="parcel"//TRIM(int_to_string(dd))//"_"//TRIM(int_to_string(orec)), &
! 	  	form="UNFORMATTED",status="REPLACE") 
! 	  !print *,threadnum+10
! 	end if 
	!$OMP END CRITICAL (output_index)
	
	WV_cont_day = 0.
	!WV_cont_day_apbl = 0.    

	! Do you want to print out statements for debugging etc?
	if (torec==0) then
	  print_test = .TRUE.
	else
	  print_test = .FALSE.
	end if

	!
    ! Determine how many parcels to release today and use precip 
    ! distribution to determine when to release parcels.
    ! The globally set nparcels is just a maximum number of parcels
    ! to release in each data timestep. We release at least one parcel
    ! per simulation timestep within a data time step when it rained.
    !
	par_release = 0
	
	!$OMP CRITICAL (par_rel_time)
	if (COUNT(MASK = precip(xx,yy,:)>0.)<(nparcels/inMM5tsteps)) then
	  npar = COUNT(MASK = precip(xx,yy,:)>0.) * inMM5tsteps
	  call parcel_release_time(precip(xx,yy,:),npar,par_release)
	else
	  npar = nparcels
	  call parcel_release_time(precip(xx,yy,:),npar,par_release)
	end if
	!$OMP END CRITICAL (par_rel_time)
	
	!!!!!!!!!!
	!testing so remove randomness
	!par_release = 0
	!par_release(1) = 1
	!!!!!!!!!!!!
	!print *,'par_release ',par_release,threadnum,yy,xx
	

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! FOR EVERY SIMULATION SUB-DAILY TIMESTEP 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	do tt = 1, daytsteps
	  if (par_release(tt)==0) then
	    CYCLE
	  end if
  
  
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! FOR EVERY PARCEL RELEASE TIME
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	  do mm = 1, par_release(tt)
	  
	    WV_cont = 0.
	    !WV_cont_apbl = 0.
	    qfac = 1.
	    x = xx
	    y = yy
	    
	  
	    !
	    !the MM5 time step before this parcel time step on the rain day
	    !
	    ttMM5day = INT((tt-1)/inMM5tsteps) + 1
	    !print *,'L3019, mm,ttMM5day=',mm,ttMM5day
	    
	    !
	    !the MM5 time step from the beginning of the loaded files
	    !
	    ttMM5 = MM5totsteps - MM5daysteps - 1 + ttMM5day
	    !print *,'L3025, mm,ttMM5=',mm,ttMM5
	  
	    !
	    !factor for linear interpolation to parcel time step
	    !
	    ttfac = MOD(tt,inMM5tsteps)*1./inMM5tsteps
	  
	    !
	    !the precip produced here at this parcel time step
	    !
	    end_precip = precip(xx,yy,ttMM5day)/inMM5tsteps
	    
	    !
	    !determine sigma level from which to release parcel
	    !
	    !$OMP CRITICAL (par_rel_height)
	    call parcel_release_height(pw(xx,yy,:,tt),par_lev)
	    !$OMP END CRITICAL (par_rel_height)
	
	

	    !testing so remove randomness
	    !par_lev = 20

	    
	    if (print_test) then
	    !print *,"day,par_lev,xx,yy,tt ",dd,par_lev,xx,yy,tt,threadnum
	    end if
	    
	    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	    
	    !release the parcel and track the back trajectory 
	    !until 99% of initial precipitable water is accounted for
	    !or the parcel leaves the domain
	    !or the user specified time to calculate the back trajectory runs out
	    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	    !
	    par_lat = lat2d(xx,yy)
	    par_lon = lon2d(xx,yy)

  
	    !$OMP CRITICAL (par_q1)
	    !calculate the parcel mixing ratio
	    par_q = lin_interp(mixtot(xx,yy,par_lev,ttMM5:ttMM5+1),ttfac) 
	    !$OMP END CRITICAL (par_q1)
	    
	   
	    !$OMP CRITICAL (pot_temp1)
	    !par_pot_temp = lin_interp(pot_temp(xx,yy,par_lev,ttMM5:ttMM5+1),ttfac)
	    !$OMP END CRITICAL (pot_temp1)
	    
	    
	    !$OMP CRITICAL (par_pres1)
	    par_pres = lin_interp(pres(xx,yy,par_lev,ttMM5:ttMM5+1),ttfac)
	    !$OMP END CRITICAL (par_pres1)
	    
	    
	     if (print_test) then
	     !print *,"par_q1",xx,yy,par_lev,par_q,ttMM5,ttfac,par_pot_temp,par_pres,threadnum
	     end if
	  
	     
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! FOR EACH PARCEL RELEASE TIME, FOR EACH SIMULATION TIME STEP IN THE 
! WHOLE BACK-TRACK PERIOD
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
	    do nn = totsteps-daytsteps+tt, 2, -1
	    !print *,'tt,(totsteps-daytsteps+tt),nn=',tt,totsteps-daytsteps+tt,nn
	      !
	      !advect the parcel back in time one step
	      !

	      
	      !current parcel stats
	      if (eachParcel) then
	        !print *,"nn ",nn,threadnum,par_lev
	        parcel_stats(1,totsteps-daytsteps+tt+1-nn) = nn*1.
			parcel_stats(2,totsteps-daytsteps+tt+1-nn) = par_lon
			parcel_stats(3,totsteps-daytsteps+tt+1-nn) = par_lat
			parcel_stats(4,totsteps-daytsteps+tt+1-nn) = par_pres
			parcel_stats(5,totsteps-daytsteps+tt+1-nn) = par_q
			parcel_stats(6,totsteps-daytsteps+tt+1-nn) = u(x,y,par_lev,ttMM5)
			parcel_stats(7,totsteps-daytsteps+tt+1-nn) = v(x,y,par_lev,ttMM5)
			parcel_stats(8,totsteps-daytsteps+tt+1-nn) = w(x,y,par_lev,ttMM5)
		end if
	    
	      !
	      !calculate the lower left location for the subsection
	      !
	      if (x+floor(ssdim/2.)>dim_j) then
	        ssx = dim_j - ssdim + 1
	      else
	        ssx = max(x-floor(ssdim/2.),1)
	      end if
	      
	      if (y+floor(ssdim/2.)>dim_i) then
	        ssy = dim_i - ssdim + 1
	      else
	        ssy = max(y-floor(ssdim/2.),1)
	      end if
	   
	      
	      ! Find where you are in the simlength timeseries
	      nnMM5 = INT(nn/inMM5tsteps) + 1
	      nnfac = MOD(nn,inMM5tsteps)*1./inMM5tsteps  
	      !print *,'L3186,nn,nnfac=',nn,nnfac
	      
	      

	      !
	      !get u,v and pot_temp for this and the previous parcel time step (just subsection)
	      !

          !$OMP CRITICAl (unow1)
	      unow(:,:,:,2) = lin_interp3D(u(ssx:ssx+ssdim-1,ssy:ssy+ssdim-1,:,nnMM5:nnMM5+1),nnfac)
	      !$OMP END CRITICAl (unow1)

	      !$OMP CRITICAl (vnow1)
	      vnow(:,:,:,2) = lin_interp3D(v(ssx:ssx+ssdim-1,ssy:ssy+ssdim-1,:,nnMM5:nnMM5+1),nnfac)
	      !$OMP END CRITICAl (vnow1)
	      
	      !$OMP CRITICAl (wnow1)
	      wnow(:,:,:,2) = lin_interp3D(w(ssx:ssx+ssdim-1,ssy:ssy+ssdim-1,:,nnMM5:nnMM5+1),nnfac)
	      !$OMP END CRITICAl (wnow1)
	      
	      !$OMP CRITICAl (tempnow1)
	      tempnow(:,:,:) = lin_interp3D(temp(ssx:ssx+ssdim-1,ssy:ssy+ssdim-1,:,nnMM5:nnMM5+1),nnfac)
	      !$OMP END CRITICAl (tempnow1)
	      
	    
		!if (nn<totsteps-daytsteps+tt) then
		!  print *,nnMM5,nnfac,unow(x,y,par_lev,2),vnow(x,y,par_lev,2)
		!end if

		  ! Find where you are in the nn timeseries
	      nnMM5 = INT((nn-1)/inMM5tsteps) + 1
	      nnfac = MOD(nn-1,inMM5tsteps)*1./inMM5tsteps
	      !print *,'L3217,nn,nnfac=',nn,nnfac
	      

          !$OMP CRITICAl (unow2)
	      unow(:,:,:,1) = lin_interp3D(u(ssx:ssx+ssdim-1,ssy:ssy+ssdim-1,:,nnMM5:nnMM5+1),nnfac)
	      !$OMP END CRITICAl (unow2)
	      
	      !$OMP CRITICAl (vnow2)
	      vnow(:,:,:,1) = lin_interp3D(v(ssx:ssx+ssdim-1,ssy:ssy+ssdim-1,:,nnMM5:nnMM5+1),nnfac)
	      !$OMP END CRITICAl (vnow2)
	      
	      !$OMP CRITICAl (wnow2)
	      wnow(:,:,:,1) = lin_interp3D(w(ssx:ssx+ssdim-1,ssy:ssy+ssdim-1,:,nnMM5:nnMM5+1),nnfac)
	      !$OMP END CRITICAl (wnow2)
	      	      
	      !$OMP CRITICAl (pot_temp2)
	      !pot_temp_then(:,:,:) = lin_interp3D(pot_temp(ssx:ssx+ssdim,ssy:ssy+ssdim,:,nnMM5:nnMM5+1),nnfac)
	      !$OMP END CRITICAl (pot_temp2)
	      
	      !$OMP CRITICAl (presthen2)
	      pres_then(:,:,:) = lin_interp3D(pres(ssx:ssx+ssdim-1,ssy:ssy+ssdim-1,:,nnMM5:nnMM5+1),nnfac)
	      !$OMP END CRITICAl (presthen2)
	      
	  
	      !$OMP CRITICAL (traj)
	      
	      !!!!!!!!!!!!!!!!!!!!!!
	      ! EXPERIMENT 0, WHERE PARCELS FOLLOW POTENTIAL TEMP??
	      
	      !call implicit_back_traj(unow,vnow,wnow,tempnow,pbl_lev(ssx:ssx+ssdim,ssy:ssy+ssdim,nnMM5), &
	      !				pot_temp_then,pres_then, &
	    	!			lon2d(ssx:ssx+ssdim,ssy:ssy+ssdim),lat2d(ssx:ssx+ssdim,ssy:ssy+ssdim), &
		!			par_lon,par_lat,par_lev,par_pot_temp,par_pres,par_q,threadnum)
	      !$OMP END CRITICAL (traj)
	      
	      !!!!!!!!!!!!!!!!!!!!!!
	      !EXPERIMENT 4?, WHERE PARCELS FOLLOW EQUIVALENT POT TEMP AND PBL IS SPLIT??
	      
	      
	      !$OMP CRITICAL (trajw)
	      !call implicit_back_traj_w(unow,vnow,wnow,tempnow,pbl_lev(ssx:ssx+ssdim-1,ssy:ssy+ssdim-1,nnMM5), &
	      !				pres_then, &
	    !				lon2d(ssx:ssx+ssdim-1,ssy:ssy+ssdim-1),lat2d(ssx:ssx+ssdim-1,ssy:ssy+ssdim-1), &
		!			par_lon,par_lat,par_lev,par_pres,par_q,threadnum)
		
		!!!!!!!!!!!!!!!!!!!!!!
		! EXPERIMENT 3, WHERE PARCELS FOLLOW TEMP AND PBL IS NOT SPLIT??
		
	      call implicit_back_traj_w(unow,vnow,wnow,tempnow, &
	      				pres_then, &
	    				lon2d(ssx:ssx+ssdim-1,ssy:ssy+ssdim-1),lat2d(ssx:ssx+ssdim-1,ssy:ssy+ssdim-1), &
					par_lon,par_lat,par_lev,par_pres,par_q,threadnum)
	      !$OMP END CRITICAL (trajw)
	      			
	
	      !
	      ! Find the grid cell nearest the new lat,lon of the parcel 
	      !
	      !$OMP CRITICAL (near)
	      call near_pt(lon2d,lat2d,par_lon,par_lat,x,y)
	      !$OMP END CRITICAL (near)
	    
	      !
	      ! Find the water mass contribution of the new grid square, at this time
	      !    
	      !$OMP CRITICAL (new_par_q1)
	      new_par_q = lin_interp(mixtot(x,y,par_lev,nnMM5:nnMM5+1),nnfac)
	      !$OMP END CRITICAL (new_par_q1)

	      !
	      !adjust the q reduction factor if we had a decrease in q
	      !so long as it isn't the first time step.
	      ! i.e. If the parcel q decreases backward in time, then the q at the current time step must not have come
	      ! from the cell evap...maybe from some other process like convection.
	      !
	      if (nn < totsteps-daytsteps+tt) then
			if (new_par_q+min_del_q < par_q) then
				qfac = MAX(qfac*(1-(par_q-new_par_q)/par_q),0.)
			end if
	      end if

	     if (print_test) then
	     !print *,"par_q",x,y,par_lev,nnMM5,nnfac,par_q,new_par_q,qfac,par_pot_temp,par_pres,threadnum
	     end if

          !was moisture contributed to the parcel?
	      !is the parcel in the pbl?
	      !$OMP CRITICAL (wv_cont1)
	      !if (par_lev >= pbl_lev(x,y,nnMM5+1)) then
		if (lin_interp(evap(x,y,nnMM5:nnMM5+1),nnfac) > 0.) then
		!print *,'L3249, evap=',lin_interp(evap(x,y,nnMM5:nnMM5+1),nnfac
		!print *,'L3250, wv_cont 1=',WV_cont(x,y)
			WV_cont(x,y) = WV_cont(x,y) + (lin_interp(evap(x,y,nnMM5:nnMM5+1),nnfac) &
					/ (inMM5tsteps*lin_interp(tpw(x,y,nnMM5:nnMM5+1),nnfac)))
		!print *,'L3383, wv_cont=',WV_cont(x,y),'evap_interp=',(lin_interp(evap(x,y,nnMM5:nnMM5+1),nnfac)),'tpw_interp=',(lin_interp(tpw(x,y,nnMM5:nnMM5+1),nnfac))
		
		end if
	      !else

	      !  if (par_q < new_par_q-min_del_q) then
	      !    WV_cont_apbl(x,y) = WV_cont_apbl(x,y) + ((new_par_q - par_q)/par_q)*qfac
		!end if
		
	      !end if
	      !$OMP END CRITICAL (wv_cont1)

	      par_q = new_par_q
	      	      
	      
	      !print *,threadnum,x,y,par_lev,SUM(WV_cont),lin_interp(tpw(x,y,nnMM5:nnMM5+1),nnfac),lin_interp(evap(x,y,nnMM5:nnMM5+1),nnfac),unow(x,y,par_lev,1),vnow(x,y,par_lev,1)
	      if (print_test) then
	      !print *,"what",x,y,par_lev,SUM(WV_cont),SUM(WV_cont_apbl),threadnum
              end if
	     
	     
	      !saving parcel stats
	      if (eachParcel) then
	        !print *,"here"
	        parcel_stats(9,totsteps-daytsteps+tt+1-nn) = evap(x,y,ttMM5)
		parcel_stats(10,totsteps-daytsteps+tt+1-nn) = WV_cont(x,y)
              end if
	      
	      
	      !
	      !if we have accounted for all the precip  then go to next parcel
	      !
	      !if (SUM(WV_cont + WV_cont_apbl)>=1.) then
	      if (SUM(WV_cont)>=1.) then
	      !print *,"all precip accounted (torec,wv_cont) ",torec,SUM(WV_cont + WV_cont_apbl)
		!if (par_lev >= pbl_lev(x,y,nnMM5+1)) then
		  WV_cont(x,y) = WV_cont(x,y) - (SUM(WV_cont) - 1)
		!WV_cont(x,y) = WV_cont(x,y) - (SUM(WV_cont+WV_cont_apbl) - 1)
		!else
		!  WV_cont_apbl(x,y) = WV_cont_apbl(x,y) - (SUM(WV_cont+WV_cont_apbl) - 1)
		!end if
	        EXIT
	      end if
        
	      !
	      !if qfac = 0 then all of the increases in the water vapor 
	      !further back along the trajectory are lost before reaching 
	      !the end precipitation point so they don't contribute and
	      !there is no need to continue
	      !
	      !the water not accounted for must have come from convection
	      !or some other process that remains unaccounted for
	      !
	      if (qfac==0) then
	        !print *,"qfac = 0",tt,torec,par_lev,xx,yy,x,y,SUM(WV_cont+WV_cont_apbl),threadnum
	        EXIT
	      end if
	      
	     !
	      !if we have left the domain then assign the remaining precip to
	      !outside and go to next parcel
	      !
	      if (x<2) then
	        !print *,"hit boundary (orec,wv_cont)",torec,SUM(WV_cont)
	        WV_cont(1,y) = 1. - SUM(WV_cont)
	        EXIT
	      else if (x>dim_j-2) then
	        !print *,"hit boundary (orec,wv_cont)",torec,SUM(WV_cont)
	        WV_cont(dim_j,y) = 1. - SUM(WV_cont)
	        EXIT
	      end if
	      if (y<2) then
	        !print *,"hit boundary (orec,wv_cont)",torec,SUM(WV_cont)
	        WV_cont(x,1) = 1. - SUM(WV_cont)
	        EXIT
	      else if (y>dim_i-2) then
	        !print *,"hit boundary (orec,wv_cont)",torec,SUM(WV_cont)
	        WV_cont(x,dim_i) = 1. - SUM(WV_cont)
	        EXIT
	      end if
	     
	      !
	      !if we have left the domain then assign the remaining precip to
	      !outside and go to next parcel
	      !
	      !if (x<2) then
	        !print *,"hit boundary (orec,wv_cont)",torec,SUM(WV_cont)
		!if (par_lev >= pbl_lev(x,y,nnMM5+1)) then
	       !   WV_cont(1,y) = 1. - SUM(WV_cont+WV_cont_apbl)
		!else
		!  WV_cont_apbl(1,y) = 1. - SUM(WV_cont+WV_cont_apbl)
		!end if
	       ! EXIT
	      !else if (x>dim_j-2) then
	        !print *,"hit boundary (orec,wv_cont)",torec,SUM(WV_cont)
		!if (par_lev >= pbl_lev(x,y,nnMM5+1)) then
	      !    WV_cont(dim_j,y) = 1. - SUM(WV_cont+WV_cont_apbl)
		!else
		!  WV_cont_apbl(dim_j,y) = 1. - SUM(WV_cont+WV_cont_apbl)
		!end if
	      !  EXIT
	     ! end if
	     ! if (y<2) then
	        !print *,"hit boundary (orec,wv_cont)",torec,SUM(WV_cont)
		!if (par_lev >= pbl_lev(x,y,nnMM5+1)) then
	      !    WV_cont(x,1) = 1. - SUM(WV_cont+WV_cont_apbl)
		!else
		!  WV_cont_apbl(x,1) = 1. - SUM(WV_cont+WV_cont_apbl)
		!end if
	      !  EXIT
	     ! else if (y>dim_i-2) then
	        !print *,"hit boundary (orec,wv_cont)",torec,SUM(WV_cont)
		!if (par_lev >= pbl_lev(x,y,nnMM5+1)) then
	      !    WV_cont(x,dim_i) = 1. - SUM(WV_cont+WV_cont_apbl)
		!else
		!  WV_cont_apbl(x,dim_i) = 1. - SUM(WV_cont+WV_cont_apbl)
		!end if
	       ! EXIT
	     ! end if
	    
	      !
	      !if we have reached here and nn=2 then we have neither left the
	      !domain nor acounted for the precip
	      !
	      if (nn==2) then
	        !print *,"couldn't go back far enough in time.(day,time,orec,lev,xx,yy,x,y,%)",dd,tt,torec,par_lev,xx,yy,x,y,SUM(WV_cont+WV_cont_apbl),threadnum
	        !print *,"L3506, couldn't go back far enough in time.(day,time,orec,lev,xx,yy,x,y,%)",dd,tt,torec,par_lev,xx,yy,x,y,SUM(WV_cont),threadnum
	      
	        if (SUM(WV_cont)<0) STOP
		!if (SUM(WV_cont+WV_cont_apbl)<0) STOP
		
	      end if
          
	    end do   !nn loop
	    
	    WV_cont_day = WV_cont_day + WV_cont/npar
	    !WV_cont_day_apbl = WV_cont_day_apbl + WV_cont_apbl/npar
	  
	    if (print_test) then
	    !print *,'end2',sum(WV_cont_day),sum(WV_cont_day_apbl),SUM(WV_cont+WV_cont_apbl)
	    end if
	  
	    if (par_lev==0) STOP
	  
	  
	    !if keeping track of each parcel
	    if (eachParcel) then
	      !print *,"output",threadnum+10
	      WRITE(threadnum+10) parcel_stats
	      CLOSE(threadnum+10)
	    end if
	  
	  end do   !mm loop
	  
        end do	  !tt loop
	  
	!
	!write output to netcdf file
	!
	!$OMP CRITICAL (output)
	status = nf90_put_var(outncid,wvcid,WV_cont_day,start=(/1,1,torec/), &
				count=(/dim_j,dim_i,1/))
	if(status /= nf90_NoErr) call handle_err(status)
	!status = nf90_put_var(outncid,wvc2id,WV_cont_day_apbl,start=(/1,1,torec/), &
	!			count=(/dim_j,dim_i,1/))
	!if(status /= nf90_NoErr) call handle_err(status)
	status = nf90_put_var(outncid,xlocid,xx,start=(/torec/))
	if(status /= nf90_NoErr) call handle_err(status)
	status = nf90_put_var(outncid,ylocid,yy,start=(/torec/))
	if(status /= nf90_NoErr) call handle_err(status)
! 	status = nf90_put_var(outncid,dayid,daylist(dd)-1,start=(/torec/))
	status = nf90_put_var(outncid,dayid,dd,start=(/torec/)) 
 	if(status /= nf90_NoErr) call handle_err(status)
	status = nf90_put_var(outncid,opreid,SUM(precip(xx,yy,:)),start=(/torec/))
	if(status /= nf90_NoErr) call handle_err(status)
	!$OMP END CRITICAL (output)

! 	else
! 	print *,'No rain in the domain on this day'
	end if
	
	
    !end do 
    
    !*********************************!!!!!!!!
    !if (torec==1) then
    !  exit
    !end if
        
    !test
    !
    !EXIT
    !
    !test
    
  end do   !xx_omp loop
  !$OMP END DO NOWAIT
  
  !$OMP END PARALLEL 
  
  status = nf90_close(outncid)
  if(status /= nf90_NoErr) call handle_err(status)

end do


END PROGRAM back_traj
