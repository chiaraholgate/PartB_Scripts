! I've run this code on storm servers by doing the following:
! (1) 
! 	export MP_SHARED_MEMORY=yes
! 	export OMP_NUM_THREADS=1
! 	ulimit -s unlimited
! 	export FC='ifort'
! 	export NETCDFF='/share/apps/netcdf-f/intel/4.4.4'
! 	export NETCDFFLIBDIR=$NETCDFF/lib
! 	export NETCDFFINCDIR=$NETCDFF/include
! 	export NETCDFC='/share/apps/netcdf-c/intel/4.4.1.1' #'/share/apps/netcdf-c/4.4.1.1'
! 	export NETCDFCLIBDIR=$NETCDFC/lib
! 	export HDF5='/share/apps/hdf5/intel/1.10.1'  #'/share/apps/hdf5/1.10.1-threadsafe'
! 	export HDF5LIBDIR=$HDF5/lib
! 	export SZIP='/share/apps/szip/2.1'
! 	export SZIPLIBDIR=$SZIP/lib
! (2) ifort -I$NETCDFFINCDIR -c Check_read_input_data.f90  -g -traceback
! (3) ifort -o Check_read_input_data Check_read_input_data.o -L$NETCDFCLIBDIR -lnetcdf -L$NETCDFFLIBDIR -lnetcdff -shared-intel -L$HDF5LIBDIR -lhdf5_hl -lhdf5
! (4) ./Check_read_input_data 3 1 1979 4 1 1979 /srv/ccrc/data03/z3131380/PartB/test_data/

! I want to check that the shape of the data that's read in matches the array slicing order we have implemented in get_data.
! >> Outcome: 8/5/19 I confirm that the data is being input as it should be. I also confirm that the raw T data (pertubation potential temperature) is being correctly converted to actual temperature for use by the QIBT program.


MODULE global_data

IMPLICIT NONE

SAVE

!
!*******user modified variables**********************
!

INTEGER :: sday,smon,syear    !start day for calculations
INTEGER :: edday,edmon,edyear !end day for calculations (Exclusive. Must be at least one day after start day)
INTEGER :: totdays
INTEGER, PARAMETER :: totbtadays = 2   !number of days of data to keep for bta; i.e. how far back in time to calc.
                                       !must be less than days you have input data for

!INTEGER, PARAMETER :: totdays = 30    !total number of days to calculate for; must be <= no. days in days_of_rain.txt
INTEGER, PARAMETER :: tstep = 10   !number of minutes for back trajectory time step (simultion time step)
                      !must divide evenly into number of minutes in day 1440 and number of minutes in MM5 time step (here 180)
INTEGER, PARAMETER :: nparcels = 100   !set the number of parcels to release if it rains
REAL, PARAMETER :: minpre = 2   !min daily precip to deal with (mm)

INTEGER, PARAMETER :: bdy = 6   !boundary layers to ignore; trajectories will be tracked to this boundary

CHARACTER(LEN=50), PARAMETER :: diri = "/srv/ccrc/data03/z3131380/PartB/test_data/"   
!CHARACTER(LEN=100), PARAMETER :: diro = "/g/data/xc0/user/Holgate/QIBT/exp02/"
CHARACTER(LEN=100) :: diro  
CHARACTER(LEN=100), PARAMETER :: dirdata_atm = "/srv/ccrc/data03/z3131380/PartB/test_data/"
CHARACTER(LEN=100), PARAMETER :: dirdata_land = "/srv/ccrc/data03/z3131380/PartB/NARCliM_postprocess/"  

INTEGER, PARAMETER :: numthreads = 1   !set the number of parallel openmp threads

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

INTEGER :: daytsteps,totsteps,indatatsteps,datadaysteps,datatotsteps
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

COMMON /global_vars/ daytsteps,totsteps,indatatsteps,datadaysteps,datatotsteps, &
		dim_i,dim_j,dim_k,sday,smon,syear,mon,year,day,dd,totpts, &
		fdim_i,fdim_j,ssdim,diro

!$OMP THREADPRIVATE(/global_vars/)

END MODULE global_data



MODULE qibt_subs

IMPLICIT NONE

CONTAINS

SUBROUTINE handle_err(status)
!---------------------------------
! handle any errors from the fortran 90 netCDF interface
!---------------------------------
  
USE netcdf
  
IMPLICIT NONE
   
integer, intent (in) :: status
  
if(status /= nf90_noerr) then
  print *, trim(nf90_strerror(status))
  stop "Stopped"
end if
  
END SUBROUTINE handle_err

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

CHARACTER(LEN=50) FUNCTION real_to_string(num)
!----------------------------------------------
! converts an real to a character string 
! here I allow a maximum of 10 decimal places
!----------------------------------------------

IMPLICIT NONE

real, intent(in) :: num

! local 
integer :: i, whole, frac
real :: fracnum


whole = AINT(num)
fracnum = num - whole

do i=1,10
  if (MOD(fracnum,1.)==0) exit
  fracnum = fracnum * 10.
end do

frac = AINT(fracnum)

real_to_string = TRIM(int_to_string(whole))//"."//TRIM(int_to_string(frac))

end FUNCTION real_to_string

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

INTEGER FUNCTION julian(year,mon,day)

! From http://aa.usno.navy.mil/faq/docs/JD_Formula.php

IMPLICIT NONE
INTEGER, INTENT(IN) :: year,mon,day

julian= day-32075+1461*(year+4800+(mon-14)/12)/4+367*(mon-2-(mon-14)/12*12)/12-3*((year+4900+(mon-14)/12)/100)/4
   
END FUNCTION julian

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
    filename_ext_atm = "wrfout_d01_"//TRIM(int_to_string(yr))//"-0"//TRIM(int_to_string(mn))//"-0"//TRIM(int_to_string(d))//"_00:00:00.nc"
    filename_ext_RAIN = "wrfhrly_d01_"//TRIM(int_to_string(yr))//"-0"//TRIM(int_to_string(mn))//"-0"//TRIM(int_to_string(d))//"_00:00:00_RAIN.nc"
    filename_ext_LH = "wrfhrly_d01_"//TRIM(int_to_string(yr))//"-0"//TRIM(int_to_string(mn))//"-0"//TRIM(int_to_string(d))//"_00:00:00_LH.nc"
    filename_ext_P = "wrfhrly_d01_"//TRIM(int_to_string(yr))//"-0"//TRIM(int_to_string(mn))//"-0"//TRIM(int_to_string(d))//"_00:00:00_PSFC.nc"
  else 
    filename_ext_atm = "wrfout_d01_"//TRIM(int_to_string(yr))//"-0"//TRIM(int_to_string(mn))//"-"//TRIM(int_to_string(d))//"_00:00:00.nc"
    filename_ext_RAIN = "wrfhrly_d01_"//TRIM(int_to_string(yr))//"-0"//TRIM(int_to_string(mn))//"-"//TRIM(int_to_string(d))//"_00:00:00_RAIN.nc"
    filename_ext_LH = "wrfhrly_d01_"//TRIM(int_to_string(yr))//"-0"//TRIM(int_to_string(mn))//"-"//TRIM(int_to_string(d))//"_00:00:00_LH.nc"
    filename_ext_P = "wrfhrly_d01_"//TRIM(int_to_string(yr))//"-0"//TRIM(int_to_string(mn))//"-"//TRIM(int_to_string(d))//"_00:00:00_PSFC.nc"
  end if
else
  if (d<10) then
    filename_ext_atm = "wrfout_d01_"//TRIM(int_to_string(yr))//"-"//TRIM(int_to_string(mn))//"-0"//TRIM(int_to_string(d))//"_00:00:00.nc"
    filename_ext_RAIN = "wrfhrly_d01_"//TRIM(int_to_string(yr))//"-"//TRIM(int_to_string(mn))//"-0"//TRIM(int_to_string(d))//"_00:00:00_RAIN.nc"
    filename_ext_LH = "wrfhrly_d01_"//TRIM(int_to_string(yr))//"-"//TRIM(int_to_string(mn))//"-0"//TRIM(int_to_string(d))//"_00:00:00_LH.nc"
    filename_ext_P = "wrfhrly_d01_"//TRIM(int_to_string(yr))//"-"//TRIM(int_to_string(mn))//"-0"//TRIM(int_to_string(d))//"_00:00:00_PSFC.nc"
  else
    filename_ext_atm = "wrfout_d01_"//TRIM(int_to_string(yr))//"-"//TRIM(int_to_string(mn))//"-"//TRIM(int_to_string(d))//"_00:00:00.nc"
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


SUBROUTINE get_data(precip,evap,u,v,w,t,q,qc,qt,pp,pb,pbl_hgt,psfc)
!-----------------------------------------------
! read in the data for the first time
!-----------------------------------------------

USE global_data
USE netcdf

IMPLICIT NONE

REAL, DIMENSION(:,:,:) :: precip,evap,pbl_hgt,psfc
REAL, DIMENSION(:,:,:,:) :: u,v,w,t,q,qc,qt,pp,pb
REAL, DIMENSION(SIZE(u,1),SIZE(u,2),SIZE(u,3),datadaysteps) :: temp

CHARACTER(LEN=100) :: filename_ext_atm,filename_ext_RAIN,filename_ext_LH,filename_ext_P

INTEGER :: ncid,prencid,lhncid,psfcncid !,uncid,vncid,wncid,tncid,qncid,ppncid,pblncid
INTEGER :: preid,lhid,uid,vid,wid,tid,qid,ppid,pbid,pblid,psfcid
INTEGER :: sind,status,sind2,i,getsteps,getsteps2
!REAL :: filedays

REAL :: dayend

INTEGER :: jd_today,jd_before,new_y,new_m,new_d

integer, dimension(nf90_max_var_dims) :: rhDimIds
integer :: numDims

integer :: dim1size, dim2size, dim3size, dim4size
real,dimension(4) :: a

call get_filename(day,mon,year,filename_ext_atm,filename_ext_RAIN,filename_ext_LH,filename_ext_P)


call open_netcdf_files(ncid,prencid,lhncid,psfcncid,preid,lhid,uid,vid,wid,tid,qid,ppid,pbid,pblid,psfcid,filename_ext_atm,filename_ext_RAIN,filename_ext_LH,filename_ext_P)

! CHECK SHAPE OF NETCDF VARIABLE AS SEEN BY FORTRAN

! Precipitation
status = nf90_inquire_variable(prencid, preid, ndims = numDims)
if(status /= nf90_NoErr) call handle_err(status)
status = nf90_inquire_variable(prencid, preid, dimids = rhDimIds(:numDims))
if(status /= nf90_NoErr) call handle_err(status)
status = nf90_inquire_dimension(prencid, rhDimIds(1), len = dim1size) 
if(status /= nf90_NoErr) call handle_err(status) 
status = nf90_inquire_dimension(prencid, rhDimIds(2), len = dim2size) 
if(status /= nf90_NoErr) call handle_err(status) 
status = nf90_inquire_dimension(prencid, rhDimIds(3), len = dim3size) 
if(status /= nf90_NoErr) call handle_err(status)
print *,'PRECIP NETCDF SHAPE [1,2,3] =',dim1size,dim2size,dim3size

! Evap
status = nf90_inquire_variable(lhncid, lhid, ndims = numDims)
if(status /= nf90_NoErr) call handle_err(status)
status = nf90_inquire_variable(lhncid, lhid, dimids = rhDimIds(:numDims))
if(status /= nf90_NoErr) call handle_err(status)
status = nf90_inquire_dimension(lhncid, rhDimIds(1), len = dim1size) 
if(status /= nf90_NoErr) call handle_err(status) 
status = nf90_inquire_dimension(lhncid, rhDimIds(2), len = dim2size) 
if(status /= nf90_NoErr) call handle_err(status) 
status = nf90_inquire_dimension(lhncid, rhDimIds(3), len = dim3size) 
if(status /= nf90_NoErr) call handle_err(status)
print *,'LH NETCDF SHAPE [1,2,3] =',dim1size,dim2size,dim3size

! U 
status = nf90_inquire_variable(ncid, uid, ndims = numDims)
if(status /= nf90_NoErr) call handle_err(status)
status = nf90_inquire_variable(ncid, uid, dimids = rhDimIds(:numDims))
if(status /= nf90_NoErr) call handle_err(status)
status = nf90_inquire_dimension(ncid, rhDimIds(1), len = dim1size) 
if(status /= nf90_NoErr) call handle_err(status) 
status = nf90_inquire_dimension(ncid, rhDimIds(2), len = dim2size) 
if(status /= nf90_NoErr) call handle_err(status) 
status = nf90_inquire_dimension(ncid, rhDimIds(3), len = dim3size) 
if(status /= nf90_NoErr) call handle_err(status)
status = nf90_inquire_dimension(ncid, rhDimIds(4), len = dim4size) 
if(status /= nf90_NoErr) call handle_err(status)
print *,'U NETCDF SHAPE [1,2,3,4] =',dim1size,dim2size,dim3size,dim4size

! T 
status = nf90_inquire_variable(ncid, tid, ndims = numDims)
if(status /= nf90_NoErr) call handle_err(status)
status = nf90_inquire_variable(ncid, tid, dimids = rhDimIds(:numDims))
if(status /= nf90_NoErr) call handle_err(status)
status = nf90_inquire_dimension(ncid, rhDimIds(1), len = dim1size) 
if(status /= nf90_NoErr) call handle_err(status) 
status = nf90_inquire_dimension(ncid, rhDimIds(2), len = dim2size) 
if(status /= nf90_NoErr) call handle_err(status) 
status = nf90_inquire_dimension(ncid, rhDimIds(3), len = dim3size) 
if(status /= nf90_NoErr) call handle_err(status)
status = nf90_inquire_dimension(ncid, rhDimIds(4), len = dim4size) 
if(status /= nf90_NoErr) call handle_err(status)
print *,'T NETCDF SHAPE [1,2,3,4] =',dim1size,dim2size,dim3size,dim4size



! Since our input files only consist of one day, open all timesteps (datadaysteps) in file (i.e. remove sind2, make it 1)

!!! SIM DAY SHOULD BE AT THE END OF THE ARRAY, DAY BEFORE JUST BEFORE THAT, ETC.
!!! LAST BACK-TRACKED DAY SHOULD BE AT THE START OF THE ARRAY.
!!! THE LAST TIME POSITION IN THE ARRAY SHOULD BE THE FIRST TIMESTEP OF SIM DAY + 1.

if (totbtadays>1) then
	! Open the first day input file
	status = nf90_get_var(prencid, preid, precip, &
	start=(/bdy,bdy,1/),count=(/dim_j,dim_i,datadaysteps/))
	if(status /= nf90_NoErr) call handle_err(status)
	
	status = nf90_get_var(lhncid, lhid, evap(:,:,(datatotsteps-datadaysteps):(datatotsteps-1)), &
	start=(/bdy,bdy,1/),count=(/dim_j,dim_i,datadaysteps/))
	if(status /= nf90_NoErr) call handle_err(status)

	status = nf90_get_var(ncid, qid, temp(:,:,:,:), &
	start=(/bdy,bdy,1,1/),count=(/dim_j,dim_i,dim_k,datadaysteps/))
	if(status /= nf90_NoErr) call handle_err(status)

	q(:,:,:,(datatotsteps-datadaysteps):(datatotsteps-1)) = temp(:,:,dim_k:1:-1,:)

	status = nf90_get_var(ncid, uid, temp(:,:,:,:), &
	start=(/bdy,bdy,1,1/),count=(/dim_j,dim_i,dim_k,datadaysteps/))
	if(status /= nf90_NoErr) call handle_err(status)
	
	u(:,:,:,(datatotsteps-datadaysteps):(datatotsteps-1)) = temp(:,:,dim_k:1:-1,:)

	status = nf90_get_var(ncid, vid, temp(:,:,:,:), &
	start=(/bdy,bdy,1,1/),count=(/dim_j,dim_i,dim_k,datadaysteps/))
	if(status /= nf90_NoErr) call handle_err(status)
	
	v(:,:,:,(datatotsteps-datadaysteps):(datatotsteps-1)) = temp(:,:,dim_k:1:-1,:)

	status = nf90_get_var(ncid, wid, temp(:,:,:,:), &
	start=(/bdy,bdy,1,1/),count=(/dim_j,dim_i,dim_k,datadaysteps/))
	if(status /= nf90_NoErr) call handle_err(status)
	
	w(:,:,:,(datatotsteps-datadaysteps):(datatotsteps-1)) = temp(:,:,dim_k:1:-1,:)

	status = nf90_get_var(ncid, tid, temp(:,:,:,:), &
	start=(/bdy,bdy,1,1/),count=(/dim_j,dim_i,dim_k,datadaysteps/))
	if(status /= nf90_NoErr) call handle_err(status)
	
	t(:,:,:,(datatotsteps-datadaysteps):(datatotsteps-1)) = temp(:,:,dim_k:1:-1,:)

	status = nf90_get_var(ncid, ppid, temp(:,:,:,:), &
	start=(/bdy,bdy,1,1/),count=(/dim_j,dim_i,dim_k,datadaysteps/))
	if(status /= nf90_NoErr) call handle_err(status)
	
	pp(:,:,:,(datatotsteps-datadaysteps):(datatotsteps-1)) = temp(:,:,dim_k:1:-1,:)
	
	status = nf90_get_var(ncid, pbid, temp(:,:,:,:), &
	start=(/bdy,bdy,1,1/),count=(/dim_j,dim_i,dim_k,datadaysteps/))
	if(status /= nf90_NoErr) call handle_err(status)
	
	pb(:,:,:,(datatotsteps-datadaysteps):(datatotsteps-1)) = temp(:,:,dim_k:1:-1,:)

	status = nf90_get_var(ncid, pblid, pbl_hgt(:,:,(datatotsteps-datadaysteps):), &
	start=(/bdy,bdy,1/),count=(/dim_j,dim_i,datadaysteps/))
	if(status /= nf90_NoErr) call handle_err(status)
	
	status = nf90_get_var(psfcncid, psfcid, psfc(:,:,(datatotsteps-datadaysteps):), &
	start=(/bdy,bdy,1/),count=(/dim_j,dim_i,datadaysteps/))
	if(status /= nf90_NoErr) call handle_err(status)

	! close the netcdf files
	status = nf90_close(ncid)
	status = nf90_close(prencid)
	status = nf90_close(lhncid)
	status = nf90_close(psfcncid)
	
	print *,'Input file of first day loaded successfully:',filename_ext_atm


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

		status = nf90_get_var(lhncid, lhid, evap(:,:,datatotsteps-(datadaysteps*(i+1)):(datatotsteps-(datadaysteps*i)-1)), &
		start=(/bdy,bdy,1/),count=(/dim_j,dim_i,datadaysteps/))
		if(status /= nf90_NoErr) call handle_err(status)
		
		status = nf90_get_var(ncid, qid, temp(:,:,:,:), &
		start=(/bdy,bdy,1,1/),count=(/dim_j,dim_i,dim_k,datadaysteps/))
		if(status /= nf90_NoErr) call handle_err(status)
		
		q(:,:,:,datatotsteps-(datadaysteps*(i+1)):(datatotsteps-(datadaysteps*i)-1)) = temp(:,:,dim_k:1:-1,:)

		status = nf90_get_var(ncid, uid, temp(:,:,:,:), &
		start=(/bdy,bdy,1,1/),count=(/dim_j,dim_i,dim_k,datadaysteps/))
		if(status /= nf90_NoErr) call handle_err(status)
		
		u(:,:,:,datatotsteps-(datadaysteps*(i+1)):(datatotsteps-(datadaysteps*i)-1)) = temp(:,:,dim_k:1:-1,:)

		status = nf90_get_var(ncid, vid, temp(:,:,:,:), &
		start=(/bdy,bdy,1,1/),count=(/dim_j,dim_i,dim_k,datadaysteps/))
		if(status /= nf90_NoErr) call handle_err(status)
		
		v(:,:,:,datatotsteps-(datadaysteps*(i+1)):(datatotsteps-(datadaysteps*i)-1)) = temp(:,:,dim_k:1:-1,:)

		status = nf90_get_var(ncid, wid, temp(:,:,:,:), &
		start=(/bdy,bdy,1,1/),count=(/dim_j,dim_i,dim_k,datadaysteps/))
		if(status /= nf90_NoErr) call handle_err(status)
		
		w(:,:,:,datatotsteps-(datadaysteps*(i+1)):(datatotsteps-(datadaysteps*i)-1)) = temp(:,:,dim_k:1:-1,:)

		status = nf90_get_var(ncid, tid, temp(:,:,:,:), &
		start=(/bdy,bdy,1,1/),count=(/dim_j,dim_i,dim_k,datadaysteps/))
		if(status /= nf90_NoErr) call handle_err(status)
		
		t(:,:,:,datatotsteps-(datadaysteps*(i+1)):(datatotsteps-(datadaysteps*i)-1)) = temp(:,:,dim_k:1:-1,:)

		status = nf90_get_var(ncid, ppid, temp(:,:,:,:), &
		start=(/bdy,bdy,1,1/),count=(/dim_j,dim_i,dim_k,datadaysteps/))
		if(status /= nf90_NoErr) call handle_err(status)
		
		pp(:,:,:,datatotsteps-(datadaysteps*(i+1)):(datatotsteps-(datadaysteps*i)-1)) = temp(:,:,dim_k:1:-1,:)
		
		status = nf90_get_var(ncid, pbid, temp(:,:,:,:), &
		start=(/bdy,bdy,1,1/),count=(/dim_j,dim_i,dim_k,datadaysteps/))
		if(status /= nf90_NoErr) call handle_err(status)
		
		pb(:,:,:,datatotsteps-(datadaysteps*(i+1)):(datatotsteps-(datadaysteps*i)-1)) = temp(:,:,dim_k:1:-1,:)

		status = nf90_get_var(ncid, pblid, pbl_hgt(:,:,datatotsteps-(datadaysteps*(i+1)):(datatotsteps-(datadaysteps*i)-1)), &
		start=(/bdy,bdy,1/),count=(/dim_j,dim_i,datadaysteps/))
		if(status /= nf90_NoErr) call handle_err(status)
		
		status = nf90_get_var(psfcncid, psfcid, psfc(:,:,datatotsteps-(datadaysteps*(i+1)):(datatotsteps-(datadaysteps*i)-1)), &
		start=(/bdy,bdy,1/),count=(/dim_j,dim_i,datadaysteps/))
		if(status /= nf90_NoErr) call handle_err(status)
		
		print *,'Input file of previous day loaded successfully:',i,filename_ext_atm
		!print *,'L964, :(datatotsteps-(datadaysteps*i)-1)=',(datatotsteps-(datadaysteps*i)-1)
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
	
	status = nf90_get_var(lhncid, lhid, evap(:,:,datatotsteps), &
	start=(/bdy,bdy,1/),count=(/dim_j,dim_i,1/))
	if(status /= nf90_NoErr) call handle_err(status)
	
	status = nf90_get_var(ncid, qid, temp(:,:,:,1), &
	start=(/bdy,bdy,1,1/),count=(/dim_j,dim_i,dim_k,1/))
	if(status /= nf90_NoErr) call handle_err(status)
	
	q(:,:,:,datatotsteps) = temp(:,:,dim_k:1:-1,1)

	status = nf90_get_var(ncid, uid, temp(:,:,:,1), &
	start=(/bdy,bdy,1,1/),count=(/dim_j,dim_i,dim_k,1/))
	if(status /= nf90_NoErr) call handle_err(status)
	
	u(:,:,:,datatotsteps) = temp(:,:,dim_k:1:-1,1)

	status = nf90_get_var(ncid, vid, temp(:,:,:,1), &
	start=(/bdy,bdy,1,1/),count=(/dim_j,dim_i,dim_k,1/))
	if(status /= nf90_NoErr) call handle_err(status)
	
	v(:,:,:,datatotsteps) = temp(:,:,dim_k:1:-1,1)

	status = nf90_get_var(ncid, wid, temp(:,:,:,1), &
	start=(/bdy,bdy,1,1/),count=(/dim_j,dim_i,dim_k,1/))
	if(status /= nf90_NoErr) call handle_err(status)
	
	w(:,:,:,datatotsteps) = temp(:,:,dim_k:1:-1,1)

	status = nf90_get_var(ncid, tid, temp(:,:,:,1), &
	start=(/bdy,bdy,1,1/),count=(/dim_j,dim_i,dim_k,1/))
	if(status /= nf90_NoErr) call handle_err(status)
	
	t(:,:,:,datatotsteps) = temp(:,:,dim_k:1:-1,1)

	status = nf90_get_var(ncid, ppid, temp(:,:,:,1), &
	start=(/bdy,bdy,1,1/),count=(/dim_j,dim_i,dim_k,1/))
	if(status /= nf90_NoErr) call handle_err(status)
	
	pp(:,:,:,datatotsteps) = temp(:,:,dim_k:1:-1,1)
	
	status = nf90_get_var(ncid, pbid, temp(:,:,:,1), &
	start=(/bdy,bdy,1,1/),count=(/dim_j,dim_i,dim_k,1/))
	if(status /= nf90_NoErr) call handle_err(status)
	
	pb(:,:,:,datatotsteps) = temp(:,:,dim_k:1:-1,1)

	status = nf90_get_var(ncid, pblid, pbl_hgt(:,:,datatotsteps), &
	start=(/bdy,bdy,1/),count=(/dim_j,dim_i,1/))
	if(status /= nf90_NoErr) call handle_err(status)
	
	status = nf90_get_var(psfcncid, psfcid, psfc(:,:,datatotsteps), &
	start=(/bdy,bdy,1/),count=(/dim_j,dim_i,1/))
	if(status /= nf90_NoErr) call handle_err(status)
	
	! close the netcdf files
	status = nf90_close(ncid)
	status = nf90_close(prencid)
	status = nf90_close(lhncid)  
	status = nf90_close(psfcncid)
	
	print *,'Input file of next day (1st time step) loaded successfully:',filename_ext_atm
	
else
	print*, 'If you only want to back-track for one day, must change how input data is retrieved.'



end if
  

!evap converted to mm > Unit conversion checked and OK 18/7/17 :)
evap = evap*(1440/datadaysteps)*60/Lv
  
qt = qt + q ! i.e. SUM(QCLD,QRAIN,QSNOW,QICE) + QVAPOUR
!qt = q
qc = qt

END SUBROUTINE get_data

SUBROUTINE calc_actual_temp(temp,pres,act_temp)
!------------------------------------------------
! calculate the actual temperature at every level and time,
! given the input perturbation potential temperature.
!-----------------------------------------------------

USE global_data

IMPLICIT NONE

REAL, INTENT(IN), DIMENSION(:,:,:,:) :: temp,pres

REAL, INTENT(OUT), DIMENSION(:,:,:,:) :: act_temp


!
! Calculate the actual temperature [K] from potential temperature.
! Need to add 300K since wrfout gives *pertubation* potential temperature as temp T.
!
act_temp = (temp + 300) * ((P0/pres)**(-Rd/Cp))

END SUBROUTINE calc_actual_temp

SUBROUTINE new_out_file(outncid,testid,wvcid,wvc2id,xlocid,ylocid,dayid,opreid,daynum,lat2d,lon2d)
!----------------------------------------------
!create the output netcdf file and prepare it to accept data
!--------------------------------------------------------

USE global_data
USE netcdf

IMPLICIT NONE

INTEGER, INTENT(INOUT) :: outncid,testid,wvcid,wvc2id,xlocid,ylocid,dayid,opreid
!REAL, INTENT(IN) :: daynum
INTEGER, INTENT(IN) :: daynum
REAL,INTENT(IN),DIMENSION(:,:) :: lat2d,lon2d

INTEGER :: status,jdimid,idimid,kdimid,tdimid,gwvcdimid,latid,lonid


!
!create the file
!
status = nf90_create(TRIM(diro)//"act_temp_check.nc",nf90_clobber,outncid)
if (status /= NF90_NOERR) call handle_err(status)


!define dimensions
!
status = nf90_def_dim(outncid,"j",dim_j,jdimid)
if (status /= NF90_NOERR) call handle_err(status)
status = nf90_def_dim(outncid,"i",dim_i,idimid)
if (status /= NF90_NOERR) call handle_err(status)
status = nf90_def_dim(outncid,"k",dim_k,kdimid)
if (status /= NF90_NOERR) call handle_err(status)
status = nf90_def_dim(outncid,"t",datatotsteps,tdimid)
if (status /= NF90_NOERR) call handle_err(status)

!
!define the variable
!
status = nf90_def_var(outncid,"act_temp",nf90_float,(/jdimid,idimid,kdimid,tdimid/),testid)
if (status /= NF90_NOERR) call handle_err(status)
! status = nf90_def_var(outncid,"x_loc",nf90_int,(/gwvcdimid/),xlocid)
! if (status /= NF90_NOERR) call handle_err(status)
! status = nf90_def_var(outncid,"y_loc",nf90_int,(/gwvcdimid/),ylocid)
! if (status /= NF90_NOERR) call handle_err(status)
! status = nf90_def_var(outncid,"day",nf90_float,(/gwvcdimid/),dayid)
! if (status /= NF90_NOERR) call handle_err(status)
! status = nf90_def_var(outncid,"pre",nf90_float,(/gwvcdimid/),opreid)
! if (status /= NF90_NOERR) call handle_err(status)
status = nf90_def_var(outncid,"latitcrs",nf90_float,(/jdimid,idimid/),latid)
if (status /= NF90_NOERR) call handle_err(status)
status = nf90_def_var(outncid,"longicrs",nf90_float,(/jdimid,idimid/),lonid)
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

END MODULE qibt_subs






PROGRAM test
USE netcdf
USE global_data
USE qibt_subs

IMPLICIT NONE

!
!netcdf id variables
!
INTEGER :: status,headncid
INTEGER :: spid,ptopid,delxid,latcrsid,loncrsid,terid,tstepid,sdayid,smonid,syearid
INTEGER :: dimjid,dimiid,sigid,dimkid
INTEGER :: outncid,testid,wvcid,wvc2id,xlocid,ylocid,dayid,opreid,latid,lonid
INTEGER :: wsncid, wsid
CHARACTER(LEN=100):: fname
CHARACTER(LEN=10):: fname_smon, fname_sday
INTEGER :: datatstep

! REAL,ALLOCATABLE,DIMENSION(:,:,:) :: precip, evap
! REAL,ALLOCATABLE,DIMENSION(:,:,:,:) :: u,temp,pp,pb,pres,act_temp
REAL,ALLOCATABLE,DIMENSION(:,:,:) :: precip
REAL,ALLOCATABLE,DIMENSION(:,:,:) :: evap,tpw,pbl_hgt,surf_pres,pstar,psfc
REAL,ALLOCATABLE,DIMENSION(:,:,:,:) :: u,v,w,temp,act_temp,mix,pp,pb,pw,mixcld,mixtot,pres
INTEGER,ALLOCATABLE,DIMENSION(:,:,:) :: pbl_lev
REAL,ALLOCATABLE,DIMENSION(:,:) :: lat2d,lon2d




!----------------------------------------------------------------
!
! Retrieve simulation start and end dates, and output directory, from the command line input
!
INTEGER :: num_args
character(len=100), dimension(:), allocatable :: args
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
call get_command_argument(7,args(7))
diro = args(7)

print *,diro

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

fname="wrfout_d01_"//TRIM(int_to_string(syear))//"-"//TRIM(fname_smon)//"-"//TRIM(fname_sday)//"_00:00:00.nc"
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
! status = nf90_get_var(headncid, ptopid, ptop)
! if(status /= nf90_NoErr) call handle_err(status)
! status = nf90_get_att(headncid, NF90_GLOBAL, "DX",delx)
! if(status /= nf90_NoErr) call handle_err(status)
status = nf90_inquire_dimension(headncid, tstepid,len = datatstep) 
datatstep=1440/datatstep ! Value must be 1440/8=180, where 8 is number of timesteps in the file (it was set up this way based on MM5 input files, where MM5 model timestep was 180mins)
print *,'L1124,datatstep(mins)=',datatstep
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

! ! Allocate the required arrays
ALLOCATE( lon2d(dim_j,dim_i),lat2d(dim_j,dim_i), STAT = status) !sigma(dim_k),pstar(dim_j,dim_i,datadaysteps), &
 
!
! Read in more variables
!
status = nf90_get_var(headncid, latcrsid, lat2d,start=(/bdy,bdy/),count=(/dim_j,dim_i/))
if(status /= nf90_NoErr) call handle_err(status)
status = nf90_get_var(headncid, loncrsid, lon2d,start=(/bdy,bdy/),count=(/dim_j,dim_i/))
if(status /= nf90_NoErr) call handle_err(status)


status = nf90_close(headncid)

!
! Calculate the number of trajectory time steps in a day and in input file time step
!
daytsteps = 1440/tstep                ! number of sub-daily simulation time steps
indatatsteps = datatstep/tstep          ! divide input file time step by the number of simulation time steps, as they may differ
totsteps = daytsteps*(totbtadays+1)   ! total number of simulation data time steps to remember

datadaysteps = 1440/datatstep           ! number of input file time steps in day
datatotsteps = (datadaysteps*(totbtadays+1)) + 1 ! total number of input file time steps over the back-track period

! Find total number of days to run simulation for, given input start and end dates
totdays=simlength(sday,smon,syear,edday,edmon,edyear)

print *,"Total days to run analysis=",totdays
print *,"Total number of back-track days=",totbtadays
print *,"Number of parcels=",nparcels
print *,"Simulation time step (mins)=",tstep

! Allocate the variable arrays
ALLOCATE( precip(dim_j,dim_i,datadaysteps), &
  evap(dim_j,dim_i,datatotsteps),u(dim_j,dim_i,dim_k,datatotsteps), &
  v(dim_j,dim_i,dim_k,datatotsteps),temp(dim_j,dim_i,dim_k,datatotsteps), &
  mix(dim_j,dim_i,dim_k,datatotsteps),pp(dim_j,dim_i,dim_k,datatotsteps),pb(dim_j,dim_i,dim_k,datatotsteps), &
  act_temp(dim_j,dim_i,dim_k,datatotsteps), &
  mixtot(dim_j,dim_i,dim_k,datatotsteps), &
  tpw(dim_j,dim_i,datatotsteps),pw(dim_j,dim_i,dim_k,daytsteps+1), &
  mixcld(dim_j,dim_i,dim_k,datatotsteps),pres(dim_j,dim_i,dim_k,datatotsteps), &
  psfc(dim_j,dim_i,datatotsteps),surf_pres(dim_j,dim_i,datatotsteps), w(dim_j,dim_i,dim_k,datatotsteps), &
  pbl_hgt(dim_j,dim_i,datatotsteps),pbl_lev(dim_j,dim_i,datatotsteps), STAT = status )

  

dd=1

call day_month_year(dd) 
print *,"day,mon,year",day,mon,year  

call get_data(precip,evap,u,v,w,temp,mix,mixcld,mixtot,pp,pb,pbl_hgt,psfc)
print *,'got data'

pres = pp + pb

call calc_actual_temp(temp,pres,act_temp)
print *,'actual temperature calculated'

call new_out_file(outncid,testid,wvcid,wvc2id,xlocid,ylocid,dayid,opreid,day,lat2d,lon2d)
print *,'created new out file'

!$OMP CRITICAL (output)
status = nf90_put_var(outncid,testid,act_temp,start=(/1,1,1,1/), count=(/dim_j,dim_i,dim_k,datatotsteps/))
if(status /= nf90_NoErr) call handle_err(status)

!$OMP END CRITICAL (output)
status = nf90_close(outncid)
if(status /= nf90_NoErr) call handle_err(status)


print *,'T[1,1,29,17]=',temp(1,1,29,17),'K'
print *,'pres[1,1,29,17]=',pres(1,1,29,17),'Pa'
print *,'manual_act_temp=',(temp(1,1,29,17)+300)*((P0/pres(1,1,29,17))**(-Rd/Cp)),'K'
print *,'manual_act_temp=',((temp(1,1,29,17)+300)*((P0/pres(1,1,29,17))**(-Rd/Cp)))-273.15,'degC'
print *,'act_temp[1,1,29,17]=',act_temp(1,1,29,17),' K'
print *,'act_temp[1,1,29,17]=',act_temp(1,1,29,17)-273.15,'degC'

END PROGRAM test
