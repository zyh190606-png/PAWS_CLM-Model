MODULE clm_time_manager
! modified by CP
! added some common storages, which are to be evaluated by PAWS
real*8 dt_rad
real*8 dt_veg
real*8 dt_vdz
integer year_PAWS, month_PAWS, day_PAWS, secs_PAWS, mcsec_Day
real*8 JD_PAWS
public:: is_end_curr_day
integer:: step_no = 0
interface get_curr_calday
MODULE PROCEDURE get_curr_calday
end interface

CONTAINS
REAL*8 function get_step_size()
get_step_size = dt_rad
!get_step_size = 1D0/24D0 * 86400D0
!get_step_size = VE(1)%dt * 1D0

END function

logical function is_end_curr_day()
use Vdata
if (abs(w%t + w%dt - ceiling(w%t - 1D-10))<w%dt/1000D0) then
  is_end_curr_day = .true.
else 
  is_end_curr_day = .false.
endif
  
end function

function get_calday(ymd, tod)
use Vdata

! Return calendar day corresponding to specified time instant.
! Calendar day 1.0 = 0Z on Jan 1.

! Arguments
   integer, intent(in) :: &
      ymd,   &! date in yearmmdd format
      tod     ! time of day (seconds past 0Z)

! Return value
   real*8 :: get_calday

! Local variables
   character(len=*), parameter :: sub = 'clm::get_calday'
   integer :: rc                 ! return code
   integer d, YJD, darray(3)
   
   call ymd2jd( (/ymd/), 3, d, YJD, darray)
   get_calday = d + tod
end function get_calday

subroutine get_prev_date(yr, mon, day, tod)
use Vdata
! Return date components valid at beginning of current timestep.

! Arguments
   integer, intent(out) ::&
      yr,    &! year
      mon,   &! month
      day,   &! day of month
      tod     ! time of day (seconds past 0Z)
      integer d3(3)
      
      d3 = datevec2(int(w%t))
      yr = d3(1)
      mon= d3(2)
      day= d3(3)
      tod= int((w%t - int(w%t))*86400D0+1D-5)

end subroutine get_prev_date

integer function get_nstep()
!get_nstep = VE(1)%dt * 1D0
USE Vdata
get_nstep = int((w%t - w%ModelStartTime)/w%dt + 1D-8) + 1
END function get_nstep

INTEGER FUNCTION get_days_per_year()
integer y
y = year_PAWS
if ((mod(y,4) == 0 .AND. mod(y,100) .ne. 0) .OR. mod(y,400) == 0) then
  get_days_per_year = 366
else
  get_days_per_year = 366
endif
END FUNCTION

REAL*8 function get_rad_step_size()
!get_rad_step_size = VE(1)%dt * 1D0
get_rad_step_size = dt_rad
END function

subroutine get_curr_date(year, month, day, secs)
integer year, month, day, secs
year = year_PAWS
month= month_PAWS
day  = day_PAWS
secs = secs_PAWS
end subroutine

function get_curr_calday(offset)
integer,OPTIONAL:: offset ! when using offset it should be seconds. to be consistent with CLM logic
real*8 get_curr_calday
if (PRESENT(offset)) then
  get_curr_calday   = JD_PAWS + DBLE(offset)/86400D0
else
  get_curr_calday   = JD_PAWS
endif

end function

logical function is_perpetual()
is_perpetual = .false.
end function

subroutine get_ref_date(yr, mon, day, tod)
integer ::&
      yr,    &! year
      mon,   &! month
      day,   &! day of month
      tod     ! time of day (seconds past 0Z)
end subroutine


logical function newyr()
! check if this is the first step of a new year
USE VData
integer d(3)
integer :: prec = 1D-8
d = datevec2(int(w%t+PREC))
newyr = .false.
if (abs(w%t - int(w%t))<PREC) then
  if ((d(2) .eq. 1) .and. (d(3) .eq. 1)) then
     newyr = .true.
  endif
endif

end function

integer function nyrSimulated()
! check if this is the first step of a new year
USE paws_varctl, Only: nyrs
nyrSimulated = nyrs
end function

END MODULE
