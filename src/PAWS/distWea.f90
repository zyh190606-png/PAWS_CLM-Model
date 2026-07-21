#include <preproc.h> 
!

      subroutine distWea(H,JD)
      !% (E) This solver normally works on hourly time step on a station level, it does:
      !% (1) distribute prcp using reference precip unit hyetograph
      !% (2) calculate net radiation, reference (FAO grass) ET
      !% (3) store some ET coefficients
      !% data structure:
      !% a struct called 'day' include data distributed for each day for each-
      !% station. day.prcp(1:24,i) is for the i-th station.
      !% Input:
      !% wid: watershed ID, H: hour of day, JD: julian day (of the year)
      !% Output:
      !% Updated: g.Wea.day, g.OVN.prcp, g.Wea.t
      use Vdata
      use polyEval, Only: linearInterp1d
      implicit none
      integer :: H, JD, k, datevec(3), yy, mm, dd, i, j
      real*8 :: gid, PREC, frac, RFRAC, loc=0D0, lastLoc=0D0, t1, p1
      real*8, dimension(:), pointer :: prcp=>NULL()
      ! two lines below added 12/13/13 for running constant forcing hypothetic case.
      ! normally, implementConst = .false. . If turned true, it will detect the YYMMDD specified in constYY, etc
      ! climate forcing everyday after constYYconstMMconstDD will be set to the same as this day.
      logical*1 ::  implementConst = .false., constantDay=.false., constantDaySWAP=.true.,constantDaySWAPini=.false.
      integer:: constYY2 = 2006, constMM2 = 4, constDD2 = 15, JD2
      integer:: constYY = 2006, constMM = 4, constDD = 15

      !global g w
      !% It can be: daily data distributed to arbitrary time steps
      !% Only at beginning of the day will re-distribution happen
      !% If this day prcp=0, then that station won't do anything throughout the
      !% day. This is achieved by making DData=[];
      !% Sta = w(wid).Stations;
      gid = w%g
      PREC = g%Wea%dt/10.0D0
      !%k = round((g(gid).Wea.t-floor(g(gid).Wea.t+PREC))/g(gid).Wea.dt)+1;
      k = H
      frac = 1.0D0
      !RFRAC = 1.08D0    !% solar radiation from SPOKAS is reduced!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      RFRAC = 1.0D0   ! after reducing Sp0 to 1320, now this is put back to 1.0
      ![yy,mm,dd]=datevec(g(gid).Wea.t);
      if (k .ne. 0) then
        datevec = datevec2(int(g%Wea%t))
        yy = datevec(1)
        mm = datevec(2)
        dd = datevec(3)
      endif
#ifndef MEX     
      if (abs(w%t-floor(w%t))<1D-10) then
        !% first step of the day, distribute the data
        !g(gid).OVN.prcp = zeros(size(g(gid).OVN.h));
        if (implementConst .AND. (.NOT. constantDay)) then
          if (yy >= constYY .and. mm >= constMM .and. dd >= constDD) then
            constantDay = .true.
          endif
        endif
        if (constantDay .and. constantDaySWAP .and. (.NOT. constantDaySWAPini)) then
            call ymd2jd((/constYY2, constMM2, constDD2/), 3, JD2)
            call makeDay(RFRAC,frac,constYY2,constMM2,JD2,g%Wea%day)
            constantDaySWAPini = .true.
        endif
        if (.not. constantDay) then
        g%OVN%prcp = 0.0D0
        !g(gid).Wea.day = makeDay(wid,RFRAC,frac,yy,mm,JD,g(gid).Wea.day);
        call makeDay(RFRAC,frac,yy,mm,JD,g%Wea%day)
        endif
      endif
#endif
      !% If this day has prcp = 0,(DData=[]) then don't do anything throughout this day
      do i = 1, size(w%Stations)
        if (associated(w%Stations(i)%gIndices)) then
            prcp => oneD_Ptr(size(g%OVN%prcp),g%OVN%prcp)
            if (k .eq. 0) then
                ! a signal that it is not called by wtrshd_day. Rather, a continuous time series
                ! just do a linear interpolation of values.
                t1 = w%t + w%dt/2D0
                call linearInterp1d(w%Stations(i)%datenums,w%Stations(i)%prcp,t1,p1,loc,lastLoc); 
                lastLoc = loc;
                prcp(w%Stations(i)%gIndices) = p1
            else    
                if (g%Wea%day%prcp(k,i)>0) then
                  prcp(w%Stations(i)%gIndices) = g%Wea%day%prcp(k,i)  
                else
                  !% no need to set to 0 as it is zero by itself, see
                  !% typeIIrainfall.mat
                endif    
            endif
        endif
      enddo

      g%Wea%t = g%Wea%t + g%Wea%dt


      end subroutine distWea

