
      subroutine makeDay(RFRAC,frac,year,mm,JD,day)
      !% (E) gather data from weather stations to create an object, 'day'
      !% Input:
      !% RFRAC: rescaling factor for radiation as the program's output is larger-
      !% than actual radiation; frac: a rescaling factor for prcp used during-
      !% testing, it should be 1 now; year: YYYY format, mm: month in MM format;-
      !% JD: julian day of the year, day: a structure template for the outgoing-
      !% object
      !% Output:
      !% day, a structure which contains weather data for the day in arrays that-
      !% are stored as fields
      use Vdata
      implicit none
      type(STA_type) :: day
      integer :: year, mm, JD, datevec(3), y, m, d, dn, i, ii, ptr, ptr1, kkk, nr, nc, nk, bk, nbb
      real*8 :: RFRAC, frac, gid, p, tmin, tmax, Tavg, dew, sat
      real*8, dimension(24):: temp, e0, hmd, hk, rad, rn, ret, rets, hb, hli, wnd
      real*8, dimension(:,:), allocatable :: nb
      integer, dimension(size(w%Stations)) :: PTRs
      integer, save:: status = 0, nss = 0

      if (.not. associated(g%Wea%month%tempM)) then
       nss = size(w%Stations)
        allocate(g%Wea%month%tempM(nss))
        allocate(g%Wea%month%tmxM(nss))
        allocate(g%Wea%month%tmnM(nss))
        allocate(g%Wea%month%radM(nss))
        allocate(g%Wea%month%wndM(nss))
        allocate(g%Wea%month%hmdM(nss))
        allocate(g%Wea%month%dewM(nss))
        status = 1
      endif



      !global g w
      gid = w%g
      !temp= zeros(24,length(w(wid).Stations));
      temp = 0.0D0
      ![y,m,d]=datevec(g(gid).Wea.t);
      datevec = datevec2(int(g%Wea%t))
      y = datevec(1)
      m = datevec(2)
      d = datevec(3)
      if (y/=year) then
        call scheduleSolar(year)
      endif

      dn = datenum2((/year, JD/))
      !PTRs = zeros(size(w(wid).Stations));
      PTRs = 0
      do i = 1, size(w%Stations)
        PTRs(i) = dn-w%Stations(i)%datenums(1)+1
      enddo

      do i = 1, size(w%Stations)
        !% precip
        !%ptr = w(wid).Stations(i).ptr;
        ptr = PTRs(i)       !% ASSUMING DAILY DATA
        p = w%Stations(i)%prcp(ptr)
        if (p>0.0D0 .and. p<100.0D0) then
            day%prcp(:,i) = w%Stations(i)%ref * p/(1000.0D0*frac)   !% Currently, divide by frac as there is no subsurface component
        elseif (p==0.0D0) then
            day%prcp(:,i) = 0.0D0
        else    !% Missing Data. Weather Generator
            !nb = w(wid).Stations(i).neighbor; kkk = 1; [nr,nc]=size(nb);
            if(allocated(nb)) deallocate(nb)
            allocate(nb(size(w%Stations(i)%neighbor,1),size(w%Stations(i)%neighbor,2)))
            nb = w%Stations(i)%neighbor
            kkk = 1_8
            nr = size(nb,1)
            nc = size(nb,2)
            do while (kkk<=nr)
                nbb = nb(kkk,1)
                !% p = w(wid).Stations(nbb).prcp(w(wid).Stations(nbb).ptr); a
                !% bug for new code!!!
                p = w%Stations(nbb)%prcp(PTRs(nbb))
                if (p<0) then
                  kkk = kkk + 1_8
                else
                  exit
                endif
            enddo
            if (kkk>nr) then
                day%prcp(:,i) = 0.0D0
            else
                day%prcp(:,i) = w%Stations(nbb)%ref * p/(1000.0D0*frac)
            endif
        endif
        !% Temperature
        Tmin = w%Stations(i)%tmin(ptr)
        Tmax = w%Stations(i)%tmax(ptr)
        if (Tmin<-90) then
            !% try to borrow
            ![Tmin,nk,bk] = borrowStations(w(wid).Stations,i,'tmin',ptr,ptr);
            call borrowStations(w%Stations,size(w%Stations),i,'tmin',ptr,ptr,Tmin,nk,bk)
            if (Tmin<-90) then
              !Tmin = g%Wea%month%tmn(mm)
              Tmin = g%Wea%month%tmnM(i)
            endif
        endif
        if (Tmax<-90) then
            !% try to borrow
            ![Tmax,nk,bk] = borrowStations(w(wid).Stations,i,'tmax',ptr,ptr);
            call borrowStations(w%Stations,size(w%Stations),i,'tmax',ptr,ptr,Tmax,nk,bk)
            if (Tmax<-90) then
                !Tmax = g%Wea%month%tmx(mm)
                Tmax = g%Wea%month%tmxM(i)
            endif
        endif
        Tavg = (Tmin + Tmax) / 2.0D0
        hk = (/(ii, ii = 1, 24)/)
        temp = Tavg + (Tmax-Tmin)*dcos(0.2618D0*(hk-15.0D0))/2.0D0      !% SWAT page 29
        !% Incoming Radiation is generated for an year at start of an year.
        ptr1 = (JD - 1) * 24 + 1
        !%ptr1 = (dd-1)*24+1;
        !try
        !    rad = w(wid).Stations(i).Rad(ptr1:ptr1+23)/RFRAC;
        !catch
        !    4;
        !end
        rad = w%Stations(i)%Rad(ptr1:ptr1+23)/RFRAC
        !% Wnd
        !%wnd = zeros(24,1) + g(gid).Wea.month.wnd;
        !wnd = zeros(24,1) + w.Stations(i).awnd(ptr);
        wnd = 0.0D0 + w%Stations(i)%awnd(ptr)
        where (wnd<-90)
            wnd = g%Wea%month%wndM(i) * 0.89D0
        endwhere
        !% TEMPORARY, ONLY USING STATION 5!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !% Hmd
        !%dew = g(gid).Wea.month.dew; %monthly averages
        !%dew = w.Stations(7).dptp(ptr); if dew<-90, dew = g(gid).Wea.month.dew; end
        dew = w%Stations(i)%dptp(ptr)
        if (dew<-90) then 
            dew = g%Wea%month%dewM(i)
        endif
        e0 = (dexp((16.78D0*temp-116.9D0)/(temp+237.3D0)))
        hmd = dexp((16.78D0*dew-116.9D0)/(dew+237.3D0))/e0
        sat = 0.98D0
        !hmd(hmd>sat)=sat;
        where(hmd>sat)
          hmd = sat
        endwhere
        !% Reference ET
        ![ret,rets,rn,hb,hli]=PenmanMonteith(rad,hmd,e0,wnd,temp,w(wid).Stations(i).XYElev(3),w(wid).Stations(i).tau(JD),0);
        call PenmanMonteith(24,rad,hmd,e0,wnd,temp,w%Stations(i)%XYElev(3),w%Stations(i)%tau(JD),0,ret,rets,rn,hb,hli)
        !% SNOW NOT INCLUDED NOW!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !% unit
        day%temp(:,i) = temp
        day%rad(:,i) = rad
        day%rn(:,i) = rn
        day%wnd(:,i) = wnd
        !day.dptp(:,i)=zeros(24,1)+dew;
        day%dptp(:,i) = 0.0D0 + dew
        day%hmd(:,i) = hmd
        day%ret(:,i) = ret/1000.0D0     !% Unit conversion from mm to m
        day%rets(:,i) = rets/1000.0D0
        day%hb(:,i) = hb
        day%hli(:,i) = hli
        day%cosz(:,i) = w%Stations(i)%cosz(ptr1:ptr1+23)
        day%radd(:,i) = w%Stations(i)%RD(ptr1:ptr1+23)  ! these are in W/m2
        day%radi(:,i) = w%Stations(i)%RI(ptr1:ptr1+23)  ! these are in W/m2
        day%decl(:,i) = w%Stations(i)%decl(ptr1:ptr1+23)
      enddo

      !% for i=1:length(w(wid).Stations)
      !%     w(wid).Stations(i).ptr = w(wid).Stations(i).ptr + 1;
      !% end
      !% NO LONGER KEEPS UPDATING EACH POINTERS!!!

      if(allocated(nb)) deallocate(nb)
      
      end subroutine makeDay
