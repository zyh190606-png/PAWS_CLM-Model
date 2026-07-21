      module wMonthlyAverages_Mod
      
      integer,parameter:: ns = 1100
      integer, public :: nSta = 0
      real*8, public :: POS(ns,2),Welev(ns),TMPMX(ns,12),TMPMN(ns,12),TMPSTDMX(ns,12),TMPSTDMN(ns,12)
      REAL*8, public :: PCPMM(ns,12),PCPSTD(ns,12),PCPSKW(ns,12),PR_W1(ns,12),PR_W2(ns,12),PCPD(ns,12)
      REAL*8, public :: RAINHHMX(ns,12),SOLARAV(ns,12),DEWPT(ns,12),WNDAV(ns,12)
      integer:: id(ns), rain_yrs
      public:: readwgn, wMonthlyAverages

      Contains

      subroutine wMonthlyAverages(M, month)
      !% extract monthly average values for month
      !% These are variables obtained from BASINS database (4112500,basin 3)
      !% from .wgn file
      use Vdata
      implicit none
      integer M, Site
      type(month_type) month
      real*8, dimension(size(w%Stations),12) :: TMX, TMN, TAVG, DEW, RH, RAD, WND, PCP
      integer, save:: status = 0, is = 0, cis = 0, nss
      integer :: cisM(200)=0
      integer I,J
      real*8 lat,long, dd,mxpct
      real*8 :: dis, disM(200)=0

      if (status .eq. 0)  then 
       call readwgn()
       status = 1
      endif

      if (nSta .eq. 0) then
       Site = 1
      else
       Site = 0 ! use wgn dataset
      endif
      nss = size(w%Stations)
      if (Site == 0) then       !% WGN
        ! to be made better later
        ! First identify a WGN station
        mxpct = 0; 

        if (is .eq. 0) then
        DO I=1,size(w%Stations)
          if (w%Stations(I)%pct(1) > mxpct) then
            is = I
            mxpct = w%Stations(I)%pct(1)
          endif
        ENDDO
        
        DO J=1,size(w%Stations)
        lat = w%Stations(J)%LatLong(1)
        long = w%Stations(J)%LatLong(2)
        dis = 1e4
        DO I=1,nSta
          dd = sqrt((pos(i,1)-lat)**2D0+(pos(i,2)-long)**2D0)
          if (dd<dis) then
            dis = dd
            disM(J) = dis
            cisM(J) = I
          endif
        ENDDO
        if (dis > 1) call display('wMonthlyAverages::Stations too far away from WGN stations')

        TMX(J,:) = TMPMX(cisM(j),:)
        TMN(J,:) = TMPMN(cisM(j),:)
        TAVG(J,:) = (TMX(J,:)+TMN(J,:))/2.0D0
        DEW(J,:) = DEWPT(cisM(J),:)      !% dew point temperature
        RH(J,:) = dexp((16.78D0*DEW(J,:)-116.9D0)/(DEW(J,:)+237.3D0))/(dexp((16.78D0*TAVG(J,:)-116.9D0)/(TAVG(J,:)+237.3D0)))
        RAD(J,:) = SOLARAV(cisM(j),:) 
        PCP(J,:) = PCPMM(cisM(j),:)
        WND(J,:) = 5D0
        ENDDO
        endif ! is_IF
        
      elseif (Site == 1) then       !% RCR
        is = 1
        DO J=1,size(w%Stations)
        TMX(J,:) = [-1.58D0, 0.14D0, 5.70D0, 13.79D0, 20.64D0, 26.01D0, 28.29D0, 27.27D0, 23.00D0, 16.52D0, 8.17D0, 1.07D0]
        TMN(J,:) = [-9.53D0, -8.60D0, -3.78D0, 2.38D0, 8.14D0, 13.61D0, 16.07D0, 15.04D0, 10.91D0,  5.05D0, -0.33D0, -6.37D0]
        TAVG(J,:) = (TMX(J,:)+TMN(J,:))/2.0D0
        DEW(J,:) = [-7.18D0, -6.92D0, -3.85D0, 1.57D0, 7.67D0, 13.36D0, 15.31D0, 15.31D0, 11.42D0, 5.57D0, -0.27D0, -5.27D0]      !% dew point temperature
        RH(J,:) = dexp((16.78D0*DEW(J,:)-116.9D0)/(DEW(J,:)+237.3D0))/(dexp((16.78D0*TAVG(J,:)-116.9D0)/(TAVG(J,:)+237.3D0)))
        RAD(J,:) = [ 5.03D0, 8.67D0, 12.78D0, 15.29D0, 19.99D0, 22.79D0, 22.79D0, 19.69D0, 15.67D0, 10.81D0,  5.70D0,  4.44D0]
        WND(J,:) = [ 5.29D0, 5.40D0,  5.77D0,  5.52D0,  4.33D0,  3.55D0,  3.16D0,  2.94D0,  3.76D0,  4.10D0,  4.92D0,  5.21D0]
        PCP(J,:) = [ 40.20D0, 38.20D0, 51.90D0, 74.20D0, 71.80D0, 87.40D0, 80.50D0, 68.60D0, 59.90D0, 50.90D0, 52.60D0, 47.80D0]
        ENDDO
      elseif (Site == 2) then
        is = 1
        DO J=1,size(w%Stations)
        TMX(J,:) = [-1.82D0, -0.29D0, 4.88D0, 13.03D0, 19.97D0, 25.43D0, 27.81D0, 26.80D0, 22.71D0, 16.28D0, 8.00D0, 0.86D0]
        TMN(J,:) = [-10.02D0, -9.54D0, -4.92D0, 1.10D0, 6.55D0, 12.09D0, 14.64D0, 13.92D0, 10.13D0, 4.59D0, -0.63D0, -6.57D0]
        TAVG(J,:) = (TMX(J,:)+TMN(J,:))/2.0D0
        DEW(J,:) = [-8.41D0, -8.22D0, -4.72D0, 1.17D0, 6.88D0, 12.42D0, 14.78D0, 14.78D0, 11.06D0, 5.78D0, -0.17D0, -5.92D0]     !% dew point temperature
        RH(J,:) = dexp((16.78D0*DEW(J,:)-116.9D0)/(DEW(J,:)+237.3D0))/(dexp((16.78D0*TAVG(J,:)-116.9D0)/(TAVG(J,:)+237.3D0)))
        RAD(J,:) = [ 5.11D0, 8.42D0, 12.70D0, 14.46D0, 20.07D0, 22.63D0, 22.54D0, 19.48D0, 14.41D0, 10.68D0, 5.61D0, 4.44D0]
        WND(J,:) = [ 4.77D0, 4.75D0, 4.65D0, 4.54D0, 4.17D0, 3.74D0, 3.41D0, 3.43D0, 3.70D0, 3.99D0, 4.46D0, 4.51D0]
        PCP(J,:) = [ 38.10D0, 33.90D0, 48.60D0, 64.70D0, 66.30D0, 77.80D0, 76.80D0, 69.60D0, 64.20D0, 61.30D0, 60.80D0, 51.10D0]
        ENDDO
      endif
      
      if (.not. associated(month%tempM)) then
        allocate(month%tempM(nss))
        allocate(month%tmxM(nss))
        allocate(month%tmnM(nss))
        allocate(month%radM(nss))
        allocate(month%wndM(nss))
        allocate(month%hmdM(nss))
        allocate(month%dewM(nss))
        allocate(month%prcpM(nss))
        status = 1
      endif

      if (M > 0 .and. M<= 12 ) then
        !m =struct('prcp',PCP(M),'temp',TAVG(M),'tmx',TMX(M),'tmn',TMN(M),'rad',RAD(M),'wnd',WND(M),'hmd',RH(M),'dew',DEW(M));month%prcp = PCP(M)
        month%tempM = TAVG(:,M)
        month%prcpM = PCP(:,M)
        month%tmxM = TMX(:,M)
        month%tmnM = TMN(:,M)
        month%radM = RAD(:,M)
        month%wndM = WND(:,M)
        month%hmdM = RH(:,M)
        month%dewM = DEW(:,M)

        month%prcp = PCP(is,M)
        month%temp = TAVG(is,M)
        month%tmx = TMX(is,M)
        month%tmn = TMN(is,M)
        month%rad = RAD(is,M)
        month%wnd = WND(is,M)
        month%hmd = RH(is,M)
        month%dew = DEW(is,M)
        !% note prcp is monthly average (month total)
      else
        !m =struct('prcp',PCP,'temp',TAVG,'tmx',TMX,'tmn',TMN,'rad',RAD,'wnd',WND,'hmd',RH,'dew',DEW);
      endif

      end subroutine wMonthlyAverages


      
          subroutine readwgn()
          ! same format as swat Statwgn database file, but remove LSTATION field!! 
          ! can't handle that field!
          integer fid, eof
          character*100 :: ST,Station,StationL1,StationL2
          character*(100) fmtstr
          fid = 50
          OPEN(UNIT = fid, FILE='Statwgn.csv',ACTION='READ', IOSTAT = eof )
          m = 0
          IF (eof .eq. 0) then
            read(fid,*,iostat = eof)
            j = 1
            do while (eof .eq. 0)
              read(fid,*,IOSTAT = eof ) ST,Station,id(j),(pos(j,k),k=1,2),welev(j),rain_yrs,(tmpmx(j,k),k=1,12),   &
         &     (tmpmn(j,k),k=1,12),(tmpstdmx(j,k),k=1,12),(tmpstdmn(j,k),k=1,12),(pcpmm(j,k),k=1,12),(pcpstd(j,k),k=1,12),  &
         &     (pcpskw(j,k),k=1,12),(pr_w1(j,k),k=1,12),(pr_w2(j,k),k=1,12),(pcpd(j,k),k=1,12),(rainhhmx(j,k),k=1,12),      &
         &     (solarav(j,k),k=1,12),(dewpt(j,k),k=1,12),(wndav(j,k),k=1,12)
              j = j + 1
            enddo
            nSta = j - 2
            CLOSE(fid)
          ENDIF

          end subroutine readwgn

      end module