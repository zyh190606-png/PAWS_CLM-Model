
    !
    subroutine scheduleSolar(Year)
    !  % (E) run through the climatic input to prepare the solar radiation for a
    !  % year using Spokas' method for every station, and record the air
    !  % transmittance, tau, for each day
    use Vdata
    use Solar_Radiation_Mod, Only: Solar_Radiation
    USE shr_orb_mod, Only: shr_orb_cosz, shr_orb_params, shr_orb_decl
    use clm_varorb, only : eccen, mvelpp, lambm0, obliqr
    implicit none
    !global w
    integer :: Year, i, days(366), tt, date(1), p, p_e, n, ii,H,K,J
    integer, dimension(:), allocatable :: idx
    !type(Stations_type),dimension(size(w%Stations)) :: Sta
    type(Stations_type),dimension(:),pointer :: Sta=>NULL() ! using it as an alias
    real*8,dimension(size(w%Stations)) ::t_noon
    real*8 :: wea(366,4)
    real*8, dimension(:), pointer :: Rs=>NULL(), tau=>NULL()
    real*8, dimension(:), pointer :: cosz=>NULL(), RD=>NULL(), RI=>NULL()
    real*8, dimension(size(w%Stations)) :: t_n
    REAL*8 decl(366*24), JDG
    real*8:: obliq, eccf,mvelp
    integer, save :: solarStatus = 0;

    !Sta = w(wid).Stations;
    Sta => w%Stations
    !%tt = floor(w(wid).t + 1e-10);
    days = [(i, i = 1, 366)]
    !try
    !    %day2=num2str(Year*1e4+101);
    date = Year*10000 + 101
    tt = datenum2(date)
    !catch
    !    str = ['01-Jan-',num2str(Year)];
    !    tt = datenum(str);
    !end
    !%t_noon=[11.6539   11.6651   11.6206 11.6262   11.6384   11.6606 11.6606  11.6206]; % annual mean value from previous run
    !% because no material difference was found using the full hourly formula
    !% for solar noon
    !t_noon = zeros(1,length(Sta));


    call shr_orb_params( year , eccen  , obliq , mvelp     ,            &
        &               obliqr   , lambm0 , mvelpp, .TRUE. )

    K = 0
    DO J = 1,366
        DO H = 1,24
            JDG      = J + (H - 0.5D0 - w%timezone)/24D0
            K        = K + 1
            CALL shr_orb_decl(JDG ,eccen ,mvelpp ,lambm0 ,obliqr ,decl(K) ,eccf)
            ! +0.5D0 is to get the JD (with fraction) at noon
        ENDDO
    ENDDO

    t_noon = 0.0D0
    t_n = 0.0D0
    call display( 't_noon is calculated inside Solar_Radiation')
    do i = 1, size(Sta)
        !if ~isfield(w.Stations(i),'dat_time') || isempty(w.Stations(i).dat_time) || w.Stations(i).dat_time(1)~=Year
        if(solarStatus .eq. 0 .or. w%Stations(i)%dat_time(1) /= Year) then
            p = tt - w%Stations(i)%datenums(1) + 1
            p_e = min(p+365,size(Sta(i)%datenums))
            n = p_e - p + 1
            !idx = p:p_e;
            !wea = [days',Sta(i).tmin(idx),Sta(i).tmax(idx),Sta(i).prcp(idx)];
            wea(1:n,1) = days
            wea(1:n,2) = Sta(i)%tmin(p:p_e)
            wea(1:n,3) = Sta(i)%tmax(p:p_e)
            wea(1:n,4) = Sta(i)%prcp(p:p_e)
            if (n .eq. 365) wea(366,:)=wea(365,:)
            !try
            !  [Rs, t_n(i), tau, cosz] = Solar_Radiation(Sta(i).LatLong(1), Sta(i).LatLong(2), Sta(i).XYElev(3),...
            !      Year, wea,t_noon(i));
            !catch
            !  4
            !end
            if(associated(Rs)) deallocate(RS)
            if(associated(RD)) deallocate(RD)
            if(associated(RI)) deallocate(RI)
            if(associated(tau)) deallocate(tau)
            if(associated(cosz)) deallocate(cosz)
            call Solar_Radiation(Sta(i)%LatLong(1),Sta(i)%LatLong(2),Sta(i)%XYElev(3), &
                &       Year,wea,t_noon(i),RD,RI,Rs,t_n(i),tau,cosz,decl,w%timezone)

            if(.not. associated(w%Stations(i)%Rad)) allocate(w%Stations(i)%Rad(size(Rs)))
            w%Stations(i)%Rad(1:size(Rs)) = Rs*(86400.0D0/1D6)      !% converting W/m2 to MJ/day/m2
            w%Stations(i)%RD(1:size(RD)) = RD
            w%Stations(i)%RI(1:size(RI)) = RI
            if(.not. associated(w%Stations(i)%tau)) allocate(w%Stations(i)%tau((size(tau)-23)/24+1))
            w%Stations(i)%tau(1:size(tau(1:size(tau)-23:24))) = tau(1:size(tau)-23:24)
            !%rad = w(wid).Stations(i).Rad(ptr1:ptr1+23);
            !%[Rs, d_T, tau] = Solar_Radiation(latitude, longitude, elevation, year, weather_data)
            if(.not. associated(w%Stations(i)%dat_time)) allocate(w%Stations(i)%dat_time(1))
            w%Stations(i)%dat_time(1) = Year
            if(.not. associated(w%Stations(i)%cosz)) allocate(w%Stations(i)%cosz(size(cosz)))
            w%Stations(i)%cosz(1:size(cosz)) = cosz
            w%Stations(i)%decl = decl
        endif
    enddo
    !4;
    solarStatus = 1
    end subroutine scheduleSolar
