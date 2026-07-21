
    !
    module degrees

    implicit none

    private

    integer, parameter :: sp = kind(1.e0), dp = kind(1.d0)
    real(sp), parameter :: pio180_sp = 4.0e0 * atan(1.e0_sp) / 180.e0_sp
    real(dp), parameter :: pio180_dp = 4.0d0 * atan(1.d0) / 180.d0

    public sind, dsind, dcosd, dasind, dacosd

    interface sind
    module procedure  sind ! Single precision
    module procedure dsind ! Double precision
    end interface

    contains

    function sind(x)
    implicit none
    real(sp) x
    real(sp) sind
    sind = sin(pio180_sp * x)
    end function sind

    function dsind(x)
    implicit none
    real(dp) x
    real(dp) dsind
    dsind = sin(pio180_dp * x)
    end function dsind

    function dcosd(x)
    implicit none
    real(dp) x
    real(dp) dcosd
    dcosd = cos(pio180_dp * x)
    end function dcosd

    function dasind(x)
    implicit none
    real(dp) x
    real(dp) dasind
    dasind = asin(x) / pio180_dp
    end function dasind

    function dacosd(x)
    implicit none
    real(dp) x
    real(dp) dacosd
    dacosd = acos(x) / pio180_dp
    end function dacosd
    end module degrees

    !=================================================================================
    MODULE Solar_Radiation_Mod


    use Vdata
    use degrees
    PUBLIC:: solar_noon,Solar_Radiation,Gregorian_time
    ! methods from Spokas and others for solar radiation generation in PAWS
    ! solar obliquity may be better described in CLM
    PRIVATE
    CONTAINS
    subroutine solar_noon(Latitude, Longitude, Date, time_rise, time_sn, time_set)
    !  % calculate the solar noon for the given latitude, longitude and date
    !  % Input arg: Latitude = North Latitude (35N = 35, 25S = -25)
    !  %            Longitude = East Longitude (75E = 75, 45W = -45)
    !  %            date = eg. ('7-sep-2002')
    !  % Output arg:   time_rise = time of sun rise
    !  %               time_sn = time of solar noon
    !  %               time_set = time of sun set
    implicit none
    real*8 :: Latitude, Longitude
    integer, dimension(3) :: Date
    real*8 :: time_rise, time_sn, time_set
    integer :: J_date, n

    real*8 :: n_star, J_star, M, C, lamda, J_transit, error, J_transit_new
    real*8 :: delta, H, J_ss, J_set, J_rise, d_set, d_rise, zd
    integer, dimension(3), parameter :: fixed_date = [2000, 1, 1]
    real temp


    !%J_date = daysact('1-Jan-2000', Date) + 2451545;   % Julian date
    J_date = datenum2(Date) - datenum2(fixed_date) + 2451545
    n_star = (dble(J_date) - 2451545.0D0 - 0.0009D0) + (Longitude/360.0D0)
    n = NINT(n_star)     !% Julian cycle
    J_star = 2451545.0D0 + 0.0009D0 - (Longitude/360.0D0) + dble(n)
    M = dmod(357.5291D0 + 0.98560028D0*(J_star - 2451545.0D0), 360.0D0)   !% mean solar anomaly
    C = (1.9148D0*dsind(M)) + (0.0200D0*dsind(2.0D0*M)) + (0.0003D0*dsind(3.0D0*M))   !% equation of center
    lamda = dmod(M+102.9372D0+C+180.0D0, 360.0D0)     !% ecliptical longitude of the sun
    J_transit = J_star + (0.0053D0*dsind(M)) - (0.0069D0*dsind(2.0D0*lamda))      !% accurate Julian date for solar noon
    error = 1.0D0
    do while (error > 1e-15)
        M = dmod(357.5291D0 + 0.98560028D0*(J_transit - 2451545), 360.0D0)      !% mean solar anomaly
        C = (1.9148D0*dsind(M)) + (0.0200D0*dsind(2.0D0*M)) + (0.0003D0*dsind(3.0D0*M))      !% equation of center
        lamda = dmod(M+102.9372D0+C+180.0D0, 360.0D0)       !% ecliptical longitude of the sun
        J_transit_new = J_star + (0.0053D0*dsind(M)) - (0.0069D0*dsind(2.0D0*lamda))        !% accurate Julian date for solar noon
        error = dabs(J_transit_new - J_transit)
        J_transit = J_transit_new
    enddo
    delta = dasind(dsind(lamda)*dsind(23.45D0))       !% declination of the sun
    H = dacosd((dsind(-0.83D0)-dsind(Latitude)*dsind(delta))/(dcosd(Latitude)*dcosd(delta)))    !% hour angle (half the arc length of the sun)
    J_ss = 2451545 + 0.0009D0 + ((H-Longitude)/360.0D0) + n
    J_set = J_ss + (0.0053D0*dsind(M)) - (0.0069D0*dsind(2.0D0*lamda))    !% Julian date of sunset on cycle n
    J_rise = J_transit - (J_set - J_transit)      !% Julian date of sunset on cycle n

    ![d_set] = Gregorian_time(J_set);
    ![d_rise] = Gregorian_time(J_rise);
    call Gregorian_time(J_set, d_set)
    call Gregorian_time(J_rise, d_rise)

    !zd = timezone(-Longitude,'degrees');
    zd = -6.0D0       ! for current model

    time_set = (d_set - dint(d_set))*24.0D0 + zd
    if (time_set < 0.0D0) then
        time_set = 24.0D0 + time_set
    endif
    time_rise = (d_rise - dint(d_rise))*24.0D0 + zd
    if (time_rise > 24) then
        time_rise = time_rise - 24.0D0
    endif
    time_sn = (time_rise + time_set) / 2.0D0      !% time of solar noon

    end subroutine solar_noon

    !=============================================================================

    subroutine Gregorian_time(J, d)
    !  % convert Julian day number to the Gregorian time
    !  % Input arg:    J = Julian day number
    !  % Output arg:    day d in month m in year j
    implicit none
    real*8 :: J, d
    real*8 :: p, s1, n, s2, i, s3, q, e, s4, m, jj

    p = dint(J + 0.5D0)
    s1 = p + 68569D0
    n = dint(4.0D0*s1/146097.0D0)
    s2 = s1 - dint((146097.0D0*n + 3.0D0)/4.0D0)
    i = dint(4000.0D0*(s2 + 1.0D0)/1461001.0D0)
    s3 = s2 - dint(1461.0D0*i/4.0D0) + 31.0D0
    q = dint(80.0D0*s3/2447.0D0)
    e = s3 - dint(2447.0D0*q/80.0D0)
    s4 = dint(q/11.0D0)
    m = q + 2.0D0 - 12.0D0*s4
    jj = 100.0D0*(n - 49D0) + i + s4
    d = e + J - p + 0.5D0

    end subroutine Gregorian_time



    !====================================================================================
    subroutine Solar_Radiation(latitude,longitude,elevation,year,weather_data,t_noon,    &
        &          Sb,Sd,Rs,t_n,tau,cosz,decl,timezone)
    ! CP: in a new version of this algorithm, radiation is always computed at half hours,-
    ! e.g., rad(1) is for 00:30, rad(2) is for 01:30, etc, to better represent a time step.-
    ! This requires that the declination angle computed outside is also computed at these times.
    ! Prior to this change, rad was calculated at the beginning of the time steps
    ! The second change is that we now use CLM's TOA radiation algorithm. This means that-
    ! solar_noon algorithm is not used. Rather, local time is turned into Greenwich JD, along-
    ! with latitude, sent into shr_orb_cosz


    !  % (E) with numerical input, calculate the solar radiation for a given year
    !  % Input arg:    latitude = North Latitude (35N = 35, 25S = -25), degree
    !  %               longitude = East Longitude (75E = 75, 45W = -45), degree
    !  %               elevation (m)
    !  %               year = eg. (2000)
    !  %               weatherdata, which has the following format
    !  % Day of year, min. air temp, max. air temp, rainfall
    !  % Each day of year is on a separate line, missing values can be skipped.
    !  % Output arg:   Rs = solar radiation (W/m2)
    !  %               d_T = T_air_max - T_air_min, extreme temperature change of the
    !  %               tau value for each day

    !  % Reference: Spokas 2006, Weed Science.
    !  % Implemented by Niu Jie, MSU
    USE shr_orb_mod, Only: shr_orb_cosz
    USE VData, Only: timer
    implicit none
    real*8 :: latitude, longitude, elevation, t_noon
    real*8, dimension(366,4) :: weather_data
    integer :: year
    real*8, dimension(:), pointer :: Rs, tau, cosz, Sd, Sb
    real*8 :: t_n
    real*8, parameter :: pi = 3.141592653589793238462D0
    real*8 :: d_T, Phi, a, Pa, Sp0, sP, cP, sDD, delta_SD, coszz, Psi, cos_is, m, alpha
    real*8, dimension(366) :: w_days, min_tmp, max_tmp, precipitation
    real*8, dimension(:), allocatable :: t_rise, t_sn, t_set
    integer :: days, n, J, today, yesterday, i, t,H
    integer, dimension(3) :: date
    character*32 :: date_char, time_char
    real*8 :: time_start, time_end
    real*8 :: decl(366*24),JDG
    real*8 :: lat, lon, timezone,tt

    !tic
    time_start = 0d0
    call timer(time_start)
    d_T = 0.0D0
    Phi = latitude/180.0D0*pi         !% convert latitude to radians
    lat = Phi; lon = longitude*pi/180D0
    tt  = pi

    a = elevation                     !% site elevation (m)
    Pa = 101.30D0*exp(-a/8200.0D0)   !% average barometric pressure (kPa)
    !Sp0 = 1367.0D0                    !% solar constant (W m^{-2})
    Sp0 = 1320.0D0                    !% solar constant (W m^{-2}) --- this is excluding ultraviolet (visible + near infrared)
    ! THIS WAS 1360, CHAOPENG CHANGED IT TO 1367
    !days = yeardays(year);
    if (eomday(year, int(2)) == 29) then
        days = 366
    else
        days = 365
    endif
    sP = dsin(Phi)     !% sP is the same for the same site, forever
    cP = dcos(Phi)
    !% calendar days have available weather data
    w_days = weather_data(:, 1)
    !% minimum air temperature
    min_tmp = weather_data(:, 2)
    !% maximum air temperature
    max_tmp = weather_data(:, 3)
    !% precipitation
    precipitation = weather_data(:, 4)

    !Rs = zeros(24*days, 1);
    !cosz = zeros(24*days, 1);
    if(associated(Rs)) deallocate(RS)
    if(associated(cosz)) deallocate(cosz)
    if(associated(tau)) deallocate(tau)
    if(associated(Sb)) deallocate(Sb)
    if(associated(Sd)) deallocate(Sd)
    allocate(Rs(24*days))
    allocate(cosz(24*days))
    allocate(tau(24*days))
    Rs = 0.0D0
    cosz = 0.0D0
    allocate(Sd(24*days))
    allocate(Sb(24*days))
    Sd = Rs
    Sb = Rs
    tau = Rs
    !t_rise = zeros(days, 1); t_sn = t_rise; t_set = t_rise;
    allocate(t_rise(days))
    allocate(t_sn(days))
    allocate(t_set(days))
    t_rise = 0.0D0
    t_sn = 0.0D0
    t_set = 0.0D0
    n = 0

    do J = 1, days
        !% J = calendar day

        !% deal with atmospheric transmittance with respect to weather
        !% check if the weather data available for today
        !today = find(w_days == J, 1);
        !if isempty(today)
        !    today = 0;
        !end
        today = 0
        do i = 1, size(w_days)
            if(int(w_days(i))==J) then
                today = i
                exit
            endif
        enddo
        !% check if the weather data available for yesterday
        !yesterday = find(w_days == (J-1), 1);
        !if isempty(yesterday)
        !    yesterday = 0;
        !end
        yesterday = 0
        do i = 1, size(w_days)
            if(int(w_days(i))==J-1) then
                yesterday = i
                exit
            endif
        enddo
        tau(1+24*(J-1) : 24*J) = 0.70D0    !% by default
        if (J == 1) then
            if (today /= 0) then
                !% No precipitation and d_T > 10C
                if ((max_tmp(today) - min_tmp(today) > 10.0D0) .and. (dabs(precipitation(today)) < 1.0D0)) then
                    tau(1+24*(J-1) : 24*J) = 0.7D0
                endif
                !% precipitation today, not not for yesterday
                if (dabs(precipitation(today)) > 1.0D0) then
                    tau(1+24*(J-1) : 24*J) = 0.4D0
                endif
            endif
        else
            !% no weather data available for both days
            if ((today==0) .and. ( (yesterday==0) .or. ( (yesterday + 1) /= today ) )) then
                !% atmospheric transmittance
                tau(1+24*(J-1) : 24*J) = 0.7D0     !% by default
            else
                !% no weather data available for today but yes for yesterday
                if ((today==0) .and. (yesterday + 1 == today)) then
                    !% precipitation for yesterday
                    if (dabs(precipitation(yesterday)) > 1.0D0) then
                        tau(1+24*(J-1) : 24*J) = 0.6D0
                    else
                        !% no precipitation for yesterday
                        tau(1+24*(J-1) : 24*J) = 0.7D0
                    endif
                else
                    !% no weather data available for yesterday but yes for today
                    if (( (yesterday==0) .or. ( (yesterday + 1) /= today ) ) .and. (today/=0)) then
                        !% No precipitation and d_T > 10C
                        if ((max_tmp(today) - min_tmp(today) > 10.0D0) .and. (dabs(precipitation(today)) < 1.0D0)) then
                            tau(1+24*(J-1) : 24*J) = 0.7D0
                        endif
                    else
                        !% make sure weather data available for yesterday
                        if (yesterday + 1 == today) then
                            !% No precipitation and d_T > 10C
                            if ((max_tmp(today) - min_tmp(today) > 10.0D0) .and. &
                                &         (dabs(precipitation(today)) < 1.0D0) .and. (dabs(precipitation(yesterday)) < 1.0D0)) then
                                tau(1+24*(J-1) : 24*J) = 0.7D0
                            endif
                            !% no precipitation today, but precipitation fell the
                            !% previous day
                            if ((dabs(precipitation(today)) < 1.0D0) .and. (dabs(precipitation(yesterday)) > 1.0D0)) then
                                tau(1+24*(J-1) : 24*J) = 0.6D0
                            endif
                            !% precipitation today, not not for yesterday
                            if ((dabs(precipitation(today)) > 1.0D0) .and. (dabs(precipitation(yesterday)) < 1.0D0)) then
                                tau(1+24*(J-1) : 24*J) = 0.4D0
                            endif
                            !% prcipitation for both days
                            if ((abs(precipitation(today)) > 1.0D0) .and. (abs(precipitation(yesterday)) > 1.0D0)) then
                                tau(1+24*(J-1) : 24*J) = 0.3D0
                            endif
                        endif
                    endif
                endif
            endif
        endif
        if ((today/=0) .and. (dabs(latitude) < 60.0D0)) then
            if (dint(max_tmp(today)) - dint(min_tmp(today)) <= 10.0D0) then
                tau(1+24*(J-1) : 24*J) = tau(1+24*(J-1) : 24*J) / dble(11 - (dint(max_tmp(today)) - dint(min_tmp(today))))
            endif
        endif

        !%plot(1+24*(J-1) : 24*J, tau(1+24*(J-1) : 24*J), 'o'); hold on
        !%xlim([1+24*(J-10), 24*(J+3)]);
        !%ylim([0 0.8]);

        !% delta_SD = solar declination angle (radians)
        !%delta_SD = asin( 0.39785 * sin(4.869 + 0.0172*J + 0.03345*sin(6.2238 + 0.0172*J)) );
        !sDD = 0.39785D0 * dsin(4.869D0 + 0.0172D0*J + 0.03345D0*dsin(6.2238D0 + 0.0172D0*J))  ! NJ's old code
        !delta_SD = dasin(sDD)   ! NJ's old code
        !delta_SD = decl(J)

        !%plot(J, delta_SD, 'o'); hold on
        !% convert calendar day to date DD-MM-YYYY
        !%date = datestr(daysadd(['01-Jan-', num2str(year)], J-1));
        !%date = datestr(datenum(['01-Jan-', num2str(year)]) + J-1);
        date = [year, 1, 1]
        !% calculate the time of solar noon
        ![t_rise(J), t_sn(J), t_set(J)] = solar_noon(latitude, longitude, date);
        !call solar_noon(latitude, longitude, date, t_rise(J), t_sn(J), t_set(J))
        !%t_sn(J) = t_noon;
        do H = 1, 24
            delta_SD = decl((J-1)*24+H) ! Declination angle is computed houly at half-hours
            JDG      = J + (H - 0.5D0 - timezone)/24D0  ! Greenwich Julian Day (1.xxx to 365.xx) of half-hours
            sDD = sin(delta_SD)
            n = n + 1
            !% Psi = zenith angle (radians)
            !% Originally:
            !% Psi = acos(sin(Phi)*sin(delta_SD) + cos(Phi)*cos(delta_SD)*cos(1/12*pi*(t-t_sn(J))));
            !% alpha = acos((sin(pi/2-Psi)*sin(Phi) - sin(delta_SD)) / (cos(pi/2-Psi)*cos(Phi)));

            !%sD = sin(delta_SD);
            !tt = sP*sDD + cP*dcos(delta_SD)*dcos(1.0D0/12.0D0*pi*(H -0.5D0-t_sn(J)))
            !coszz = sP*sDD + cP*dcos(delta_SD)*dcos(1.0D0/12.0D0*pi*(H -0.5D0-t_sn(J)))
            !CP comments: the solar_noon algorithm is kind of strange, how do you know the time zone? H can be arbitary (e.g. summer time)
            ! need to check to source literature
            coszz = shr_orb_cosz(JDG,lat,lon,delta_SD)
            !coszz = tt
            cosz(n) = coszz;
            Psi = dacos(coszz)
            !% alpha = solar azimuth angle
            !%alpha = acos((sin(pi/2-Psi)*sP - sD) / (cos(pi/2-Psi)*cP));
            !alpha = dacos((coszz*sP - sDD) / (dcos(pi/2.0D0-Psi)*cP))
            !if (t < 12) then
            !    alpha = - alpha
            !endif
            !%plot(n, Psi,'o'); hold on
            !%plot(n, cos(Psi)*100,'*'); hold on
            !% m = optical air mass number
            m = Pa / (101.3D0*coszz)
            !%subplot(1,2,1)
            !%plot(n, m,'o'); hold on
            !if ((t+1 > dint(t_rise(J))) .and. (t < dint(t_set(J)))) then
            if (coszz > 0D0) then
                !% beam radiation
                if (Psi*180.0D0/pi < 90.0D0) then
                    !% diffuse radiation
                    !%Sd(n) = 0.3*(1-tau(n)^m)*Sp0*cos(Psi);
                    Sd(n) = 0.30D0*(1.0D0-tau(n)**m)*Sp0*coszz
                    !%cos_is = sin(pi/2-Psi)*cos(0) + cos(pi/2-Psi)*cos(pi-alpha)*sin(0);
                    cos_is = coszz      !% + cos(pi/2-Psi)*cos(pi-alpha)*sin(0);
                    Sb(n) = Sp0*tau(n)**m*cos_is
                else
                    Sb(n) = 0.0D0
                    Sd(n) = 0.0D0
                endif
                !% total solar radiation
                Rs(n) = Sb(n) + Sd(n)
            else
                Sb(n) = 0.0D0
                Sd(n) = 0.0D0
                Rs(n) = 0.0D0
            endif
            !write(333,*) Sb(n),Sd(n)
            !%subplot(2,1,1)
            !%plot(n, Sb(n), 'r+', n, Sd(n), 'b*', n, Rs(n), 'ko'); hold on
            !%xlim([1+24*(J-3), 24*(J)]);
            !%drawnow
        enddo
        !%subplot(2,1,1)
        !%plot(J, max(Rs(1+24*(J-1) : 24*(J))), '.-'); hold on;
        !%subplot(2,1,2)
        !%plot(J, t_rise(J), 'b*', J, t_sn(J), 'ro', J, t_set(J), 'gs', J, t_set(J)-t_rise(J), 'k^'); hold on;
        !%drawnow

    enddo
    !toc;
    time_end = 0d0
    call timer(time_end)
    !write(*,'(A16,F9.6,A9)') 'Elapsed time is ', time_end - time_start, ' seconds.'
    call display('Elapsed time is ', time_end - time_start, ' seconds.')
    !t_n = mean(t_sn);
    t_n = sum(t_sn)/size(t_sn)

    !if(allocated(Sd)) deallocate(Sd)
    !if(allocated(Sb)) deallocate(Sb)
    if(allocated(t_rise)) deallocate(t_rise)
    if(allocated(t_sn)) deallocate(t_sn)
    if(allocated(t_set)) deallocate(t_set)

    end subroutine Solar_Radiation

    END MODULE Solar_Radiation_Mod