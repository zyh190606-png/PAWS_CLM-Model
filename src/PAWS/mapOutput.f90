    !
#include <preproc.h> 

    subroutine mapOutput(cmd, output)
    !% (E) accumulatively output map averages
    !% This produces a long term spatial trend of variables, states and fluxes
    use Vdata
    implicit none
    integer :: cmd
    integer :: nr , nc, nz
    type num_type
        real*8,dimension(:,:),allocatable::Prcp,Inf,PET,ET,EvapG,Dperc,Dperc2,SatE
        real*8,dimension(:,:),allocatable::InfE,botT,Qoc,Qgc,SUBM,SMELT,GW,ovnH
        integer :: n = 0
    end type num_type
    type(num_type), save :: num
    type(maps_type) :: output

    nr = size(g%VDZ%h,1)
    nc = size(g%VDZ%h,2)
    nz = size(g%VDZ%h,3)
    if(.not. allocated(num%Prcp)) then
        allocate(num%Prcp(nr,nc))
        num%Prcp = 0.0D0
    endif
    if(.not. allocated(num%Inf)) then
        allocate(num%Inf(nr,nc))
        num%Inf = 0.0D0
    endif
    if(.not. allocated(num%PET)) then
        allocate(num%PET(nr,nc))
        num%PET = 0.0D0
    endif
    if(.not. allocated(num%ET)) then
        allocate(num%ET(nr,nc))
        num%ET = 0.0D0
    endif
    if(.not. allocated(num%EvapG)) then
        allocate(num%EvapG(nr,nc))
        num%EvapG = 0.0D0
    endif
    if(.not. allocated(num%Dperc)) then
        allocate(num%Dperc(nr,nc))
        num%Dperc = 0.0D0
    endif
    if(.not. allocated(num%Dperc2)) then
        allocate(num%Dperc2(nr,nc))
        num%Dperc2 = 0.0D0
    endif
    if(.not. allocated(num%SatE)) then
        allocate(num%SatE(nr,nc))
        num%SatE = 0.0D0
    endif
    if(.not. allocated(num%InfE)) then
        allocate(num%InfE(nr,nc))
        num%InfE = 0.0D0
    endif
    if(.not. allocated(num%botT)) then
        allocate(num%botT(nr,nc))
        num%botT = 0.0D0
    endif
    if(.not. allocated(num%Qoc)) then
        allocate(num%Qoc(nr,nc))
        num%Qoc = 0.0D0
    endif
    if(.not. allocated(num%Qgc)) then
        allocate(num%Qgc(nr,nc))
        num%Qgc = 0.0D0
    endif
    if(.not. allocated(num%SUBM)) then
        allocate(num%SUBM(nr,nc))
        num%SUBM = 0.0D0
    endif
    if(.not. allocated(num%SMELT)) then
        allocate(num%SMELT(nr,nc))
        num%SMELT = 0.0D0
    endif
    if(.not. allocated(num%GW)) then
        allocate(num%GW(nr,nc))
        num%GW = 0.0D0
    endif
    if(.not. allocated(num%ovnH)) then
        allocate(num%ovnH(nr,nc))
        num%ovnH = 0.0D0
    endif


    if (cmd==0) then
        !% recharge, evap, PET, ET, Inf, Dperc
        num%Prcp = (num%Prcp+g%OVN%prcp)
        num%Inf = (num%Inf+g%VDZ%DF)
        num%PET = (num%PET+g%Veg%PET)
        num%ET = (num%ET+g%Veg%ET)
        num%EvapG = (num%EvapG+g%Veg%EvapG)
        num%Dperc = (num%Dperc+g%VDZ%Dperc)
        num%Dperc2 = (num%Dperc2+g%GW(1)%Qperc)
        num%SatE = num%SatE + g%VDZ%R(:,:,2) * g%OVN%Rf
        num%InfE = num%InfE + g%OVN%Rf * (1 - g%VDZ%R(:,:,2))
        num%botT = num%botT + g%VDZ%R(:,:,1)
        num%SUBM = num%SUBM + g%Veg%SNOM
        num%SMELT = num%SMELT + g%VDZ%SMR
        num%GW = num%GW + g%GW(1)%h
        num%ovnH = num%ovnH + g%OVN%h
        !output = 0.0D0
        num%n = num%n + 1   !% number of times it is added
    else    !% output and clear data
        !output = num
        output%Prcp = num%Prcp
        output%Inf = num%Inf
        output%PET = num%PET
        output%ET = num%ET
        output%EvapG = num%EvapG
        output%Dperc = num%Dperc
        output%Dperc2 = num%Dperc2
        output%SatE = num%SatE
        output%InfE = num%InfE
        output%botT = num%botT
        output%Qoc = num%Qoc
        output%Qgc = num%Qgc
        output%SUBM = num%SUBM
        output%SMELT = num%SMELT
        output%GW = num%GW
        output%ovnH = num%ovnH
        output%n = num%n
        !clear num
        if(allocated(num%Prcp)) deallocate(num%Prcp)
        if(allocated(num%Inf)) deallocate(num%Inf)
        if(allocated(num%PET)) deallocate(num%PET)
        if(allocated(num%ET)) deallocate(num%ET)
        if(allocated(num%EvapG)) deallocate(num%EvapG)
        if(allocated(num%Dperc)) deallocate(num%Dperc)
        if(allocated(num%Dperc2)) deallocate(num%Dperc2)
        if(allocated(num%SatE)) deallocate(num%SatE)
        if(allocated(num%InfE)) deallocate(num%InfE)
        if(allocated(num%botT)) deallocate(num%botT)
        if(allocated(num%Qoc)) deallocate(num%Qoc)
        if(allocated(num%Qgc)) deallocate(num%Qgc)
        if(allocated(num%SUBM)) deallocate(num%SUBM)
        if(allocated(num%SMELT)) deallocate(num%SMELT)
        if(allocated(num%GW)) deallocate(num%GW)
        if(allocated(num%ovnH)) deallocate(num%ovnH)
        num%n = 0
    endif

    end subroutine mapOutput



    subroutine mapOutput2(cmd, num)
    !% (E) accumulatively output map averages
    !% This produces a long term spatial trend of variables, states and fluxes
    use Vdata
    use memoryChunk, Only: tAVG, tempStorage
    use pawsControl, Only: CaliStartTime
#if (defined CLM)
#if (defined CLM45)
    use PhotosynthesisType     , only : photosyns_type
    use CNCarbonFluxType       , only : carbonflux_type
    use EnergyFluxType         , only : energyflux_type
    use SolarAbsorbedType      , only : solarabs_type
    use TemperatureType        , only : temperature_type
    use WaterstateType         , only : waterstate_type
    use CanopyStateType        , only : canopystate_type
#else
    USE clmtype
    use subgridAveMod
    use filtermod
#endif       
#endif      
    implicit none
    real*8, pointer :: ptrp(:)         ! pointer to input pft array
    real*8, pointer :: ptrc(:)         ! pointer to output column array
    real*8 tempdata(size(g%VDZ%h,1), size(g%VDZ%h,2))
    real*8 tempdata2(size(g%VDZ%h,1), size(g%VDZ%h,2))
    integer :: cmd, i, j, rid, nni
    integer :: nr , nc, nz, nzz,nv, ndays = 121 ! From September 1 to end of year is 121 days
    type(maps_type):: num
    real*8, pointer:: p2d(:,:), p3d(:,:,:),temphead(:,:),Qgc1D(:),temphead1d(:)
    logical varianceMap
    integer ds(1),de(1)
    integer, dimension(:), pointer:: idx

    varianceMap = .true.
    ! Start counting only after ndays from simulation
    !IF (w%t - w%ModelStartTime > ndays) THEN
    !ds = int(20060101); de = int(20091231);
    !IF (w%t >= datenum2(ds) .and. w%t <= datenum2(de) ) THEN
    !IF (w%t >= (w%ModelStartTime + ndays)) THEN
    !IF (.TRUE.) THEN
    IF (cmd<10) THEN ! recording operation

        IF (w%t >= CaliStartTime) THEN

            nr = size(g%VDZ%h,1)
            nc = size(g%VDZ%h,2)
            nz = size(g%VDZ%h,3)
            nv = size(g%Veg%RPT,3)

            !% recharge, evap, PET, ET, Inf, Dperc
            num%Prcp = (num%Prcp+g%Veg%Gprcp)
            num%Inf = (num%Inf+g%VDZ%DF)
            num%PET = (num%PET+g%Veg%PET)
            num%ET = (num%ET+g%Veg%ET)
            num%EvapG = (num%EvapG+g%Veg%EvapG)
            num%Dperc = (num%Dperc+g%VDZ%Dperc)
            num%Dperc2 = (num%Dperc2+g%GW(1)%Qperc)
            num%SatE = num%SatE + g%VDZ%R(:,:,2) * g%OVN%Rf
            num%InfE = num%InfE + g%OVN%Rf * (1 - g%VDZ%R(:,:,2))
            num%botT = num%botT + g%VDZ%R(:,:,1)
            num%SUBM = num%SUBM + g%Veg%SNOM
            num%SMELT = num%SMELT + g%VDZ%SMR
            num%GW = num%GW + g%GW(1)%h
            num%ovnH = num%ovnH + g%OVN%h

            !call tempStorage('Qgc',Qgc,1)
            Qgc1D => OneD_ptr(size(g%OVN%h),num%Qgc)
            DO i=1,size(w%CR)
                rid = w%CR(i)
                idx => g%Riv%Lidx(i)%ptr
                Qgc1D(idx) = Qgc1D(idx) + r(rid)%Riv%fG
            ENDDO
            !%Qgc = num%Qgc + Qgc
            !Qgc = 0D0
            !output = 0.0D0
            num%n = num%n + 1   !% number of times it is added

#if (defined CLM) 
#if (defined CLM45)
            i = 1;
            p2d => twoD_ptr(nr,nc,carbonflux_vars%npp_col)
            num%clmmap(:,:,i)= num%clmmap(:,:,i)+p2d;      i = i + 1;  ! column average npp (gC/m2/s)
            p2d => twoD_ptr(nr,nc,carbonflux_vars%gpp_col)
            num%clmmap(:,:,i)= num%clmmap(:,:,i)+p2d;      i = i + 1;  ! column average gpp
            num%clmmap(:,:,i)= num%clmmap(:,:,i)+g%Veg%TP; i = i + 1;   ! column average TP
            num%clmmap(:,:,i)= num%clmmap(:,:,i)+g%OVN%Ex; i = i + 1;   ! column average Exfiltration

            ! i = 5
            ptrp => energyflux_vars%eflx_lh_tot_patch ! latent heat flux (W/m**2)
            ptrc => oneD_ptr(nr*nc,tempdata)
            ptrc = 0D0
            call p2c (filter(1)%num_nolakec, filter(1)%nolakec, ptrp, ptrc) ! no lake either
            num%clmmap(:,:,i) = num%clmmap(:,:,i)+ tempdata;

            i = i + 1; ! i = 6
            ptrp => energyflux_vars%eflx_sh_tot_patch ! sensible heat flux (W/m**2)
            ptrc => oneD_ptr(nr*nc,tempdata)
            ptrc = 0D0
            call p2c (filter(1)%num_nolakec, filter(1)%nolakec, ptrp, ptrc) ! no lake either
            num%clmmap(:,:,i) = num%clmmap(:,:,i)+ tempdata

            i = i + 1; ! i = 7
            ptrp => solarabs_vars%parsun_z_patch !average absorbed PAR for sunlit leaves (W/m**2)
            ptrc => oneD_ptr(nr*nc,tempdata)
            ptrc = 0D0
            call p2c (filter(1)%num_nolakec, filter(1)%nolakec, ptrp, ptrc) ! no lake either
            num%clmmap(:,:,i) = num%clmmap(:,:,i)+ tempdata

            i = i + 1; ! i = 8
            ptrp => solarabs_vars%parsha_z_patch !average absorbed PAR for shaded leaves (W/m**2)
            ptrc => oneD_ptr(nr*nc,tempdata)
            ptrc = 0D0
            call p2c (filter(1)%num_nolakec, filter(1)%nolakec, ptrp, ptrc) ! no lake either
            num%clmmap(:,:,i) = num%clmmap(:,:,i)+ tempdata

            i = i + 1; ! i = 9
            ptrp => photosyns_vars%fpsn_patch   !photosynthesis (umol CO2 /m**2 /s)
            ptrc => oneD_ptr(nr*nc,tempdata)
            ptrc = 0D0
            call p2c (filter(1)%num_nolakec, filter(1)%nolakec, ptrp, ptrc) ! no lake either
            num%clmmap(:,:,i) = num%clmmap(:,:,i)+ tempdata

            !j = i + 1
            i = i + 1; ! i = 10
            ptrp => temperature_vars%t_grnd_col !ground temperature (Kelvin)
            ptrc => oneD_ptr(nr*nc,tempdata)
            ptrc = ptrp
            num%clmmap(:,:,i) = num%clmmap(:,:,i)+ tempdata

            !j = i + 1
            i = i + 1; ! i = 11
            ptrp => temperature_vars%t_soi10cm_col !ground temperature (Kelvin)
            ptrc => oneD_ptr(nr*nc,tempdata)
            ptrc = ptrp
            num%clmmap(:,:,i) = num%clmmap(:,:,i)+ tempdata - 273.15D0

            !j = i + 1
            i = i + 1; ! i = 12
            ptrp => waterstate_vars%h2osoi_liqice_10cm_col !liquid water + ice lens in top 10cm of soil (kg/m2)
            ptrc => oneD_ptr(nr*nc,tempdata)
            ptrc = ptrp
            num%clmmap(:,:,i) = num%clmmap(:,:,i)+ tempdata

            !j = i + 1
            i = i + 1; ! i = 13
            ptrp => waterstate_vars%h2osno_col ! snow water (mm H2O)
            ptrc => oneD_ptr(nr*nc,tempdata)
            ptrc = ptrp
            num%clmmap(:,:,i) = num%clmmap(:,:,i)+ tempdata

            !j = i + 1
            i = i + 1; ! i = 14
            ptrp => carbonflux_vars%er_col ! (gC/m2/s) total ecosystem respiration, autotrophic + heterotrophic
            ptrc => oneD_ptr(nr*nc,tempdata)
            ptrc = ptrp
            num%clmmap(:,:,i) = num%clmmap(:,:,i)+ tempdata

            !j = i + 1
            i = i + 1; ! i = 15
            ptrp => carbonflux_vars%hr_col ! (gC/m2/s) total heterotrophic respiration
            ptrc => oneD_ptr(nr*nc,tempdata)
            ptrc = ptrp
            num%clmmap(:,:,i) = num%clmmap(:,:,i)+ tempdata

            !j = i + 1
            i = i + 1; ! i = 16 top 10cm soil h2o content [-] including liquid and ice
            call tempStorage('hcAfterLand',temphead,1)
            ptrp => waterstate_vars%h2osoi_liqice_10cm_col
            ptrc => oneD_ptr(nr*nc,tempdata)
            ptrc = ptrp
            temphead1d => oneD_ptr(nr*nc,temphead)
            num%clmmap(:,:,i) = num%clmmap(:,:,i)+ (tempdata/1000D0 - temphead)/0.1D0

            i = i + 1; ! i = 17
            ptrp => canopystate_vars%tlai_patch ! total leaf area index
            ptrc => oneD_ptr(nr*nc,tempdata)
            ptrc = 0D0
            call p2c (filter(1)%num_nolakec, filter(1)%nolakec, ptrp, ptrc) ! no lake either
            num%clmmap(:,:,i) = num%clmmap(:,:,i)+ tempdata
            if (any( isnan(ptrp) )) then
                print *, 'LAI nan found'
            endif

            i = i + 1; ! i = 18
            ptrp => canopystate_vars%tsai_patch ! total stem area index
            ptrc => oneD_ptr(nr*nc,tempdata)
            ptrc = 0D0
            call p2c (filter(1)%num_nolakec, filter(1)%nolakec, ptrp, ptrc) ! no lake either
            num%clmmap(:,:,i) = num%clmmap(:,:,i)+ tempdata

            i = i + 1; ! i = 19 - 40
            num%clmmap(:,:,i:i+nz-1) = num%clmmap(:,:,i:i+nz-1) + g%VDZ%THE; i = i + nz - 1 ! soil moisture layer-by-layer, set to end of layer, if nz=22,

            i = i + 1  ! i = 41 - 43
            p3d => threeD_ptr((/nr,nc,3/),carbonflux_vars%npp_patch)  ! pft-level npp
            num%clmmap(:,:,i:i+3-1)= num%clmmap(:,:,i:i+3-1)+p3d; i = i + 3 - 1

            i = i + 1; ! i = 44
            p2d => twoD_ptr(nr,nc,carbonflux_vars%nee_col) ! net ecosystem exchange
            num%clmmap(:,:,i) = num%clmmap(:,:,i) + p2d;

            i = i + 1; ! i = 45
            num%clmmap(:,:,i) = num%clmmap(:,:,i) + g%GW(1)%DR; ! lateral inflow (m/day)

            i = i + 1; ! i = 46
            ptrp => canopystate_vars%elai_patch  ! effective leaf area index
            ptrc => oneD_ptr(nr*nc,tempdata)
            ptrc = 0D0
            call p2c (filter(1)%num_nolakec, filter(1)%nolakec, ptrp, ptrc) ! no lake either
            num%clmmap(:,:,i) = num%clmmap(:,:,i)+ tempdata

            i = i + 1; ! i = 47
            ptrp => energyflux_vars%btran_patch  !transpiration wetness factor (0 to 1)
            ptrc => oneD_ptr(nr*nc,tempdata)
            ptrc = 0D0
            call p2c (filter(1)%num_nolakec, filter(1)%nolakec, ptrp, ptrc) ! no lake either
            num%clmmap(:,:,i) = num%clmmap(:,:,i)+ tempdata

            !j = i + 1
            i = i + 1; ! i = 48
            p2d => twoD_ptr(nr,nc,carbonflux_vars%ar_col) ! (gC/m2/s) autotrophic respiration (MR + GR)
            num%clmmap(:,:,i) = num%clmmap(:,:,i)+ p2d;

            !j = i + 1
            i = i + 1; ! i = 49
            num%clmmap(:,:,i) = num%clmmap(:,:,i)+ g%VDZ%R(:,:,2); ! saturation state (==1 for saturated, ==0 for unsaturated). therefore, this field * dt =  total time saturated

            !j = i + 1
            i = i + 1; ! i = 50
            num%clmmap(:,:,i) = num%clmmap(:,:,i)+ g%OVN%hd; ! depression storage depth (m)

            !j = i + 1
            call tempStorage('EvapSoil',temphead,1)
            i = i + 1; ! i = 51
            num%clmmap(:,:,i) = num%clmmap(:,:,i)+ temphead; ! soil evapotration (m)

            !j = i + 1
            call tempStorage('irrigation',temphead,1)  ! irrigation (m)
            i = i + 1; ! i = 52
            num%clmmap(:,:,i) = num%clmmap(:,:,i)+ temphead;

            !j = i + 1
            i = i + 1; ! i = 53
            num%clmmap(:,:,i) = num%clmmap(:,:,i)+ g%GW(1)%DR;  ! groundwater lateral inflow (m/day). looks like repetitive with 45

            i = i + 1  ! i = 54-58 : save five layers of ice content
            nni = 5
            num%clmmap(:,:,i:i+nni-1) = num%clmmap(:,:,i:i+nni-1) + g%VDZ%ICE(:,:,2:nni+1); i = i + nni - 1 !


            ! to calculate temporal variation:
            ! we do sth adhoc here, assuming num%clmmap(:,:,end) stores the temporal mean after simulation
            ! we compute again to get the variance. If num%clmmap(:,:,end) = 0 this value has no meaning
            ! now just sum up the squares
            !IF (w%t - w%ModelStartTime > 2*ndays .and. varianceMap) THEN ! TEMPORAL variance measured after a longer time
            i = i + 1; ! i = 59
            nzz = size(num%clmmap,3)   ! GW temporal sum of squares accumulation -- need to load previous mean into clmmap(:,:,nzz)
            num%clmmap(:,:,i) = num%clmmap(:,:,i)+ (g%GW(1)%h - num%clmmap(:,:,nzz) )**2D0

            i = i + 1; ! i = 60 - 81  ! pressure head in all layers
            !num%clmmap(:,:,i:i+nz-1) = num%clmmap(:,:,i:i+nz-1)+ (g%VDZ%THE - num%clmmap(:,:,nzz - nz : nzz -1) )**2D0;  i = i+nz-1
            num%clmmap(:,:,i:i+nz-1) = num%clmmap(:,:,i:i+nz-1)+ g%VDZ%h;  i = i+nz-1

            i = i + 1; ! i = 82
            ptrp => temperature_vars%t_grnd_col ! temporal sum of squares of temperature
            ptrc => oneD_ptr(nr*nc,tempdata)
            ptrc = ptrp
            call tAVG(tempdata,tempdata2,'dailyAVGt_grnd',24) ! XY: may need change too.
            !num%clmmap(:,:,i) = num%clmmap(:,:,i) + tempdata2
            num%clmmap(:,:,i) = num%clmmap(:,:,i)+ (tempData - tempdata2 )**2D0;

            i = i + 1  ! i = 83
            num%clmmap(:,:,i) = num%clmmap(:,:,i)+ (tempData - num%clmmap(:,:,nzz - nz -1) )**2D0;
            !ENDIF

            i = i + 1  ! i = 84
            call tempStorage('rchg2',temphead,1)   !
            call depthRecharge(temphead)
            num%clmmap(:,:,i) = num%clmmap(:,:,i)+ temphead; ! accumulated m/hr

            i = i + 1  ! i = 85
            call tempStorage('satWetland0',temphead,1)  ! wetland saturation time
            call satWetland(temphead,0D0)
            num%clmmap(:,:,i) = num%clmmap(:,:,i)+ temphead; ! accumulated m/hr

            i = i + 1  ! i = 86
            call tempStorage('satWetland30',temphead,1) ! wetland saturation time with GW 30cm below definition
            call satWetland(temphead,0.3D0)
            num%clmmap(:,:,i) = num%clmmap(:,:,i)+ temphead; ! accumulated m/hr
#else
            i = 1;
            p2d => twoD_ptr(nr,nc,clm3%g%l%c%ccf%pcf_a%npp)
            num%clmmap(:,:,i)= num%clmmap(:,:,i)+p2d;      i = i + 1;  ! column average npp (gC/m2/s)
            p2d => twoD_ptr(nr,nc,clm3%g%l%c%ccf%pcf_a%gpp)
            num%clmmap(:,:,i)= num%clmmap(:,:,i)+p2d;      i = i + 1;  ! column average gpp
            num%clmmap(:,:,i)= num%clmmap(:,:,i)+g%Veg%TP; i = i + 1;   ! column average TP
            num%clmmap(:,:,i)= num%clmmap(:,:,i)+g%OVN%Ex; i = i + 1;   ! column average Exfiltration

            ! i = 5
            ptrp => clm3%g%l%c%p%pef%eflx_lh_tot ! latent heat flux (W/m**2)
            ptrc => oneD_ptr(nr*nc,tempdata)
            ptrc = 0D0
            call p2c (filter(1)%num_nolakec, filter(1)%nolakec, ptrp, ptrc) ! no lake either
            num%clmmap(:,:,i) = num%clmmap(:,:,i)+ tempdata;

            i = i + 1; ! i = 6
            ptrp => clm3%g%l%c%p%pef%eflx_sh_tot ! sensible heat flux (W/m**2)
            ptrc => oneD_ptr(nr*nc,tempdata)
            ptrc = 0D0
            call p2c (filter(1)%num_nolakec, filter(1)%nolakec, ptrp, ptrc) ! no lake either
            num%clmmap(:,:,i) = num%clmmap(:,:,i)+ tempdata

            i = i + 1; ! i = 7
            ptrp => clm3%g%l%c%p%pef%parsun !average absorbed PAR for sunlit leaves (W/m**2)
            ptrc => oneD_ptr(nr*nc,tempdata)
            ptrc = 0D0
            call p2c (filter(1)%num_nolakec, filter(1)%nolakec, ptrp, ptrc) ! no lake either
            num%clmmap(:,:,i) = num%clmmap(:,:,i)+ tempdata

            i = i + 1; ! i = 8
            ptrp => clm3%g%l%c%p%pef%parsha !average absorbed PAR for shaded leaves (W/m**2)
            ptrc => oneD_ptr(nr*nc,tempdata)
            ptrc = 0D0
            call p2c (filter(1)%num_nolakec, filter(1)%nolakec, ptrp, ptrc) ! no lake either
            num%clmmap(:,:,i) = num%clmmap(:,:,i)+ tempdata

            i = i + 1; ! i = 9
            ptrp => clm3%g%l%c%p%pcf%fpsn   !photosynthesis (umol CO2 /m**2 /s)
            ptrc => oneD_ptr(nr*nc,tempdata)
            ptrc = 0D0
            call p2c (filter(1)%num_nolakec, filter(1)%nolakec, ptrp, ptrc) ! no lake either
            num%clmmap(:,:,i) = num%clmmap(:,:,i)+ tempdata

            !j = i + 1
            i = i + 1; ! i = 10
            ptrp => clm3%g%l%c%ces%t_grnd !ground temperature (Kelvin)
            ptrc => oneD_ptr(nr*nc,tempdata)
            ptrc = ptrp
            num%clmmap(:,:,i) = num%clmmap(:,:,i)+ tempdata

            !j = i + 1
            i = i + 1; ! i = 11
            ptrp => clm3%g%l%c%ces%t_soi_10cm !ground temperature (Kelvin)
            ptrc => oneD_ptr(nr*nc,tempdata)
            ptrc = ptrp
            num%clmmap(:,:,i) = num%clmmap(:,:,i)+ tempdata - 273.15D0

            !j = i + 1
            i = i + 1; ! i = 12
            ptrp => clm3%g%l%c%cws%h2osoi_liqice_10cm !liquid water + ice lens in top 10cm of soil (kg/m2)
            ptrc => oneD_ptr(nr*nc,tempdata)
            ptrc = ptrp
            num%clmmap(:,:,i) = num%clmmap(:,:,i)+ tempdata

            !j = i + 1
            i = i + 1; ! i = 13
            ptrp => clm3%g%l%c%cws%h2osno ! snow water (mm H2O)
            ptrc => oneD_ptr(nr*nc,tempdata)
            ptrc = ptrp
            num%clmmap(:,:,i) = num%clmmap(:,:,i)+ tempdata

            !j = i + 1
            i = i + 1; ! i = 14
            ptrp => clm3%g%l%c%ccf%er ! (gC/m2/s) total ecosystem respiration, autotrophic + heterotrophic
            ptrc => oneD_ptr(nr*nc,tempdata)
            ptrc = ptrp
            num%clmmap(:,:,i) = num%clmmap(:,:,i)+ tempdata

            !j = i + 1
            i = i + 1; ! i = 15
            ptrp => clm3%g%l%c%ccf%hr ! (gC/m2/s) total heterotrophic respiration
            ptrc => oneD_ptr(nr*nc,tempdata)
            ptrc = ptrp
            num%clmmap(:,:,i) = num%clmmap(:,:,i)+ tempdata

            !j = i + 1
            i = i + 1; ! i = 16 top 10cm soil h2o content [-] including liquid and ice
            call tempStorage('hcAfterLand',temphead,1)
            ptrp => clm3%g%l%c%cws%h2osoi_liqice_10cm
            ptrc => oneD_ptr(nr*nc,tempdata)
            ptrc = ptrp
            temphead1d => oneD_ptr(nr*nc,temphead)
            num%clmmap(:,:,i) = num%clmmap(:,:,i)+ (tempdata/1000D0 - temphead)/0.1D0

            i = i + 1; ! i = 17
            ptrp => clm3%g%l%c%p%pps%tlai ! total leaf area index
            ptrc => oneD_ptr(nr*nc,tempdata)
            ptrc = 0D0
            call p2c (filter(1)%num_nolakec, filter(1)%nolakec, ptrp, ptrc) ! no lake either
            num%clmmap(:,:,i) = num%clmmap(:,:,i)+ tempdata
            if (any( isnan(ptrp) )) then
                print *, 'LAI nan found'
            endif

            i = i + 1; ! i = 18
            ptrp => clm3%g%l%c%p%pps%tsai ! total stem area index
            ptrc => oneD_ptr(nr*nc,tempdata)
            ptrc = 0D0
            call p2c (filter(1)%num_nolakec, filter(1)%nolakec, ptrp, ptrc) ! no lake either
            num%clmmap(:,:,i) = num%clmmap(:,:,i)+ tempdata

            i = i + 1; ! i = 19 - 40
            num%clmmap(:,:,i:i+nz-1) = num%clmmap(:,:,i:i+nz-1) + g%VDZ%THE; i = i + nz - 1 ! soil moisture layer-by-layer, set to end of layer, if nz=22,

            i = i + 1  ! i = 41 - 43
            p3d => threeD_ptr((/nr,nc,3/),clm3%g%l%c%p%pcf%npp)  ! pft-level npp
            num%clmmap(:,:,i:i+3-1)= num%clmmap(:,:,i:i+3-1)+p3d; i = i + 3 - 1

            i = i + 1; ! i = 44
            p2d => twoD_ptr(nr,nc,clm3%g%l%c%ccf%nee) ! net ecosystem exchange
            num%clmmap(:,:,i) = num%clmmap(:,:,i) + p2d;

            i = i + 1; ! i = 45
            num%clmmap(:,:,i) = num%clmmap(:,:,i) + g%GW(1)%DR; ! lateral inflow (m/day)

            i = i + 1; ! i = 46
            ptrp => clm3%g%l%c%p%pps%elai  ! effective leaf area index
            ptrc => oneD_ptr(nr*nc,tempdata)
            ptrc = 0D0
            call p2c (filter(1)%num_nolakec, filter(1)%nolakec, ptrp, ptrc) ! no lake either
            num%clmmap(:,:,i) = num%clmmap(:,:,i)+ tempdata

            i = i + 1; ! i = 47
            ptrp => clm3%g%l%c%p%pps%btran  !transpiration wetness factor (0 to 1)
            ptrc => oneD_ptr(nr*nc,tempdata)
            ptrc = 0D0
            call p2c (filter(1)%num_nolakec, filter(1)%nolakec, ptrp, ptrc) ! no lake either
            num%clmmap(:,:,i) = num%clmmap(:,:,i)+ tempdata

            !j = i + 1
            i = i + 1; ! i = 48
            p2d => twoD_ptr(nr,nc,clm3%g%l%c%ccf%pcf_a%ar) ! (gC/m2/s) autotrophic respiration (MR + GR)
            num%clmmap(:,:,i) = num%clmmap(:,:,i)+ p2d;

            !j = i + 1
            i = i + 1; ! i = 49
            num%clmmap(:,:,i) = num%clmmap(:,:,i)+ g%VDZ%R(:,:,2); ! saturation state (==1 for saturated, ==0 for unsaturated). therefore, this field * dt =  total time saturated

            !j = i + 1
            i = i + 1; ! i = 50
            num%clmmap(:,:,i) = num%clmmap(:,:,i)+ g%OVN%hd; ! depression storage depth (m)

            !j = i + 1
            call tempStorage('EvapSoil',temphead,1)
            i = i + 1; ! i = 51
            num%clmmap(:,:,i) = num%clmmap(:,:,i)+ temphead; ! soil evapotration (m)

            !j = i + 1
            call tempStorage('irrigation',temphead,1)  ! irrigation (m)
            i = i + 1; ! i = 52
            num%clmmap(:,:,i) = num%clmmap(:,:,i)+ temphead;

            !j = i + 1
            i = i + 1; ! i = 53
            num%clmmap(:,:,i) = num%clmmap(:,:,i)+ g%GW(1)%DR;  ! groundwater lateral inflow (m/day). looks like repetitive with 45

            i = i + 1  ! i = 54-58 : save five layers of ice content
            nni = 5
            num%clmmap(:,:,i:i+nni-1) = num%clmmap(:,:,i:i+nni-1) + g%VDZ%ICE(:,:,2:nni+1); i = i + nni - 1 !


            ! to calculate temporal variation:
            ! we do sth adhoc here, assuming num%clmmap(:,:,end) stores the temporal mean after simulation
            ! we compute again to get the variance. If num%clmmap(:,:,end) = 0 this value has no meaning
            ! now just sum up the squares
            !IF (w%t - w%ModelStartTime > 2*ndays .and. varianceMap) THEN ! TEMPORAL variance measured after a longer time
            i = i + 1; ! i = 59
            nzz = size(num%clmmap,3)   ! GW temporal sum of squares accumulation -- need to load previous mean into clmmap(:,:,nzz)
            num%clmmap(:,:,i) = num%clmmap(:,:,i)+ (g%GW(1)%h - num%clmmap(:,:,nzz) )**2D0

            i = i + 1; ! i = 60 - 81  ! pressure head in all layers
            !num%clmmap(:,:,i:i+nz-1) = num%clmmap(:,:,i:i+nz-1)+ (g%VDZ%THE - num%clmmap(:,:,nzz - nz : nzz -1) )**2D0;  i = i+nz-1
            num%clmmap(:,:,i:i+nz-1) = num%clmmap(:,:,i:i+nz-1)+ g%VDZ%h;  i = i+nz-1

            i = i + 1; ! i = 82
            ptrp => clm3%g%l%c%ces%t_grnd ! temporal sum of squares of temperature
            ptrc => oneD_ptr(nr*nc,tempdata)
            ptrc = ptrp
            call tAVG(tempdata,tempdata2,'dailyAVGt_grnd',24)
            !num%clmmap(:,:,i) = num%clmmap(:,:,i) + tempdata2
            num%clmmap(:,:,i) = num%clmmap(:,:,i)+ (tempData - tempdata2 )**2D0;

            i = i + 1  ! i = 83
            num%clmmap(:,:,i) = num%clmmap(:,:,i)+ (tempData - num%clmmap(:,:,nzz - nz -1) )**2D0;
            !ENDIF

            i = i + 1  ! i = 84
            call tempStorage('rchg2',temphead,1)   !
            call depthRecharge(temphead)
            num%clmmap(:,:,i) = num%clmmap(:,:,i)+ temphead; ! accumulated m/hr

            i = i + 1  ! i = 85
            call tempStorage('satWetland0',temphead,1)  ! wetland saturation time
            call satWetland(temphead,0D0)
            num%clmmap(:,:,i) = num%clmmap(:,:,i)+ temphead; ! accumulated m/hr

            i = i + 1  ! i = 86
            call tempStorage('satWetland30',temphead,1) ! wetland saturation time with GW 30cm below definition
            call satWetland(temphead,0.3D0)
            num%clmmap(:,:,i) = num%clmmap(:,:,i)+ temphead; ! accumulated m/hr
#endif
#endif
        ENDIF ! time IF

    ELSEIF (cmd .eq. 10) THEN ! wiping operation
        num%Prcp   = 0D0
        num%Inf    = 0D0
        num%PET    = 0D0
        num%ET     = 0D0
        num%EvapG  = 0D0
        num%Dperc  = 0D0
        num%Dperc2 = 0D0
        num%SatE   = 0D0
        num%InfE   = 0D0
        num%botT   = 0D0
        num%SUBM   = 0D0
        num%SMELT  = 0D0
        num%GW     = 0D0
        num%ovnH   = 0D0
        num%Qgc    = 0D0
        !output = 0.0D0
        num%clmmap = 0D0
        num%n      = 0
        do i=1,size(w%CR)
            rid = w%CR(i)
            r(rid)%Riv%vrv(16:,:)=0D0
        enddo

    ENDIF

    end subroutine mapOutput2


    subroutine timeSeriesMapControl()
    ! this subroutine checks if a file named TimeSeriesMap.dat exists
    ! if it does, we try to read the time control points in the file
    ! at every control point, write out the mapOutput and clean mapOutput2 memory
    ! time stamp is daily. called after syncTime
    use vdata
    USE clmtype
    use memoryChunk, Only: tAVG, tempStorage
    integer i,eof
    integer, save:: k = 0, uid=0, m = 1, j = 0, dsig = 0, nr,nc,nz,nv
    real*8,dimension(:),allocatable, save :: times
    real*8 d
    character*280 :: fname
    character*30 :: varn
    real*8, pointer :: ptrp(:)         ! pointer to input pft array
    real*8, pointer :: ptrc(:),temphead1d(:),temphead(:,:)         ! pointer to output column array
    real*8 tempdata(size(g%VDZ%h,1), size(g%VDZ%h,2))
    real*8 tempdata2(size(g%VDZ%h,1), size(g%VDZ%h,2))
    integer:: tsMapOption = 0

    if (size(w%control)>2) then
        if (w%control(3)>0) then
            tsMapOption = int(w%control(3))
        endif
    endif

    if (k .eq. 0) then
        OPEN(unit = uid, file = 'TimeSeriesMap.dat', ACTION='READ',IOSTAT = eof )
        do while (eof .eq. 0)
            read(uid,*,iostat = eof)
            if (eof .eq. 0) j = j + 1
        enddo
        if (j>0) then
            allocate(times(j))
            rewind(uid)
            do i=1,j
                read(uid,*,iostat = eof) times(i)
            enddo
        endif
        close(uid)
        do i=1,j
            if (times(i)>=w%t-1D-13) then
                m = i
                EXIT
            endif
        enddo
    endif

    nr = size(g%VDZ%h,1)
    nc = size(g%VDZ%h,2)
    nz = size(g%VDZ%h,3)
    nv = size(g%Veg%RPT,3)

    k = k + 1
    if (j > 0 .and. m <= j) then
        if (abs(w%t - times(m))<w%dt/10D0) then
            ! the first one is normally the starting point
            if (tsMapOption .eq. 0) then
                call saveW([trim('tsMap_')//trim(num2str(m-1))//'.mat']) !
            elseif (tsMapOption .eq. 1) then
                call saveW([trim('tsMap_')//trim(num2str(m-1))//'.mat'],'r') ! an optional input
            else
                dsig = 0; call saveVar('the',size(g%VDZ%the),g%VDZ%the,dsig,m)
                call tempStorage('hcAfterLand',temphead,1)
#if (defined CLM45)
                ptrp => waterstate_vars%h2osoi_liqice_10cm_col
#else
                ptrp => clm3%g%l%c%cws%h2osoi_liqice_10cm
#endif             
                ptrc => oneD_ptr(nr*nc,tempdata)
                ptrc = ptrp
                temphead1d => oneD_ptr(nr*nc,temphead)
                tempData2 = (tempdata/1000D0 - temphead)/0.1D0
                dsig = 1; call saveVar('the10',size(tempData2),tempData2,dsig,m)
                dsig = 1; call saveVar('ET',size(tempData2),w%tData%maps(1)%ET,dsig,m)
                dsig = 1; call saveVar('LAI',size(tempData2),w%tData%maps(1)%clmmap(:,:,17),dsig,m)
                dsig = 2; call saveVar('prcp',size(tempData2),w%tData%maps(1)%prcp,dsig,m)
            endif

            call mapOutput2(10,w%tData%maps(1)) !clear
            m = m + 1
        endif
    endif


    end subroutine



    subroutine depthRecharge(Rchg)
    ! return recharge option #2. same unit as K
    use vdata
    real*8 Rchg(size(g%VDZ%h,1), size(g%VDZ%h,2))
    ! logic: take the minimum of (WTL - 2) and 5m depth
    ! take the maximum of the above result and 2. Calculate the flux going through the bottom of this layer
    integer i,j,k,km
    integer, save :: status = 0
    integer, allocatable, save :: k5m(:,:)
    real*8 KZ,dt,dh
    real*8, pointer:: h(:),kk(:),the(:)
    integer, pointer:: wtl(:,:)

    if (status .eq. 0) then
        allocate(k5m(g%VDZ%m(1),g%VDZ%m(2)) )
        DO I = 1, g%VDZ%m(1)
            DO J = 1,g%VDZ%m(2)
                IF (w%DM%mask(i,j)) then
                    k5m(i,j) = g%VDZ%m(3)-2
                    DO k = 1,g%VDZ%m(3)
                        IF (g%VDZ%EB(I,J,1)-g%VDZ%EB(I,J,k) >= 5D0) then
                            k5m(i,j) = k
                            EXIT
                        endif
                    ENDDO
                endif
            enddo
        enddo
        status = 1
    endif

    Rchg = 0D0
    wtl => g%VDZ%WTL
    DO I = 1, g%DM%msize(1)
        DO J = 1,g%DM%msize(2)
            IF (w%DM%mask(i,j)) then
                h => g%VDZ%h(i,j,:)
                the => g%VDZ%the(i,j,:)
                kk => g%VDZ%K(i,j,:)
                km = min(WTL(I,J)-2, k5m(i,j))
                km = max(km, 2)
                if (g%VDZ%h(i,j,km)>0 .AND. g%VDZ%h(i,j,km+1)>0) then
                    KZ = (g%VDZ%KS(i,j,km)+g%VDZ%KS(i,j,km+1))/2D0
                else
                    KZ = (g%VDZ%K(i,j,km)+g%VDZ%K(i,j,km+1))/2D0
                    !dt = g%VDZ%dtP(3)
                endif
                dh = (g%VDZ%h(i,j,km) - g%VDZ%h(i,j,km+1))/g%VDZ%DDZ(i,j,km)
                dh = min(max(dh,-5D0),5D0)
                rchg(i,j) = (dH + 1D0) *KZ*g%VDZ%dtP(3)
            endif
        enddo
    enddo
    end subroutine

    subroutine satWetland(satTime,d,jdRanges)
    ! decide the amount of time wetland is saturated.
    ! wetland is decided saturated when groundwater table is within d distance from the surface
    use vdata
    use memoryChunk, Only: tAVG, tempStorage
    !real*8 Rchg(size(g%VDZ%h,1), size(g%VDZ%h,2))
    real*8, dimension(size(g%VDZ%h,1), size(g%VDZ%h,2)) :: satTime
    real*8 d
    real*8, dimension(:,:),pointer:: dormant
    integer jDRanges(2)
    integer i,j,k,l
    call tempStorage('dormant',dormant,1)

    if (any(dormant>0.5D0)) then
        satTime = 0D0
        return
    endif
    satTime(1,1) = w%dt ! I know this cell is not in the mask. Use it to save total growing season time
    DO I = 1, g%DM%msize(1)
        DO J = 1,g%DM%msize(2)
            if (w%DM%mask(i,j)) then
                if (g%GW(1)%h(i,j)>= (g%OVN%dP(i,j,11)-d)) then
                    satTime(i,j) = w%dt
                else
                    satTime(i,j) = 0D0
                endif
            else
                satTime(i,j) = 0D0
            endif
        ENDDO
    ENDDO

    end subroutine