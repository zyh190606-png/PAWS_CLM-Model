#include <misc.h>
#include <preproc.h>
    ! Bridging functionalities between PAWS and CLM, this is the only place where both PAWS data structures and CLM data structures are loaded
    ! currently there is an exception in mkarbinitMod, where VDATA is used. consider removing that later
    MODULE PAWS_CLM_MOD
#if (defined CLM)
    USE clmtype
    USE clmtypeInitMod
    use clm_varpar  , only : maxpatch_pft, nlevsno, nlevgrnd, numrad, nlevlak, &
        numpft, ndst, nvoc, nlevurb, nlevsoi
    USE VDATA
    USE memoryChunk, only: mn,nn
    PUBLIC :: PAWS_CLM, initialSTemp,passV,SURFACE_URB

    PRIVATE
    LOGICAL, PUBLIC:: WRITEOUT
    INTEGER obsc(1000),nobsc,obsp(1000),nobsp
    integer, dimension(:),pointer:: filterDBG
    integer K
    !#######################################################################################
    CONTAINS

    subroutine PAWS_CLM(dir,state,H)
    ! main gateway subroutine for the linkage between PAWS and CLM
    ! These tasks are being performed in this subroutine
    ! I: At initialization:
    ! 1: Memory Allocation (Assignment of pointer to central memory storage)
    ! 2: CLM interna physical Initialization (Assign physically meaningful default values to CLM variables) including
    !    --> 2.0. set signaling NaN for some strctures (cps,pps,pepv)
    !    --> 2.1. obtain site information (from PAWS): Geography, Topography, Time
    !    --> 2.2. gridcell,landunit,column,pft relationships, set up filters
    !    --> 2.3. read parameter files (pftcon,urban,crop,etc)
    !    --> 2.4. default state/parameter values (from mkarinit, iniTimeConst,urbanInit)
    ! 3. Use PAWS data (states) to evaluate some CLM data (one-time) (use Pointbacks -- pointers in CLM that points to PAWS existing data)
    !    --> 3.2. Soil water retention parameters
    !    --> 3.3. Initial states: Carbon, Nitrogen (In place of EcosystemDynini)
    ! II. At Run-time:
    ! 0. If in Matlab-mex mode, pointers are assigned to central memory storage
    ! 1. Land-based observations are cast into CLM forcing
    ! 2. Water state updates coming from PAWS
    ! 3. CLM computed fluxes/states passing to PAWS
    ! LOOK AT subroutine LAND in prism.F90 for state changes -- that is more up to date than the description here!!!!!!!!!!!!!!!!!!
    ! state = 1, during simulation.
    ! state = 2, cold initialization. memory needs to be allocated and initialized.
    ! state = 3, restart run, in this run, model is resumed from a previous saved state.
    ! state = 3, reading a mat file, but values need to be initiated (only that no allocation of memory is needed), after initialization, state should be set to 1
    !       memory chunks are load and structure pointers are assigned.

    ! dir = 1, run-time linking from PAWS to CLM  -- will do passV(1) and driver
    ! dir = 2, run-time linking from CLM to PAWS  -- will do passV(2) only
    ! dir = 3, run-time linking from PAWS to CLM, with pointer assignment (Matlab Mex mode) -- will do linkCLMtype, passV(1) and driver
    ! dir = 4, finish-up procedures, (steps done after soil moisture is updated by PAWS), including EcoSystemDyn and MassBalanceCheck  -- will do driverPart2, ecosystem
    ! support for dir =3 may not work currently
    USE memoryChunk, Only: R3,R4,L3,L4,i3,i4,j3,j4,mm,VBB,VBCD
    USE filterMod, Only : filter
    USE abortutils, Only : endrun
    use clm_varctl, only : fpftcon
    implicit none
    integer H
    integer dir,state
    integer nc,m(4)
    integer,dimension(2):: mv2
    integer T4(4)
    type(gridcell_type),pointer::gg
    type(GlobalG),pointer::gp
    real*8,save:: savett
    logical :: dbgsave = .false.

    gp      => g
    gg => clm3%g
    if (state .ne. 1) then
        ! Initialization
        m = (/g%VDZ%m(1),g%VDZ%m(2),g%VDZ%m(3),g%Veg%nt/)
        mm = m                               ! in memoryChunk
        mn = m(1)*m(2); nn = m(1)*m(2)*m(3)  ! in memoryChunk
        i3 =1200; i4 = 550  ! CHaopeng, counted 10/23/10: i3=1200;i4=476; j3=47;j4=7. with mm=[64 70 22 3], this is about 100 MB memory
        j3 =  50; j4 =  20
        if (state .eq. 0) then
            ALLOCATE(R3(m(1),m(2),i3))
            ALLOCATE(R4(m(1),m(2),m(4),i4))
            ALLOCATE(L3(m(1),m(2),j3))
            ALLOCATE(L4(m(1),m(2),m(4),j4))
            R3 = 0D0; R4 = 0D0;
            L3 = 0D0; L4 = 0D0;
        else
            ! already linked in from mat reading
        endif

        !ALLOCATE(VBB(17,60))  ![ntype nproperty], I have not counted the number of properties. This is in memoryChunk for pftcon
        ! For crop.dat
        mv2 = (/200,38/) ! haven't counted how many crops are here
        ALLOCATE(VBCD(mv2(1),mv2(2)))

        call linkCLMtype                                         ! I.1 & II.0
        ! This is to create some mutual pointing of data structures with shared memory
        ! It handles the input from PAWS to CLM including:
        ! frac,lai,itype,htop
        call initCLM(m,state)                                          ! I.2
        if (state .eq. 3 .and. clm3%g%l%c%cws%h2osoi_liq(filter(1)%soilc(1),1) .eq. 0D0) then
            call ch_assert('PAWS_CLM::Renew job but initialized as 0 water content--Forgot to use correct *_CLM.mat??')
        endif

        state = 1

        !!!!!!!!!!! TEMPORARY TEST  !!!!!!!!!!!!!!!!!!!!!
        !     pftcon%mp = pftcon%mp * 5D0
        ! THIS IS TO ADJUST ET!!!!!!!!!!!!!!!!!



    endif

    if (dir .eq. 1 .OR. dir .eq. 3) then
        call clmForcing(g%Wea%day,H,m,g%Veg%SMAP)                ! II.1
        call clm_driver(1,H)       ! passV(2) must be used in initialization
        call passV(1)            ! From CLM to PAWS, then outside, do unsat

    elseif (dir .eq. 2) then
        ! just after doing UNSAT

        call passV(2)                                            ! II.3

        ! from PAWS to CLM
        call clm_driver(2,H)
        ! With the current CLM code, there is no need for further variable passing from CLM to PAWS
        ! part 2 of clm_driver does Hydrology2(part2), SnowAge, CNEcosystemDyn,CNAnnualUpdate,SurfaceAlbedo(for next time step)
        ! which will update vegetation variables and snow variables
        ! LAI,hc is directly linked via memory
        savett = w%t
        call writeCLMrec(nobsc,obsc,15,nobsp,obsp,16)
    else
        call endrun('incorrect direction of simulation')
    endif

    end subroutine

    !###########################################################################################
    subroutine initCLM(m,state)
    ! I: At initialization:
    ! 2: CLM internal physical Initialization (Assign physically meaningful default values to CLM variables)
    ! Emulating Initialize1() and Initialize2() in clm_initializeMod, and init3 in clm_comp
    ! Following their orders of execution
    USE mkarbinitMod, Only: mkarbinit
    USE filterMod, Only: allocFilters,setFilters,filter
    USE pftvarcon, Only: pftconrd
    USE initSurfAlbMod, Only: initSurfAlb
    USE clm_time_manager, Only: get_curr_calday
    USE shr_orb_mod, Only: shr_orb_params, shr_orb_decl
    use clm_varorb, only : eccen, mvelpp, lambm0, obliqr
    USE clm_varpar, Only: clm_varpar_init,MAXPATCH_PFT
    use decompMod, only : get_proc_clumps, get_clump_bounds
    USE CNEcosystemDynMod, Only : CNEcosystemDynInit
    USE accFldsMod,  Only: initAccFlds, initAccClmtype
    implicit none
    integer m(4),T4(4),state
    type(gridcell_type),pointer::gg
    real*8 calday,caldaym1,declin,declinm1
    real*8:: obliq, eccf,mvelp
    integer nc, begg, endg, begl, endl, begc, endc, begp, endp

    gg => clm3%g
    MAXPATCH_PFT = g%Veg%nt

    if (state .ne. 1) then
        !    --> 2.0. set signaling NaN for some strctures (cps,pps,pepv)
        !if (state .eq. 2) call nanFlagging

        !    --> 2.1. obtain site information (from PAWS): Geography, Topography, Fractions, Time
        call clm_varpar_init
        call paws_par_transfer()
        call fileInit
        call mapGeography(g%Wea%Day)
        T4 = set_time_PAWS()
        call bridgeVegPointers(m)    ! Fractions
        call pftconrd()

        !    --> 2.2. gridcell,landunit,column,pft relationships, set up filters
        call buildRelations(m,w%DM%mapL,g%Veg%frac) ! land use classification: derive l%itype, c%itype
        call allocFilters()
        call setFilters(1)
        !call fixFilters()
        call iniTimeConst(state) ! modified from CLM. cps,pps

        call get_clump_bounds(nc, begg, endg, begl, endl, begc, endc, begp, endp)
        call initAccFlds()    ! allocate memory (later make a link (for restart))
        call initAccClmtype() ! evaluate (from maybe restart)

        call CNEcosystemDynInit(begc, endc, begp, endp )
    endif

    if (state .ne. 1 .and. state .ne. 3) then
        ! these operations will not be done during

        !    --> 2.3. read parameter files (pftcon,urban,crop,etc)
        !    --> 2.4. default state/parameter values (from mkarinit, iniTimeConst,urbanInit)
        nc = 1
        call CNiniTimeVar()
        call waterPar       !THES, THER, FC, WP, water retention/hydraulic parameters
        CALL PAWS_Veg_Initial(filter(nc)%num_soilp,filter(nc)%soilp) ! use PAWS LAI to evaluate some carbon state variables. This need to be done after iniTimeConst
        call mkarbinit      ! default arbitrary settings for water state/energy state
        call col_state_init() ! new column water/energy state initializations, more meaningful
        call passV(2)       ! initialize water states from PAWS

        ! init3
        call shr_orb_params( T4(1) , eccen  , obliq , mvelp     ,            &
            &               obliqr   , lambm0 , mvelpp, .TRUE. )
        calday = get_curr_calday()
        call shr_orb_decl( calday, eccen, mvelpp, lambm0, obliqr, declin, eccf )
        caldaym1 = get_curr_calday(offset=-int(w%dt*86400D0))
        call shr_orb_decl( caldaym1, eccen, mvelpp, lambm0, obliqr, declinm1, eccf )
        call initSurfAlb( calday, declin, declinm1 )
    endif
    end subroutine

    !#######################################################################################
    function set_time_PAWS
    USE clm_time_manager, Only: year_PAWS, month_PAWS, day_PAWS, secs_PAWS, JD_PAWS
    implicit none
    REAL*8 tt
    INTEGER C3(3)
    INTEGER d,YJD,darray(3)
    integer set_time_PAWS(4)
    real*8 offset,GT
    ! the solar declination algorithm of CLM (I think) uses Greenwich time over the globe
    ! It is necessary to convert local time to Greenwich Time (GT)

    offset = w%timezone/24D0
    GT = w%t - offset + w%dt ! apparently, these is used as end points in CLM
    C3 = datevec2(int(GT  + 1D-8))
    year_PAWS = C3(1)
    month_PAWS = C3(2)
    day_PAWS = C3(3)
    secs_PAWS = floor((GT - int(GT  + 1D-8))*86400D0 + 0.5D0 - 1D-5 )

    call ymd2jd([year_PAWS, month_PAWS, day_PAWS], 3, d, YJD, darray)
    JD_PAWS = DBLE(d) + GT - int(GT  + 1D-8)
    set_time_PAWS(1:3)=C3
    set_time_PAWS(4)  =secs_PAWS
    end function

    !#######################################################################################

    subroutine clm_driver(dir,H)
    use shr_kind_mod        , only : r8 => shr_kind_r8
    use clm_varctl          , only : wrtdia, fpftdyn
    use clm_varctl          , only : iulog
    use spmdMod             , only : masterproc,mpicom
    use decompMod           , only : get_proc_clumps, get_clump_bounds
    use filterMod           , only : filter, setFilters
#if (defined CNDV)
    use CNDVMod             , only : dv, histCNDV
    use pftdynMod           , only : pftwt_interp
#endif
    !use pftdynMod           , only : pftdyn_interp, pftdyn_wbal_init, pftdyn_wbal
#if (defined CN)
    !use pftdynMod           , only : pftdyn_cnbal
    !use clm_varctl          , only : fndepdyn, use_ndepstream
#endif
    !use dynlandMod          , only : dynland_hwcontent
    use clm_varcon          , only : zlnd, isturb
    use clm_time_manager    , only : get_step_size, get_curr_calday, mcsec_Day, &
        get_curr_date, get_ref_date, get_nstep, is_perpetual
    !use histFileMod         , only : hist_update_hbuf, hist_htapes_wrapup
    !use restFileMod         , only : restFile_write, restFile_write_binary, restFile_filename
    !use inicFileMod         , only : inicfile_perp
    use accFldsMod          , only : updateAccFlds
    use clm_driverInitMod   , only : clm_driverInit
    use BalanceCheckMod     , only : BeginWaterBalance, BalanceCheck
    use SurfaceRadiationMod , only : SurfaceRadiation
    use Hydrology1Mod       , only : Hydrology1
    use Hydrology2Mod       , only : Hydrology2
    use HydrologyLakeMod    , only : HydrologyLake
    use Biogeophysics1Mod   , only : Biogeophysics1
    use BareGroundFluxesMod , only : BareGroundFluxes
    use CanopyFluxesMod     , only : CanopyFluxes
    use Biogeophysics2Mod   , only : Biogeophysics2
    use BiogeophysicsLakeMod, only : BiogeophysicsLake
    use SurfaceAlbedoMod    , only : SurfaceAlbedo
    use pft2colMod          , only : pft2col
#if (defined CN)
    use CNSetValueMod       , only : CNZeroFluxes_dwt
    use CNEcosystemDynMod   , only : CNEcosystemDyn
    use CNAnnualUpdateMod   , only : CNAnnualUpdate
    use CNBalanceCheckMod   , only : BeginCBalance, BeginNBalance, &
        CBalanceCheck, NBalanceCheck
    !use ndepFileMod         , only : ndepdyn_interp
    !use ndepStreamMod       , only : ndep_interp, &
    !                           stream_year_first_ndep, stream_year_last_ndep
#else
    use STATICEcosysDynMod  , only : EcosystemDyn
#endif
#if (defined DUST)
    use DUSTMod             , only : DustDryDep, DustEmission
#endif
    !use VOCEmissionMod      , only : VOCEmission
    !use seq_drydep_mod      , only : n_drydep, drydep_method, DD_XLND
    !use STATICEcosysDynMod  , only : interpMonthlyVeg
    !use DryDepVelocity      , only : depvel_compute
#if (defined CASA)
    use CASAMod             , only : casa_ecosystemDyn
#endif
#if (defined RTM)
    use RtmMod              , only : Rtmriverflux
#endif
    use abortutils          , only : endrun
    !use UrbanMod            , only : UrbanAlbedo, UrbanRadiation, UrbanFluxes
    !use perf_mod
    use SNICARMod           , only : SnowAge_grain
    !use aerdepMod           , only : interpMonthlyAerdep
    !use clm_varctl          , only : set_caerdep_from_file, set_dustdep_from_file
    use clm_varorb, only : eccen, mvelpp, lambm0, obliqr
    USE shr_orb_mod, Only: shr_orb_cosz, shr_orb_params, shr_orb_decl
    USE shr_const_mod, Only: SHR_CONST_PI
    USE clm_atmlnd

    implicit none
    integer dir,H
    integer begg, endg, begl, endl, begc, endc, begp, endp
    logical,save:: doalb
    type(column_type), pointer ::cptr
    type(gridcell_type),pointer::gg
    real*8,save :: nextsw_cday, declinp1, calday, declin, eccf,cosz,d2pi,tt
    integer T4(4),nc
    integer I,C,L
    integer, save:: Kalb=0, Light=0,Light_Old=0
    type(Sta_type),pointer::day
    logical :: dbgsave = .false.

    gg => clm3%g
    cptr => gg%l%c
    day  => g%wea%day
    nc = 1
    d2pi = SHR_CONST_PI/180D0
    call get_clump_bounds(nc, begg, endg, begl, endl, begc, endc, begp, endp)

    if (dir .eq. 1) then
        ! 11111111111111111111111111111111111111111111111111111111111111
        ! some filters changed from original driver callings
        T4 = set_time_PAWS()
        calday = get_curr_calday()
        nextsw_cday = get_curr_calday(int(w%dt*86400))
        call shr_orb_decl( calday     , eccen, mvelpp, lambm0, obliqr, declin  , eccf )
        call shr_orb_decl( nextsw_cday, eccen, mvelpp, lambm0, obliqr, declinp1, eccf )
        ! decide if doalb:
        ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! THIS IS THE ONLY PIECE IN THIS SUBROUTINE THAT USES IMPLICIT DATA FROM PAWS
        Light_Old = Light
        Light = 0
        DO L = 1, size(day%lat,1)
            cosz = shr_orb_cosz (nextsw_cday, day%lat(L)*d2pi, day%lon(L)*d2pi, declinp1) ! what is saved in day is radians
            if (cosz>0D0) then
                light = 1; EXIT
            endif
        ENDDO
        IF (light .eq. 1) then
            !Kalb = Kalb + 1
            !if (mod(Kalb,2) .eq. 1) then
            ! currently, the station method may be off as many stations are quite far apart, just use albedo wheneve there is light next time step
            doalb = .true.
            !else
            !  doalb = .false.
            !endif
            if (Light_Old .eq. 0) then ! first sight of light
                tt = w%t - w%dt
                mcsec_Day = T4(4)
            else
                mcsec_Day = -999
            endif
        else
            if (Light_Old .eq. 1) then ! so the darkness begins...
                Kalb = 0
                doalb = .true. ! do last time alb today (to clear out matrices)
            else
                doalb = .false. ! darkness continues ...
            endif
        ENDIF
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if (defined CLM_Only)
        call BeginWaterBalance(begc, endc, begp, endp, &
            filter(nc)%num_nolakec, filter(nc)%nolakec, filter(nc)%num_lakec, filter(nc)%lakec, &
            filter(nc)%num_hydrologyc, filter(nc)%hydrologyc)
#endif          

        call BeginCBalance(begc, endc, filter(nc)%num_soilc, filter(nc)%soilc)
        call BeginNBalance(begc, endc, filter(nc)%num_soilc, filter(nc)%soilc)


        call clm_driverInit(nc)

        call Hydrology1(begc, endc, begp, endp, &
            filter(nc)%num_nolakec, filter(nc)%nolakec, &
            filter(nc)%num_nolakep, filter(nc)%nolakep)



        call SurfaceRadiation(begp, endp, &
            filter(nc)%num_nourbanp, filter(nc)%nourbanp)
#if (defined Urban)                     
        call UrbanRadiation(nc, begl, endl, begc, endc, begp, endp, &
            filter(nc)%num_nourbanl, filter(nc)%nourbanl, &
            filter(nc)%num_urbanl, filter(nc)%urbanl, &
            filter(nc)%num_urbanc, filter(nc)%urbanc, &
            filter(nc)%num_urbanp, filter(nc)%urbanp)
#endif
        call Biogeophysics1(begg, endg, begc, endc, begp, endp, &
            filter(nc)%num_nolakec, filter(nc)%nolakec, &
            filter(nc)%num_nolakep, filter(nc)%nolakep)

        !! all PFTs that is num_nolakep and frac_veg_nosno(p) == 0 is computed inside
        call BareGroundFluxes(begp, endp, &
            filter(nc)%num_nolakeurbanp, filter(nc)%nolakeurbanp)

        !! Urban water fluxes are calculated inside BareGroundFluxes but are zeroed out here. I wonder if I should actually keep the energy fluxes
        call fixUrban(filter(nc)%num_nolakep,filter(nc)%nolakep)
#if (defined Urban)  
        call UrbanFluxes(nc, begp, endp, begl, endl, begc, endc, &
            filter(nc)%num_nourbanl, filter(nc)%nourbanl, &
            filter(nc)%num_urbanl, filter(nc)%urbanl, &
            filter(nc)%num_urbanc, filter(nc)%urbanc, &
            filter(nc)%num_urbanp, filter(nc)%urbanp)
#endif


        call CanopyFluxes(begg, endg, begc, endc, begp, endp, &
            filter(nc)%num_nolakep, filter(nc)%nolakep)

#if (defined Lake)                 
        call BiogeophysicsLake(begc, endc, begp, endp, &
            filter(nc)%num_lakec, filter(nc)%lakec, &
            filter(nc)%num_lakep, filter(nc)%lakep)
#endif

#if (defined DUST)
        ! Dust mobilization (C. Zender's modified codes)
        call DustEmission(begp, endp, begc, endc, begl, endl, &
            filter(nc)%num_nolakep, filter(nc)%nolakep)

        ! Dust dry deposition (C. Zender's modified codes)
        call DustDryDep(begp, endp)
#endif

        call Biogeophysics2(begl, endl, begc, endc, begp, endp, &
            filter(nc)%num_urbanl,  filter(nc)%urbanl, &
            filter(nc)%num_nolakec, filter(nc)%nolakec, &
            filter(nc)%num_nolakep, filter(nc)%nolakep)



        call pft2col(begc, endc, filter(nc)%num_nolakec, filter(nc)%nolakec)

        call Hydrology2(begc, endc, begp, endp, &
            filter(nc)%num_nolakec, filter(nc)%nolakec, &
            filter(nc)%num_hydrologyc, filter(nc)%hydrologyc, &
            filter(nc)%num_urbanc, filter(nc)%urbanc, &
            filter(nc)%num_snowc, filter(nc)%snowc, &
            filter(nc)%num_nosnowc, filter(nc)%nosnowc,1)


    elseif (dir .eq. 2) then

        call Hydrology2(begc, endc, begp, endp, &
            filter(nc)%num_nolakec, filter(nc)%nolakec, &
            filter(nc)%num_hydrologyc, filter(nc)%hydrologyc, &
            filter(nc)%num_urbanc, filter(nc)%urbanc, &
            filter(nc)%num_snowc, filter(nc)%snowc, &
            filter(nc)%num_nosnowc, filter(nc)%nosnowc,2)

        do i = 1,filter(nc)%num_simc
            c = filter(nc)%simc(i)
            l = clm3%g%l%c%landunit(c)
            if (clm3%g%l%itype(l) == isturb) then
                ! Urban landunit use Bonan 1996 (LSM Technical Note)
                cptr%cps%frac_sno(c) = min( cptr%cps%snowdp(c)/0.05_r8, 1._r8)
            else
                ! snow cover fraction in Niu et al. 2007
                cptr%cps%frac_sno(c) = 0.0_r8
                if(cptr%cps%snowdp(c) .gt. 0.0_r8)  then
                    cptr%cps%frac_sno(c) = tanh(cptr%cps%snowdp(c)/(2.5_r8*zlnd* &
                        (min(800._r8,cptr%cws%h2osno(c)/cptr%cps%snowdp(c))/100._r8)**1._r8) )
                endif
            end if
        end do

        call SnowAge_grain(begc, endc, &
            filter(nc)%num_snowc, filter(nc)%snowc, &
            filter(nc)%num_nosnowc, filter(nc)%nosnowc)

        ! 44444444444444444444444444444444444444444444444444444444444444444444444
#if (defined CN)
        ! fully prognostic canopy structure and C-N biogeochemistry
        ! - CNDV defined: prognostic biogeography; else prescribed
        call CNZeroFluxes_dwt(filter(nc)%num_soilc, filter(nc)%soilc,   &
            filter(nc)%num_soilp, filter(nc)%soilp) !CP added to keep mass balance with crop

        call CNEcosystemDyn(begc,endc,begp,endp,filter(nc)%num_soilc,&
            filter(nc)%soilc, filter(nc)%num_soilp, &
            filter(nc)%soilp, filter(nc)%num_pcropp, &
            filter(nc)%pcropp, doalb)
        call CNAnnualUpdate(filter(nc)%num_soilc,&
            filter(nc)%soilc, filter(nc)%num_soilp, &
            filter(nc)%soilp)


        call CBalanceCheck(begc, endc, filter(nc)%num_soilc, filter(nc)%soilc)
        call NBalanceCheck(begc, endc, filter(nc)%num_soilc, filter(nc)%soilc)
#endif

#if (defined CLM_Only)
        call BalanceCheck(begp, endp, begc, endc, begl, endl, begg, endg)
#endif

        if (doalb) then
            !call t_startf('surfalb')

            ! Albedos for non-urban columns
            ! nextsw_cday, declinp1!!!!!!!!!!!!!
            call SurfaceAlbedo(begg, endg, begc, endc, begp, endp, &
                filter(nc)%num_simc, filter(nc)%simc, &
                filter(nc)%num_simp, filter(nc)%simp, &
                nextsw_cday, declinp1)


            !call t_stopf('surfalb')

            ! Albedos for urban columns

            !call t_startf('urbsurfalb')

            !        if (filter(nc)%num_urbanl > 0) then
            !           call UrbanAlbedo(nc, begl, endl, begc, endc, begp, endp,   &
            !                            filter(nc)%num_urbanl, filter(nc)%urbanl, &
            !                            filter(nc)%num_urbanc, filter(nc)%urbanc, &
            !                            filter(nc)%num_urbanp, filter(nc)%urbanp)
            !        end if

            !call t_stopf('urbsurfalb')
        end if

        call updateAccFlds()

    else
        call endrun('clm_driver direction error')

    endif

    end subroutine
    ! ##########################################################################################
    subroutine fixUrban(num_nolakep,nolakep)
    ! Urban treatment:
    ! fixFilters is not correct, because heat fluxes need to be computed on these urban land uses
    ! before a better approach is used, the urban is treated as bareground wrt heat fluxes
    ! however, water fluxes are computed at the same time, so, these fluxes are set to zero here,
    ! as they will be calculated separately in SURFACE_URB in Flow.F90
    ! This function is called before pft2col, in clm_driver

    USE memoryChunk, Only: TEMPZ,ITZ,TEMPL,ITL,FINDEMPTY,mm
    USE clm_varpar, Only: numpft, numveg
    implicit none

    integer num_nolakep
    integer nolakep(num_nolakep)
    integer i,j,k,p,c,nn
    integer,pointer :: RPT(:)=>NULL(), purb(:)=>NULL()
    INTEGER,SAVE:: TNI = 0, num_urb
    INTEGER :: isurb
    real*8,dimension(:),pointer::qflx_evap_soi=>NULL(),qflx_evap_tot=>NULL(),qflx_rain_grnd=>NULL()
    real*8,dimension(:),pointer::qflx_dew_grnd=>NULL(),qflx_evap_grnd=>NULL()

    qflx_evap_soi => clm3%g%l%c%p%pwf%qflx_evap_soi
    qflx_evap_tot => clm3%g%l%c%p%pwf%qflx_evap_tot
    qflx_rain_grnd=> clm3%g%l%c%p%pwf%qflx_rain_grnd
    qflx_dew_grnd => clm3%g%l%c%p%pwf%qflx_dew_grnd
    qflx_evap_grnd => clm3%g%l%c%p%pwf%qflx_evap_grnd
    !% LUidx is a special array that records the codes of
    !% [urban wetland 0  0  0]
    !% in the landuse datasets
    !% This array is needed because different LULC codes may be used, yet some
    ! % codes inside fortran needs to discern which one is urban,wetland, ...etc
    !nn = numveg

    isurb = int(g%Veg%LUidx(1))
    RPT => oneD_ptrI(nn,g%Veg%RPT)
    IF (TNI .eq. 0) THEN
        TNI = FINDEMPTY(ITL,1,size(ITL))
        purb => oneD_ptrI(nn,TEMPL(:,:,TNI:TNI+mm(4)-1));    ITL(TNI:TNI+mm(4)-1) = 1;
        K = 0
        DO I = 1,num_nolakep
            p = nolakep(I)
            if (RPT(p) .eq. isurb) then
                K = K + 1
                purb(K) = p
            endif
        ENDDO
        num_urb = K
    ELSE
        purb => oneD_ptrI(nn,TEMPL(:,:,TNI:TNI+mm(4)-1));
    ENDIF

    DO I = 1,num_urb
        p = purb(I)
        qflx_evap_tot(p) = 0D0 ! will be done in SURFACE_URB
        qflx_evap_soi(p) = 0D0
        qflx_rain_grnd(p)= 0D0
        qflx_dew_grnd(p) = 0D0
        qflx_evap_grnd(p) = 0D0
        clm3%g%l%c%p%pwf%qflx_efpot(p) = 0D0
        clm3%g%l%c%p%pef%cgrndl(p) = 0D0
        !clm3%g%l%c%p%pwf%qflx_snow_grnd(p) = 0D0 ! will be used in balance check because intercepted=forcing-through.
    ENDDO

    end subroutine

    ! ##########################################################################################
    subroutine nanFlagging
    ! apply NaN to: cps, pps, pepv
    USE memoryChunk, only: R3,R4,R3C,R4C
    implicit none
    type(column_pstate_type),pointer::cps
    type(pft_pstate_type),pointer::pps
    type(pft_epv_type),pointer::pepv
    integer b,e
    integer,dimension(2)::I

    cps => clm3%g%l%c%cps
    pps => clm3%g%l%c%p%pps
    pepv=> clm3%g%l%c%p%pepv !CP just give me a handle to see what's going on inside.

    I = R3C%cps
    R3(:,:,I(1):I(2)) = NaN

    I = R4C%pps
    R4(:,:,:,I(1):I(2)) = NaN

    I = R4C%pepv
    R4(:,:,:,I(1):I(2)) = NaN
    end subroutine

    ! ##########################################################################################

    subroutine fixFilters
    ! fix some filter issues
    ! already fixed in filterMod is the MASK is taken into account
    ! here we resolve some issues that result from some minor incompatibilities, e.g. urban is resolved using PAWS itself/bare ground has another code
    USE filtermod, Only: filter
    USE clm_varpar, Only: numpft, numveg
    ! soilp, nolakep, simp
    implicit none
    integer i,j,k,nc,p,c,ivt
    type(pft_type),pointer::pp

    pp => clm3%g%l%c%p
    nc = 1

    k = 0
    do i=1,filter(nc)%num_soilp
        p = filter(nc)%soilp(i)
        ivt = pp%itype(p)
        if (ivt>0 .and. ivt<=numveg) then
            k = k + 1
            filter(nc)%soilp(k) = p
        endif
    enddo
    filter(nc)%num_soilp = k

    ! not doing urban using their code right now
    k = 0
    do i=1,filter(nc)%num_nolakep
        p = filter(nc)%nolakep(i)
        ivt = pp%itype(p)
        if (ivt>0 .and. ivt<=numveg) then
            k = k + 1
            filter(nc)%nolakep(k) = p
        endif
    enddo
    filter(nc)%num_nolakep = k


    k = 0
    do i=1,filter(nc)%num_simp
        p = filter(nc)%simp(i)
        ivt = pp%itype(p)
        if (ivt>0 .and. ivt<=numveg) then
            k = k + 1
            filter(nc)%simp(k) = p
        endif
    enddo
    filter(nc)%num_simp = k

    end subroutine

    !#######################################################################################
    subroutine buildRelations(m,mapp,frac)
    ! Establish geographic relations between grid,lu,col and pft
    ! in PAWS-CLM linking, 1:1:1 mapping is assumed between gridcell, land unit and column
    ! l%itype and c%itype need to be derived
    ! Urban is not currently activated, since SoilTemperature is applied on nolakec, the main impact-
    ! would be that some urban heat fluxes are not considered. The impervious cover impact on hydrology is separately-
    ! considered in PAWS
    USE clm_varcon,  Only: istsoil,icol_road_imperv,istcrop
    USE pftvarcon, Only: npcropmin, numpft
    implicit none
    integer m(4)
    integer mapp(m(1),m(2)),map3(m(1),m(2),m(4)),smap(m(1),m(2))
    real*8 frac(m(1),m(2),m(4))
    type(pft_type),pointer::p
    type(column_type),pointer::c
    type(gridcell_type),pointer::gg
    type(landunit_type),pointer::l
    integer i,j,k
    integer kk,kp
    integer, pointer :: iP(:),iI(:),iJ(:),itt(:)
    REAL*8,DIMENSION(:),POINTER::tt
    INTEGER,DIMENSION(:),POINTER:: MAP1D,ivt
    gg => clm3%g
    l =>    gg%l
    c =>    l%c
    p =>    c%p


    do i= 1,m(4)
        map3(:,:,i)=mapp
    enddo

    MAP1D => oneD_ptrI(nn,map3)
    ivt   => oneD_ptrI(nn,g%Veg%RPT)

    l%lakpoi    = .TRUE.
    l%urbpoi    = .FALSE.
    l%glcmecpoi = .FALSE.
    l%ifspecial = .FALSE.


    !where (MAP1D>0.5D0)
    l%lakpoi = .FALSE.
    !endwhere
    ! CURRENTLY, USE LAKPOI as a way to tell CLM not to work on these columns. Later if lake is activated, then need to change this!!!!
    do i = 1,m(1)
        do j = 1,m(2)


            kk = (j-1)*m(1)+i ! 1D index
            ! g, g:lu:c = 1:1:1
            gg%luni(kk) = kk
            gg%lunf(kk) = kk
            gg%nlandunits(kk) = 1
            gg%coli(kk) = kk
            gg%colf(kk) = kk
            gg%ncolumns(kk) = 1
            gg%pfti(kk) = kk
            gg%pftf(kk) = kk + m(4)*mn - 1
            gg%npfts(kk) = m(4)

            ! landunit
            l%gridcell(kk) = kk
            l%wtgcell(kk) = 1D0
            l%coli(kk) = kk
            l%colf(kk) = kk
            l%ncolumns(kk) = 1
            l%pfti(kk) = kk
            l%pftf(kk) = kk + m(4)*mn - 1
            l%npfts(kk) = m(4)

            ! column
            c%gridcell(kk) = kk
            c%landunit(kk) = kk
            c%wtlunit(kk) = 1D0
            c%wtgcell(kk) = 1D0
            c%pfti(kk) = kk
            c%pftf(kk) = kk + m(4)*mn - 1
            c%npfts(kk) = m(4)


            ! pft
            do k = 1,m(4)
                kp = kk + (k-1)*mn ! 1D index of pft
                p%column(kp) = kk
                p%landunit(kp) = kk
                p%gridcell(kp) = kk
                !p%wtcol(kp) = frac(i,j,k) !! Frac
                !p%wtlunit(kp) = frac(i,j,k)
                !p%wtgcell(kp) = frac(i,j,k)
            enddo

            if (mapp(i,j)>0.5D0) then
                if (g%Veg%RPT(i,j,1) >= npcropmin .and. g%Veg%RPT(i,j,1)<numpft) then
                    l%itype(kk)  = istcrop
                    c%itype(kk)  = istcrop
                    l%lakpoi(kk) = .false.
                else
                    l%itype(kk) = istsoil ! istsoil in clm_varcon
                    c%itype(kk) = istsoil ! urban will be handled later
                    l%lakpoi(kk) = .false.
                endif
            else
                l%itype(kk) = -999 ! istsoil in clm_varcon
                c%itype(kk) = -999 ! urban will be handled later
                l%lakpoi(kk) = .true. !invalid cells are now lakes
            endif
        enddo
    enddo

    end subroutine buildRelations

    ! ##########################################################################################
    subroutine mapGeography(DAY)
    ! extract geography/elevation/time info from PAWS data to CLM
    ! I don't have large domains currently, therefore, I'm using the lat/lon
    !   from stations as the lat/lon values.
    USE clm_time_manager
    USE shr_const_mod, Only: SHR_CONST_PI
    USE domainMod, Only: domain_init, ldomain
    implicit none
    type(sta_type)::day
    integer I,J,K,L,P,T4(4),size1,size2
    integer,dimension(:),pointer::SMAP=>NULL()
    real*8,dimension(:),pointer::ELEV=>NULL()
    type(gridcell_type),pointer::gg=>NULL()
    Type(Vege_type),POINTER::V


    V    => g%Veg
    SMAP => oneD_ptrI(mn,V%SMAP)
    gg   => clm3%g

    DO I = 1,mn
        L = SMAP(I)
        gg%lat(I) = DAY%lat(L)  * (SHR_CONST_PI/180D0)
        gg%lon(I) = DAY%lon(L)  * (SHR_CONST_PI/180D0)
        gg%latdeg(I) = DAY%lat(L)
        gg%londeg(I) = DAY%lon(L)
        gg%lat_a(I) = gg%lat(I)
        gg%lon_a(I) = gg%lon(I)
        ! elevation?
    ENDDO

    ! evaluate time steps in clm_time_manager, modified by CP
    dt_rad = g%Veg%dt * 86400D0
    dt_vdz = g%VDZ%dtP(3) * 86400D0
    dt_veg = g%Veg%dt * 86400D0

    call domain_init(ldomain,V%ny,V%nx,1,mn,'PAWS')
    ELEV => oneD_ptr(mn,g%Topo%E)
    size1 = size(ldomain%topo)
    size2 = size(ELEV)
    ldomain%topo = ELEV

    clm3%ngridcells = mn
    clm3%area       = mn*g%DM%d(1)*g%DM%d(2)/1D6 !km2

    end subroutine

    !###########################################################################################
    subroutine fileInit
    USE clm_varctl, Only: iulog,fpftcon,fsnowaging,fsnowoptics
    !USE SNICARMod,  Only: SnowAge_init,SnowOptics_init
    USE SNICARMod
    use matTable
    integer i,j,k,mode,id, uid, eof
    character(len=300)::line

    OPEN(iulog,file = 'CLM4PAWS.log')

    fsnowoptics= 'snicar_optics_5bnd_c090915.mat'
    fsnowaging = 'snicar_drdt_bst_fit_60_c070416.mat'
    fpftcon  = 'pftcon.mat'

    call buildMatTable(fsnowoptics)
    call buildMatTable(fsnowaging)
    call buildMatTable(fpftcon)

    call SnowOptics_init()
    call SnowAge_init()

    nobsc = 0
    obsc  = 0
    nobsp = 0
    obsp  = 0
    eof   = 0
    uid   = 9
    OPEN(unit = uid, file = 'CLMpoints.dat', ACTION='READ',IOSTAT = eof )
    !read(uid,'(A200)',iostat = eof) line
    !read(line,*) mode, id
    do while (eof .eq. 0)
        read(uid,*,iostat = eof) mode, id
        if (eof .eq. 0) then
            if (mode .eq. 1) THEN
                ! C
                nobsc = nobsc + 1
                obsc(nobsc) = id
            elseif (mode .eq. 2) THEN
                ! P
                nobsp = nobsp + 1
                obsp(nobsp) = id
            endif
        endif
    enddo
    CLOSE(3)

    OPEN(unit = 15,file = 'cRec.txt')
    OPEN(unit = 16,file = 'pRec.txt')

    end subroutine


    ! ##########################################################################################
    subroutine clmForcing(DAY,H,m,SMAP)
    ! From Land-based observation stations to clm forcing variables
    ! so-called data atmospheric model
    !TEMPORARY: potential temperature?? snow???
    ! in Day:
    ! RADD : shortwave direct (beam) radiation (MJ/day)
    ! RADI : diffuse radiation (MJ/day)
    ! WND  : wind speed (m/s)
    ! TEMP : temperature (Celsius)
    ! PA   : atmospheric pressure (kPa)
    ! prcp : precipitation (m/day)
    ! HMD  : relative humidity (volume fraction)
    USE clm_atmlnd
    USE shr_const_mod
    use VDB, Only : VDBCON_RAD2PAR, VDBCON_CP
    USE filterMod, Only: filter
    USE QSatMod, Only :QSat
    USE SNOWUEB_MOD, Only: PARTSNOW
    implicit none
    TYPE(Sta_type)::DAY
    INTEGER H,m(4)
    INTEGER SMAP(m(1),m(2))
    type(atm2lnd_type),pointer::a2l
    type(column_pstate_type),pointer::cps

    REAL*8,DIMENSION(:),POINTER ::forc_t=>NULL(),forc_u=>NULL(),forc_v=>NULL(), PET=>null()
    REAL*8,DIMENSION(:),POINTER ::forc_wind=>NULL(),forc_q=>NULL(),forc_rh=>NULL(),forc_hgt=>NULL()
    INTEGER,DIMENSION(:),POINTER::S=>NULL()
    INTEGER :: I,J,K,L,C
    REAL*8 hsq2,dptp,rr,r1,r2,t
    REAL*8 esdT,qsdT
    real*8,dimension(size(DAY%lat))::es,snow,rain,qs

    !rr   = 1320D0/1367D0 ! 1367 is all radiation (which was used in Spokas), 1320 is excluding ultra-violet
    rr = 1D0 ! should be 1 if rescaling is applied while generating the weather file
    S => oneD_ptrI(mn,SMAP)
    hsq2 = sqrt(2D0)/2D0
    r1 = rr * VDBCON_RAD2PAR
    r2 = rr * (1D0 - VDBCON_RAD2PAR)
    a2l => clm_a2l
    cps => clm3%g%l%c%cps
    PET => oneD_ptr(mn,g%Veg%PET)

    DO L = 1,size(DAY%lat)
        DPTP = min(DAY%DPTP(H,L) + 273.15D0,DAY%TEMP(H,L) + 273.15D0)
        CALL QSat (DPTP, DAY%PA(H,L) * 1D3 , es(L), esdT, qs(L), qsdT)
        SNOW(L) = PARTSNOW(DAY%prcp(H,L) /86.4D0,DAY%TEMP(H,L),3D0,-1D0)
        RAIN(L) = DAY%prcp(H,L) /86.4D0 - SNOW(L)
        DAY%SNOWF(H,L) = SNOW(L)*86.4D0
    ENDDO
    ! question is does it make sense to force potential temperature/specific humidity in this fashion
    DO I = 1,filter(1)%num_soilc
        C = filter(1)%soilc(I)
        L = S(C)
        IF (L>0) THEN
            ! P235 of CLM4 technical doc
            a2l%forc_t(C)        = DAY%TEMP(H,L) + 273.15D0  ! to Kelvin!
            a2l%forc_u(C)        = DAY%WND(H,L) * hsq2
            a2l%forc_v(C)        = DAY%WND(H,L) * hsq2
            !a2l%forc_wind
            !a2l%forc_rh
            a2l%forc_hgt(C)      = 2D0
            a2l%forc_hgt_u(C)    = 2D0 ! what is the impact except for resistance values???? ? ? ? ????
            a2l%forc_hgt_t(C)    = 2D0
            a2l%forc_hgt_q(C)    = 2D0
            a2l%forc_pbot(C)     = DAY%PA(H,L) * 1D3 ! originally in kPa! may need to downscale to location (correct for elevation!)
            a2l%forc_th(C)       = a2l%forc_t(C) ! P235 of CLM4 tech !a2l%forc_t(C)*(DAY%PA(H,L)/101.3D0)**(SHR_CONST_RGAS/VDBCON_CP)  ! potential temperature,K
            !a2l%forc_rho(C)      = 1.292D0 - (0.00428D0 * a2l%forc_t(C))  ! THIS IS A BIG BUG CAUSING 1 WEEK OF DEBUGGING EFFORT!!!!
            a2l%forc_rho(C)      = 1.292D0 - (0.00428D0 * (a2l%forc_t(C)-273.15D0))
            !a2l%forc_q(C)        = a2l%forc_pbot(C)/a2l%forc_t(C)*((1D0-DAY%HMD(H,L))/SHR_CONST_RDAIR+DAY%HMD(H,L)/SHR_CONST_RWV)   ! specific humidity, kg/kg
            a2l%forc_q(C)        = qs(L) !0.622D0*es(I)/(a2l%forc_pbot(C)-0.378D0*es(I))  !specific humidity, kg/kg
            a2l%forc_psrf(C)     = a2l%forc_pbot(C) ! What is the difference?
            a2l%forc_pco2(C)     = DAY%CO2(H,L)/1D6 * a2l%forc_pbot(C) ! DAY%CO2 is ppmv
            a2l%forc_solar(C)    = DAY%RAD(H,L) *(1D6/86400D0) ! incident solar radiation, including shortwave and near-infrared
            a2l%forc_lwrad(C)    = DAY%HLi(H,L) *(1D6/86400D0) ! also from MJ/m2/day to W/m2
            a2l%forc_solad(C,1)  = r1 * DAY%RADD(H,L)
            a2l%forc_solad(C,2)  = r2 * DAY%RADD(H,L)
            a2l%forc_solai(C,1)  = r1 * DAY%RADI(H,L)
            a2l%forc_solai(C,2)  = r2 * DAY%RADI(H,L)

            a2l%forc_rain(C)     = RAIN(L) !DAY%prcp(H,L) /86.4D0 ! from m/day to mm/s!
            a2l%forc_snow(C)     = SNOW(L) !0D0 ! snow handled differently!
            a2l%forc_ndep(C)     = 1.5D0/(86400D0*365D0) ! Annual net N deposition rate (current) from REAY et al, Nature Geoscience, Vol 1, July 2008
            a2l%rainf(C)         = DAY%prcp(H,L) /86.4D0 ! from m/day to mm/s!
#ifdef C13
            !a2l%forc_pc13o2
#endif
            a2l%forc_po2(C)      = DAY%PA(H,L)* 210D0 ! need to convert from kPa to Pa
            !a2l%forc_aer => twoD_ptr(mn,1,R3(:,:,i3)); i3=i3+1 ! not using snowhydrology or anything related to aerosol
            cps%decl(c)          = DAY%decl(H,L)
            PET(C)               = day%rets(H,L)*g%Veg%dt ! m/day bare soil potential ET according to P-M. convert it to m by * dt
            !if (RAIN(L) >0D0) THEN
            !PRINT*,'clmForcing'
            !endif
            !  write(c,FMT='(100E12.5)') a2l%forc_th(c),a2l%forc_lwrad(c),a2l%forc_q(c),	a2l%forc_u(c),	a2l%forc_pbot(c),	a2l%forc_solad(c,1),	a2l%forc_solad(c,2),	a2l%forc_solai(c,1),	a2l%forc_solai(c,2)

        ENDIF
    ENDDO

    end subroutine
    ! ##########################################################################################
    subroutine waterPar
    !watsat
    !bsw
    !sucsat
    !bsw2
    !psisat
    !vwcsat
    !xksat
    !hksat
    !hksat
    !tkmg
    !tksatu
    !tkdry
    !csol
    !watdry
    !watopt
    !watfc
    ! many of them not used if soilHydrology is turned off (using PAWS)

    end subroutine

    ! ##########################################################################################
    subroutine dailyVegUpdateNew
    ! this subroutine replaces the old dailyVegUpdate and does two things:
    ! 1. use the new LAI to compute ho, hod
    ! 2. manage agricultural operations, including, planting and harvesting (automatically scheduled by heat units)




    end subroutine

    ! ##########################################################################################
    subroutine passV(dir)
    ! True Physical Linking variable passing
    ! water state may be updated in PAWS, this is the linking/converting done at every time step
    USE filterMod, only: filter
    use clm_varpar  , only : maxpatch_pft, nlevsno, nlevgrnd, numrad, nlevlak, &
        numpft, ndst, nvoc, nlevurb, nlevsoi
    USE pftvarcon, Only: npcropmin, numpft
    use clm_varcon  , only : denh2o, denice, wimp, ssi
    USE abortutils, Only : endrun
    USE memoryChunk, Only: TEMP,IT,FINDEMPTY,mm,tempStorage
    USE paws_varctl, Only: dbgflg,writeout
    USE clm_atmlnd, Only: clm_a2l
    implicit none
    integer dir
    type(GlobalG),pointer::gp
    type(vege_type),pointer::V
    type(vdz_type),pointer:: VD
    type(pft_type),pointer ::pp
    type(column_type),pointer ::cc
    type(pft_wflux_type),pointer ::pwf_a
    real*8,dimension(:),pointer::Gprcp=>NULL(),EvapG=>NULL(),TP=>NULL(),ET=>NULL()
    real*8,dimension(:),pointer::EvapCS=>NULL(),Evap=>NULL(),CS=>NULL(),Dperc=>NULL()
    real*8,dimension(:),pointer::Rf=>NULL(),DEW=>NULL(),SWE=>NULL(),SNOM=>NULL()
    real*8,dimension(:),pointer::CSNEW=>NULL(),PET=>NULL(),Df=>NULL(),Imp=>NULL()
    real*8,dimension(:,:),pointer::THE=>NULL(),DZ=>NULL(),ICE=>NULL(),SRC=>NULL()
    integer,dimension(:),pointer:: WTL=>NULL(),NZC=>NULL()
    REAL*8,DIMENSION(:),POINTER::VS=>NULL(),h1=>NULL(),h2=>NULL(),h3=>NULL(),h4=>NULL()

    real*8,dimension(:),pointer::h2ocan=>NULL(),qflx_drain=>NULL(),qflx_surf=>NULL()
    real*8,dimension(:),pointer::h2osno=>NULL(),qflx_sub_snow=>NULL(),qflx_dew_snow=>NULL()
    real*8,dimension(:),pointer::qflx_infl=>NULL()
    real*8,dimension(:,:),pointer::rootr_col=>NULL(),h2osoi_liq=>NULL(),h2osoi_ice=>NULL(), dormant_flag=>null()
    real*8,dimension(:,:),pointer::h2osoi_vol=>NULL(),t_soisno,rresis=>NULL()
    real*8,dimension(:,:),pointer::qflx_tran_veg_pft=>NULL(),wtcol=>NULL()
    real*8,dimension(:,:),pointer::qflx_prec_intr=>NULL(),qflx_candrip=>NULL()
    real*8,dimension(:,:,:),pointer::rootr_pft=>NULL(),rootfr_pft=>NULL() ! calling rootr(P,k,z) is the same as calling (I,J,K,Z)?
    integer , pointer :: snl(:)=>NULL(),PFT(:,:)=>NULL(),ivt(:,:)=>null()           ! snl: number of snow layers
    integer,dimension(:),pointer::RDFN=>NULL()
    real*8 dt,conv,conv2,conv3,conv4,wtcolVeg,liqChange,iceChange,totChange,dhSurf,theIceBuffer,prcpCan,nonUrbSnow,the2Temp2,deniceFactor
    REAL*8 SS,m2,output(100),t1,t2,ST,THRSD,GAM,sumthe,sumthez,t3,t4,t5,t6,Qout
    REAL*8 snoLiqEvap,nonSnow2air,dtime
    REAL*8 THERR
    REAL*8 :: THER_MIN = 0.002D0 ! should be consistent with subroutine vdz1c
    REAL*8,DIMENSION(:,:),POINTER::THER=>NULL(),THES=>NULL(),N=>NULL(),ALPHA=>NULL(),iceMove=>null()
    REAL*8,DIMENSION(:,:),POINTER::HH=>NULL(),FC=>NULL(),WP=>NULL()
    REAL*8,DIMENSION(:),POINTER:: iceMovePTR=>null()
    REAL*8,POINTER::THE2tempPtr(:)=>null(),THE2temp(:,:)=>NULL(),liq1TEMP(:,:)=>null(),liq1TEMPptr(:)=>NULL(),lastTHE(:,:,:)=>null(),lastForc(:,:)=>null()
    REAL*8,POINTER::liq0TEMP(:,:)=>null(),liq0TEMPptr(:)=>NULL(),ice0TEMP(:,:)=>null(),ice0TEMPptr(:)=>NULL(),lastTHEptr(:,:)=>null(),lastForcPtr(:)=>null()
    REAL*8,POINTER::ice1TEMP(:,:)=>null(),ice1TEMPptr(:)=>NULL(),Combo(:,:)=>null(),FrozenInf(:,:)=>null(),Soil2AirDirectPTR(:)=>null()
    REAL*8,POINTER::irrigationPTR(:)=>null()
    REAL*8,POINTER::ReserveH1(:,:)=>null(),ReserveH1ptr(:)=>null(),Soil2AirDirect(:,:)=>null(),irrigation(:,:)=>null(),dormant(:,:)=>null(),dormantPtr(:)=>null()
    REAL*8,POINTER::ComboPtr(:)=>null(),FrozenInfPtr(:)=>null(),ptrp=>null(),ptrc=>null()
    REAL*8,POINTER::qflx_snow_grnd(:,:)=>null(),qflx_dew_grnd(:)=>null(),qflx_evap_grnd(:)=>null()
    integer P,C,m(3),I,J,K,nt,ntf
    integer,pointer:: snl_old(:), snl2d(:,:) ! to save the info if snow layer existed
    integer,save :: snl_state = 0
    REAL*8,DIMENSION(:),POINTER::icesum0=>null(), the_sum0=>null(), h2osum0=>null(), cs0=>null(),swe0=>null()
    REAL*8,DIMENSION(:),POINTER::the_sum=>null(), icesum=>null(),h2osum=>null()
    REAL*8,DIMENSION(:,:),POINTER::icesum02D=>null(), the_sum02D=>null(), h2osum02D=>null(), cs02D=>null(),swe02D=>null()
    REAL*8,DIMENSION(:,:),POINTER::the_sum2D=>null(), icesum2D=>null(),h2osum2D=>null()
    integer,save:: TNI,TNI2
    logical :: passVdbg = .true.
    ! dir = 1, from PAWS to CLM
    ! dir = 2, from CLM to PAWS
    ! Think -- coupled water flow/energy problem? small cells stability? soil evaporation??
    gp      => g
    V       => g%Veg
    VD      => g%VDZ
    DZ      => twoD_ptr(mn,VD%m(3),VD%DZ)
    THE     => twoD_ptr(mn,VD%m(3),VD%THE)
    hh      => twoD_ptr(mn,VD%m(3),VD%h)
    THES    => twoD_ptr(mn,VD%m(3),VD%THES)
    THER    => twoD_ptr(mn,VD%m(3),VD%THER)
    ALPHA   => twoD_ptr(mn,VD%m(3),VD%ALPHA)
    N       => twoD_ptr(mn,VD%m(3),VD%N)
    ICE     => twoD_ptr(mn,VD%m(3),VD%ICE)          ! ice kg/m2
    FC      => twoD_ptr(mn,VD%m(3),VD%FC)
    WP      => twoD_ptr(mn,VD%m(3),VD%WP)
    PFT     => twoD_ptrI(mn,VD%m(3),V%RPT)
    Rf      => oneD_ptr(mn,g%OVN%Rf)
    Df      => oneD_ptr(mn,g%VDZ%Df)
    Imp     => oneD_ptr(mn,g%Veg%Imp)
    nt      = V%nt

    cc      => clm3%g%l%c
    pp      => cc%p
    pwf_a   => cc%cwf%pwf_a
    h2osoi_liq    => cc%cws%h2osoi_liq
    h2osoi_ice    => cc%cws%h2osoi_ice
    h2osoi_vol    => cc%cws%h2osoi_vol
    qflx_surf     => clm3%g%l%c%cwf%qflx_surf
    qflx_infl     => clm3%g%l%c%cwf%qflx_infl
    wtcol         => twoD_ptr(mn,V%nt,clm3%g%l%c%p%wtcol)
    qflx_snow_grnd=> twoD_ptr(mn,V%nt,clm3%g%l%c%p%pwf%qflx_snow_grnd)
    rootr_pft     => threeD_ptr(m,pp%pps%rootr)            ! reshaped/re-dimensioned. now layer numbering is the same as PAWS. root(c,k,z) is for k-th PTF, z-th layer
    rootfr_pft    => threeD_ptr(m,pp%pps%rootfr)
    h2osno        => clm3%g%l%c%cws%h2osno
    t_soisno      => clm3%g%l%c%ces%t_soisno
    rresis        => clm3%g%l%c%p%pps%rresis
    GPrcp => oneD_ptr(mn,V%GPrcp)
    PET   => oneD_ptr(mn,V%PET)

    conv    =  V%dt * 86400D0/1000D0 ! from mm/s (CLM) to m/day (PAWS) and then times it by dt
    conv2   =  86400D0/1000D0        ! from mm/s (CLM) to m/day (PAWS)
    dt      =  V%dt
    deniceFactor = denh2o/denice
    dtime   = dt * 86400D0

    IF (TNI .eq. 0) THEN
        TNI = FINDEMPTY(IT,4,size(IT)); IT(TNI:TNI+3) = 1;
    ENDIF
    h1 => oneD_ptr(mn,TEMP(:,:,TNI))
    h2 => oneD_ptr(mn,TEMP(:,:,TNI+1))
    h3 => oneD_ptr(mn,TEMP(:,:,TNI+2))
    h4 => oneD_ptr(mn,TEMP(:,:,TNI+3))

    call tempStorage('THE2temp',THE2temp,1)
    call tempStorage('liq1TEMP',liq1TEMP,1)
    call tempStorage('liq0TEMP',liq0TEMP,1)
    call tempStorage('ice1TEMP',ice1TEMP,1)
    call tempStorage('ice0TEMP',ice0TEMP,1)
    call tempStorage('FrozenInf',FrozenInf,1)
    call tempStorage('Combo',Combo,1)
    call tempStorage('ReserveH1',ReserveH1,1)
    call tempStorage('snl2d',snl2d,1)
    call tempStorage('Soil2AirDirect',Soil2AirDirect,1)
    call tempStorage('iceMove',iceMove,1)
    call tempStorage('irrigation',irrigation,1)
    call tempStorage('lastTHE',lastTHE,VD%m(3))
    call tempStorage('dormant',dormant,1)
    !call tempStorage('lastTHE',lastForc,1)
    irrigationPtr => oneD_Ptr(mn,irrigation)
    THE2tempPtr => oneD_Ptr(mn,THE2temp)
    liq1TEMPPtr => oneD_Ptr(mn,liq1TEMP)
    liq0TEMPPtr => oneD_Ptr(mn,liq0TEMP)
    ice1TEMPPtr => oneD_Ptr(mn,ice1TEMP)
    ice0TEMPPtr => oneD_Ptr(mn,ice0TEMP)
    dormantPtr => oneD_Ptr(mn,dormant)
    FrozenInfPtr => oneD_Ptr(mn,FrozenInf)
    ComboPtr => oneD_Ptr(mn, Combo)
    ReserveH1Ptr => oneD_Ptr(mn, ReserveH1)
    Soil2AirDirectPtr => oneD_Ptr(mn, Soil2AirDirect)
    iceMovePtr => oneD_Ptr(mn, iceMove)
    snl_old => oneD_PtrI(mn, snl2d)
    lastTHEptr   => twoD_ptr(mn,VD%m(3),lastTHE)
    !lastForcPtr => oneD_Ptr(mn,lastForcPtr)
    theIceBuffer = 0.005D0 ! don't get too close to THER

    IF (debugg) then
        call tempStorage('swe0',SWE02D,1)
        call tempStorage('cs0',cs02D,1)
        call tempStorage('the_sum0',the_sum02D,1)
        call tempStorage('h2osum0',h2osum02D,1)
        call tempStorage('icesum0',icesum02D,1)
        call tempStorage('the_sum',the_sum2D,1)
        call tempStorage('h2osum',h2osum2D,1)
        call tempStorage('icesum',icesum2D,1)

        SWE0 => oneD_Ptr(mn, SWE02D)
        cs0 => oneD_Ptr(mn, cs02D)
        the_sum0 => oneD_Ptr(mn, the_sum02D)
        h2osum0 => oneD_Ptr(mn, h2osum02D)
        icesum0 => oneD_Ptr(mn, icesum02D)
        the_sum => oneD_Ptr(mn, the_sum2D)
        h2osum => oneD_Ptr(mn, h2osum2D)
        icesum => oneD_Ptr(mn, icesum2D)
    ENDIF



    h3 => oneD_ptr(mn,g%VDZ%R(:,:,3))  ! Previously, I was using a local save variable TEMP for h3, but I found that this information needs to be saved
    ! therefore, I am using g%VDZ%R(:,:,3) as a place to save this value

    IF (DIR .eq. 2) then
        !! 111111111111111111111111111111111111111111111111111111111111111111111111
        ! From PAWS to CLM
        ! h2o2soi_vol, h2osoi_liq
        ! psi, qflx_drain (for NLeaching)
        qflx_drain    => clm3%g%l%c%cwf%qflx_drain
        GAM = VD%VParm(1)

        Dperc         => oneD_ptr(mn,VD%Dperc)
        WTL           => oneD_ptrI(mn,VD%WTL)


        DO I = 1,filter(1)%num_nolakec
            C = filter(1)%nolakec(I)
            ! these conversions below don't actually matter except for budgeting. They ought to be cleaned by the end of the time step
            qflx_drain(C) = Dperc(C) / conv
            qflx_surf(C)  = Rf(C)    / conv
            ! layer numbering
            ! (in PAWS we don't modify h2osoi_ice, only use it to rescale conductivity)
            !CP: I found that if simply using THE(2) to compute h2osoi_liq(1), ET is so much less than PET. I suspect the ponding layer was not properly used.
            !the strategy is to use both THE(2) and h(1) to compute h2osoi_liq, and then recover them after computing inside CLM
            ! However, must be careful. h2osoi_liq may be modified inside CLM due to: snowmelt, freeze/thaw; also, h(1) may become slightly negative in VDZ1C calculations
            ! This combination is done: only h(1)>0,

            h1(C) = hh(C,1)
            h2(C) = THE(C,2)
            J = 1
            t1    = max(THES(C,2)-THER(C,2)-ICE(C,2)/1000D0/DZ(C,2)*deniceFactor,0D0)
            ReserveH1Ptr(C)= max(0D0,hh(C,J)-t1*DZ(C,2))     ! This portion of h1 cannot be combined with THE(C,2) because there is no more pore space
            if (hh(C,1) > 0D0) THEN
                h2osoi_liq(C,J) = ((THE(C,J+1)-THER(C,J+1))*DZ(C,J+1) + hh(C,J)-ReserveH1Ptr(C))*1000D0 ! m3/m2 -> 1000kg/m2
                h2osoi_vol(C,J) = THE(C,J+1) - THER(C,J+1)
                h2osoi_ice(C,J) = ICE(C,J+1)
                h3(C)           = hh(C,J)    ! h3 stores the amount of water on the surface ponding layer that is counted toward h2osoi_liq (because CLM requires so). For mass balance, this needs to be subtracted from h2osoi_liq while converting back to THE and h(1)
                !hh(C,J)         = 0D0
            ELSE
                h2osoi_vol(C,J) = THE(C,J+1) - THER(C,J+1)  ! Verified that it is indeed making the right index on 10/25/10
                h2osoi_liq(C,J) = (THE(C,J+1)-THER(C,J+1))*DZ(C,J+1)*1000D0 ! m3/m2 -> 1000kg/m2
                h2osoi_ice(C,J) = ICE(C,J+1)
                h3(C)           = 0D0
            ENDIF
            THE2tempPtr(C) = THE(C,2);
            liq1tempPtr(C) = h2osoi_liq(C,1);
            liq0tempPtr(C) = h2osoi_liq(C,0);
            ice1tempPtr(C) = h2osoi_ice(C,1);
            ice0tempPtr(C) = h2osoi_ice(C,0);

            DO J = 2,nlevgrnd
                !DO J = 1,nlevgrnd
                ! Because in CLM soil parameterization, there is no THER, however, without THER, these plants may suck too much water, therefore
                ! we take THER out of available soil water storage, and add it back when passing variables out
                h2osoi_vol(C,J) = THE(C,J+1) - THER(C,J+1)  ! Verified that it is indeed making the right index on 10/25/10
                h2osoi_liq(C,J) = (THE(C,J+1)-THER(C,J+1))*DZ(C,J+1)*1000D0 ! m3/m2 -> 1000kg/m2
                h2osoi_ice(C,J) = ICE(C,J+1) ! a fix-up in vdz1c may modify its content
                ! Before adding THER:
                !h2osoi_vol(C,J) = THE(C,J+1)  ! Verified that it is indeed making the right index on 10/25/10
                !h2osoi_liq(C,J) = THE(C,J+1)*DZ(C,J+1)*1000D0 ! m3/m2 -> 1000kg/m2

                ! root water uptake function:
                ! use LK00 paper
                ! this resistance currently does not include compensation effect. can be implemented in canopyFluxes


            ENDDO

            DO J = 1,nlevgrnd
                THRSD = .05D0*(FC(C,J+1)-WP(C,J+1))+WP(C,J+1) ! safety threshold not to drain water
                ntf   = C + mn*(nt-1)
                IF (THE(C,J+1)<=THRSD) THEN
                    rresis(C:ntf:mn,J) = 0D0
                ELSE
                    ST = (THE(C,J+1)-WP(C,J+1))/THES(C,J+1)
                    if (ST > 1D0) then
                        call display('PAWS_CLM_MOD large ST?')
                    endif
                    rresis(C:ntf:mn,J) = min(ST**(GAM/(THE(C,J+1)-WP(C,J+1))),1D0)

                    t3 = Gprcp(C) - PET(C)
                    if (THE(C,J+1) <= WP(C,J+1)+.05D0*(THES(C,J+1)-WP(C,J+1)) .AND. (THE(C,J+1)<lastTHEptr(C,J+1) .or. t3<0D0)) then ! be careful now
                        ! assuming ET this time step is close to last time step's, predict where theta will be
                        !t1 = max(2D0*THE(C,J+1)-lastTHEptr(C,J+1),WP(C,J+1)+1D-10) ! This doesn't work because using THE^(n-1) is not reliable
                        t1 = THE(C,J+1)
                        do k = 1, V%nt
                            if (wtcol(C,K) > 0D0) then
                                t1 = t1 + t3 * wtcol(C,K)*rootfr_pft(C,K,J)/DZ(C,J+1) ! note here we use rootfr_pft (root distribution) instead of rootr_pft (root-soil distribution) to keep stability because the later is changing
                            endif
                        enddo
                        t1 = max(min(t1,2D0*THE(C,J+1)-lastTHEptr(C,J+1)), WP(C,J+1)+1D-20)
                        ST = (t1-WP(C,J+1))/THES(C,J+1)
                        t2 = ST**(GAM/(t1-WP(C,J+1)))
                        !rresis(C:ntf:mn,J) = (rresis(C:ntf:mn,J)+t2)/2D0
                        rresis(C:ntf:mn,J) = t2
                    endif
                ENDIF


            ENDDO
            !lastForcPtr(C) = Gprcp(C) - PET(C) ! [m] save this for better determining the rresis next time step

            ! also first layer need to be treated differently
            ! in CLM, there is no ponding layer, all water goes to first layer of soil (or runoff/evap)
            ! in PAWS, we have a ponding layer which can store some water
            !!!! CURRENTLY: FROZEN SOIL AND SOIL PONDING NOT WELL SUPPORTED. NOW JUST GET IT RUNNING

            ! print debug
            if (debugg) then
                the_sum0(C) = sum(THE(C,2:w%m(3))*DZ(C,2:w%m(3)))
                h2osum0(C)  = sum(h2osoi_liq(C,1:nlevgrnd))
                icesum0(C)  = sum(h2osoi_ice(C,1:nlevgrnd))
                swe0(C)  = h2osno(C)/1000D0
            endif
        ENDDO
        !! 111111111111111111111111111111111111111111111111111111111111111111111111
    elseif (DIR .eq. 1) THEN
        !! 222222222222222222222222222222222222222222222222222222222222222222222222
        ! From CLM to PAWS
        ! GPrcp,TP,ET,EvapG,EvapCS,CS,RDFN
        ! Snow variables
        ! h2osoi_ice/h2osoi_liq
        TP    => oneD_ptr(mn,V%TP)
        ET    => oneD_ptr(mn,V%ET)
        EvapG => oneD_ptr(mn,V%EvapG)
        EvapCS=> oneD_ptr(mn,V%EvapCS)
        Evap  => oneD_ptr(mn,V%Evap)
        CS    => oneD_ptr(mn,V%CS)
        DEW   => oneD_ptr(mn,V%DEW)
        SWE   => oneD_ptr(mn,V%SWE)
        CSNEW => oneD_ptr(mn,V%CSNEW)
        SNOM  => oneD_ptr(mn,V%SNOM)  ! sublimation - dew_on_snow
        SRC   => twoD_ptr(mn,VD%m(3),VD%SRC)
        NZC   => oneD_ptrI(mn,VD%NZC)

        h2ocan            => cc%cws%pws_a%h2ocan
        m                 = (/VD%m(1)*VD%m(2),V%nt,nlevgrnd/)
        qflx_tran_veg_pft => twoD_ptr(mn,V%nt,pp%pwf%qflx_tran_veg) ! qflx(c,k) is for k-th PTF
        !qflx_candrip      => twoD_ptr(mn,V%nt,pp%pwf%qflx_candrip) ! qflx(c,k) is for k-th PTF
        qflx_prec_intr    => twoD_ptr(mn,V%nt,pp%pwf%qflx_prec_intr) ! qflx(c,k) is for k-th PTF
        snl               => clm3%g%l%c%cps%snl
        qflx_dew_snow     => cc%cwf%pwf_a%qflx_dew_snow
        qflx_sub_snow     => cc%cwf%pwf_a%qflx_sub_snow
        qflx_dew_grnd     => cc%cwf%pwf_a%qflx_dew_grnd
        qflx_evap_grnd     => cc%cwf%pwf_a%qflx_evap_grnd
        dormant_flag       => twoD_ptr(mn,V%nt,cc%p%pepv%dormant_flag)
        ivt               => twoD_ptrI(mn,V%nt,cc%p%itype)

        WHERE (snl_old .ne. 0) snl_old = 0D0

        DO I = 1,filter(1)%num_snowc
            ! ################################################################
            ! the reason to do this is that in CombineSnowLayers, the snl(c)<0 info is lost. This info exists in the filter
            ! in Hydrology2 step1::snowWater, filter(1)%snowc is what is used to apply qflx_rain_grnd(c) + qflx_dew_grnd(c) - qflx_evap_grnd(c)
            ! Otherwise, these fluxes are applied to ground (qflx_rain_grnd included in qflx_top_soil; qflx_dew_grnd - qflx_evap_grnd is qflx_evap_grnd)
            C = filter(1)%snowc(I)
            snl_old(C) = -1
        ENDDO

        DO I = 1,filter(1)%num_nolakec
            C = filter(1)%nolakec(I)
            ! 3 evaporative sources: CanopyStorage, SoilSurface, VegetatedSurface
            ! in CLM qflx_evap_veg is  (TP + EvapCS)
            !if (C .eq. 2695) then
            !write(*,*) 'paws_clm::debug'
            !endif


            ! CLM surface water logic:
            ! nosnowc: qflx_top_soil = qflx_rain_grnd + qflx_snomelt (SnowWater); qflx_dew_grnd was later added to h2osoi_liq in Drainage in SoilHydrology, which we skipped
            ! therefore, qflx_dew_grnd ought to be added to the flux term
            ! snowc  : qflx_top_soil = qout(c) / dtime (SnowWater); However, (qflx_rain_grnd(c) + qflx_dew_grnd(c) - qflx_evap_grnd(c))*dtime were added to h2osoi_liq
            ! as a result, they should be removed from the source term, but should be in Evap for mass balance
            ! UPDATE: 10/11/2013: This above thought is proven to be wrong, as snl(c)<0 for snowc, and therefore never touches h2osoi_liq(1)
            ! Any water to reaches ground for snowc comes from Qout
            FrozenInfPtr(C) = 0D0
            Soil2AirDirectPTR(C)=0D0
            iceMovePTR(C)   = 0D0

            Gprcp(C)  = cc%cwf%qflx_top_soil(C) * conv  ! Hydrology1
            TP(C)     = pwf_a%qflx_tran_veg(C)  * conv  ! CanopyFluxes
            EvapCS(C) = pwf_a%qflx_evap_can(C)  * conv  ! BG2 (qflx_evap_can(p)  = qflx_evap_veg(p) - qflx_tran_veg(p) -- a summation is added in pft2colMod)

            ! EvapG and Gprcp will be source term for ponding layer. Others are for budgets
            if (snl_old(c) < 0) then
                ! BG2. Currently, qflx_evap_grnd is (qflx_evap_soi - evap_over_snow), the total soil evaporation, and has not been distributed to soil layers (applied as flux)
                ! snowc qflx_dew_grnd(c)
                EvapG(C)  = 0D0 !pwf_a%qflx_evap_grnd(C) * conv
                SNOM(C)  = (qflx_sub_snow(C)-qflx_dew_snow(C)) * conv
                snoLiqEvap = (pwf_a%qflx_evap_grnd(C)-pwf_a%qflx_dew_grnd(c)) * conv
                Evap(C)   = EvapG(C) + EvapCS(C) + snoLiqEvap    ! BG2 (qflx_evap_tot(p) = qflx_evap_veg(p) + qflx_evap_soi(p))
            else
                ! nosnowc
                EvapG(C)  = (pwf_a%qflx_evap_grnd(C)-pwf_a%qflx_dew_grnd(C)) * conv  ! Hydrology1
                !EvapG(C)  = qflx_evap_grnd(C)*conv
                SNOM(C) = 0D0
                Evap(C)   = EvapG(C) + EvapCS(C)            ! BG2 (qflx_evap_tot(p) = qflx_evap_veg(p) + qflx_evap_soi(p))
                snoLiqEvap = 0D0
                IF (abs(qflx_sub_snow(C)-qflx_dew_snow(C))>1D-10) THEN
                    ! These fluxes may exist when h2osoi_ice(c,1)>0
                    ! If there is no snow, these fluxes were applied in Drainage() in SoilHydrology. However, we don't do that. so do it here
                    t1 = qflx_dew_snow(c) * dtime; h2osoi_ice(c,1) = h2osoi_ice(c,1) + t1
                    t2 = min(h2osoi_liq(c,1), qflx_sub_snow(C)*dtime); h2osoi_ice(c,1) = h2osoi_ice(c,1) - t2
                    Soil2AirDirectPTR(C) = (t2 - t1)/1000D0
                ENDIF

            endif ! snl IF

            ET(C)     = TP(C) + Evap(C) ! Check if this is pwf_a%qflx_evap_grnd + pwf_a%qflx_evap_veg and pwf_a%qflx_evap_tot
            CS(C)     = h2ocan(C)/1000D0
            DEW(C)    = (pwf_a%qflx_dew_grnd(C) - min(0D0, pwf_a%qflx_evap_can(C))) * conv ! Dew for the mass balance check
            ! Source term input into vdz1c: (impervious?)
            VS        => SRC(C,:)
            WHERE (VS .ne. 0D0) VS = 0D0
            irrigationPTR(C) = clm3%g%l%c%cwf%qflx_irrig(C)  * conv
            VS(1)     = (Gprcp(C) - EvapG(C))/dt ! This is hoping the flux boundary condition can work

            IF (pwf_a%qflx_tran_veg(C) > 0D0) THEN
                DO K = 1,V%nt
                    if (wtcol(C,K) > 0D0) then
                        conv3 = conv2 * wtcol(C,K) * qflx_tran_veg_pft(C,K)
                        !IF (rootr_pft(C,K,1)<=1D-5 .AND. qflx_tran_veg_pft(C,K)>0D0) THEN ! bug
                        IF (rootfr_pft(C,K,1)<=1D-5 .AND. qflx_tran_veg_pft(C,K)>0D0) THEN
                            ! must be special types, i.e., wetland/water/urban
                            VS(1) = VS(1)- conv3
                        ELSEIF (conv3 > 0D0) THEN
                            DO J = 2,nlevsoi+1
                                ! Remember the 1st layer in CLM is soil, but in PAWS it is ground surface
                                IF (J < NZC(C)-1) THEN
                                    VS(J)   = VS(J) - conv3 * rootr_pft(C,K,J-1) /DZ(C,J)  !rootr_pft is corrected
                                ELSEIF (J .eq. NZC(C)-1) THEN
                                    VS(J)   = VS(J) - conv3 * sum(rootr_pft(C,K,J-1:nlevgrnd))/DZ(C,J)
                                ENDIF
                            ENDDO
                            ! the CLM root can go very very deep (unncessary)
                            ! Full form: VS(J)   = VS(J) - qflx_tran_veg_pft(C,K) * conv3 * rootr_pft(C,K,J-1) /DZ(C,J)
                        ENDIF
                    endif
                ENDDO
            ENDIF

            SWE(C)  = h2osno(C)/1000D0  ! This in fact includes both liq and
            !PET(C)  = 0D0 !min(pwf_a%qflx_efpot(C) * conv,0.5D0)

            dormantPtr(C) = 0D0
            DO J=1,V%nt
                if (dormant_flag(C,J) > 0D0 .AND. pftcon%season_decid(ivt(C,J)) .eq. 1D0) then
                    dormantPtr(C)=1
                endif
            ENDDO

            !      CSNEW(C) = 0D0
            !      DO K = 1, V%nt
            !         CSNEW(C) = CSNEW(C)+ wtcol(C,K)*(qflx_prec_intr(C,K))*conv
            !      ENDDO

            ! Although the linked CLM components do not compute soil water flow
            ! Phase may change that require modification
            the2Temp2 = THE(C,2)
            DO J = 2,nlevgrnd
                ICE(C,J)  = h2osoi_ice(C,J-1) ! layer numbering
                !THE(C,J)  = h2osoi_liq(C,J-1)/DZ(C,J)/1000D0
                THE(C,J)  = h2osoi_liq(C,J-1)/DZ(C,J)/1000D0 + THER(C,J) ! add THER back
            ENDDO
            !THE2tempPtr(C) = THE(C,2);
            !liq1tempPtr(C) = h2osoi_liq(C,1);
            !ice2tempPtr(C) = h2osoi_ice(C,2);
            ! now THE(C,2) has all of THE(C,2) and hh(C,1)
            liqChange = h2osoi_liq(C,1) - liq1TEMPptr(C)
            iceChange = h2osoi_ice(C,1) - ice1TEMPptr(C) + Soil2AirDirectPTR(C)*1000D0
            !ComboPtr(C) = 0D0 ! flux that comes from combineSnowLayer
            IF (abs(liqChange) < 1D-14 .AND. abs(iceChange)<1D-14) THEN ! not changed in clm, therefore (THE2(C,2) - THE2tempPtr(C) )*DZ(C,2)=h3(C)
                !hh(C,1)   = hh(C,1) - h3(C)
                !hh(C,1)  = hh(C,1) + (THE2(C,2) - THE2tempPtr(C) )*DZ(C,2) ! no change needs to be made to hh(C,1)
                THE(C,2) = THE2tempPtr(C)
            ELSEIF (abs(liqChange + iceChange)<1D-11) THEN
                t1 = iceChange/1000D0/DZ(C,2)  ! frozen in the time step
                t2 = THE2tempPtr(C) - t1       ! THE(C,2)^{star}. a test value for next THE(C,2)
                dhSurf = hh(C,1) - h3(C)
                t4 = max(THER(C,2)+theIceBuffer,THES(C,2) - ICE(C,2)/1000D0/DZ(C,2)*deniceFactor)  ! available pore space for water
                t5 = min(0D0, max(0D0, the2Temp2 - THER(C,2)) )  ! how close do we allow THE to be to THER
                !hh(C,1)  = hh(C,1) - h3(C) + h2osoi_liq(C,1) -  nothing needs to happen to h(1) unless rescue is needed
                IF (t2 < THER(C,2)+t5) THEN ! too much gets frozen. Shouldn't get too close to THER
                    ! FrozenInf: some surface water is absorbed into the soil to satisfy the freezing demand
                    t3 = (THER(C,2) + t5 - t2)*DZ(C,2) ! deficit. need some fix
                    FrozenInfPtr(C) = t3
                    hh(C,1) = hh(C,1) - t3
                    THE(C,2)= THER(C,2)+t5
                    ! It would be bad if: h(1)
                ELSEIF (t2 > t4+1D-10) THEN
                    t3 = (t2 - t4)*DZ(C,2)
                    FrozenInfPtr(C) = - t3
                    hh(C,1) = hh(C,1) + t3
                    THE(C,2)= t4
                ELSE
                    ! no fix is needed
                    THE(C,2)= t2
                ENDIF
            ELSEIF ((snl(C) .eq. 0 .and. snl_old(C) <0 .and. clm3%g%l%c%cps%snowdp(c)<0.01D0) .OR. ComboPTR(C)>0D0) THEN ! CombineSnowLayers
                !ELSEIF (ComboPTR(C)>0D0) then
                IF (ComboPtr(C) .eq. 0D0) ComboPtr(C) = (liqChange+iceChange)/1000D0  ! This quantifies extra input into ponding. If freezing happens (iceChange>0), there must in fact be more input. If melt happened (iceChange<0). Then there has in fact not been so much extra input (+ a negative iceChange)
                hh(C,1) = hh(C,1) + ComboPtr(C)
                THE(C,2)= THE2tempPtr(C)-iceChange/1000D0/DZ(C,2)
                t4 = max(THER(C,2)+theIceBuffer,THES(C,2) - ICE(C,2)/1000D0/DZ(C,2)*deniceFactor)   ! available pore space for water
                t5 = min(0D0, max(0D0, the2Temp2 - THER(C,2)) )  ! distance from THER
                IF (THE(C,2)< THER(C,2)+t5) THEN
                    t1 = (THER(C,2)+t5)-THE(C,2)
                    FrozenInfPtr(C) = t1*DZ(C,2)
                    !THE(C,2)= THE(C,2)+t1
                    THE(C,2)= THE(C,2)+t1
                    hh(C,1) = hh(C,1) - FrozenInfPtr(C)
                ELSEIF (THE(C,2)> t4) THEN
                    t1 = THE(C,2) - t4
                    THE(C,2) = THE(C,2)-t1
                    FrozenInfPtr(C) = - t1*DZ(C,2) ! Thaw Exfiltration
                    hh(C,1) = hh(C,1) - FrozenInfPtr(C)
                ENDIF
            ELSEIF (snl(C) .eq. 0 .and. snl_old(C) <0 .and. liqChange+iceChange>0D0 .and.    &
                &     (liq0tempPtr(C)+ice0tempPtr(C))>0D0 .and. (liq0tempPtr(C)+ice0tempPtr(C))<50D0) THEN
                ! in snowWater, the 0-th layer is directly absorbed into the 1st
                ! just remove the water that was from h(1)
                iceMovePTR(C) = (iceChange+liqChange)/1000D0 ! both ice and water
                THE(C,2) = THE(C,2) - (h3(C)-ReserveH1PTR(C)) /DZ(C,2)
            ELSE
                call display("passV(1): assumption of totChange=0 broken--C: ", C)
                !THE(C,2) = THE2tempPtr(C)
                !call ch_assert("passV(1): assumption of totChange=0 broken")
            ENDIF
            lastTHEptr(C,:) = THE(C,:) ! save this for better determining the rresis next time step

            ! check
            if (debugg .AND. passVdbg) then

                !    ptrp => clm3%g%l%c%p%pwf%qflx_evap_can
                !    ptrc => clm3%g%l%c%cwf%pwf_a%qflx_evap_can
                !    call p2c (num_nolakec, filter_nolakec, ptrp, ptrc)

                the_sum(C) = sum(THE(C,2:w%m(3))*DZ(C,2:w%m(3)))
                h2osum(C)  = sum(h2osoi_liq(C,1:nlevgrnd))
                icesum(C)  = sum(h2osoi_ice(C,1:nlevgrnd))
                if (abs(h2osum(C)+icesum(C) -(h2osum0(C)+icesum0(C)))>1D-10) then
                    !call display('pass dir==1:', h2osum(C)+icesum(C) -(h2osum0(C)+icesum0(C)))
                endif
                t1 = sum(VS(2:nlevsoi+1)*DZ(C,2:nlevsoi+1))*dt
                t1 = t1+VS(1)*dt
                t2 = abs(t1 - (Gprcp(C)- ET(C) + snoLiqEvap + EvapCS(C)) )
                wtcolVeg = 0D0
                nonUrbSnow = 0 ! will save snowgrnd
                do K= 1,nt
                    if (PFT(C,K)<=numpft) wtcolVeg = wtcolVeg+wtcol(C,K)
                    if (PFT(C,K) .ne. g%Veg%LUidx(1)) nonUrbSnow= nonUrbSnow + qflx_snow_grnd(C,K) * wtCol(C,K)
                enddo
                prcpCan=(clm_a2l%forc_snow(C)+clm_a2l%forc_rain(C))*(1D0-Imp(C)) - cc%cwf%pwf_a%qflx_rain_grnd(C) - nonUrbSnow
                t3 = (h2ocan(C) - cs0(C))/1000D0 - (prcpCan*3.6D0 -EvapCS(C) )


                if (snl_old(c)<0) then
                    Qout = cc%cwf%qflx_top_soil(C)
                    t4   = clm3%g%l%c%cwf%pwf_a%qflx_rain_grnd(C)
                else
                    Qout = cc%cwf%qflx_snomelt(C)
                    t4   = 0D0  ! If no snow, the snoLiqEvap is in fact applied as a source term in EvapG. It has been
                endif

                !t4 = (h2ocan(C) - cs0(C))/1000D0 - (clm_a2l%rainf(C)*(1D0-Imp(C))*conv -Gprcp(C)-EvapCS(C) )
                t5 = (SWE(C) - swe0(C)) - (clm3%g%l%c%cwf%pwf_a%qflx_snow_grnd(C)+t4- Qout)*conv + SNOM(C) + snoLiqEvap + iceMovePTR(C) + ComboPTR(C)
                if (abs(t2)>1D-10 .or. abs(t3)>1D-10 .or. abs(t5)>1D-10) then
                    call display('pass dir==1: mass balance error')
                endif


                t1 = hh(C,1) + THE(C,2)*DZ(C,2) - (h2osoi_liq(C,1)/1000D0 + ReserveH1PTR(C) + THER(C,2)*DZ(C,2) + dhSurf)
                IF (abs(t1)>1D-12) then
                    call display("passV(1): MBE during fix")
                ENDIF
                IF (THE(C,2)<THER(C,2)) THEN
                    call display("passV(1): negative THE(2) generated")
                ENDIF

                cs0(C)      = cc%cws%pws_a%h2ocan(C)
                swe0(C)     = h2osno(C)/1000D0
            endif

            !      IF (THE(C,2)>THES(C,2)) THEN
            !        hh(C,1) = hh(C,1)+ (THE(C,2)-THES(C,2))*DZ(C,2)
            !        Gprcp(C)= Gprcp(C)+ (THE(C,2)-THES(C,2))*DZ(C,2)
            !        THE(C,2)= THES(C,2)
            !      ENDIF
            if (hh(C,1)<0) then
                if (hh(C,1)>-1D10) then
                    hh(C,1)=0D0 ! round off
                else
                    call display( 'pass dir==1: how could this be...?')
                endif
            endif

#if (defined CLM_Only)      
            ! then runoff and infiltration is also computed by CLM
            Rf(C) = qflx_surf(C) * conv
            Df(C) = qflx_infl(c) * conv
            Dperc(C) = qflx_drain(C) * conv
#endif

            if (any(isnan(VS))) THEN
                call display('column: ',C)
                call CH_ASSERT('NaN VS in passV')
            ENDif
            ! test output region
            !      Output(1:15) = (/w%t,THE(C,2),hh(C,1),cc%p%pps%laisun(C),pp%pps%tlai(C),      &
            !     &   pp%pcs%leafc(C),pp%pcf%psnsun(C),cc%cws%pws_a%h2ocan(C),cc%cws%h2osoi_ice(C,1),pp%pns%leafn(C),  &
            !     &   pp%pes%t_veg(C),pp%pps%btran(C),cc%ces%t_grnd(C),cc%ces%t_soisno(C,3),pp%pps%rssun(C)/)
            !       Output(1:9) = (/w%t, ET(C), PET(C), clm_a2l%forc_lwrad(C),clm_a2l%forc_solad(C,1), &
            !      & clm_a2l%forc_solai(C,1),clm_a2l%forc_solad(C,2),clm_a2l%forc_solai(C,2),clm3%g%l%c%ces%forc_t(c)/)
            !       WRITE(222,FMT='(15E16.8)') Output(1:9)
            if (writeout .AND. abs(w%t - int(w%t+0.5))<100 ) then
                WRITE(int(C+10000),FMT='(16E16.8)') w%t,cc%ces%t_soisno(C,1:15)
            endif
        ENDDO

        DO I = 1,filter(1)%num_snowc
            C = filter(1)%snowc(I)
            snl_old(C) = 0 ! put it back
        ENDDO
        !22222222222222222222222222222222222222222222222222222222222222222222222
    else
        call endrun('incorrect passing direction')
    endif


    if ( any( isnan(THE) ) .or. any( isnan(VS) ) ) then

        return
    endif

    end subroutine

    ! ##########################################################################################
    subroutine writeCLMrec(nc,filter_c,cid,np,filter_p,pid)
    ! in this subroutine, a prescribed (in here) variables are written into the file with Unit id cid/pid.
    ! each row contains the record for each index in the list (filter_c, filter_p)
    ! the first item is always time, the second row is the index, followed by output items
    ! this NOT ONLY output clm variables, it may also output PAWS variables for convenience
    use clmtype
    use vdata      , Only: w, g
    use memoryChunk, Only: mn, nn
    implicit none
    integer nc,np,cid,pid
    integer filter_c(nc), filter_p(np),i,j
    integer, parameter :: nitems = 33, nitemsp = 12
    real*8 output(nitems+2)
    real*8 outputP(nitemsp+3)
    type(vege_type),pointer::V
    type(vdz_type),pointer:: VD
    type(pft_type),pointer ::pp
    type(column_type),pointer ::cc
    real*8, pointer,dimension(:) :: c1=>NULL(),c2=>NULL(),c3=>NULL(),c4=>NULL()
    real*8, pointer,dimension(:) :: c5=>NULL(),c6=>NULL(),c7=>NULL(),c8=>NULL()
    real*8, pointer,dimension(:,:) :: p1=>NULL(),p2=>NULL(),p3=>NULL(),p4=>NULL()
    real*8, pointer,dimension(:,:) :: p5=>NULL(),dz=>NULL(),E=>NULL(),THE=>NULL()
    integer n,c,p,k1,k2
    cc      => clm3%g%l%c
    pp      => cc%p
    V       => g%Veg
    VD      => g%VDZ

    c1      => oneD_ptr(mn,g%GW(1)%h)
    c2      => oneD_ptr(mn,g%GW(2)%h)
    c3      => oneD_ptr(mn,g%VDZ%DF)
    c4      => oneD_ptr(mn,g%VDZ%DPerc)
    c5      => oneD_ptr(mn,g%GW(1)%DR)
    c6      => oneD_ptr(mn,g%Veg%Gprcp)
    c7      => oneD_ptr(mn,g%OVN%Rf)
    c8      => oneD_ptr(mn,g%Veg%Tp)
    dz      => twoD_ptr(mn,VD%m(3),g%VDZ%DZ)
    E       => twoD_ptr(mn,VD%m(3),g%VDZ%E)
    THE     => twoD_ptr(mn,VD%m(3),g%VDZ%THE)

    n = nitems + 2
    do i = 1,nc
        c = filter_c(i)
        output(1)  = w%t
        output(2)  = c
        k1 = 1; k2 = 1
        if (E(c,1)-E(c,2)>0.25D0) then
        else
            do j=2,VD%m(3)-1
                if (E(c,1)-E(c,j)<=0.15D0 .and. E(c,1)-E(c,j+1)>0.15D0) then
                    k1 = j
                endif
                if (E(c,1)-E(c,j)<=0.45D0 .and. E(c,1)-E(c,j+1)>0.45D0) then
                    k2 = j
                    exit
                endif
            enddo
        endif
        output(3:n) = (/cc%ces%t_soi_10cm(c)-273.15D0,(cc%cws%h2osoi_liqice_10cm(c)/1000D0)/0.1D0,THE(C,k1),  & ! 3 - 5
            &              THE(C,k2),c1(c),c2(c),c3(c),c4(c),             & ! 6 - 10
            &              c5(c),c6(c),c7(c),c8(c),pp%pps%btran(c),       & ! 11 - 15
            &              THE(C,2:21)/) ! note btran is define on pft level ! 16 - 35
        WRITE(cid,FMT='(100E16.8)') Output(1:n)
    enddo

    n = nitemsP + 3
    do i = 1,np
        p = filter_p(i)
        c = pp%column(p)
        outputP(1)  = w%t
        outputP(2)  = p
        outputP(3)  = pp%wtcol(P)
        OutputP(4:n) = (/cc%p%pps%laisun(p),pp%pps%tlai(p),      &     ! 1 --> 5 (count includes %t and c)
            &   pp%pcf%npp(p),pp%pcs%leafc(p),pp%pns%leafn(p),pp%pcf%fpsn(p),pp%pcs%frootc(P),  &   ! 6 --> 10
            &   pp%pps%btran(p),pp%pef%fsa(p),pp%pef%fsr(p),pp%pes%t_veg(p),pp%pps%rssun(p)/)   ! 11 --> 15
        WRITE(pid,FMT='(100E16.8)') OutputP(1:n)
    enddo



    end subroutine
    ! ##########################################################################################
    subroutine paws_par_transfer
    use clm_varpar, Only: paws_clm_par
    if (size(g%VDZ%VParm)>=14) paws_clm_par(1) = paws_clm_par(1) + g%VDZ%VParm(14) ! multiplier for N fixation
    if (size(g%VDZ%VParm)>=15) paws_clm_par(2) = paws_clm_par(2) + g%VDZ%VParm(15)
    if (size(g%VDZ%VParm)>=16) paws_clm_par(3) = paws_clm_par(3) + g%VDZ%VParm(16)
    if (size(g%VDZ%VParm)>=17) paws_clm_par(4) = paws_clm_par(4) + g%VDZ%VParm(17)

    end subroutine
    ! ##########################################################################################
    subroutine bridgeVegPointers(m)
    ! Establish some interrelationships between VEGE and CLM datatypes
    ! link a handful of Veg new pointers to CLM data chunks
    ! This subroutine only links, NO allocation!
    !USE VData
    USE clm_atmlnd
    use clm_varpar  , only : maxpatch_pft, nlevsno, nlevgrnd, numrad, nlevlak, &
        numpft, ndst, nvoc, nlevurb, nlevsoi
    implicit none
    Type(Vege_type),POINTER::V
    type(pft_type),pointer::pp
    type(column_type),pointer::cc
    type(gridcell_type),pointer::gg
    type(landunit_type),pointer::ll
    integer i,j,k,m(4),P,L,m2(3)
    integer kk,kp
    integer, pointer :: iP(:)=>NULL(),iI(:)=>NULL(),iJ(:)=>NULL()
    m2 = (/m(1), m(2), m(4)/)
    V  =>    g%Veg
    gg =>    clm3%g
    ll =>    gg%l
    cc =>    ll%c
    pp =>    cc%p

    ! pointbacks (CLM sharing PAWS data)
    iP          => oneD_ptrI(nn,V%RPT)
    pp%itype    = iP ! changed. no longer pointing, because I need to change some values
    do i=1,nn
        if (pp%itype(i) > numpft) then
            pp%itype(i) = 0  ! urban heat fluxes currently modeled as bare ground
        endif
    enddo
    pp%pps%tlai => oneD_ptr(nn,V%LAI)
    !pp%pps%hbot => oneD_ptr(nn,V%Root) ! hbot is not root
    pp%pps%htop => oneD_ptr(nn,V%hc)
    pp%wtgcell   => oneD_ptr(nn,V%FRAC)
    pp%wtlunit   => oneD_ptr(nn,V%FRAC)
    pp%wtcol     => oneD_ptr(nn,V%FRAC)

    end subroutine

    ! ##########################################################################################
    SUBROUTINE PAWS_Veg_Initial(n,Lidx)
    ! This subroutine evaluates the initial state of vegetation according to prescribed LAI,canopy height, etc
    ! It also evaluates flags
    ! wrappers include cps, pps, pepv
    USE pftvarcon
    USE shr_const_mod
    USE clm_time_manager, Only: JD_PAWS
    implicit none
    INTEGER MONTH, DATE, n
    INTEGER Lidx(n)
    INTEGER JD
    INTEGER I,J,K,P,L,C,ID,T
    type(pft_epv_type),POINTER::pepv
    type(pft_pstate_type),POINTER::pps
    type(pft_type),POINTER::pp
    type(pft_epc_type),pointer::pc
    type(pft_cstate_type),pointer::pcs
    type(pft_nstate_type),pointer::pns
    TYPE(VEGE_type),POINTER::V
    TYPE(VDZ_type),POINTER::VD
    REAL*8,DIMENSION(:),POINTER::LAI=>NULL(),hc=>NULL(),RMX=>NULL(),annavg_t2m=>NULL()
    INTEGER,DIMENSION(:),POINTER::ivt=>NULL(),smap=>NULL()
    REAL*8 taper,stocking,m1,m2,TMAX,TMIN,TM
    REAL*8,DIMENSION(mn):: annavg
    REAL*8,DIMENSION(:,:),POINTER :: EB=>NULL()
    INTEGER METHOD

    JD = JD_PAWS

    pp   => clm3%g%l%c%p
    pepv => pp%pepv
    pps  => pp%pps
    pcs  => pp%pcs
    pns  => pp%pns
    pc   => pftcon

    V    => g%Veg
    VD   => g%VDZ
    EB   => twoD_ptr(mn,VD%m(3),VD%EB)
    ivt  => pp%itype
    LAI  => pp%pps%tlai
    annavg_t2m   => pp%pepv%annavg_t2m
    smap => oneD_ptrI(mn,g%Veg%SMAP)
    pepv%dormant_flag = 0D0

    DO I = 1,g%wea%day%ns
        annavg(I) = 0D0
        K = 0
        J = 0
        DO WHILE (K < 365)
            J = J + 1
            TMAX = w%Stations(I)%TMAX(J)
            TMIN = w%Stations(I)%TMIN(J)
            IF (TMAX>-90D0 .AND. TMAX<80D0 .AND. TMIN>-90D0 .AND. TMIN<80D0) THEN
                annavg(I) = annavg(I) + (TMAX + TMIN)*(0.5D0/365D0)
                K = K + 1
            ENDIF
        ENDDO
    ENDDO
    ! constant allometric parameters
    taper = 200._r8
    stocking = 1000._r8

    ! convert from stems/ha -> stems/m^2
    stocking = stocking / 10000._r8


    DO I = 1,n
        P = Lidx(I)
        T = ivt(P)
        C = pp%column(P)
        L = smap(C)

        pcs%leafc(P)           = 0D0
        pcs%livestemc(p)       = 0D0
        pcs%livestemc_xfer(p)  = 0D0
        pcs%deadstemc(p)       = 0D0
        pcs%deadstemc_xfer(p)  = 0D0
        pcs%livecrootc(p)      = 0D0
        pcs%livecrootc_xfer(p) = 0D0
        pcs%deadcrootc(p)      = 0D0
        pcs%deadcrootc_xfer(p) = 0D0
        pcs%xsmrpool(p)        = 0D0
        pcs%cpool(p)           = 0D0
        clm3%g%l%c%p%pps%forc_hgt_u_pft(P) = 30D0  ! Is this really used in the growth model?

        METHOD  = 1
        IF (METHOD .eq. 1) THEN
            !111111111111111111111111111111111111111111111111111
            if (T>0 .AND. T<numpft) then
                ! Vegetation
                ! evaluation: leafc,deadstemc,livestemc,livecrootc,deadcrootc,frootc,leafcmax (...storage)

                IF (JD<60) THEN
                    pepv%dormant_flag(P) = 1
                    pepv%onset_flag(P)   = 0D0
                    pepv%offset_flag(P)  = 0D0
                ELSEIF ((JD>300) .AND. (LAI(p) >0D0)) THEN
                    pepv%dormant_flag(P) = 0D0
                    pepv%onset_flag(P)   = 0D0
                    pepv%offset_flag(P)  = 0D0
                ELSEIF (LAI(P) <= 1D-10) THEN
                    ! if still no leaf, should grow!
                    ! since at this point CNVegStructUpdate has not been applied, LAI may still be 0
                    ! this causes a series of chain reactions which made it impossible for this species to grow
                    pepv%dormant_flag(P) = 0D0
                    pepv%onset_flag(P)   = 0D0
                    !pepv%onset_flag(P)   = 1D0 ! this is a bug
                    pepv%offset_flag(P)  = 0D0
                ELSE
                    pepv%dormant_flag(P) = 0D0
                    pepv%onset_flag(P)   = 0D0
                    pepv%offset_flag(P)  = 0D0
                ENDIF
                ! Initialize leafc and htop/hbot from prescribed values
                IF (pc%dsladlai(T) > 0D0) THEN
                    pcs%leafc(P) = log(LAI(P)*pc%dsladlai(T)/pc%slatop(T)+1D0)/pc%dsladlai(T)
                ELSE
                    pcs%leafc(P) = LAI(P)/pc%slatop(T)
                ENDIF

                if (pc%woody(T) == 1._r8) then
                    if (ivt(p) >= nbrdlf_evr_shrub .and. ivt(p) <= nbrdlf_dcd_brl_shrub) then
                        taper = 10._r8
                        ! otherwise have a tall taper
                    else
                        taper = 200._r8
                    endif
                    !htop(p) = ((3._r8 * deadstemc(p) * taper * taper)/               &
                    !&     (SHR_CONST_PI * stocking * dwood(ivt(p))))**(1._r8/3._r8)
                    ! original equation CLM is suspected to use tooooo large taper
                    ! for a 6m tree, the deadstemc is only 144 gC/m2! -- 1.44 kgC dry weight per stem
                    pcs%deadstemc(P) = pc%dwood(T)*pps%htop(P)**3D0*SHR_CONST_PI*      &
                        &        stocking/(3D0*taper * taper)
                    !CP, the below evaluation data are from Baskerville, Autumn 1965, Ecology
                    ! livecrootc + deadcrootc + frootc = 23%
                    ! livestemc  + deadstemc  + leafc  = 77%
                    ! also, assume livestem:deadsteam = 1:1
                    ! livecroot:deadcroot:froot = 1:1:1,
                    ! current:storage = 1:1
                    pcs%livestemc(P) = pcs%deadstemc(P) * (1D0/5D0) ! assuming 1/5 livewood:deadwood ratio
                    m1 = pcs%livestemc(P) + pcs%deadstemc(P) + pcs%leafc(P)
                    m2 = m1 * (0.23D0/0.77D0)
                    pcs%livecrootc(P) = m2*(1D0/3D0)
                    pcs%deadcrootc(P) = m2*(1D0/3D0)
                    pcs%frootc(P)     = m2*(1D0/3D0)

                    pcs%leafc_storage(P)     = pcs%leafc(P)
                    pcs%livestemc_storage(P) = pcs%livestemc(P)
                    pcs%deadstemc_storage(P) = pcs%deadstemc(P)
                    pcs%frootc_storage(P)    = pcs%frootc(P)
                    pcs%livecrootc_storage(P) = pcs%livecrootc(P)
                    pcs%deadcrootc_storage(P) = pcs%deadcrootc(P)

                else
                    pcs%leafc_storage(P) = pcs%leafc(P)
                    pcs%frootc(P) = pcs%leafc(P)
                    pcs%frootc_storage(P) = pcs%frootc(P)
                endif

                pns%leafn(P) = pcs%leafc(P) / pc%leafcn(T)
                pns%leafn_storage(P) = pcs%leafc_storage(P) / pc%leafcn(T)
                pns%frootn(P) = pcs%frootc(P) / pc%frootcn(T)
                pns%frootn_storage(P) = pcs%frootc_storage(P) / pc%frootcn(T)
                pns%livestemn(P) = pcs%livestemc(P) / pc%livewdcn(T)
                pns%livestemn_storage(P) = pcs%livestemc_storage(P) / pc%livewdcn(T)
                pns%livecrootn(P) = pcs%livecrootc(P) / pc%livewdcn(T)
                pns%livecrootn_storage(P) = pcs%livecrootc_storage(P) / pc%livewdcn(T)
                pns%deadstemn(P) = pcs%deadstemc(P) / pc%deadwdcn(T)
                pns%deadstemn_storage(P) = pcs%deadstemc_storage(P) / pc%deadwdcn(T)
                pns%deadcrootn(P) = pcs%deadcrootc(P) / pc%deadwdcn(T)
                pns%deadcrootn_storage(P) = pcs%deadcrootc_storage(P) / pc%deadwdcn(T)

                if (isnan(pns%leafn(p))) pns%leafn(p) = 0D0
                if (isnan(pns%leafn_storage(p))) pns%leafn_storage(p) = 0D0
                if (isnan(pns%leafn_xfer(p))) pns%leafn_xfer(p) = 0D0
                if (isnan(pns%frootn(p))) pns%frootn(p) = 0D0
                if (isnan(pns%frootn_storage(p))) pns%frootn_storage(p) = 0D0
                if (isnan(pns%frootn_xfer(p))) pns%frootn_xfer(p) = 0D0
                if (isnan(pns%livestemn(p))) pns%livestemn(p) = 0D0
                if (isnan(pns%livestemn_storage(p))) pns%livestemn_storage(p) = 0D0
                if (isnan(pns%livestemn_xfer(p))) pns%livestemn_xfer(p) = 0D0
                if (isnan(pns%deadstemn(p))) pns%deadstemn(p) = 0D0
                if (isnan(pns%deadstemn_storage(p))) pns%deadstemn_storage(p) = 0D0
                if (isnan(pns%deadstemn_xfer(p))) pns%deadstemn_xfer(p) = 0D0
                if (isnan(pns%livecrootn(p))) pns%livecrootn(p) = 0D0
                if (isnan(pns%livecrootn_storage(p))) pns%livecrootn_storage(p) = 0D0
                if (isnan(pns%livecrootn_xfer(p))) pns%livecrootn_xfer(p) = 0D0
                if (isnan(pns%deadcrootn(p))) pns%deadcrootn(p) = 0D0
                if (isnan(pns%deadcrootn_storage(p))) pns%deadcrootn_storage(p) = 0D0
                if (isnan(pns%deadcrootn_xfer(p))) pns%deadcrootn_xfer(p) = 0D0
                if (isnan(pns%retransn(p))) pns%retransn(p) = 0D0
                if (isnan(pns%npool(p))) pns%npool(p) = 0D0
                if (isnan(pns%pft_ntrunc(p))) pns%pft_ntrunc(p) = 0D0
                if (isnan(pns%dispvegn(p))) pns%dispvegn(p) = 0D0
                if (isnan(pns%storvegn(p))) pns%storvegn(p) = 0D0
                if (isnan(pns%totvegn(p))) pns%totvegn(p) = 0D0
                if (isnan(pns%totpftn(p))) pns%totpftn(p) = 0D0
            else
            endif
        ELSE
            !111111111111111111111111111111111111111111111111111
            ! Tape Reading method to initialize carbon pools
            ! In this mode, carbon pools are recorded on Julian days throughout a year on a file, one record per day
            ! on the columns are the pools, first line is the header
            ! We are simple reading into this tape


        ENDIF
        !annavg_t2m(P) = annavg(L) + 273.15D0
        DO K = 2,min(VD%m(3),nlevsoi)
            IF (EB(C,1)-EB(C,K) > 0.1D0) EXIT
        ENDDO
        pp%pepv%spheno_lyr(P) = max(1,K - 1 - 1) ! e.g. if the 3rd layer is the one to use, it will exit at layer K = 4, correspond to layer 2 in clm
        ! 0.1D0 in the criteria comes from orignal CLM layer 3 setting. it may change to a more sensible setting later
    ENDDO

    END SUBROUTINE
    !#######################################################################################

    subroutine col_state_init
    ! after arbitrary initialization of temperature (mkarbinit), do a meaningful initialization of water/energy states
    ! use current air temperature as surface temperature, linearly interpolate to the bottom, which is set at 10 C
    USE shr_const_mod, Only: tkfrz => SHR_CONST_TKFRZ
    USE clm_time_manager, Only: JD_PAWS
    USE filtermod, Only: filter
    USE CNSetValueMod, Only: CNSetCcs,CNSetCns
    use memoryChunk, Only: cLevelMap
    implicit none
    integer i,j,k,c,L,JJ,kk,eof,nl,ncLevel,nk
    real*8,dimension(:,:),pointer:: sminn2=>NULL(),t_soisno2=>NULL()
    real*8,dimension(:,:,:),pointer:: t_soisno=>NULL()
    type(column_type),pointer::cc
    type(sta_type),pointer:: day
    type(vege_type),pointer:: V
    type(vdz_type),pointer::VD
    type(pft_type),pointer::pp
    type(pft_cstate_type),pointer::pcs
    type(pft_nstate_type),pointer::pns
    type(column_cstate_type),pointer::ccs
    type(column_nstate_type),pointer::cns
    type(pft_epv_type),pointer:: pepv

    real*8 dTemp,dE,dTdE, rN
    real*8:: STEMP(20),Z(20),ZZ(100)
    real*8,dimension(:,:),pointer:: z2=>NULL()
    integer strategy,uid,nc,p,nz,nq
    character(256) inifile
    integer,pointer :: ivt(:)=>NULL()
    real*8,dimension(:),pointer::wtcol=>NULL()
    integer,dimension(:),pointer::cLevelMapL=>NULL()

    real*8,dimension(0:numpft):: leafc,leafc_storage,leafc_xfer,frootc,frootc_storage,frootc_xfer, &
        & livestemc,livestemc_storage,livestemc_xfer,deadstemc,deadstemc_storage,deadstemc_xfer,&
        &  livecrootc,livecrootc_storage,livecrootc_xfer,deadcrootc,deadcrootc_storage,  &
        &  deadcrootc_xfer,gresp_storage,gresp_xfer,cpool,xsmrpool,pft_ctrunc,dispvegc, &
        &  storvegc,totvegc,totpftc,leafcmax,leafn,leafn_storage,leafn_xfer,frootn, &
        &  frootn_storage,frootn_xfer,livestemn,livestemn_storage,livestemn_xfer,deadstemn, &
        &  deadstemn_storage,deadstemn_xfer,livecrootn,livecrootn_storage, livecrootn_xfer,deadcrootn,&
        &  deadcrootn_storage,deadcrootn_xfer,retransn,npool,pft_ntrunc,dispvegn,storvegn,totvegn,&
        & totpftn,soilc,cwdc,litr1c,litr2c,litr3c,soil1c,soil2c,soil3c,soil4c,seedc,&
        & col_ctrunc,prod10c,prod100c,totprodc,totlitc,totsomc,totecosysc,totcolc, cwdn,&
        & litr1n,litr2n,litr3n,soil1n,soil2n,soil3n,soil4n,sminn,col_ntrunc,seedn,prod10n,&
        & prod100n,totprodn,totlitn,totsomn,totecosysn,totcoln, &
        & annsum_npp,annsum_potential_gpp,annavg_t2m,annmax_retransn



    strategy = 1
    nc   = 1

    cc   => clm3%g%l%c
    pp   => cc%p
    day  => g%Wea%Day
    V    => g%Veg
    VD   => g%VDZ
    sminn2    => twoD_ptr(w%m(1),w%m(2),cc%cns%sminn)
    t_soisno => threeD_ptr(w%m,cc%ces%t_soisno)
    ivt      => cc%p%itype
    wtcol    => pp%wtcol
    nq = min(int(JD_PAWS) - 1,364)

    if (strategy .eq. 0) then
        do i = 1,w%m(1)
            do j = 1,w%m(2)
                L = V%SMAP(I,J)
                ! soil temperature
                if (w%DM%mask(I,J)) then
                    t_soisno(i,j,1) = day%temp(8,L)+tkfrz
                    t_soisno(i,j,nlevgrnd) = 10D0 + tkfrz
                    dTemp = t_soisno(i,j,nlevgrnd) - t_soisno(i,j,1)
                    dE    = VD%E(i,j,nlevgrnd) - VD%E(i,j,1)
                    dTdE  = dTemp/dE
                    do k = 2,nlevgrnd-1
                        t_soisno(i,j,k) = (VD%E(i,j,k)-VD%E(i,j,1))*dTdE + t_soisno(i,j,1)
                    enddo

                    sminn2(i,j) = 300D0
                endif
            enddo
        enddo

    elseif (strategy .eq. 1) then
        ! tape reading initialization!
        ! run the model for representative column for long time, after it converges, get the average state variables in a year
        ! they are written in a file, this file is read here and values assigned to state variables in the model according to the start julian day
        inifile = 'initialCStates.dat'

        pcs => pp%pcs
        pns => pp%pns
        pepv=> pp%pepv
        ccs => cc%ccs
        cns => cc%cns

        uid = 55
        open(uid, file = inifile)

        ! count how many lines there are
        eof = 0
        i = 0
        do while(eof.ge.0)
            i=i+1
            read(uid,*,iostat=eof)
        end do
        close(uid)
        nl = i-1

        rN = (nl - 1)/365D0/numpft
        if (abs(rN - int(rN))>1D-10) call ch_assert('incorrect number of lines in initialCStates.dat')
        ncLevel = rN
        cLevelMapL=> oneD_ptrI(mn, cLevelMap)

        open(uid, file = inifile)
        read(uid, *) ! title line

        DO k = 1,ncLevel ! ncLevel loop

            nk = 0
            DO I = 1, nq
                DO J = 1, numpft
                    read(uid,*) ! skip
                ENDDO
            ENDDO

            ! all carbon pools
            DO J = 1, numpft
                read(uid,*) leafc(J),leafc_storage(J),leafc_xfer(J),frootc(J),frootc_storage(J),frootc_xfer(J), &
                    & livestemc(J),livestemc_storage(J),livestemc_xfer(J),deadstemc(J),deadstemc_storage(J),deadstemc_xfer(J),&
                    &  livecrootc(J),livecrootc_storage(J),livecrootc_xfer(J),deadcrootc(J),deadcrootc_storage(J),  &
                    &  deadcrootc_xfer(J),gresp_storage(J),gresp_xfer(J),cpool(J),xsmrpool(J),pft_ctrunc(J),dispvegc(J), &
                    &  storvegc(J),totvegc(J),totpftc(J),leafcmax(J),leafn(J),leafn_storage(J),leafn_xfer(J),frootn(J), &
                    &  frootn_storage(J),frootn_xfer(J),livestemn(J),livestemn_storage(J),livestemn_xfer(J),deadstemn(J), &
                    &  deadstemn_storage(J),deadstemn_xfer(J),livecrootn(J),livecrootn_storage(J), livecrootn_xfer(J),deadcrootn(J),&
                    &  deadcrootn_storage(J),deadcrootn_xfer(J),retransn(J),npool(J),pft_ntrunc(J),dispvegn(J),storvegn(J),totvegn(J),&
                    & totpftn(J),soilc(J),cwdc(J),litr1c(J),litr2c(J),litr3c(J),soil1c(J),soil2c(J),soil3c(J),soil4c(J),seedc(J),&
                    & col_ctrunc(J),prod10c(J),prod100c(J),totprodc(J),totlitc(J),totsomc(J),totecosysc(J),totcolc(J), cwdn(J),&
                    & litr1n(J),litr2n(J),litr3n(J),soil1n(J),soil2n(J),soil3n(J),soil4n(J),sminn(J),col_ntrunc(J),seedn(J),prod10n(J),&
                    & prod100n(J),totprodn(J),totlitn(J),totsomn(J),totecosysn(J),totcoln(J),&
                    & annsum_npp(J),annsum_potential_gpp(J),annavg_t2m(J),annmax_retransn(J)
            ENDDO
            call CNSetCcs(filter(nc)%num_simc,filter(nc)%simc, 0D0, cc%ccs)
            call CNSetCns(filter(nc)%num_simc,filter(nc)%simc, 0D0, cc%cns)

            DO I = 1,filter(nc)%num_soilp
                p = filter(nc)%soilp(I)
                c = pp%column(p)

                if (cLevelMapL(c) .eq. k) then

                    nk = nk + 1
                    J = ivt(p)
                    !Chaopeng: these codes below were generated by my smart program-writing program !CP
                    pcs%leafc(p)=leafc(J)
                    pcs%leafc_storage(p)=leafc_storage(J)
                    pcs%leafc_xfer(p)=leafc_xfer(J)
                    pcs%frootc(p)=frootc(J)
                    pcs%frootc_storage(p)=frootc_storage(J)
                    pcs%frootc_xfer(p)=frootc_xfer(J)
                    pcs%livestemc(p)=livestemc(J)
                    pcs%livestemc_storage(p)=livestemc_storage(J)
                    pcs%livestemc_xfer(p)=livestemc_xfer(J)
                    pcs%deadstemc(p)=deadstemc(J)
                    pcs%deadstemc_storage(p)=deadstemc_storage(J)
                    pcs%deadstemc_xfer(p)=deadstemc_xfer(J)
                    pcs%livecrootc(p)=livecrootc(J)
                    pcs%livecrootc_storage(p)=livecrootc_storage(J)
                    pcs%livecrootc_xfer(p)=livecrootc_xfer(J)
                    pcs%deadcrootc(p)=deadcrootc(J)
                    pcs%deadcrootc_storage(p)=deadcrootc_storage(J)
                    pcs%deadcrootc_xfer(p)=deadcrootc_xfer(J)
                    pcs%gresp_storage(p)=gresp_storage(J)
                    pcs%gresp_xfer(p)=gresp_xfer(J)
                    pcs%cpool(p)=cpool(J)
                    pcs%xsmrpool(p)=xsmrpool(J)
                    pcs%pft_ctrunc(p)=pft_ctrunc(J)
                    pcs%dispvegc(p)=dispvegc(J)
                    pcs%storvegc(p)=storvegc(J)
                    pcs%totvegc(p)=totvegc(J)
                    pcs%totpftc(p)=totpftc(J)
                    pcs%leafcmax(p)=leafcmax(J)


                    pns%leafn(p)=leafn(J)
                    pns%leafn_storage(p)=leafn_storage(J)
                    pns%leafn_xfer(p)=leafn_xfer(J)
                    pns%frootn(p)=frootn(J)
                    pns%frootn_storage(p)=frootn_storage(J)
                    pns%frootn_xfer(p)=frootn_xfer(J)
                    pns%livestemn(p)=livestemn(J)
                    pns%livestemn_storage(p)=livestemn_storage(J)
                    pns%livestemn_xfer(p)=livestemn_xfer(J)
                    pns%deadstemn(p)=deadstemn(J)
                    pns%deadstemn_storage(p)=deadstemn_storage(J)
                    pns%deadstemn_xfer(p)=deadstemn_xfer(J)
                    pns%livecrootn(p)=livecrootn(J)
                    pns%livecrootn_storage(p)=livecrootn_storage(J)
                    pns%livecrootn_xfer(p)=livecrootn_xfer(J)
                    pns%deadcrootn(p)=deadcrootn(J)
                    pns%deadcrootn_storage(p)=deadcrootn_storage(J)
                    pns%deadcrootn_xfer(p)=deadcrootn_xfer(J)
                    pns%retransn(p)=retransn(J)
                    pns%npool(p)=npool(J)
                    pns%pft_ntrunc(p)=pft_ntrunc(J)
                    pns%dispvegn(p)=dispvegn(J)
                    pns%storvegn(p)=storvegn(J)
                    pns%totvegn(p)=totvegn(J)
                    pns%totpftn(p)=totpftn(J)

                    ccs%soilc(c)=ccs%soilc(c)+soilc(J)*wtcol(p)
                    ccs%cwdc(c)=ccs%cwdc(c)+cwdc(J)*wtcol(p)
                    ccs%litr1c(c)=ccs%litr1c(c)+litr1c(J)*wtcol(p)
                    ccs%litr2c(c)=ccs%litr2c(c)+litr2c(J)*wtcol(p)
                    ccs%litr3c(c)=ccs%litr3c(c)+litr3c(J)*wtcol(p)
                    ccs%soil1c(c)=ccs%soil1c(c)+soil1c(J)*wtcol(p)
                    ccs%soil2c(c)=ccs%soil2c(c)+soil2c(J)*wtcol(p)
                    ccs%soil3c(c)=ccs%soil3c(c)+soil3c(J)*wtcol(p)
                    ccs%soil4c(c)=ccs%soil4c(c)+soil4c(J)*wtcol(p)
                    ccs%seedc(c)=ccs%seedc(c)+seedc(J)*wtcol(p)
                    ccs%col_ctrunc(c)=ccs%col_ctrunc(c)+col_ctrunc(J)*wtcol(p)
                    ccs%prod10c(c)=ccs%prod10c(c)+prod10c(J)*wtcol(p)
                    ccs%prod100c(c)=ccs%prod100c(c)+prod100c(J)*wtcol(p)
                    ccs%totprodc(c)=ccs%totprodc(c)+totprodc(J)*wtcol(p)
                    ccs%totlitc(c)=ccs%totlitc(c)+totlitc(J)*wtcol(p)
                    ccs%totsomc(c)=ccs%totsomc(c)+totsomc(J)*wtcol(p)
                    ccs%totecosysc(c)=ccs%totecosysc(c)+totecosysc(J)*wtcol(p)
                    ccs%totcolc(c)=ccs%totcolc(c)+totcolc(J)*wtcol(p)


                    cns%cwdn(c)= cns%cwdn(c)+cwdn(J)*wtcol(p)
                    cns%litr1n(c)= cns%litr1n(c)+litr1n(J)*wtcol(p)
                    cns%litr2n(c)= cns%litr2n(c)+litr2n(J)*wtcol(p)
                    cns%litr3n(c)= cns%litr3n(c)+litr3n(J)*wtcol(p)
                    cns%soil1n(c)= cns%soil1n(c)+soil1n(J)*wtcol(p)
                    cns%soil2n(c)= cns%soil2n(c)+soil2n(J)*wtcol(p)
                    cns%soil3n(c)= cns%soil3n(c)+soil3n(J)*wtcol(p)
                    cns%soil4n(c)= cns%soil4n(c)+soil4n(J)*wtcol(p)
                    cns%sminn(c)= cns%sminn(c)+sminn(J)*wtcol(p)
                    cns%col_ntrunc(c)= cns%col_ntrunc(c)+col_ntrunc(J)*wtcol(p)
                    cns%seedn(c)= cns%seedn(c)+seedn(J)*wtcol(p)
                    cns%prod10n(c)= cns%prod10n(c)+prod10n(J)*wtcol(p)
                    cns%prod100n(c)= cns%prod100n(c)+prod100n(J)*wtcol(p)
                    cns%totprodn(c)= cns%totprodn(c)+totprodn(J)*wtcol(p)
                    cns%totlitn(c)= cns%totlitn(c)+totlitn(J)*wtcol(p)
                    cns%totsomn(c)= cns%totsomn(c)+totsomn(J)*wtcol(p)
                    cns%totecosysn(c)= cns%totecosysn(c)+totecosysn(J)*wtcol(p)
                    cns%totcoln(c)= cns%totcoln(c)+totcoln(J)*wtcol(p)

                    !CP: added 09/29/12 to suppres initial hump
                    pepv%annsum_npp(p) = annsum_npp(J)
                    pepv%annsum_potential_gpp(p) = annsum_potential_gpp(J)
                    pepv%annavg_t2m(p) = annavg_t2m(J)
                    pepv%annmax_retransn(p) = annmax_retransn(J)
                    if (p <= mn) then
                        cc%cps%cannsum_npp(c) = annsum_npp(J)*wtcol(p)
                    else
                        cc%cps%cannsum_npp(c) = cc%cps%cannsum_npp(c) + annsum_npp(J)*wtcol(p)
                    endif

                    !WRITE(777,*) p,pepv%annsum_npp(p),pepv%annsum_potential_gpp(p),pepv%annavg_t2m(p),pepv%annmax_retransn(p),cc%cps%cannsum_npp(c)

                endif
            ENDDO

            ! skipping rest of the lines in the ncLevel
            DO I = 1, 365 - (nq+1)
                DO J = 1, numpft
                    read(uid,*) ! skip
                ENDDO
            ENDDO


            IF (k > 1) THEN
                call display( 'multiple cLevel: ',k, ' # pfts: ' , nk)
            ENDIF

        ENDDO ! ncLevel loop
        close(uid)

        call initialSTemp(nq)
    endif
    end subroutine

    SUBROUTINE initialSTemp(nq)
    USE shr_const_mod, Only: tkfrz => SHR_CONST_TKFRZ
    USE clm_time_manager, Only: JD_PAWS
    USE filtermod, Only: filter
    implicit none
    integer i,j,k,c,L,JJ,kk,eof,nl,ncLevel,nk
    real*8,dimension(:,:),pointer:: sminn2=>NULL(),t_soisno2=>NULL()
    real*8,dimension(:,:,:),pointer:: t_soisno=>NULL()
    type(column_type),pointer::cc
    type(sta_type),pointer:: day
    type(vege_type),pointer:: V
    type(vdz_type),pointer::VD
    type(pft_type),pointer::pp
    type(pft_cstate_type),pointer::pcs
    type(pft_nstate_type),pointer::pns
    type(column_cstate_type),pointer::ccs
    type(column_nstate_type),pointer::cns

    real*8 dTemp,dE,dTdE, rN
    real*8:: STEMP(20),Z(20),ZZ(100)
    real*8,dimension(:,:),pointer:: z2=>NULL()
    integer strategy,uid,nc,p,nz,nq
    character(256) inifile
    integer,pointer :: ivt(:)=>NULL()
    real*8,dimension(:),pointer::wtcol=>NULL()

    VD=> g%VDZ
    nc   = 1

    cc   => clm3%g%l%c
    pp   => cc%p
    open(uid,file = 'initialSTemp.dat')
    read(uid,*) (z(i), i = 1,nlevgrnd)  ! z = 0 at soil surface
    DO I = 1, nq
        read(uid,*) ! skip
    ENDDO
    read(uid,*) (STEMP(i), i = 1,nlevgrnd)
    z2  => twoD_ptr(mn,VD%m(3),VD%E)
    nz  = VD%m(3)
    t_soisno2 => cc%ces%t_soisno
    DO I = 1,filter(nc)%num_soilc
        C = filter(nc)%soilc(I)
        !       DO J = 1,nz
        !         ZZ(J) = z2(C,1)-z2(C,J)
        !       END
        K = 1
        DO J = 2,nlevgrnd+1
            ! interpolate temperature for cell J (in PAWS), therefore cell J-1 in CLM
            JJ = J - 1
            K = insertLoc(nlevgrnd,z,z2(C,J)-z2(C,1),K,1)
            KK = max(1,min(K,nlevgrnd))
            IF (z2(C,J)-z2(C,1)>=z(K) .OR. K>=nlevgrnd .OR. K<=0) THEN
                ! above highest point
                t_soisno2(C,JJ) = STEMP(KK)
            ELSE
                dTemp = STEMP(K+1) - STEMP(K)
                dE    = z(K) - z(K+1)
                dTdE  = dTemp/dE
                t_soisno2(C,JJ) = (z(K)- (z2(C,J) - z2(C,1)))*dTdE + STEMP(KK)
            ENDIF
        ENDDO

    ENDDO

    CLOSE(uid)
    END SUBROUTINE


    !=====================================================================================

    SUBROUTINE SURFACE_URB(m,n,nz,nt,ns,Veg,V,H,dt,OVH,OVF)
    ! To use with CLM, this subroutine clears some flux terms, and compute the evaporation from urban as well as wetland (impervious) covers
    ! This is because currently we have not used CLM urban code, as it is not purely compatible with our discretization (urban is defined on column level)

    USE Vdata
    USE clmtype
    use clm_varcon        , only : hvap
    use memoryChunk, Only: tempStorage
    IMPLICIT NONE
    integer m,n,nt,ns,nz,H ! H is hour of day
    integer,DIMENSION(:,:),POINTER:: nRPT=>NULL()
    integer,dimension(:,:),pointer:: SMAP=>NULL()
    integer,DIMENSION(:,:,:),POINTER::RPT=>NULL()
    REAL*8 dt
    REAL*8,DIMENSION(:,:),POINTER::ET=>NULL(),Cpc=>NULL(),hh=>NULL(),SRC=>NULL()
    REAL*8,DIMENSION(:,:,:),POINTER::FRAC=>NULL(),LAI=>NULL(),ROOT=>NULL(),KC=>NULL()
    REAL*8,DIMENSION(:,:),POINTER::MAPP=>NULL(),Imp=>NULL(),hI=>NULL(),EvapG=>NULL()
    REAL*8,DIMENSION(:,:),POINTER::Evap=>NULL(),Gprcp=>NULL(),PET=>NULL()
    REAL*8,DIMENSION(:,:),POINTER::prcp=>NULL(),rad=>NULL(),wnd=>NULL(),hmd=>NULL()
    REAL*8,DIMENSION(:,:),POINTER::ret=>NULL(),rets=>NULL(),tempe=>NULL(),snowf=>NULL()
    REAL*8,DIMENSION(:,:,:),POINTER::EB=>NULL()
    REAL*8,DIMENSION(:),POINTER::EE=>NULL(),THE=>NULL(),THES=>NULL(),THER=>NULL()
    REAL*8,DIMENSION(:),POINTER::VS=>NULL(),FC=>NULL(),WP=>NULL(),hhh=>NULL()
    REAL*8,DIMENSION(m,n)::OVH,OVF
    TYPE(Vege_type)::Veg
    !TYPE(STA)::DAY(2)
    TYPE(VDZ_type)::V
    TYPE(THERMAL_type),POINTER:: th=>NULL()
    integer I,J,K,L,P,RDN,EDN,ETN,II,RMX
    REAL*8 RPET,PREC,NewS,RPrcp,FPET,EvapF
    REAL*8 PT(nt),RDF(nz),RFrac(nt)
    !REAL*8,DIMENSION(nz):: PET,TP,EVAP,EDF,F
    REAL*8 THRSD, SS, RATIO, TEMP, temp2,temp3,temp4,temp5,temp6
    REAL*8,DIMENSION(:),POINTER:: TEMPP=>NULL(), tempp2=>NULL(),tempp3=>NULL(),tempp4=>NULL()
    integer NR,NCODE(3)
    REAL*8::EU=0D0,EL=0D0,RSnow=0D0,RS=0D0,RadS,Sa,TT,albedo,SNOMP
    REAL*8 ETS,ST,ETA,test,CPET,VA,GAM,OM,test2,OCOV,Newsnow,tt2, conv

    REAL*8,DIMENSION(7) :: inpt
    REAL*8,DIMENSION(:),POINTER:: sitev=>NULL(),statev=>NULL(),outv=>NULL()
    REAL*8 cump, cume ,cummr, input(7), fa,errmbal,w1,evap1
    integer iflag(5),nstepday ! consider change later
    INTEGER, SAVE:: state
    REAL*8, DIMENSION(:,:,:), pointer:: ptrp=>NULL()
    REAL*8, DIMENSION(:,:), pointer:: PrcpI=>NULL(), EvapI=>NULL()


    NR = 2; NCODE = 0; nstepday = 24;
    fa = 1000D0/24D0
    NCODE(1) =1; ! impervious treatment
    NCODE(2) =1; ! snowmelt

    th => g%VDZ%th
    RPT => Veg%RPT
    KC => Veg%Kc
    !CS => Veg%CS
    !CSEvap => Veg%CSEvap
    Cpc => Veg%Cpc
    !DELTA => Veg%Delta
    LAI => Veg%LAI
    Imp => Veg%Imp
    PREC = 1e-9
    !RDFN => Veg%RDFN
    FRAC => Veg%Frac
    hI   => Veg%hI
    Evap   => Veg%Evap
    EvapG  => Veg%EvapG
    ET  => Veg%ET
    !BSoil => Veg%BSoil
    hh => V%h(:,:,1)
    SRC => V%SRC(:,:,1)
    !ETC => Veg%ETC
    SMap => Veg%SMAP
    Gprcp  => Veg%GPrcp
    PET  => Veg%PET

    prcp => g%Wea%Day%prcp
    snowf=> g%Wea%Day%snowf
    ret => g%Wea%Day%ret
    tempe => g%Wea%Day%temp
    hmd => g%Wea%Day%hmd
    rad => g%Wea%Day%rad
    GAM = V%VParm(1)


    iflag = (/0, 0, 3, 1, 4/) ! iflag(1)=0 means inpt(5)=qsi, inpt(6)=qli
    cump = 0D0; cume = 0D0; cummr = 0D0

    EB => V%EB
    MAPP=> V%MAPP
    ETS = 0D0
    !OPEN(UNIT=3,file='rec.dat')

    ! Here is the logic:
    ! In CLM, we don't do either evap or prcp on the impervious covers
    ! hI is updated here. dHI/dt = Prcp - Evap - Runoff (no infiltration, Runoff is later computed in subroutine Runoff)
    ptrp => threeD_ptr( (/m,n,nt/),clm3%g%l%c%p%pef%eflx_lh_tot)  ! latent heat flux for mapOutput
    conv    =  Veg%dt * 86400D0/1000D0 ! from mm/s (CLM) to m/day (PAWS) and then times it by dt

    call tempStorage('PrcpI',PrcpI,1)
    call tempStorage('EvapI',EvapI,1)

    DO I = 1,m
        DO J = 1,n
            IF (MAPP(I,J)>PREC) THEN
                ! ##11111111111111111111111111111111111111111111111111111111111111111111111111
                !        if (i .eq. 26 .and. j .eq. 3) then
                !        write(*,*) 'surface_urb debug'
                !        endif

                L = SMap(I,J)
                Rprcp = (prcp(H,L)-snowf(H,L))*dt
                RPET=ret(H,L)*dt
                !Veg%PET(I,J)=RPET
                ! SNOW?

                ! now treat impervious surfaces
                ! it has its own storage: hI
                ! rule of applying ET: each class working on its own portion,
                ! do not use others, except for shared Canopy Storage (for now)
                VA = 1D0 - Veg%Imp(I,J);tt2=0D0
                PrcpI(I,J) = 0D0
                EvapI(I,J) = 0D0
                IF (Imp(I,J)>0D0) THEN
                    DO K = 1,nt
                        IF (Veg%RPT(I,J,K) .eq. int(Veg%LUidx(1))) THEN ! so that we have K, and an index to ptrp
                            tt2 = RPrcp*Imp(I,J)
                            if (tt2 > 0D0) then
                                hI(I,J) = hI(I,J) + tt2
                                CPET = RPET*Imp(I,J)
                                !SRC(I,J) = SRC(I,J) - tt2
                                Gprcp(I,J) = Gprcp(I,J) + tt2
                                PrcpI(I,J) = PrcpI(I,J) + tt2
                            endif
                            Evap1= max(min(CPET,hI(I,J)),0D0)
                            EvapI(I,J) = EvapI(I,J) + Evap1
                            hI(I,J) = hI(I,J) - Evap1
                            EvapG(I,J) = EvapG(I,J) + Evap1
                            Evap(I,J)  = Evap(I,J)  + Evap1
                            ET(I,J)    = ET(I,J) + Evap1
                            PET(I,J)   = PET(I,J) + RPET*Imp(I,J)

                            ! not necessary:
                            RPET = RPET - CPET
                            Rprcp = Rprcp - tt2
                            ! now both RPET and Rprcp are for vegetated portion
                            IF (Evap1 > 0D0) then
                                ptrp(I,J,K) = hvap * Evap1/conv
                            else
                                ptrp(I,J,K) = 0D0
                            endif

                        ENDIF
                    ENDDO
                ENDIF

                ! CONSIDER THIS:
                ! HOW CAN WE IMPLEMENT THIS WITH CLM? (RESCALING)
                ! now also evaporate from flow domain
                !         if (OVH(I,J)>0D0) then
                !           ! RPET in the flow domain preferrably evaporate the flow water
                !           ! if extra, will come back join rpet
                !           FPET = RPET*OVF(I,J)
                !           EvapF= max(min(OVH(I,J),FPET),0D0)
                !           OVH(I,J)=OVH(I,J)-EvapF
                !           RPET = RPET - FPET
                !           Evap(1)=Evap(1)+EvapF
                !         else
                !           FPET = 0D0
                !         endif
            ENDIF
        ENDDO
    ENDDO

    END SUBROUTINE SURFACE_URB

#endif
    END MODULE
