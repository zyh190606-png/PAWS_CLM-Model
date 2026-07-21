#include <misc.h>
#include <preproc.h>

module clmtypeInitMod
! This module is modified by Chaopeng
! This module is turned into a linking code in stead of a allocation code,-
! in order to use a same patch of memory between mex calls.
! commented (in stead of deleted) allocation statements are those not initiated to NaN
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: clmtypeInitMod
!
! !DESCRIPTION:
! Allocate clmtype components and initialize them to signaling NaN.
!
! !USES:
  
  use shr_kind_mod, only : r8 => shr_kind_r8
  use nanMod      , only : nan, bigint
  use clmtype
  use clm_varpar  , only : maxpatch_pft, nlevsno, nlevgrnd, numrad, nlevlak, &
                           numpft, ndst, nvoc, nlevurb, nlevsoi
  USE memoryChunk , only : R3,R4,L3,L4,mn,nn,VBB,i3,i4,j3,j4,iv,R3C,R4C,L3C,L4C
  USE VData       , only : oneD_Ptr,oneD_PtrI,twoD_Ptr,twoD_ptrI,twoD_ptrF,    &
                          & threeD_ptr,threeD_ptrI,num2str,str2num,             &
                          & oneD_PtrO,oneD_PtrIO,twoD_PtrO,twoD_ptrIO,         &
                          & threeD_ptrO,threeD_ptrIO
  USE abortutils  , only : endrun
!
! !PUBLIC TYPES:
  implicit none
  
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: initClmtype,linkCLMtype
!
! !REVISION HISTORY:
! Created by Peter Thornton and Mariana Vertenstein
! Modified by Colette L. Heald (05/06) for VOC emission factors
! 3/17/08 David Lawrence, changed nlevsoi to nlevgrnd where appropriate
!
! !PRIVATE MEMBER FUNCTIONS:
  private :: init_pft_type
  private :: init_column_type
  private :: init_landunit_type
  private :: init_gridcell_type
  private :: init_energy_balance_type
  private :: init_water_balance_type
  private :: init_pft_ecophys_constants
#if (defined CNDV)
  private :: init_pft_DGVMecophys_constants
#endif
  private :: init_pft_pstate_type
  private :: init_pft_epv_type
#if (defined CNDV)
  private :: init_pft_pdgvstate_type
#endif
  private :: init_pft_vstate_type
  private :: init_pft_estate_type
  private :: init_pft_wstate_type
  private :: init_pft_cstate_type
  private :: init_pft_nstate_type
  private :: init_pft_eflux_type
  private :: init_pft_mflux_type
  private :: init_pft_wflux_type
  private :: init_pft_cflux_type
  private :: init_pft_nflux_type
  private :: init_pft_vflux_type
  private :: init_pft_dflux_type
  private :: init_pft_depvd_type
  private :: init_column_pstate_type
  private :: init_column_estate_type
  private :: init_column_wstate_type
  private :: init_column_cstate_type
  private :: init_column_nstate_type
  private :: init_column_eflux_type
  private :: init_column_wflux_type
  private :: init_column_cflux_type
  private :: init_column_nflux_type
  private :: init_landunit_pstate_type
  private :: init_landunit_eflux_type
  private :: init_gridcell_pstate_type
  private :: init_gridcell_efstate_type
  private :: init_gridcell_wflux_type
!EOP
!----------------------------------------------------
contains
  ! ##########################################################################################
  ! THIS ONE IS PUBLIC
  subroutine linkCLMtype
  USE clmtype, Only: clm3
  USE memoryChunk, Only: mm,i3,i4,j3,j4,R3C,R4C
  USE clm_atmlnd,  Only: clm_a2l
  USE decompMod,   Only: get_proc_bounds
  ! Pointer assignment (to central memroy storage)
  IMPLICIT NONE
  integer m(4) ! [ny,nx,nz,nRPT]
  integer n(4) ! [nd3,nd4,nl3,nl4]
  integer begp,begc,endp,endc,begl,endl,begg,endg
  integer state
  integer,pointer :: ii3,ii4,jj3,jj4
  m  = mm
  mn = m(1)*m(2); nn = m(1)*m(2)*m(4)
  I3=1; I4=1;
  J3=1; J4=1;
  iv = 1;
  call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)
  ii3 => I3
  ii4 => I4
  jj3 => J3
  jj4 => J4
  
  ! These steps are done only once in Fortran code, before simulation starts
  ! If in Matlab mode, these steps are done every time it comes from Matlab
  ! Those codes in section 'else (state = 0)' will be done only once in either mode
  call init_pft_type     (begp, endp, clm3%g%l%c%p)
  call init_column_type  (begc, endc, clm3%g%l%c)
  call init_landunit_type(begl, endl, clm3%g%l)
  call init_gridcell_type(begg, endg, clm3%g)
  
  call init_column_cstate_type(begc, endc, clm3%g%l%c%ccs)
  call init_column_cflux_type(begc, endc, clm3%g%l%c%ccf)
  call init_column_nstate_type(begc, endc, clm3%g%l%c%cns)
  call init_column_nflux_type(begc, endc, clm3%g%l%c%cnf)
  call init_column_pstate_type(begc, endc, clm3%g%l%c%cps) ! This is a memory guzzler
  call init_column_estate_type(begc, endc, clm3%g%l%c%ces)
  call init_column_eflux_type(begc, endc, clm3%g%l%c%cef)
  call init_column_wstate_type(begc, endc, clm3%g%l%c%cws)
  call init_column_wflux_type(begc, endc, clm3%g%l%c%cwf)
    
  call init_pft_nstate_type(begp, endp, clm3%g%l%c%p%pns)
  call init_pft_nflux_type(begp, endp, clm3%g%l%c%p%pnf)
  call init_pft_cstate_type(begp, endp, clm3%g%l%c%p%pcs)
  call init_pft_cflux_type(begp, endp, clm3%g%l%c%p%pcf)
  call init_pft_wstate_type(begp, endp, clm3%g%l%c%p%pws)
  call init_pft_wflux_type(begp, endp, clm3%g%l%c%p%pwf)
  call init_pft_estate_type(begp, endp, clm3%g%l%c%p%pes)
  call init_pft_eflux_type(begp, endp, clm3%g%l%c%p%pef)
  call init_pft_mflux_type(begp, endp, clm3%g%l%c%p%pmf)
  call init_pft_pstate_type(begp, endp, clm3%g%l%c%p%pps)
  call init_pft_epv_type(begp, endp, clm3%g%l%c%p%pepv)
  
  call init_atm2lnd_type(begp, endp, clm_a2l)
  
  call init_pft_ecophys_constants()
  ! verified that ii3,ii4,ii5,ii6 sums up correctly as number of items in each subroutine on Oct 3 2010
  
  call init_col_avg(begc, endc, clm3%g%l%c%ccf%pcf_a,clm3%g%l%c%cnf%pnf_a,clm3%g%l%c%ccs%pcs_a, &
               clm3%g%l%c%cns%pns_a,clm3%g%l%c%cwf%pwf_a,clm3%g%l%c%cws%pws_a)
               
  call init_misc()
  
  ! check size
  if (i3 >= size(R3,dim=3) .OR. i4 >= size(R4,dim=4)    &
 &  .OR. j3 >= size(L3,dim=3) .OR. j4 >= size(L4,dim=4)) then
  
  call endrun('allocated sizes for either R3, R4, L3 or L4 not sufficient')
  endif
  end subroutine linkCLMtype
  ! ##########################################################################################
    SUBROUTINE init_misc()
    ! initialize misc. items that I don't want to spend that much memory on
    type(gridcell_type),pointer::gg
    type(pft_type),pointer::pp
    type(column_type),pointer::cc
    
    gg              => clm3%g
    cc              => gg%l%c
    pp              => cc%p
    
    pp%pebal%errsoi => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4+1
    pp%pebal%errsol => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4+1
    pp%pebal%errseb => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4+1
    pp%pebal%errlon => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4+1
    
    cc%cebal%errsoi => oneD_ptr(mn,R3(:,:,i3));  i3 = i3 + 1
    cc%cwbal%endwb  => oneD_ptr(mn,R3(:,:,i3));  i3 = i3 + 1
    cc%cwbal%begwb  => oneD_ptr(mn,R3(:,:,i3));  i3 = i3 + 1
    cc%cwbal%errh2o => oneD_ptr(mn,R3(:,:,i3));  i3 = i3 + 1
    cc%cebal%errsoi => oneD_ptr(mn,R3(:,:,i3));  i3 = i3 + 1
    cc%ccbal%begcb  => oneD_ptr(mn,R3(:,:,i3));  i3 = i3 + 1
    cc%ccbal%endcb  => oneD_ptr(mn,R3(:,:,i3));  i3 = i3 + 1
    cc%ccbal%errcb  => oneD_ptr(mn,R3(:,:,i3));  i3 = i3 + 1
    
    cc%cnbal%begnb  => oneD_ptr(mn,R3(:,:,i3));  i3 = i3 + 1
    cc%cnbal%endnb  => oneD_ptr(mn,R3(:,:,i3));  i3 = i3 + 1
    cc%cnbal%errnb  => oneD_ptr(mn,R3(:,:,i3));  i3 = i3 + 1
    
    
    END SUBROUTINE
  ! ##########################################################################################
    SUBROUTINE init_col_avg(begp, endp, pcf_a,pnf_a,pcs_a,pns_a,pwf_a,pws_a)
    ! there is only a handful of data pieces used and I really don't want to spend that much memory storage on these two
    USE memoryChunk, Only: R3,R4,i3,i4
    implicit none
    integer begp,endp
    type (pft_cflux_type), intent(inout) :: pcf_a
    type (pft_nflux_type), intent(inout) :: pnf_a
    type (pft_cstate_type), intent(inout) :: pcs_a
    type (pft_nstate_type), intent(inout) :: pns_a
    type (pft_wflux_type), intent(inout) :: pwf_a
    type (pft_wstate_type), intent(inout) :: pws_a
    
    pcf_a%hrv_xsmrpool_to_atm  => oneD_ptr(mn,R3(:,:,i3));  i3 = i3 + 1
    pcf_a%ar                   => oneD_ptr(mn,R3(:,:,i3));  i3 = i3 + 1
    pcf_a%gpp                  => oneD_ptr(mn,R3(:,:,i3));  i3 = i3 + 1
    pcf_a%npp                  => oneD_ptr(mn,R3(:,:,i3));  i3 = i3 + 1
    pcf_a%pft_fire_closs       => oneD_ptr(mn,R3(:,:,i3));  i3 = i3 + 1
    pcf_a%litfall              => oneD_ptr(mn,R3(:,:,i3));  i3 = i3 + 1
    pcf_a%rr                   => oneD_ptr(mn,R3(:,:,i3));  i3 = i3 + 1
    pcf_a%vegfire              => oneD_ptr(mn,R3(:,:,i3));  i3 = i3 + 1
    pcf_a%wood_harvestc        => oneD_ptr(mn,R3(:,:,i3));  i3 = i3 + 1
    
    pnf_a%pft_fire_nloss       => oneD_ptr(mn,R3(:,:,i3));  i3 = i3 + 1
    pnf_a%wood_harvestn        => oneD_ptr(mn,R3(:,:,i3));  i3 = i3 + 1
    
    pcs_a%totpftc              => oneD_ptr(mn,R3(:,:,i3));  i3 = i3 + 1
    pcs_a%totvegc              => oneD_ptr(mn,R3(:,:,i3));  i3 = i3 + 1
    
    pns_a%totpftn              => oneD_ptr(mn,R3(:,:,i3));  i3 = i3 + 1
    pns_a%totvegn              => oneD_ptr(mn,R3(:,:,i3));  i3 = i3 + 1
    
    
    pwf_a%qflx_prec_intr           => oneD_ptr(mn,R3(:,:,i3));  i3 = i3 + 1
    pwf_a%qflx_prec_grnd       => oneD_ptr(mn,R3(:,:,i3));  i3 = i3 + 1
    pwf_a%qflx_rain_grnd       => oneD_ptr(mn,R3(:,:,i3));  i3 = i3 + 1
    pwf_a%qflx_snow_grnd       => oneD_ptr(mn,R3(:,:,i3));  i3 = i3 + 1
    pwf_a%qflx_snwcp_liq       => oneD_ptr(mn,R3(:,:,i3));  i3 = i3 + 1
    pwf_a%qflx_snwcp_ice       => oneD_ptr(mn,R3(:,:,i3));  i3 = i3 + 1
    pwf_a%qflx_evap_veg        => oneD_ptr(mn,R3(:,:,i3));  i3 = i3 + 1
    pwf_a%qflx_tran_veg        => oneD_ptr(mn,R3(:,:,i3));  i3 = i3 + 1
    pwf_a%qflx_evap_can        => oneD_ptr(mn,R3(:,:,i3));  i3 = i3 + 1
    pwf_a%qflx_evap_soi        => oneD_ptr(mn,R3(:,:,i3));  i3 = i3 + 1
    pwf_a%qflx_evap_tot        => oneD_ptr(mn,R3(:,:,i3));  i3 = i3 + 1
    pwf_a%qflx_evap_grnd       => oneD_ptr(mn,R3(:,:,i3));  i3 = i3 + 1
    pwf_a%qflx_dew_grnd        => oneD_ptr(mn,R3(:,:,i3));  i3 = i3 + 1
    pwf_a%qflx_sub_snow        => oneD_ptr(mn,R3(:,:,i3));  i3 = i3 + 1
    pwf_a%qflx_dew_snow        => oneD_ptr(mn,R3(:,:,i3));  i3 = i3 + 1
    pwf_a%qflx_efpot           => oneD_ptr(mn,R3(:,:,i3));  i3 = i3 + 1
    
    pws_a%h2ocan               => oneD_ptr(mn,R3(:,:,i3));  i3 = i3 + 1
    
    end subroutine
  ! ###################################################################################
  
  
  subroutine init_atm2lnd_type(beg, end, a2l)
  use memoryChunk !CP
  USE clm_atmlnd, Only: atm2lnd_type
!
! !DESCRIPTION:
! Initialize atmospheric variables required by the land
!
! !ARGUMENTS:
  implicit none
  integer, intent(in) :: beg, end
  type (atm2lnd_type), intent(inout):: a2l
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
! Modified by T Craig, 11/01/05 for finemesh project
!
!
! !LOCAL VARIABLES:
!EOP
  real(r8) :: ival   ! initial value
!------------------------------------------------------------------------

  ! ival = nan      ! causes core dump in map_maparray, tcx fix
  ival = 0.0_r8

  a2l%forc_t => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1
  a2l%forc_u => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1
  a2l%forc_v => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1
  a2l%forc_wind => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1
  a2l%forc_q => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1
  a2l%forc_rh => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1
  a2l%forc_hgt => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1
  a2l%forc_hgt_u => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1
  a2l%forc_hgt_t => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1
  a2l%forc_hgt_q => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1
  a2l%forc_pbot => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1
  a2l%forc_th => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1
  a2l%forc_vp => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1
  a2l%forc_rho => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1
  a2l%forc_psrf => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1
  a2l%forc_pco2 => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1
  a2l%forc_lwrad => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1
  a2l%forc_solad => twoD_ptr(mn,numrad,R3(:,:,i3)); i3=i3+numrad
  a2l%forc_solai => twoD_ptr(mn,numrad,R3(:,:,i3)); i3=i3+numrad
  a2l%forc_solar => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1
  a2l%forc_rain => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1
  a2l%forc_snow => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1
  a2l%forc_ndep => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1
  a2l%rainf => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1
#ifdef C13
  a2l%forc_pc13o2 => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1
#endif
  a2l%forc_po2 => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1
  a2l%forc_aer => twoD_ptr(mn,14,R3(:,:,i3)); i3=i3+14

end subroutine init_atm2lnd_type
  ! ###################################################################################
  
  ! END CP added    
  
  subroutine initClmtype()
!
! !DESCRIPTION:
! Initialize clmtype components to signaling nan
! The following clmtype components should NOT be initialized here
! since they are set in routine clm_map which is called before this
! routine is invoked
!    *%area, *%wt, *%wtlnd, *%wtxy, *%ixy, *%jxy, *%mxy, %snindex
!    *%ifspecial, *%ityplun, *%itype
!    *%pfti, *%pftf, *%pftn
!    *%coli, *%colf, *%coln
!    *%luni, *%lunf, *%lunn
!
! !USES:
    !use decompMod , only : get_proc_bounds, get_proc_global
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! LOCAL VARAIBLES:
    integer :: begp, endp   ! per-proc beginning and ending pft indices
    integer :: begc, endc   ! per-proc beginning and ending column indices
    integer :: begl, endl   ! per-proc beginning and ending landunit indices
    integer :: begg, endg   ! per-proc gridcell ending gridcell indices
    integer :: numg         ! total number of gridcells across all processors
    integer :: numl         ! total number of landunits across all processors
    integer :: numc         ! total number of columns across all processors
    integer :: nump         ! total number of pfts across all processors
!------------------------------------------------------------------------

    ! Determine necessary indices

    !call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp) !CP
   ! call get_proc_global(numg, numl, numc, nump) !CP

    call init_pft_type     (begp, endp, clm3%g%l%c%p)
    call init_column_type  (begc, endc, clm3%g%l%c)
    call init_landunit_type(begl, endl, clm3%g%l)
    call init_gridcell_type(begg, endg, clm3%g)

    ! pft ecophysiological constants

    call init_pft_ecophys_constants()

#if (defined CNDV)
    ! pft DGVM-specific ecophysiological constants

    call init_pft_DGVMecophys_constants()
#endif

    ! pft physical state variables at pft level and averaged to the column

    call init_pft_pstate_type(begp, endp, clm3%g%l%c%p%pps)
    call init_pft_pstate_type(begc, endc, clm3%g%l%c%cps%pps_a)

    ! pft ecophysiological variables (only at the pft level for now)
    call init_pft_epv_type(begp, endp, clm3%g%l%c%p%pepv)

#if (defined CNDV)
    ! pft DGVM state variables at pft level and averaged to column

    call init_pft_pdgvstate_type(begp, endp, clm3%g%l%c%p%pdgvs)
#endif
#if (defined CNDV)
    call init_pft_pdgvstate_type(begc, endc, clm3%g%l%c%cdgvs%pdgvs_a)
#endif
    call init_pft_vstate_type(begp, endp, clm3%g%l%c%p%pvs)

    ! pft energy state variables at the pft level and averaged to the column

    call init_pft_estate_type(begp, endp, clm3%g%l%c%p%pes)
    call init_pft_estate_type(begc, endc, clm3%g%l%c%ces%pes_a)

    ! pft water state variables at the pft level and averaged to the column

    call init_pft_wstate_type(begp, endp, clm3%g%l%c%p%pws)
    call init_pft_wstate_type(begc, endc, clm3%g%l%c%cws%pws_a)

    ! pft carbon state variables at the pft level and averaged to the column

    call init_pft_cstate_type(begp, endp, clm3%g%l%c%p%pcs)
    call init_pft_cstate_type(begc, endc, clm3%g%l%c%ccs%pcs_a)
#if (defined C13)
    ! 4/14/05: PET
    ! Adding isotope code
    call init_pft_cstate_type(begp, endp, clm3%g%l%c%p%pc13s)
    call init_pft_cstate_type(begc, endc, clm3%g%l%c%cc13s%pcs_a)
#endif

    ! pft nitrogen state variables at the pft level and averaged to the column

    call init_pft_nstate_type(begp, endp, clm3%g%l%c%p%pns)
    call init_pft_nstate_type(begc, endc, clm3%g%l%c%cns%pns_a)

    ! pft energy flux variables at pft level and averaged to column

    call init_pft_eflux_type(begp, endp, clm3%g%l%c%p%pef)
    call init_pft_eflux_type(begc, endc, clm3%g%l%c%cef%pef_a)

    ! pft momentum flux variables at pft level and averaged to the column

    call init_pft_mflux_type(begp, endp, clm3%g%l%c%p%pmf)
    call init_pft_mflux_type(begc, endc, clm3%g%l%c%cmf%pmf_a)

    ! pft water flux variables

    call init_pft_wflux_type(begp, endp, clm3%g%l%c%p%pwf)
    call init_pft_wflux_type(begc, endc, clm3%g%l%c%cwf%pwf_a)

    ! pft carbon flux variables at pft level and averaged to column

    call init_pft_cflux_type(begp, endp, clm3%g%l%c%p%pcf)
    call init_pft_cflux_type(begc, endc, clm3%g%l%c%ccf%pcf_a)
#if (defined C13)
    ! 4/14/05: PET
    ! Adding isotope code
    call init_pft_cflux_type(begp, endp, clm3%g%l%c%p%pc13f)
    call init_pft_cflux_type(begc, endc, clm3%g%l%c%cc13f%pcf_a)
#endif

    ! pft nitrogen flux variables at pft level and averaged to column

    call init_pft_nflux_type(begp, endp, clm3%g%l%c%p%pnf)
    call init_pft_nflux_type(begc, endc, clm3%g%l%c%cnf%pnf_a)

    ! pft VOC flux variables at pft level and averaged to column

    call init_pft_vflux_type(begp, endp, clm3%g%l%c%p%pvf)
    call init_pft_vflux_type(begc, endc, clm3%g%l%c%cvf%pvf_a)

    ! gridcell VOC emission factors (heald, 05/06)

    call init_gridcell_efstate_type(begg, endg, clm3%g%gve)

    ! pft dust flux variables at pft level and averaged to column

    call init_pft_dflux_type(begp, endp, clm3%g%l%c%p%pdf)
    call init_pft_dflux_type(begc, endc, clm3%g%l%c%cdf%pdf_a)

    ! pft dry dep velocity variables at pft level and averaged to column

    call init_pft_depvd_type(begp, endp, clm3%g%l%c%p%pdd)

    ! column physical state variables at column level and averaged to
    ! the landunit and gridcell and model

    call init_column_pstate_type(begc, endc, clm3%g%l%c%cps)
    call init_column_pstate_type(begl, endl, clm3%g%l%lps%cps_a)
    call init_column_pstate_type(begg, endg, clm3%g%gps%cps_a)
    call init_column_pstate_type(1,       1, clm3%mps%cps_a)

    ! column energy state variables at column level and averaged to
    ! the landunit and gridcell and model

    call init_column_estate_type(begc, endc, clm3%g%l%c%ces)
    call init_column_estate_type(begl, endl, clm3%g%l%les%ces_a)
    call init_column_estate_type(begg, endg, clm3%g%ges%ces_a)
    call init_column_estate_type(1,       1, clm3%mes%ces_a)

    ! column water state variables at column level and averaged to
    ! the landunit and gridcell and model

    call init_column_wstate_type(begc, endc, clm3%g%l%c%cws)
    call init_column_wstate_type(begl, endl, clm3%g%l%lws%cws_a)
    call init_column_wstate_type(begg, endg, clm3%g%gws%cws_a)
    call init_column_wstate_type(1,       1, clm3%mws%cws_a)

    ! column carbon state variables at column level and averaged to
    ! the landunit and gridcell and model

    call init_column_cstate_type(begc, endc, clm3%g%l%c%ccs)
    call init_column_cstate_type(begl, endl, clm3%g%l%lcs%ccs_a)
    call init_column_cstate_type(begg, endg, clm3%g%gcs%ccs_a)
    call init_column_cstate_type(1,       1, clm3%mcs%ccs_a)
#if (defined C13)
    ! 4/14/05: PET
    ! Adding isotope code
    call init_column_cstate_type(begc, endc, clm3%g%l%c%cc13s)
#endif

    ! column nitrogen state variables at column level and averaged to
    ! the landunit and gridcell and model

    call init_column_nstate_type(begc, endc, clm3%g%l%c%cns)
    call init_column_nstate_type(begl, endl, clm3%g%l%lns%cns_a)
    call init_column_nstate_type(begg, endg, clm3%g%gns%cns_a)
    call init_column_nstate_type(1,       1, clm3%mns%cns_a)

    ! column energy flux variables at column level and averaged to
    ! the landunit and gridcell and model

    call init_column_eflux_type(begc, endc, clm3%g%l%c%cef)
    call init_column_eflux_type(begl, endl, clm3%g%l%lef%cef_a)
    call init_column_eflux_type(begg, endg, clm3%g%gef%cef_a)
    call init_column_eflux_type(1,       1, clm3%mef%cef_a)

    ! column water flux variables at column level and averaged to
    ! landunit, gridcell and model level

    call init_column_wflux_type(begc, endc, clm3%g%l%c%cwf)
    call init_column_wflux_type(begl, endl, clm3%g%l%lwf%cwf_a)
    call init_column_wflux_type(begg, endg, clm3%g%gwf%cwf_a)
    call init_column_wflux_type(1,       1, clm3%mwf%cwf_a)

    ! column carbon flux variables at column level

    call init_column_cflux_type(begc, endc, clm3%g%l%c%ccf)
#if (defined C13)
    ! 4/14/05: PET
    ! Adding isotope code
    call init_column_cflux_type(begc, endc, clm3%g%l%c%cc13f)
#endif

    ! column nitrogen flux variables at column level

    call init_column_nflux_type(begc, endc, clm3%g%l%c%cnf)

    ! land unit physical state variables

    call init_landunit_pstate_type(begl, endl, clm3%g%l%lps)

    ! land unit energy flux variables 

    call init_landunit_eflux_type(begl, endl, clm3%g%l%lef)

#if (defined CNDV)
    ! gridcell DGVM variables

    call init_gridcell_dgvstate_type(begg, endg, clm3%g%gdgvs)
#endif

    ! gridcell physical state variables

    call init_gridcell_pstate_type(begg, endg, clm3%g%gps)

    ! gridcell: water flux variables

    call init_gridcell_wflux_type(begg, endg, clm3%g%gwf)

    ! gridcell: energy flux variables

    call init_gridcell_eflux_type(begg, endg, clm3%g%gef)

    ! gridcell: water state variables

    call init_gridcell_wstate_type(begg, endg, clm3%g%gws)

    ! gridcell: energy state variables

    call init_gridcell_estate_type(begg, endg, clm3%g%ges)

  end subroutine initClmtype

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_type
!
! !INTERFACE:
  subroutine init_pft_type (beg, end, p)
!
! !DESCRIPTION:
! Initialize components of pft_type structure
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type(pft_type), intent(inout):: p
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------
    p%gridcell => oneD_ptrI(nn,L4(:,:,:,j4));j4=j4+1
    p%landunit => oneD_ptrI(nn,L4(:,:,:,j4));j4=j4+1
    p%column => oneD_ptrI(nn,L4(:,:,:,j4));j4=j4+1
    p%itype => oneD_ptrI(nn,L4(:,:,:,j4));j4=j4+1 
    p%mxy => oneD_ptrI(nn,L4(:,:,:,j4)); j4 = j4+1
    
    !p%wtgcell => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4+1  !LINKED
    !p%wtlunit => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4+1  !LINKED
    !p%wtcol => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4+1  !LINKED

  end subroutine init_pft_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_column_type
!
! !INTERFACE:
  subroutine init_column_type (beg, end, c)
!
! !DESCRIPTION:
! Initialize components of column_type structure
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type(column_type), intent(inout):: c
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------
   c%gridcell => oneD_ptrI(mn,L3(:,:,j3));j3=j3+1
   c%landunit => oneD_ptrI(mn,L3(:,:,j3));j3=j3+1
   c%pfti => oneD_ptrI(mn,L3(:,:,j3));j3=j3+1
   c%pftf => oneD_ptrI(mn,L3(:,:,j3));j3=j3+1
   c%npfts => oneD_ptrI(mn,L3(:,:,j3));j3=j3+1
   c%itype => oneD_ptrI(mn,L3(:,:,j3));j3=j3+1
   
   c%wtgcell => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
   c%wtlunit => oneD_ptr(mn,R3(:,:,i3));i3=i3+1

  end subroutine init_column_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_landunit_type
!
! !INTERFACE:
  subroutine init_landunit_type (beg, end,l)
!
! !DESCRIPTION:
! Initialize components of landunit_type structure
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type(landunit_type), intent(inout):: l
    integer,pointer:: ii3
    ii3 =>i3
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

   l%gridcell => oneD_ptrI(mn,L3(:,:,j3));j3=j3+1
   l%coli => oneD_ptrI(mn,L3(:,:,j3));j3=j3+1
   l%colf => oneD_ptrI(mn,L3(:,:,j3));j3=j3+1
   l%ncolumns => oneD_ptrI(mn,L3(:,:,j3));j3=j3+1
   l%pfti => oneD_ptrI(mn,L3(:,:,j3));j3=j3+1
   l%pftf => oneD_ptrI(mn,L3(:,:,j3));j3=j3+1
   l%npfts    => oneD_ptrI(mn,L3(:,:,j3));j3=j3+1

   l%itype => oneD_ptrI(mn,L3(:,:,j3));j3=j3+1
   !l%ifspecial => oneD_ptrI(mn,L3(:,:,j3));j3=j3+1
   !l%lakpoi => oneD_ptrI(mn,L3(:,:,j3));j3=j3+1
   !l%urbpoi => oneD_ptrI(mn,L3(:,:,j3));j3=j3+1
   !l%glcmecpoi => oneD_ptrI(mn,L3(:,:,j3));j3=j3+1
   ! current don't do logical, think about it later!
   ! TEMPORARY!!!!
   allocate(l%ifspecial(beg:end))
   allocate(l%lakpoi(beg:end))
   allocate(l%urbpoi(beg:end))
   allocate(l%glcmecpoi(beg:end))


   ! MV - these should be moved to landunit physical state -MV
   l%wtgcell => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
   l%wtgcell = 1D0
   l%canyon_hwr => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
   l%wtroad_perv => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
   l%ht_roof => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
   l%wtlunit_roof => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
   l%wind_hgt_canyon => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
   l%z_0_town => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
   l%z_d_town => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
   
   !l%glcmecpoi(beg:end) = .false.

  end subroutine init_landunit_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_gridcell_type
!
! !INTERFACE:
  subroutine init_gridcell_type (beg, end,g)
!
! !DESCRIPTION:
! Initialize components of gridcell_type structure
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type(gridcell_type), intent(inout):: g
    integer,pointer:: ii3
    ii3 =>i3
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

   g%luni => oneD_ptrI(mn,L3(:,:,j3));j3=j3+1
   g%lunf => oneD_ptrI(mn,L3(:,:,j3));j3=j3+1
   g%nlandunits => oneD_ptrI(mn,L3(:,:,j3));j3=j3+1
   g%coli => oneD_ptrI(mn,L3(:,:,j3));j3=j3+1
   g%colf => oneD_ptrI(mn,L3(:,:,j3));j3=j3+1
   g%ncolumns   => oneD_ptrI(mn,L3(:,:,j3));j3=j3+1
   g%pfti => oneD_ptrI(mn,L3(:,:,j3));j3=j3+1
   g%pftf => oneD_ptrI(mn,L3(:,:,j3));j3=j3+1
   g%npfts => oneD_ptrI(mn,L3(:,:,j3));j3=j3+1
   g%gindex_a => oneD_ptrI(mn,L3(:,:,j3));j3=j3+1

   g%gindex => oneD_ptrI(mn,L3(:,:,j3));j3=j3+1
   g%area => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
   g%lat => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
   g%lon => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
   g%latdeg => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
   g%londeg => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
   g%lat_a => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
   g%lon_a => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
   g%latdeg_a => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
   g%londeg_a => oneD_ptr(mn,R3(:,:,i3));i3=i3+1

   g%gris_mask => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
   g%gris_area => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
   g%aais_mask => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
   g%aais_area => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
  end subroutine init_gridcell_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_energy_balance_type
!
! !INTERFACE:
  subroutine init_energy_balance_type(beg, end, ebal)
!
! !DESCRIPTION:
! Initialize energy balance variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type(energy_balance_type), intent(inout):: ebal
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(ebal%errsoi(beg:end))
    allocate(ebal%errseb(beg:end))
    allocate(ebal%errsol(beg:end))
    allocate(ebal%errlon(beg:end))

    ebal%errsoi(beg:end) = nan
    ebal%errseb(beg:end) = nan
    ebal%errsol(beg:end) = nan
    ebal%errlon(beg:end) = nan

  end subroutine init_energy_balance_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_water_balance_type
!
! !INTERFACE:
  subroutine init_water_balance_type(beg, end, wbal)
!
! !DESCRIPTION:
! Initialize water balance variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type(water_balance_type), intent(inout):: wbal
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(wbal%begwb(beg:end))
    allocate(wbal%endwb(beg:end))
    allocate(wbal%errh2o(beg:end))

    wbal%begwb(beg:end) = nan
    wbal%endwb(beg:end) = nan
    wbal%errh2o(beg:end) = nan

  end subroutine init_water_balance_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_carbon_balance_type
!
! !INTERFACE:
  subroutine init_carbon_balance_type(beg, end, cbal)
!
! !DESCRIPTION:
! Initialize carbon balance variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type(carbon_balance_type), intent(inout):: cbal
!
! !REVISION HISTORY:
! Created by Peter Thornton, 12/11/2003
!
!EOP
!------------------------------------------------------------------------

    allocate(cbal%begcb(beg:end))
    allocate(cbal%endcb(beg:end))
    allocate(cbal%errcb(beg:end))

    cbal%begcb(beg:end) = nan
    cbal%endcb(beg:end) = nan
    cbal%errcb(beg:end) = nan

  end subroutine init_carbon_balance_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_nitrogen_balance_type
!
! !INTERFACE:
  subroutine init_nitrogen_balance_type(beg, end, nbal)
!
! !DESCRIPTION:
! Initialize nitrogen balance variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type(nitrogen_balance_type), intent(inout):: nbal
!
! !REVISION HISTORY:
! Created by Peter Thornton, 12/11/2003
!
!EOP
!------------------------------------------------------------------------

    allocate(nbal%begnb(beg:end))
    allocate(nbal%endnb(beg:end))
    allocate(nbal%errnb(beg:end))

    nbal%begnb(beg:end) = nan
    nbal%endnb(beg:end) = nan
    nbal%errnb(beg:end) = nan

  end subroutine init_nitrogen_balance_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_ecophys_constants
!
! !INTERFACE:
  subroutine init_pft_ecophys_constants()
!
! !DESCRIPTION:
! Initialize pft physical state
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(pftcon%noveg(0:numpft))
    allocate(pftcon%tree(0:numpft))
    allocate(pftcon%smpso(0:numpft)) 
    allocate(pftcon%smpsc(0:numpft)) 
    allocate(pftcon%fnitr(0:numpft))
    allocate(pftcon%foln(0:numpft))
    allocate(pftcon%dleaf(0:numpft))
    allocate(pftcon%c3psn(0:numpft))
    allocate(pftcon%mp(0:numpft))
    allocate(pftcon%qe25(0:numpft))
    allocate(pftcon%xl(0:numpft))
    allocate(pftcon%rhol(0:numpft,numrad))
    allocate(pftcon%rhos(0:numpft,numrad))
    allocate(pftcon%taul(0:numpft,numrad))
    allocate(pftcon%taus(0:numpft,numrad))
    allocate(pftcon%z0mr(0:numpft))
    allocate(pftcon%displar(0:numpft))
    allocate(pftcon%roota_par(0:numpft))
    allocate(pftcon%rootb_par(0:numpft))
    allocate(pftcon%sla(0:numpft))
    allocate(pftcon%slatop(0:numpft))
    allocate(pftcon%dsladlai(0:numpft))
    allocate(pftcon%leafcn(0:numpft))
    allocate(pftcon%flnr(0:numpft))
    allocate(pftcon%woody(0:numpft))
    allocate(pftcon%lflitcn(0:numpft))
    allocate(pftcon%frootcn(0:numpft))
    allocate(pftcon%livewdcn(0:numpft))
    allocate(pftcon%deadwdcn(0:numpft))
    allocate(pftcon%graincn(0:numpft))
    allocate(pftcon%froot_leaf(0:numpft))
    allocate(pftcon%stem_leaf(0:numpft))
    allocate(pftcon%croot_stem(0:numpft))
    allocate(pftcon%flivewd(0:numpft))
    allocate(pftcon%fcur(0:numpft))
    allocate(pftcon%lf_flab(0:numpft))
    allocate(pftcon%lf_fcel(0:numpft))
    allocate(pftcon%lf_flig(0:numpft))
    allocate(pftcon%fr_flab(0:numpft))
    allocate(pftcon%fr_fcel(0:numpft))
    allocate(pftcon%fr_flig(0:numpft))
    allocate(pftcon%leaf_long(0:numpft))
    allocate(pftcon%evergreen(0:numpft))
    allocate(pftcon%stress_decid(0:numpft))
    allocate(pftcon%season_decid(0:numpft))
    allocate(pftcon%resist(0:numpft))
    allocate(pftcon%dwood(0:numpft))

    pftcon%noveg(:) = bigint
    pftcon%tree(:) = bigint
    pftcon%smpso(:) = nan
    pftcon%smpsc(:) = nan
    pftcon%fnitr(:) = nan
    pftcon%foln(:) = nan
    pftcon%dleaf(:) = nan
    pftcon%c3psn(:) = nan
    pftcon%mp(:) = nan
    pftcon%qe25(:) = nan
    pftcon%xl(:) = nan
    pftcon%rhol(:,:numrad) = nan
    pftcon%rhos(:,:numrad) = nan
    pftcon%taul(:,:numrad) = nan
    pftcon%taus(:,:numrad) = nan
    pftcon%z0mr(:) = nan
    pftcon%displar(:) = nan
    pftcon%roota_par(:) = nan
    pftcon%rootb_par(:) = nan
    pftcon%sla(:) = nan
    pftcon%slatop(:) = nan
    pftcon%dsladlai(:) = nan
    pftcon%leafcn(:) = nan
    pftcon%flnr(:) = nan
    pftcon%woody(:) = nan
    pftcon%lflitcn(:) = nan
    pftcon%frootcn(:) = nan
    pftcon%livewdcn(:) = nan
    pftcon%deadwdcn(:) = nan
    pftcon%graincn(:) = nan
    pftcon%froot_leaf(:) = nan
    pftcon%stem_leaf(:) = nan
    pftcon%croot_stem(:) = nan
    pftcon%flivewd(:) = nan
    pftcon%fcur(:) = nan
    pftcon%lf_flab(:) = nan
    pftcon%lf_fcel(:) = nan
    pftcon%lf_flig(:) = nan
    pftcon%fr_flab(:) = nan
    pftcon%fr_fcel(:) = nan
    pftcon%fr_flig(:) = nan
    pftcon%leaf_long(:) = nan
    pftcon%evergreen(:) = nan
    pftcon%stress_decid(:) = nan
    pftcon%season_decid(:) = nan
    pftcon%resist(:) = nan
    pftcon%dwood(:) = nan

  end subroutine init_pft_ecophys_constants

#if (defined CNDV)
!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_DGVMecophys_constants
!
! !INTERFACE:
  subroutine init_pft_DGVMecophys_constants()
!
! !DESCRIPTION:
! Initialize pft physical state
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(dgv_pftcon%crownarea_max(0:numpft))
    allocate(dgv_pftcon%tcmin(0:numpft))
    allocate(dgv_pftcon%tcmax(0:numpft))
    allocate(dgv_pftcon%gddmin(0:numpft))
    allocate(dgv_pftcon%twmax(0:numpft))
    allocate(dgv_pftcon%reinickerp(0:numpft))
    allocate(dgv_pftcon%allom1(0:numpft))
    allocate(dgv_pftcon%allom2(0:numpft))
    allocate(dgv_pftcon%allom3(0:numpft))

    dgv_pftcon%crownarea_max(:) = nan
    dgv_pftcon%tcmin(:) = nan
    dgv_pftcon%tcmax(:) = nan
    dgv_pftcon%gddmin(:) = nan
    dgv_pftcon%twmax(:) = nan
    dgv_pftcon%reinickerp(:) = nan
    dgv_pftcon%allom1(:) = nan
    dgv_pftcon%allom2(:) = nan
    dgv_pftcon%allom3(:) = nan

  end subroutine init_pft_DGVMecophys_constants
#endif

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_pstate_type
!
! !INTERFACE:
  subroutine init_pft_pstate_type(beg, end, pps)
!
! !DESCRIPTION:
! Initialize pft physical state
!
! !USES:
    use clm_varcon, only : spval
    use surfrdMod , only : crop_prog
#if (defined CASA)
    use CASAMod   , only : npools, nresp_pools, nlive, npool_types
#endif
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    integer,pointer::ii4
    type (pft_pstate_type), intent(inout):: pps
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------
    ! 4/14/05: PET
    ! Adding isotope code
    ! heald: added from CASA definition
    !allocate(pps%sandfrac(beg:end))
    !allocate(pps%clayfrac(beg:end))
    !pps%sandfrac(beg:end) = nan
    !pps%clayfrac(beg:end) = nan
    !allocate(pps%mlaidiff(beg:end))
    !allocate(pps%rb1(beg:end))
    !allocate(pps%annlai(12,beg:end))
    !pps%mlaidiff(beg:end) = nan
    !pps%rb1(beg:end) = nan
    !pps%annlai(:,:) = nan
    
    ii4 => i4
    R4C%pps(1) = i4 ! counter
    !pps%frac_veg_nosno(beg:end) = bigint
    !pps%frac_veg_nosno_alb(beg:end) = 0
    
    pps%frac_veg_nosno => oneD_ptrI(nn,L4(:,:,:,j4));j4=j4+1
    pps%frac_veg_nosno_alb => oneD_ptrI(nn,L4(:,:,:,j4));j4=j4+1
    pps%emv => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pps%z0mv => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pps%z0hv => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pps%z0qv => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    
    !pps%rootfr => twoD_ptr(nn,nlevgrnd,R4(:,:,:,i4:i4+nlevgrnd-1)); i4 = i4+nlevgrnd
    !pps%rootr => twoD_ptr(nn,nlevgrnd,R4(:,:,:,i4:i4+nlevgrnd-1)); i4 = i4+nlevgrnd
    !pps%rresis => twoD_ptr(nn,nlevgrnd,R4(:,:,:,i4:i4+nlevgrnd-1)); i4 = i4+nlevgrnd
    pps%rootfr => twoD_ptr(nn,nlevgrnd,R4(:,:,:,i4)); i4 = i4+nlevgrnd
    pps%rootr => twoD_ptr(nn,nlevgrnd,R4(:,:,:,i4)); i4 = i4+nlevgrnd
    pps%rresis => twoD_ptr(nn,nlevgrnd,R4(:,:,:,i4)); i4 = i4+nlevgrnd
    !pps%rootfr(beg:end,:nlevgrnd) = spval !CP except for CNMResp, nothing really use this
    !pps%rootr (beg:end,:nlevgrnd) = spval !CP geophysic uses this. Not useful for geochem
    !pps%rresis(beg:end,:nlevgrnd) = spval !CP geophysic uses this. Not useful for geochem
    pps%dewmx => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pps%rssun => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pps%rssha => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pps%laisun => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pps%laisha => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pps%btran => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    !pps%fsun(beg:end) = spval
    !pps%tlai(beg:end) = 0._r8
    !pps%tsai(beg:end) = 0._r8
    !pps%elai(beg:end) = 0._r8
    !pps%esai(beg:end) = 0._r8
    pps%fsun => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    !pps%tlai => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1  !LINKED
    pps%tsai => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pps%elai => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pps%esai => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1 ! state variables, initialize
    pps%fwet => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pps%fdry => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pps%dt_veg => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    !pps%htop(beg:end) = 0._r8
    !pps%hbot(beg:end) = 0._r8
    !pps%htop => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1  !LINKED
    pps%hbot => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pps%z0m => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pps%displa => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pps%albd => twoD_ptr(nn,numrad,R4(:,:,:,i4)); i4 = i4 + numrad
    pps%albi => twoD_ptr(nn,numrad,R4(:,:,:,i4)); i4 = i4 + numrad
    pps%fabd => twoD_ptr(nn,numrad,R4(:,:,:,i4)); i4 = i4 + numrad
    pps%fabi => twoD_ptr(nn,numrad,R4(:,:,:,i4)); i4 = i4 + numrad
    pps%ftdd => twoD_ptr(nn,numrad,R4(:,:,:,i4)); i4 = i4 + numrad
    pps%ftid => twoD_ptr(nn,numrad,R4(:,:,:,i4)); i4 = i4 + numrad
    pps%ftii => twoD_ptr(nn,numrad,R4(:,:,:,i4)); i4 = i4 + numrad
    pps%u10 => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pps%u10_clm => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pps%va => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    if ( crop_prog )then
        pps%hdidx => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
        pps%cumvd => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
        pps%htmx => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
        pps%vf => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
        pps%gddmaturity => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
        pps%gdd0 => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
        pps%gdd8 => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
        pps%gdd10 => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
        pps%gdd020 => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
        pps%gdd820 => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
        pps%gdd1020 => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
        pps%gddplant => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
        pps%gddtsoi => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
        pps%huileaf => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
        pps%huigrain => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
        pps%aleafi => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
        pps%astemi => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
        pps%aleaf => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
        pps%astem => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
        
       allocate(pps%croplive(beg:end))
       allocate(pps%cropplant(beg:end)) !,numpft)) ! make 2-D if using
       allocate(pps%harvdate(beg:end))  !,numpft)) ! crop rotation
        
       pps%cropplant = .false.
       pps%croplive  = .false.
       pps%harvdate  = 999
       
        pps%idop => oneD_ptrI(nn,L4(:,:,:,j4));j4=j4+1
        pps%peaklai => oneD_ptrI(nn,L4(:,:,:,j4));j4=j4+1
    end if
    pps%fv => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pps%ram1 => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pps%vds => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pps%slasun => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pps%slasha => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pps%lncsun => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pps%lncsha => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pps%vcmxsun => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pps%vcmxsha => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pps%gdir => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pps%omega => twoD_ptr(nn,numrad,R4(:,:,:,i4)); i4 = i4 + numrad
    pps%eff_kid => twoD_ptr(nn,numrad,R4(:,:,:,i4)); i4 = i4 + numrad
    pps%eff_kii => twoD_ptr(nn,numrad,R4(:,:,:,i4)); i4 = i4 + numrad
    pps%sun_faid => twoD_ptr(nn,numrad,R4(:,:,:,i4)); i4 = i4 + numrad
    pps%sun_faii => twoD_ptr(nn,numrad,R4(:,:,:,i4)); i4 = i4 + numrad
    pps%sha_faid => twoD_ptr(nn,numrad,R4(:,:,:,i4)); i4 = i4 + numrad
    pps%sha_faii => twoD_ptr(nn,numrad,R4(:,:,:,i4)); i4 = i4 + numrad
    pps%forc_hgt_u_pft => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pps%forc_hgt_t_pft => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pps%forc_hgt_q_pft => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    ! 4/14/05: PET
    ! Adding isotope code
    pps%cisun => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pps%cisha => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    
    pps%sandfrac => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pps%clayfrac => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pps%mlaidiff => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pps%rb1 => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    !allocate(pps%annlai(12,beg:end))
    pps%annlai => twoD_ptr(12,nn,R4(:,:,:,i4)); i4 = i4 + 12
    
    R4C%pps(2) = i4-1 ! counter
  end subroutine init_pft_pstate_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_epv_type
!
! !INTERFACE:
  subroutine init_pft_epv_type(beg, end, pepv)
!
! !DESCRIPTION:
! Initialize pft ecophysiological variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (pft_epv_type), intent(inout):: pepv
!
! !REVISION HISTORY:
! Created by Peter Thornton
!
!EOP
!------------------------------------------------------------------------
    R4C%pepv(1) = i4 ! counter
    pepv%dormant_flag => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pepv%days_active => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pepv%onset_flag => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pepv%onset_counter => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pepv%onset_gddflag => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pepv%onset_fdd => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pepv%onset_gdd => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pepv%onset_swi => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pepv%offset_flag => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    
    pepv%spheno_lyr => oneD_ptrI(nn,L4(:,:,:,j4)); j4=j4+1 !CP added field
    
    pepv%offset_counter => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pepv%offset_fdd => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pepv%offset_swi => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pepv%lgsf => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pepv%bglfr => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pepv%bgtr => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pepv%dayl => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pepv%prev_dayl => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pepv%annavg_t2m => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pepv%tempavg_t2m => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pepv%gpp => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pepv%availc => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pepv%xsmrpool_recover => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
#if (defined C13)
    pepv%xsmrpool_c13ratio => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
#endif
    pepv%alloc_pnow => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pepv%c_allometry => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pepv%n_allometry => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pepv%plant_ndemand => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pepv%tempsum_potential_gpp => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pepv%annsum_potential_gpp => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pepv%tempmax_retransn => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pepv%annmax_retransn => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pepv%avail_retransn => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pepv%plant_nalloc => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pepv%plant_calloc => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pepv%excess_cflux => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pepv%downreg => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pepv%prev_leafc_to_litter => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pepv%prev_frootc_to_litter => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pepv%tempsum_npp => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pepv%annsum_npp => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
#if (defined CNDV)
    pepv%tempsum_litfall => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pepv%annsum_litfall => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
#endif
#if (defined C13)
    ! 4/21/05, PET
    ! Adding isotope code
    pepv%rc13_canair => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pepv%rc13_psnsun => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pepv%rc13_psnsha => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
#endif
    R4C%pepv(2) = i4 - 1
  end subroutine init_pft_epv_type

#if (defined CNDV)
!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_pdgvstate_type
!
! !INTERFACE:
  subroutine init_pft_pdgvstate_type(beg, end, pdgvs)
!
! !DESCRIPTION:
! Initialize pft DGVM state variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (pft_dgvstate_type), intent(inout):: pdgvs
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(pdgvs%agddtw(beg:end))
    allocate(pdgvs%agdd(beg:end))
    allocate(pdgvs%t_mo(beg:end))
    allocate(pdgvs%t_mo_min(beg:end))
    allocate(pdgvs%prec365(beg:end))
    allocate(pdgvs%present(beg:end))
    allocate(pdgvs%pftmayexist(beg:end))
    allocate(pdgvs%nind(beg:end))
    allocate(pdgvs%lm_ind(beg:end))
    allocate(pdgvs%lai_ind(beg:end))
    allocate(pdgvs%fpcinc(beg:end))
    allocate(pdgvs%fpcgrid(beg:end))
    allocate(pdgvs%fpcgridold(beg:end))
    allocate(pdgvs%crownarea(beg:end))
    allocate(pdgvs%greffic(beg:end))
    allocate(pdgvs%heatstress(beg:end))

    pdgvs%agddtw(beg:end)           = nan
    pdgvs%agdd(beg:end)             = nan
    pdgvs%t_mo(beg:end)             = nan
    pdgvs%t_mo_min(beg:end)         = nan
    pdgvs%prec365(beg:end)          = nan
    pdgvs%present(beg:end)          = .false.
    pdgvs%pftmayexist(beg:end)      = .true.
    pdgvs%nind(beg:end)             = nan
    pdgvs%lm_ind(beg:end)           = nan
    pdgvs%lai_ind(beg:end)          = nan
    pdgvs%fpcinc(beg:end)           = nan
    pdgvs%fpcgrid(beg:end)          = nan
    pdgvs%fpcgridold(beg:end)       = nan
    pdgvs%crownarea(beg:end)        = nan
    pdgvs%greffic(beg:end)          = nan
    pdgvs%heatstress(beg:end)       = nan

  end subroutine init_pft_pdgvstate_type
#endif

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_vstate_type
!
! !INTERFACE:
  subroutine init_pft_vstate_type(beg, end, pvs)
!
! !DESCRIPTION:
! Initialize pft VOC variables
!
! !USES:
    use clm_varcon, only : spval
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (pft_vstate_type), intent(inout):: pvs
!
! !REVISION HISTORY:
! Created by Erik Kluzek
!
!EOP
!------------------------------------------------------------------------

    allocate(pvs%t_veg24 (beg:end))
    allocate(pvs%t_veg240(beg:end))
    allocate(pvs%fsd24   (beg:end))
    allocate(pvs%fsd240  (beg:end))
    allocate(pvs%fsi24   (beg:end))
    allocate(pvs%fsi240  (beg:end))
    allocate(pvs%fsun24  (beg:end))
    allocate(pvs%fsun240 (beg:end))
    allocate(pvs%elai_p  (beg:end))

    pvs%t_veg24 (beg:end)   = spval
    pvs%t_veg240(beg:end)   = spval
    pvs%fsd24   (beg:end)   = spval
    pvs%fsd240  (beg:end)   = spval
    pvs%fsi24   (beg:end)   = spval
    pvs%fsi240  (beg:end)   = spval
    pvs%fsun24  (beg:end)   = spval
    pvs%fsun240 (beg:end)   = spval
    pvs%elai_p  (beg:end)   = spval
  end subroutine init_pft_vstate_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_estate_type
!
! !INTERFACE:
  subroutine init_pft_estate_type(beg, end, pes)
!
! !DESCRIPTION:
! Initialize pft energy state
!
! !USES:
    use clm_varcon, only : spval
    use surfrdMod, only : crop_prog
! !AGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (pft_estate_type), intent(inout):: pes
    integer, pointer::ii4
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------
    ii4 => i4
    pes%t_ref2m => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1
    pes%t_ref2m_min => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1
    pes%t_ref2m_max => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1
    pes%t_ref2m_min_inst => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1
    pes%t_ref2m_max_inst => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1
    pes%q_ref2m => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1
    pes%t_ref2m_u => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1
    pes%t_ref2m_r => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1
    pes%t_ref2m_min_u => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1
    pes%t_ref2m_min_r => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1
    pes%t_ref2m_max_u => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1
    pes%t_ref2m_max_r => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1
    pes%t_ref2m_min_inst_u => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1
    pes%t_ref2m_min_inst_r => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1
    pes%t_ref2m_max_inst_u => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1
    pes%t_ref2m_max_inst_r => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1
#if (defined CNDV) || (defined CROP)
    pes%t10=>oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1
#endif
    if ( crop_prog )then
       pes%a10tmin =>oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1
       pes%a5tmin  =>oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1
    end if
    pes%rh_ref2m => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1
    pes%rh_ref2m_u => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1
    pes%rh_ref2m_r => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1
    pes%t_veg => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1
    pes%thm => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1

  end subroutine init_pft_estate_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_wstate_type
!
! !INTERFACE:
  subroutine init_pft_wstate_type(beg, end, pws)
!
! !DESCRIPTION:
! Initialize pft water state
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (pft_wstate_type), intent(inout):: pws !pft water state
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    pws%h2ocan => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1

  end subroutine init_pft_wstate_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_cstate_type
!
! !INTERFACE:
  subroutine init_pft_cstate_type(beg, end, pcs)
!
! !DESCRIPTION:
! Initialize pft carbon state
!
! !USES:
    use surfrdMod, only : crop_prog
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (pft_cstate_type), intent(inout):: pcs !pft carbon state
!
! !REVISION HISTORY:
! Created by Peter Thornton
!
!EOP
!------------------------------------------------------------------------
    pcs%leafc => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcs%leafc_storage => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcs%leafc_xfer => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcs%frootc => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcs%frootc_storage => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcs%frootc_xfer => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcs%livestemc => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcs%livestemc_storage => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcs%livestemc_xfer => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcs%deadstemc => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcs%deadstemc_storage => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcs%deadstemc_xfer => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcs%livecrootc => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcs%livecrootc_storage => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcs%livecrootc_xfer => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcs%deadcrootc => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcs%deadcrootc_storage => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcs%deadcrootc_xfer => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcs%gresp_storage => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcs%gresp_xfer => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcs%cpool => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcs%xsmrpool => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcs%pft_ctrunc => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcs%dispvegc => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcs%storvegc => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcs%totvegc => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcs%totpftc => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcs%leafcmax => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    if ( crop_prog )then
       pcs%grainc => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
       pcs%grainc_storage => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
       pcs%grainc_xfer => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    end if
#ifdef CN
    pcs%woodc => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
#endif

  end subroutine init_pft_cstate_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_nstate_type
!
! !INTERFACE:
  subroutine init_pft_nstate_type(beg, end, pns)
!
! !DESCRIPTION:
! Initialize pft nitrogen state
!
! !USES:
    use surfrdMod, only : crop_prog
    implicit none
    integer, intent(in) :: beg, end
    type (pft_nstate_type), intent(inout):: pns !pft nitrogen state
!
! !REVISION HISTORY:
! Created by Peter Thornton
!
!EOP
!------------------------------------------------------------------------
    if ( crop_prog )then
       pns%grainn => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
       pns%grainn_storage => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
       pns%grainn_xfer => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    end if
    pns%leafn => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pns%leafn_storage => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pns%leafn_xfer => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pns%frootn => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pns%frootn_storage => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pns%frootn_xfer => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pns%livestemn => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pns%livestemn_storage => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pns%livestemn_xfer => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pns%deadstemn => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pns%deadstemn_storage => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pns%deadstemn_xfer => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pns%livecrootn => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pns%livecrootn_storage => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pns%livecrootn_xfer => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pns%deadcrootn => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pns%deadcrootn_storage => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pns%deadcrootn_xfer => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pns%retransn => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pns%npool => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pns%pft_ntrunc => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pns%dispvegn => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pns%storvegn => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pns%totvegn => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pns%totpftn => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1

  end subroutine init_pft_nstate_type
!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_eflux_type
!
! !INTERFACE:
  subroutine init_pft_eflux_type(beg, end, pef)
!
! !DESCRIPTION:
! Initialize pft energy flux variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (pft_eflux_type), intent(inout):: pef
    integer, pointer :: ii4
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------
    ii4 => i4
    pef%sabg => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1
    pef%sabv => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1
    pef%fsa => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1
    pef%fsa_u => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1
    pef%fsa_r => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1
    pef%fsr => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1
    pef%parsun => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1
    pef%parsha => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1
    pef%dlrad => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1
    pef%ulrad => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1
    pef%eflx_lh_tot => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1
    pef%eflx_lh_tot_u => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1
    pef%eflx_lh_tot_r => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1
    pef%eflx_lh_grnd => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1
    pef%eflx_soil_grnd => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1
    pef%eflx_soil_grnd_u => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1
    pef%eflx_soil_grnd_r => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1
    pef%eflx_sh_tot => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1
    pef%eflx_sh_tot_u => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1
    pef%eflx_sh_tot_r => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1
    pef%eflx_sh_grnd => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1
    pef%eflx_sh_veg => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1
    pef%eflx_lh_vege => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1
    pef%eflx_lh_vegt => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1
    pef%eflx_wasteheat_pft => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1
    pef%eflx_heat_from_ac_pft => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1
    pef%eflx_traffic_pft => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1
    pef%eflx_anthro => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1
    pef%cgrnd => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1
    pef%cgrndl => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1
    pef%cgrnds => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1
    pef%eflx_gnet => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1
    pef%dgnetdT => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1
    pef%eflx_lwrad_out => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1
    pef%eflx_lwrad_net => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1
    pef%eflx_lwrad_net_u => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1
    pef%eflx_lwrad_net_r => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1
    pef%netrad => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1
    pef%fsds_vis_d => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1
    pef%fsds_nir_d => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1
    pef%fsds_vis_i => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1
    pef%fsds_nir_i => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1
    pef%fsr_vis_d => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1
    pef%fsr_nir_d => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1
    pef%fsr_vis_i => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1
    pef%fsr_nir_i => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1
    pef%fsds_vis_d_ln => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1
    pef%fsds_nir_d_ln => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1
    pef%fsr_vis_d_ln => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1
    pef%fsr_nir_d_ln => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1
    pef%sun_add => twoD_ptr(nn,numrad,R4(:,:,:,i4)); i4 = i4 + numrad
    pef%tot_aid => twoD_ptr(nn,numrad,R4(:,:,:,i4)); i4 = i4 + numrad
    pef%sun_aid => twoD_ptr(nn,numrad,R4(:,:,:,i4)); i4 = i4 + numrad
    pef%sun_aii => twoD_ptr(nn,numrad,R4(:,:,:,i4)); i4 = i4 + numrad
    pef%sha_aid => twoD_ptr(nn,numrad,R4(:,:,:,i4)); i4 = i4 + numrad
    pef%sha_aii => twoD_ptr(nn,numrad,R4(:,:,:,i4)); i4 = i4 + numrad
    pef%sun_atot => twoD_ptr(nn,numrad,R4(:,:,:,i4)); i4 = i4 + numrad
    pef%sha_atot => twoD_ptr(nn,numrad,R4(:,:,:,i4)); i4 = i4 + numrad
    pef%sun_alf => twoD_ptr(nn,numrad,R4(:,:,:,i4)); i4 = i4 + numrad
    pef%sha_alf => twoD_ptr(nn,numrad,R4(:,:,:,i4)); i4 = i4 + numrad
    pef%sun_aperlai => twoD_ptr(nn,numrad,R4(:,:,:,i4)); i4 = i4 + numrad
    pef%sha_aperlai => twoD_ptr(nn,numrad,R4(:,:,:,i4)); i4 = i4 + numrad
    pef%sabg_lyr => twoD_ptrO(nn,numrad,R4(:,:,:,i4),-nlevsno+1); i4 = i4 + nlevsno+1
    !pef%sabg_lyr(beg:end,-nlevsno+1:1) = nan
    pef%sfc_frc_aer => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1
    pef%sfc_frc_bc => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1
    pef%sfc_frc_oc => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1
    pef%sfc_frc_dst => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1
    pef%sfc_frc_aer_sno => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1
    pef%sfc_frc_bc_sno => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1
    pef%sfc_frc_oc_sno => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1
    pef%sfc_frc_dst_sno => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1
    pef%fsr_sno_vd => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1
    pef%fsr_sno_nd => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1
    pef%fsr_sno_vi => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1
    pef%fsr_sno_ni => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1
    pef%fsds_sno_vd => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1
    pef%fsds_sno_nd => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1
    pef%fsds_sno_vi => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1
    pef%fsds_sno_ni => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1
  end subroutine init_pft_eflux_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_mflux_type
!
! !INTERFACE:
  subroutine init_pft_mflux_type(beg, end, pmf)
!
! !DESCRIPTION:
! Initialize pft momentum flux variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (pft_mflux_type), intent(inout) :: pmf
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------
    pmf%taux => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1
    pmf%tauy => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1

  end subroutine init_pft_mflux_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_wflux_type
!
! !INTERFACE:
  subroutine init_pft_wflux_type(beg, end, pwf)
!
! !DESCRIPTION:
! Initialize pft water flux variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (pft_wflux_type), intent(inout) :: pwf
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    pwf%qflx_prec_intr => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1
    pwf%qflx_prec_grnd => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1
    pwf%qflx_rain_grnd => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1
    pwf%qflx_snow_grnd => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1
    pwf%qflx_snwcp_liq => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1
    pwf%qflx_snwcp_ice => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1
    pwf%qflx_evap_veg => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1
    pwf%qflx_tran_veg => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1
    pwf%qflx_evap_can => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1
    pwf%qflx_evap_soi => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1
    pwf%qflx_evap_tot => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1
    pwf%qflx_evap_grnd => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1
    pwf%qflx_dew_grnd => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1
    pwf%qflx_sub_snow => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1
    pwf%qflx_dew_snow => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1
    pwf%qflx_efpot    => oneD_ptr(nn,R4(:,:,:,i4)); i4 = i4 + 1

  end subroutine init_pft_wflux_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_cflux_type
!
! !INTERFACE:
  subroutine init_pft_cflux_type(beg, end, pcf)
!
! !DESCRIPTION:
! Initialize pft carbon flux variables
!
! !USES:
    use surfrdMod , only : crop_prog
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (pft_cflux_type), intent(inout) :: pcf
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------
    
    pcf%psnsun => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcf%psnsha => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcf%fpsn => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcf%fco2 => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    ! removed all m_*** (mortality) and hrv_*** (harvest) pointers
    ! they now all point to the same patch, which will remain 0. the reduction will be taken in SUBROUTINE ****   
    pcf%m_leafc_to_litter => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1
    pcf%m_frootc_to_litter => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1
    pcf%m_leafc_storage_to_litter => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1
    pcf%m_frootc_storage_to_litter => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1
    pcf%m_livestemc_storage_to_litter => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1
    pcf%m_deadstemc_storage_to_litter => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1
    pcf%m_livecrootc_storage_to_litter => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1
    pcf%m_deadcrootc_storage_to_litter => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1
    pcf%m_leafc_xfer_to_litter => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1
    pcf%m_frootc_xfer_to_litter => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1
    pcf%m_livestemc_xfer_to_litter => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1
    pcf%m_deadstemc_xfer_to_litter => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1
    pcf%m_livecrootc_xfer_to_litter => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1
    pcf%m_deadcrootc_xfer_to_litter => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1
    pcf%m_livestemc_to_litter => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1
    pcf%m_deadstemc_to_litter => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1
    pcf%m_livecrootc_to_litter => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1
    pcf%m_deadcrootc_to_litter => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1
    pcf%m_gresp_storage_to_litter => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1
    pcf%m_gresp_xfer_to_litter => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1
    pcf%hrv_leafc_to_litter => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1             
    pcf%hrv_leafc_storage_to_litter => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1     
    pcf%hrv_leafc_xfer_to_litter => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1        
    pcf%hrv_frootc_to_litter => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1            
    pcf%hrv_frootc_storage_to_litter => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1    
    pcf%hrv_frootc_xfer_to_litter => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1       
    pcf%hrv_livestemc_to_litter => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1         
    pcf%hrv_livestemc_storage_to_litter => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1 
    pcf%hrv_livestemc_xfer_to_litter => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1    
    pcf%hrv_deadstemc_to_prod10c => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1        
    pcf%hrv_deadstemc_to_prod100c => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1       
    pcf%hrv_deadstemc_storage_to_litter => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1 
    pcf%hrv_deadstemc_xfer_to_litter => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1    
    pcf%hrv_livecrootc_to_litter => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1        
    pcf%hrv_livecrootc_storage_to_litter => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1
    pcf%hrv_livecrootc_xfer_to_litter => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1   
    pcf%hrv_deadcrootc_to_litter => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1        
    pcf%hrv_deadcrootc_storage_to_litter => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1
    pcf%hrv_deadcrootc_xfer_to_litter => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1   
    pcf%hrv_gresp_storage_to_litter => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1     
    pcf%hrv_gresp_xfer_to_litter => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1        
    pcf%hrv_xsmrpool_to_atm => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1                 
    pcf%m_leafc_to_fire => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1
    pcf%m_frootc_to_fire => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1
    pcf%m_leafc_storage_to_fire => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1
    pcf%m_frootc_storage_to_fire => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1
    pcf%m_livestemc_storage_to_fire => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1
    pcf%m_deadstemc_storage_to_fire => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1
    pcf%m_livecrootc_storage_to_fire => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1
    pcf%m_deadcrootc_storage_to_fire => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1
    pcf%m_leafc_xfer_to_fire => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1
    pcf%m_frootc_xfer_to_fire => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1
    pcf%m_livestemc_xfer_to_fire => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1
    pcf%m_deadstemc_xfer_to_fire => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1
    pcf%m_livecrootc_xfer_to_fire => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1
    pcf%m_deadcrootc_xfer_to_fire => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1
    pcf%m_livestemc_to_fire => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1
    pcf%m_deadstemc_to_fire => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1
    pcf%m_deadstemc_to_litter_fire => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1
    pcf%m_livecrootc_to_fire => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1
    pcf%m_deadcrootc_to_fire => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1
    pcf%m_deadcrootc_to_litter_fire => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1
    pcf%m_gresp_storage_to_fire => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1
    pcf%m_gresp_xfer_to_fire => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1
    i4 = i4 + 1;
    ! end multi-pointing
    
    pcf%leafc_xfer_to_leafc => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcf%frootc_xfer_to_frootc => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcf%livestemc_xfer_to_livestemc => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcf%deadstemc_xfer_to_deadstemc => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcf%livecrootc_xfer_to_livecrootc => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcf%deadcrootc_xfer_to_deadcrootc => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcf%leafc_to_litter => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcf%frootc_to_litter => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcf%leaf_mr => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcf%froot_mr => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcf%livestem_mr => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcf%livecroot_mr => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcf%leaf_curmr => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcf%froot_curmr => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcf%livestem_curmr => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcf%livecroot_curmr => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcf%leaf_xsmr => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcf%froot_xsmr => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcf%livestem_xsmr => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcf%livecroot_xsmr => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcf%psnsun_to_cpool => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcf%psnshade_to_cpool => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcf%cpool_to_xsmrpool => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcf%cpool_to_leafc => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcf%cpool_to_leafc_storage => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcf%cpool_to_frootc => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcf%cpool_to_frootc_storage => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcf%cpool_to_livestemc => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcf%cpool_to_livestemc_storage => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcf%cpool_to_deadstemc => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcf%cpool_to_deadstemc_storage => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcf%cpool_to_livecrootc => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcf%cpool_to_livecrootc_storage => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcf%cpool_to_deadcrootc => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcf%cpool_to_deadcrootc_storage => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcf%cpool_to_gresp_storage => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcf%cpool_leaf_gr => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcf%cpool_leaf_storage_gr => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcf%transfer_leaf_gr => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcf%cpool_froot_gr => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcf%cpool_froot_storage_gr => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcf%transfer_froot_gr => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcf%cpool_livestem_gr => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcf%cpool_livestem_storage_gr => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcf%transfer_livestem_gr => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcf%cpool_deadstem_gr => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcf%cpool_deadstem_storage_gr => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcf%transfer_deadstem_gr => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcf%cpool_livecroot_gr => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcf%cpool_livecroot_storage_gr => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcf%transfer_livecroot_gr => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcf%cpool_deadcroot_gr => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcf%cpool_deadcroot_storage_gr => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcf%transfer_deadcroot_gr => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcf%leafc_storage_to_xfer => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcf%frootc_storage_to_xfer => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcf%livestemc_storage_to_xfer => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcf%deadstemc_storage_to_xfer => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcf%livecrootc_storage_to_xfer => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcf%deadcrootc_storage_to_xfer => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcf%gresp_storage_to_xfer => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcf%livestemc_to_deadstemc => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcf%livecrootc_to_deadcrootc => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcf%gpp => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcf%mr => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcf%current_gr => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcf%transfer_gr => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcf%storage_gr => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcf%gr => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcf%ar => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcf%rr => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcf%npp => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcf%agnpp => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcf%bgnpp => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcf%litfall => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcf%vegfire => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcf%wood_harvestc => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcf%pft_cinputs => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcf%pft_coutputs => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcf%pft_fire_closs => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    if ( crop_prog )then
       pcf%xsmrpool_to_atm => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
       pcf%grainc_xfer_to_grainc => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
       pcf%livestemc_to_litter => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
       pcf%grainc_to_food => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
       pcf%cpool_to_grainc => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
       pcf%cpool_to_grainc_storage => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
       pcf%cpool_grain_gr => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
       pcf%cpool_grain_storage_gr => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
       pcf%transfer_grain_gr => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
       pcf%grainc_storage_to_xfer => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    end if
#ifdef CN
    pcf%frootc_alloc=> oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcf%frootc_loss => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcf%leafc_alloc => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcf%leafc_loss  => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcf%woodc_alloc => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pcf%woodc_loss  => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
#endif

  end subroutine init_pft_cflux_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_nflux_type
!
! !INTERFACE:
  subroutine init_pft_nflux_type(beg, end, pnf)
!
! !DESCRIPTION:
! Initialize pft nitrogen flux variables
!
! !USES:
    use surfrdMod , only : crop_prog
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (pft_nflux_type), intent(inout) :: pnf
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!-----------------------------------------------------------------------
! point all m_***, hrv_*** to the same chunk, which will always remain 0. 
! the actual reduction due to mortality and harvesting will be done in SUBROUTINE ****
    pnf%m_leafn_to_litter => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1
    pnf%m_frootn_to_litter => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1
    pnf%m_leafn_storage_to_litter => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1
    pnf%m_frootn_storage_to_litter => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1
    pnf%m_livestemn_storage_to_litter => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1
    pnf%m_deadstemn_storage_to_litter => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1
    pnf%m_livecrootn_storage_to_litter => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1
    pnf%m_deadcrootn_storage_to_litter => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1
    pnf%m_leafn_xfer_to_litter => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1
    pnf%m_frootn_xfer_to_litter => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1
    pnf%m_livestemn_xfer_to_litter => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1
    pnf%m_deadstemn_xfer_to_litter => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1
    pnf%m_livecrootn_xfer_to_litter => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1
    pnf%m_deadcrootn_xfer_to_litter => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1
    pnf%m_livestemn_to_litter => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1
    pnf%m_deadstemn_to_litter => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1
    pnf%m_livecrootn_to_litter => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1
    pnf%m_deadcrootn_to_litter => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1
    pnf%m_retransn_to_litter => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1
    pnf%hrv_leafn_to_litter => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1             
    pnf%hrv_frootn_to_litter => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1            
    pnf%hrv_leafn_storage_to_litter => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1     
    pnf%hrv_frootn_storage_to_litter => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1    
    pnf%hrv_livestemn_storage_to_litter => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1 
    pnf%hrv_deadstemn_storage_to_litter => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1 
    pnf%hrv_livecrootn_storage_to_litter => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1
    pnf%hrv_deadcrootn_storage_to_litter => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1
    pnf%hrv_leafn_xfer_to_litter => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1        
    pnf%hrv_frootn_xfer_to_litter => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1       
    pnf%hrv_livestemn_xfer_to_litter => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1    
    pnf%hrv_deadstemn_xfer_to_litter => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1    
    pnf%hrv_livecrootn_xfer_to_litter => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1   
    pnf%hrv_deadcrootn_xfer_to_litter => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1   
    pnf%hrv_livestemn_to_litter => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1         
    pnf%hrv_deadstemn_to_prod10n => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1        
    pnf%hrv_deadstemn_to_prod100n => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1       
    pnf%hrv_livecrootn_to_litter => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1        
    pnf%hrv_deadcrootn_to_litter => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1        
    pnf%hrv_retransn_to_litter => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1           
    pnf%m_leafn_to_fire => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1
    pnf%m_frootn_to_fire => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1
    pnf%m_leafn_storage_to_fire => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1
    pnf%m_frootn_storage_to_fire => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1
    pnf%m_livestemn_storage_to_fire => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1
    pnf%m_deadstemn_storage_to_fire => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1
    pnf%m_livecrootn_storage_to_fire => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1
    pnf%m_deadcrootn_storage_to_fire => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1
    pnf%m_leafn_xfer_to_fire => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1
    pnf%m_frootn_xfer_to_fire => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1
    pnf%m_livestemn_xfer_to_fire => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1
    pnf%m_deadstemn_xfer_to_fire => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1
    pnf%m_livecrootn_xfer_to_fire => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1
    pnf%m_deadcrootn_xfer_to_fire => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1
    pnf%m_livestemn_to_fire => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1
    pnf%m_deadstemn_to_fire => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1
    pnf%m_deadstemn_to_litter_fire => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1
    pnf%m_livecrootn_to_fire => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1
    pnf%m_deadcrootn_to_fire => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1
    pnf%m_deadcrootn_to_litter_fire => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1
    pnf%m_retransn_to_fire => oneD_ptr(nn,R4(:,:,:,i4)); i4=i4+1
    i4 = i4 + 1
    ! done multi-pointing
    
    pnf%leafn_xfer_to_leafn => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pnf%frootn_xfer_to_frootn => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pnf%livestemn_xfer_to_livestemn => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pnf%deadstemn_xfer_to_deadstemn => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pnf%livecrootn_xfer_to_livecrootn => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pnf%deadcrootn_xfer_to_deadcrootn => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pnf%leafn_to_litter => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pnf%leafn_to_retransn => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pnf%frootn_to_litter => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pnf%retransn_to_npool => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pnf%sminn_to_npool => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pnf%npool_to_leafn => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pnf%npool_to_leafn_storage => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pnf%npool_to_frootn => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pnf%npool_to_frootn_storage => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pnf%npool_to_livestemn => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pnf%npool_to_livestemn_storage => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pnf%npool_to_deadstemn => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pnf%npool_to_deadstemn_storage => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pnf%npool_to_livecrootn => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pnf%npool_to_livecrootn_storage => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pnf%npool_to_deadcrootn => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pnf%npool_to_deadcrootn_storage => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pnf%leafn_storage_to_xfer => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pnf%frootn_storage_to_xfer => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pnf%livestemn_storage_to_xfer => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pnf%deadstemn_storage_to_xfer => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pnf%livecrootn_storage_to_xfer => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pnf%deadcrootn_storage_to_xfer => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pnf%livestemn_to_deadstemn => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pnf%livestemn_to_retransn => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pnf%livecrootn_to_deadcrootn => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pnf%livecrootn_to_retransn => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pnf%ndeploy => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pnf%pft_ninputs => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pnf%pft_noutputs => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pnf%wood_harvestn => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    pnf%pft_fire_nloss => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    if ( crop_prog )then
       pnf%grainn_xfer_to_grainn => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
       pnf%livestemn_to_litter => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
       pnf%grainn_to_food => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
       pnf%npool_to_grainn => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
       pnf%npool_to_grainn_storage => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
       pnf%grainn_storage_to_xfer => oneD_ptr(nn,R4(:,:,:,i4));i4=i4+1
    end if
  end subroutine init_pft_nflux_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_vflux_type
!
! !INTERFACE:
  subroutine init_pft_vflux_type(beg, end, pvf)
!
! !DESCRIPTION:
! Initialize pft VOC flux variables
!
    use clm_varcon, only : spval
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (pft_vflux_type), intent(inout) :: pvf
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
! (heald, 08/06)
!
!EOP
!------------------------------------------------------------------------

    allocate(pvf%vocflx_tot(beg:end))
    allocate(pvf%vocflx(beg:end,1:nvoc))
    allocate(pvf%vocflx_1(beg:end))
    allocate(pvf%vocflx_2(beg:end))
    allocate(pvf%vocflx_3(beg:end))
    allocate(pvf%vocflx_4(beg:end))
    allocate(pvf%vocflx_5(beg:end))
    allocate(pvf%Eopt_out(beg:end))
    allocate(pvf%topt_out(beg:end))
    allocate(pvf%alpha_out(beg:end))
    allocate(pvf%cp_out(beg:end))
    allocate(pvf%para_out(beg:end))
    allocate(pvf%par24a_out(beg:end))
    allocate(pvf%par240a_out(beg:end))
    allocate(pvf%paru_out(beg:end))
    allocate(pvf%par24u_out(beg:end))
    allocate(pvf%par240u_out(beg:end))
    allocate(pvf%gamma_out(beg:end))
    allocate(pvf%gammaL_out(beg:end))
    allocate(pvf%gammaT_out(beg:end))
    allocate(pvf%gammaP_out(beg:end))
    allocate(pvf%gammaA_out(beg:end))
    allocate(pvf%gammaS_out(beg:end))

    pvf%vocflx_tot(beg:end) = spval
    pvf%vocflx(beg:end,1:nvoc) = spval
    pvf%vocflx_1(beg:end) = spval
    pvf%vocflx_2(beg:end) = spval
    pvf%vocflx_3(beg:end) = spval
    pvf%vocflx_4(beg:end) = spval
    pvf%vocflx_5(beg:end) = spval
    pvf%Eopt_out(beg:end) = nan
    pvf%topt_out(beg:end) = nan
    pvf%alpha_out(beg:end) = nan
    pvf%cp_out(beg:end) = nan
    pvf%para_out(beg:end) = nan
    pvf%par24a_out(beg:end) = nan
    pvf%par240a_out(beg:end) = nan
    pvf%paru_out(beg:end) = nan
    pvf%par24u_out(beg:end) = nan
    pvf%par240u_out(beg:end) = nan
    pvf%gamma_out(beg:end) = nan
    pvf%gammaL_out(beg:end) = nan
    pvf%gammaT_out(beg:end) = nan
    pvf%gammaP_out(beg:end) = nan
    pvf%gammaA_out(beg:end) = nan
    pvf%gammaS_out(beg:end) = nan

  end subroutine init_pft_vflux_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_dflux_type
!
! !INTERFACE:
  subroutine init_pft_dflux_type(beg, end, pdf)
!
! !DESCRIPTION:
! Initialize pft dust flux variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (pft_dflux_type), intent(inout):: pdf
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(pdf%flx_mss_vrt_dst(beg:end,1:ndst))
    allocate(pdf%flx_mss_vrt_dst_tot(beg:end))
    allocate(pdf%vlc_trb(beg:end,1:ndst))
    allocate(pdf%vlc_trb_1(beg:end))
    allocate(pdf%vlc_trb_2(beg:end))
    allocate(pdf%vlc_trb_3(beg:end))
    allocate(pdf%vlc_trb_4(beg:end))

    pdf%flx_mss_vrt_dst(beg:end,1:ndst) = nan
    pdf%flx_mss_vrt_dst_tot(beg:end) = nan
    pdf%vlc_trb(beg:end,1:ndst) = nan
    pdf%vlc_trb_1(beg:end) = nan
    pdf%vlc_trb_2(beg:end) = nan
    pdf%vlc_trb_3(beg:end) = nan
    pdf%vlc_trb_4(beg:end) = nan

  end subroutine init_pft_dflux_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_depvd_type
!
! !INTERFACE:
  subroutine init_pft_depvd_type(beg, end, pdd)
! !DESCRIPTION:
! Initialize pft dep velocity variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (pft_depvd_type), intent(inout):: pdd
    integer :: i
!
#if (defined WTF)
    use seq_drydep_mod, only:  n_drydep, drydep_method, DD_XLND
!
! !REVISION HISTORY:
! Created by James Sulzman 541-929-6183
!
!EOP
!------------------------------------------------------------------------

    if ( n_drydep > 0 .and. drydep_method == DD_XLND )then
       allocate(pdd%drydepvel(beg:end,n_drydep))
       pdd%drydepvel = nan
    end if
#endif
  end subroutine init_pft_depvd_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_column_pstate_type
!
! !INTERFACE:
  subroutine init_column_pstate_type(beg, end, cps)
!
! !DESCRIPTION:
! Initialize column physical state variables
!
! !USES:
    use clm_varcon, only : spval
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    integer nlev
    type (column_pstate_type), intent(inout):: cps
    integer,pointer:: ii3,jj3
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------
!CP Too much toll on memory, this is not a complete list, some have been removed (fields starting with mss...this can be added later)
    ii3 => i3; jj3=>j3
    R3C%cps(1) = i3
    nlev = nlevsno + nlevgrnd
    cps%snl => oneD_ptrI(mn,L3(:,:,j3)); j3 = j3 + 1 
    cps%isoicol => oneD_ptrI(mn,L3(:,:,j3)); j3 = j3 + 1 
    cps%bsw => twoD_ptr(mn,nlevgrnd,R3(:,:,i3)); i3 = i3 + nlevgrnd 
    cps%watsat => twoD_ptr(mn,nlevgrnd,R3(:,:,i3)); i3 = i3 + nlevgrnd 
    cps%watfc => twoD_ptr(mn,nlevgrnd,R3(:,:,i3)); i3 = i3 + nlevgrnd 
    cps%watdry => twoD_ptr(mn,nlevgrnd,R3(:,:,i3)); i3 = i3 + nlevgrnd 
    cps%watopt => twoD_ptr(mn,nlevgrnd,R3(:,:,i3)); i3 = i3 + nlevgrnd 
    cps%hksat => twoD_ptr(mn,nlevgrnd,R3(:,:,i3)); i3 = i3 + nlevgrnd 
    cps%sucsat => twoD_ptr(mn,nlevgrnd,R3(:,:,i3)); i3 = i3 + nlevgrnd 
    cps%csol => twoD_ptr(mn,nlevgrnd,R3(:,:,i3)); i3 = i3 + nlevgrnd 
    cps%tkmg => twoD_ptr(mn,nlevgrnd,R3(:,:,i3)); i3 = i3 + nlevgrnd 
    cps%tkdry => twoD_ptr(mn,nlevgrnd,R3(:,:,i3)); i3 = i3 + nlevgrnd 
    cps%tksatu => twoD_ptr(mn,nlevgrnd,R3(:,:,i3)); i3 = i3 + nlevgrnd 
    cps%smpmin => oneD_ptr(mn,R3(:,:,i3)); i3 = i3 + 1 
    cps%hkdepth => oneD_ptr(mn,R3(:,:,i3)); i3 = i3 + 1 
    cps%wtfact => oneD_ptr(mn,R3(:,:,i3)); i3 = i3 + 1 
    cps%fracice => twoD_ptr(mn,nlevgrnd,R3(:,:,i3)); i3 = i3 + nlevgrnd 
    cps%gwc_thr => oneD_ptr(mn,R3(:,:,i3)); i3 = i3 + 1 
    cps%mss_frc_cly_vld => oneD_ptr(mn,R3(:,:,i3)); i3 = i3 + 1 
    cps%mbl_bsn_fct => oneD_ptr(mn,R3(:,:,i3)); i3 = i3 + 1 
    
    allocate(cps%do_capsnow(beg:end))
    allocate(cps%snow_zeroed(beg:end)) !CP added field
    cps%do_capsnow = .false.
    cps%snow_zeroed = .false.
    
    cps%snowdp => oneD_ptr(mn,R3(:,:,i3)); i3 = i3 + 1 
    cps%frac_sno => oneD_ptr(mn,R3(:,:,i3)); i3 = i3 + 1 
    cps%zi => twoD_ptrO(mn,nlevgrnd,R3(:,:,i3),-nlevsno); i3 = i3 + nlev + 1 
    cps%dz => twoD_ptrO(mn,nlevgrnd,R3(:,:,i3),-nlevsno+1); i3 = i3 + nlev 
    cps%z  => twoD_ptrO(mn,nlevgrnd,R3(:,:,i3),-nlevsno+1); i3 = i3 + nlev 
    cps%frac_iceold => twoD_ptrO(mn,nlevgrnd,R3(:,:,i3),-nlevsno+1); i3 = i3 + nlev 
    
    cps%imelt => twoD_ptrIO(mn,nlevgrnd,L3(:,:,j3),-nlevsno+1); j3 = j3 + nlev 
    
    cps%eff_porosity => twoD_ptr(mn,nlevgrnd,R3(:,:,i3)); i3 = i3 + nlevgrnd 
    cps%emg => oneD_ptr(mn,R3(:,:,i3)); i3 = i3 + 1 
    cps%z0mg => oneD_ptr(mn,R3(:,:,i3)); i3 = i3 + 1 
    cps%z0hg => oneD_ptr(mn,R3(:,:,i3)); i3 = i3 + 1 
    cps%z0qg => oneD_ptr(mn,R3(:,:,i3)); i3 = i3 + 1 
    cps%htvp => oneD_ptr(mn,R3(:,:,i3)); i3 = i3 + 1 
    cps%beta => oneD_ptr(mn,R3(:,:,i3)); i3 = i3 + 1 
    cps%zii => oneD_ptr(mn,R3(:,:,i3)); i3 = i3 + 1 
    cps%albgrd => twoD_ptr(mn,numrad,R3(:,:,i3)); i3 = i3 + numrad
    cps%albgri => twoD_ptr(mn,numrad,R3(:,:,i3)); i3 = i3 + numrad
    cps%rootr_column => twoD_ptr(mn,nlevgrnd,R3(:,:,i3)); i3 = i3 + nlevgrnd
    cps%rootfr_road_perv => twoD_ptr(mn,nlevurb,R3(:,:,i3)); i3 = i3 + nlevurb
    cps%rootr_road_perv => twoD_ptr(mn,nlevurb,R3(:,:,i3)); i3 = i3 + nlevurb
    cps%wf => oneD_ptr(mn,R3(:,:,i3)); i3 = i3 + 1 
!   cps%xirrig(beg:end) = 0._r8
    cps%max_dayl => oneD_ptr(mn,R3(:,:,i3)); i3 = i3 + 1 
    cps%bsw2 => twoD_ptr(mn,nlevgrnd,R3(:,:,i3)); i3 = i3 + nlevgrnd 
    cps%psisat => twoD_ptr(mn,nlevgrnd,R3(:,:,i3)); i3 = i3 + nlevgrnd 
    cps%vwcsat => twoD_ptr(mn,nlevgrnd,R3(:,:,i3)); i3 = i3 + nlevgrnd 
    cps%soilpsi => twoD_ptr(mn,nlevgrnd,R3(:,:,i3)); i3 = i3 + nlevgrnd 
    cps%decl => oneD_ptr(mn,R3(:,:,i3)); i3 = i3 + 1 
    cps%coszen => oneD_ptr(mn,R3(:,:,i3)); i3 = i3 + 1 
    cps%fpi => oneD_ptr(mn,R3(:,:,i3)); i3 = i3 + 1 
    cps%fpg => oneD_ptr(mn,R3(:,:,i3)); i3 = i3 + 1 
    cps%annsum_counter => oneD_ptr(mn,R3(:,:,i3)); i3 = i3 + 1 
    cps%cannsum_npp => oneD_ptr(mn,R3(:,:,i3)); i3 = i3 + 1 
    cps%cannavg_t2m => oneD_ptr(mn,R3(:,:,i3)); i3 = i3 + 1 
    cps%me => oneD_ptr(mn,R3(:,:,i3)); i3 = i3 + 1 
    cps%fire_prob => oneD_ptr(mn,R3(:,:,i3)); i3 = i3 + 1 
    cps%mean_fire_prob => oneD_ptr(mn,R3(:,:,i3)); i3 = i3 + 1 
    cps%fireseasonl => oneD_ptr(mn,R3(:,:,i3)); i3 = i3 + 1 
    cps%farea_burned => oneD_ptr(mn,R3(:,:,i3)); i3 = i3 + 1 
    cps%ann_farea_burned => oneD_ptr(mn,R3(:,:,i3)); i3 = i3 + 1 
    cps%albsnd_hst => twoD_ptr(mn,numrad,R3(:,:,i3)); i3 = i3 + numrad
    cps%albsni_hst => twoD_ptr(mn,numrad,R3(:,:,i3)); i3 = i3 + numrad
    cps%albsod => twoD_ptr(mn,numrad,R3(:,:,i3)); i3 = i3 + numrad
    cps%albsoi => twoD_ptr(mn,numrad,R3(:,:,i3)); i3 = i3 + numrad
    
    cps%flx_absdv => twoD_ptrO(mn,1,R3(:,:,i3),-nlevsno+1); i3 = i3 + nlevsno + 1
    cps%flx_absdn => twoD_ptrO(mn,1,R3(:,:,i3),-nlevsno+1); i3 = i3 + nlevsno + 1
    cps%flx_absin => twoD_ptrO(mn,1,R3(:,:,i3),-nlevsno+1); i3 = i3 + nlevsno + 1
    cps%flx_absiv => twoD_ptrO(mn,1,R3(:,:,i3),-nlevsno+1); i3 = i3 + nlevsno + 1
    cps%snw_rds   => twoD_ptrO(mn,0,R3(:,:,i3),-nlevsno+1); i3 = i3 + nlevsno 
    !cps%flx_absin(beg:end,-nlevsno+1:1) 
    !cps%snw_rds(beg:end,-nlevsno+1:0) 
    
    cps%snw_rds_top => oneD_ptr(mn,R3(:,:,i3)); i3 = i3 + 1 
    cps%sno_liq_top => oneD_ptr(mn,R3(:,:,i3)); i3 = i3 + 1 
    cps%h2osno_top  => oneD_ptr(mn,R3(:,:,i3)); i3 = i3 + 1 
    
    !CP NOTICE, these fields do not have i3+1 (or i3+nlevsno) behind them, so they just share one piece of data, because they are not currently used.
    cps%mss_bcpho=> twoD_ptrO(mn,0,R3(:,:,i3),-nlevsno+1); i3 = i3 + nlevsno 
    cps%mss_bcphi=> twoD_ptrO(mn,0,R3(:,:,i3),-nlevsno+1); i3 = i3 + nlevsno 
    cps%mss_bctot=> twoD_ptrO(mn,0,R3(:,:,i3),-nlevsno+1); i3 = i3 + nlevsno 
    cps%mss_bc_col=> oneD_ptr(mn,R3(:,:,i3));   i3 = i3 + 1
    cps%mss_bc_top=> oneD_ptr(mn,R3(:,:,i3));   i3 = i3 + 1
    cps%mss_ocpho=> twoD_ptrO(mn,0,R3(:,:,i3),-nlevsno+1); i3 = i3 + nlevsno 
    cps%mss_ocphi=> twoD_ptrO(mn,0,R3(:,:,i3),-nlevsno+1); i3 = i3 + nlevsno 
    cps%mss_octot=> twoD_ptrO(mn,0,R3(:,:,i3),-nlevsno+1); i3 = i3 + nlevsno 
    cps%mss_oc_col=> oneD_ptr(mn,R3(:,:,i3));   i3 = i3 + 1
    cps%mss_oc_top=> oneD_ptr(mn,R3(:,:,i3));   i3 = i3 + 1
    cps%mss_dst1=> twoD_ptrO(mn,0,R3(:,:,i3),-nlevsno+1); i3 = i3 + nlevsno 
    cps%mss_dst2=> twoD_ptrO(mn,0,R3(:,:,i3),-nlevsno+1); i3 = i3 + nlevsno 
    cps%mss_dst3=> twoD_ptrO(mn,0,R3(:,:,i3),-nlevsno+1); i3 = i3 + nlevsno 
    cps%mss_dst4=> twoD_ptrO(mn,0,R3(:,:,i3),-nlevsno+1); i3 = i3 + nlevsno 
    cps%mss_dsttot=> twoD_ptrO(mn,0,R3(:,:,i3),-nlevsno+1); i3 = i3 + nlevsno 
    cps%mss_dst_col=> oneD_ptr(mn,R3(:,:,i3));   i3 = i3 + 1
    cps%mss_dst_top=> oneD_ptr(mn,R3(:,:,i3));   i3 = i3 + 1
    cps%h2osno_top=> oneD_ptr(mn,R3(:,:,i3));   i3 = i3 + 1
    cps%mss_cnc_bcphi=> twoD_ptrO(mn,0,R3(:,:,i3),-nlevsno+1); i3 = i3 + nlevsno 
    cps%mss_cnc_bcpho=> twoD_ptrO(mn,0,R3(:,:,i3),-nlevsno+1); i3 = i3 + nlevsno 
    cps%mss_cnc_ocphi=> twoD_ptrO(mn,0,R3(:,:,i3),-nlevsno+1); i3 = i3 + nlevsno 
    cps%mss_cnc_ocpho=> twoD_ptrO(mn,0,R3(:,:,i3),-nlevsno+1); i3 = i3 + nlevsno 
    cps%mss_cnc_dst1=> twoD_ptrO(mn,0,R3(:,:,i3),-nlevsno+1); i3 = i3 + nlevsno 
    cps%mss_cnc_dst2=> twoD_ptrO(mn,0,R3(:,:,i3),-nlevsno+1); i3 = i3 + nlevsno 
    cps%mss_cnc_dst3=> twoD_ptrO(mn,0,R3(:,:,i3),-nlevsno+1); i3 = i3 + nlevsno 
    cps%mss_cnc_dst4=> twoD_ptrO(mn,0,R3(:,:,i3),-nlevsno+1); i3 = i3 + nlevsno 
    
    cps%albgrd_pur => twoD_ptr(mn,numrad,R3(:,:,i3)); i3 = i3 + numrad
    cps%albgri_pur => twoD_ptr(mn,numrad,R3(:,:,i3)); i3 = i3 + numrad
    cps%albgrd_bc => twoD_ptr(mn,numrad,R3(:,:,i3)); i3 = i3 + numrad
    cps%albgri_bc => twoD_ptr(mn,numrad,R3(:,:,i3)); i3 = i3 + numrad
    cps%albgrd_oc => twoD_ptr(mn,numrad,R3(:,:,i3)); i3 = i3 + numrad
    cps%albgri_oc => twoD_ptr(mn,numrad,R3(:,:,i3)); i3 = i3 + numrad 
    cps%albgrd_dst => twoD_ptr(mn,numrad,R3(:,:,i3)); i3 = i3 + numrad
    cps%albgri_dst => twoD_ptr(mn,numrad,R3(:,:,i3)); i3 = i3 + numrad
    cps%dTdz_top => oneD_ptr(mn,R3(:,:,i3)); i3 = i3 + 1
    cps%snot_top => oneD_ptr(mn,R3(:,:,i3)); i3 = i3 + 1
    cps%irrig_rate => oneD_ptr(mn,R3(:,:,i3)); i3 = i3 + 1
    cps%n_irrig_steps_left => oneD_ptrI(mn,L3(:,:,j3)); j3 = j3 + 1   
    cps%forc_pbot => oneD_ptr(mn,R3(:,:,i3)); i3 = i3 + 1
    cps%forc_rho => oneD_ptr(mn,R3(:,:,i3)); i3 = i3 + 1
    cps%glc_topo => oneD_ptr(mn,R3(:,:,i3)); i3 = i3 + 1
    
    R3C%cps(2) = i3 - 1
    
  end subroutine init_column_pstate_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_column_estate_type
!
! !INTERFACE:
  subroutine init_column_estate_type(beg, end, ces)
!
! !DESCRIPTION:
! Initialize column energy state variables
!
! !USES:
    use clm_varcon, only : spval
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (column_estate_type), intent(inout):: ces
    integer nlev
    
    nlev = nlevsno + nlevgrnd
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    ces%t_grnd => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1    
    ces%t_grnd_u => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1  
    ces%t_grnd_r => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1  
    ces%dt_grnd => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1   
    
    ces%t_soisno => twoD_ptrO(mn,nlevgrnd,R3(:,:,i3),-nlevsno+1); i3 = i3 + nlev 
    !ces%t_soisno(beg:end,-nlevsno+1:nlevgrnd) = spval
    
    ces%t_soi_10cm => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1
    ces%t_lake => twoD_ptr(mn,nlevlak,R3(:,:,i3)); i3=i3+nlevlak
    !ces%t_lake(beg:end,1:nlevlak)            
    ces%tssbef => twoD_ptrO(mn,nlevgrnd,R3(:,:,i3),-nlevsno+1); i3 = i3 + nlev 
    ces%thv => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1       
    ces%hc_soi => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1    
    ces%hc_soisno => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1 
    ces%forc_t => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1 
    ces%forc_th => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1 

  end subroutine init_column_estate_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_column_wstate_type
!
! !INTERFACE:
  subroutine init_column_wstate_type(beg, end, cws)
!
! !DESCRIPTION:
! Initialize column water state variables
!
! !USES:
    use clm_varcon, only : spval
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (column_wstate_type), intent(inout):: cws !column water state
    integer nlev
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------
    nlev = nlevgrnd + nlevsno
    
    cws%h2osno => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1 
    
    cws%h2osoi_liq => twoD_ptrO(mn,nlevgrnd,R3(:,:,i3),-nlevsno+1); i3 = i3 + nlev 
    cws%h2osoi_ice => twoD_ptrO(mn,nlevgrnd,R3(:,:,i3),-nlevsno+1); i3 = i3 + nlev 
    
    cws%h2osoi_liqice_10cm => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1
    cws%h2osoi_vol => twoD_ptr(mn,nlevgrnd,R3(:,:,i3)); i3 = i3+nlevgrnd
    cws%h2osno_old => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1 
    cws%qg => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1 
    cws%dqgdT => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1 
    cws%snowice => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1 
    cws%snowliq => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1 
    cws%soilalpha => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1 
    cws%soilbeta => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1 
    cws%soilalpha_u => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1 
    cws%zwt => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1 
    cws%fcov => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1 
    cws%fsat => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1 
    cws%wa => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1 
    cws%wt => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1 
    cws%qcharge => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1 
    cws%smp_l => twoD_ptr(mn,nlevgrnd,R3(:,:,i3)); i3 = i3+nlevgrnd
    cws%hk_l => twoD_ptr(mn,nlevgrnd,R3(:,:,i3)); i3 = i3+nlevgrnd
    cws%forc_q => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1 

  end subroutine init_column_wstate_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_column_cstate_type
!
! !INTERFACE:
  subroutine init_column_cstate_type(beg, end, ccs)
!
! !DESCRIPTION:
! Initialize column carbon state variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (column_cstate_type), intent(inout):: ccs
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------
    ccs%soilc => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
    ccs%cwdc => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
    ccs%litr1c => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
    ccs%litr2c => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
    ccs%litr3c => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
    ccs%soil1c => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
    ccs%soil2c => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
    ccs%soil3c => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
    ccs%soil4c => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
    ccs%seedc => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
    ccs%col_ctrunc => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
    ccs%prod10c => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
    ccs%prod100c => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
    ccs%totprodc => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
    ccs%totlitc => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
    ccs%totsomc => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
    ccs%totecosysc => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
    ccs%totcolc => oneD_ptr(mn,R3(:,:,i3));i3=i3+1

  end subroutine init_column_cstate_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_column_nstate_type
!
! !INTERFACE:
  subroutine init_column_nstate_type(beg, end, cns)
!
! !DESCRIPTION:
! Initialize column nitrogen state variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (column_nstate_type), intent(inout):: cns
!
! !REVISION HISTORY:
! Created by Peter Thornton
!
!EOP
!------------------------------------------------------------------------
    cns%cwdn => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
    cns%litr1n => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
    cns%litr2n => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
    cns%litr3n => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
    cns%soil1n => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
    cns%soil2n => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
    cns%soil3n => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
    cns%soil4n => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
    cns%sminn => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
    cns%col_ntrunc => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
    cns%seedn => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
    cns%prod10n => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
    cns%prod100n => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
    cns%totprodn => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
    cns%totlitn => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
    cns%totsomn => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
    cns%totecosysn => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
    cns%totcoln => oneD_ptr(mn,R3(:,:,i3));i3=i3+1

  end subroutine init_column_nstate_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_column_eflux_type
!
! !INTERFACE:
  subroutine init_column_eflux_type(beg, end, cef)
!
! !DESCRIPTION:
! Initialize column energy flux variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (column_eflux_type), intent(inout):: cef
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    cef%eflx_snomelt => oneD_ptr(mn,R3(:,:,i3)); i3 = i3 + 1      
    cef%eflx_snomelt_u => oneD_ptr(mn,R3(:,:,i3)); i3 = i3 + 1      
    cef%eflx_snomelt_r => oneD_ptr(mn,R3(:,:,i3)); i3 = i3 + 1      
    cef%eflx_impsoil => oneD_ptr(mn,R3(:,:,i3)); i3 = i3 + 1      
    cef%eflx_fgr12 => oneD_ptr(mn,R3(:,:,i3)); i3 = i3 + 1        
    cef%eflx_building_heat => oneD_ptr(mn,R3(:,:,i3)); i3 = i3 + 1
    cef%eflx_urban_ac => oneD_ptr(mn,R3(:,:,i3)); i3 = i3 + 1
    cef%eflx_urban_heat => oneD_ptr(mn,R3(:,:,i3)); i3 = i3 + 1
    cef%eflx_bot => oneD_ptr(mn,R3(:,:,i3)); i3 = i3 + 1

  end subroutine init_column_eflux_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_column_wflux_type
!
! !INTERFACE:
  subroutine init_column_wflux_type(beg, end, cwf)
!
! !DESCRIPTION:
! Initialize column water flux variables
!
! !USES:
    use clm_varcon, only : spval
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (column_wflux_type), intent(inout):: cwf
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    cwf%qflx_infl => oneD_ptr(mn,R3(:,:,i3)); i3 = i3 + 1
    cwf%qflx_surf => oneD_ptr(mn,R3(:,:,i3)); i3 = i3 + 1
    cwf%qflx_drain => oneD_ptr(mn,R3(:,:,i3)); i3 = i3 + 1
    cwf%qflx_top_soil => oneD_ptr(mn,R3(:,:,i3)); i3 = i3 + 1
    cwf%qflx_snomelt => oneD_ptr(mn,R3(:,:,i3)); i3 = i3 + 1
    cwf%qflx_qrgwl => oneD_ptr(mn,R3(:,:,i3)); i3 = i3 + 1
    cwf%qflx_runoff => oneD_ptr(mn,R3(:,:,i3)); i3 = i3 + 1
    cwf%qflx_runoff_u => oneD_ptr(mn,R3(:,:,i3)); i3 = i3 + 1
    cwf%qflx_runoff_r => oneD_ptr(mn,R3(:,:,i3)); i3 = i3 + 1
    cwf%qmelt => oneD_ptr(mn,R3(:,:,i3)); i3 = i3 + 1
    cwf%h2ocan_loss => oneD_ptr(mn,R3(:,:,i3)); i3 = i3 + 1
    cwf%qflx_rsub_sat => oneD_ptr(mn,R3(:,:,i3)); i3 = i3 + 1
    cwf%flx_bc_dep_dry => oneD_ptr(mn,R3(:,:,i3)); i3 = i3 + 1
    cwf%flx_bc_dep_wet => oneD_ptr(mn,R3(:,:,i3)); i3 = i3 + 1
    cwf%flx_bc_dep_pho => oneD_ptr(mn,R3(:,:,i3)); i3 = i3 + 1
    cwf%flx_bc_dep_phi => oneD_ptr(mn,R3(:,:,i3)); i3 = i3 + 1
    cwf%flx_bc_dep => oneD_ptr(mn,R3(:,:,i3)); i3 = i3 + 1
    cwf%flx_oc_dep_dry => oneD_ptr(mn,R3(:,:,i3)); i3 = i3 + 1
    cwf%flx_oc_dep_wet => oneD_ptr(mn,R3(:,:,i3)); i3 = i3 + 1
    cwf%flx_oc_dep_pho => oneD_ptr(mn,R3(:,:,i3)); i3 = i3 + 1
    cwf%flx_oc_dep_phi => oneD_ptr(mn,R3(:,:,i3)); i3 = i3 + 1
    cwf%flx_oc_dep => oneD_ptr(mn,R3(:,:,i3)); i3 = i3 + 1
    cwf%flx_dst_dep_dry1 => oneD_ptr(mn,R3(:,:,i3)); i3 = i3 + 1
    cwf%flx_dst_dep_wet1 => oneD_ptr(mn,R3(:,:,i3)); i3 = i3 + 1
    cwf%flx_dst_dep_dry2 => oneD_ptr(mn,R3(:,:,i3)); i3 = i3 + 1
    cwf%flx_dst_dep_wet2 => oneD_ptr(mn,R3(:,:,i3)); i3 = i3 + 1
    cwf%flx_dst_dep_dry3 => oneD_ptr(mn,R3(:,:,i3)); i3 = i3 + 1
    cwf%flx_dst_dep_wet3 => oneD_ptr(mn,R3(:,:,i3)); i3 = i3 + 1
    cwf%flx_dst_dep_dry4 => oneD_ptr(mn,R3(:,:,i3)); i3 = i3 + 1
    cwf%flx_dst_dep_wet4 => oneD_ptr(mn,R3(:,:,i3)); i3 = i3 + 1
    cwf%flx_dst_dep => oneD_ptr(mn,R3(:,:,i3)); i3 = i3 + 1
    !cwf%qflx_snofrz_lyr(beg:end,-nlevsno+1:0) = spval
    cwf%qflx_snofrz_lyr => twoD_ptrO(mn,0,R3(:,:,i3),-nlevsno+1); i3 = i3 + nlevsno
    cwf%qflx_glcice => oneD_ptr(mn,R3(:,:,i3)); i3 = i3 + 1
    cwf%qflx_irrig => oneD_ptr(mn,R3(:,:,i3)); i3 = i3 + 1
    cwf%glc_rofi => oneD_ptr(mn,R3(:,:,i3)); i3 = i3 + 1
    cwf%glc_rofl => oneD_ptr(mn,R3(:,:,i3)); i3 = i3 + 1

  end subroutine init_column_wflux_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_column_cflux_type
!
! !INTERFACE:
  subroutine init_column_cflux_type(beg, end, ccf)
!
! !DESCRIPTION:
! Initialize column carbon flux variables
!
! !USES:
    use surfrdMod , only : crop_prog
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (column_cflux_type), intent(inout):: ccf
!
! !REVISION HISTORY:
! Created by Peter Thornton
!
!EOP
!------------------------------------------------------------------------
    ccf%m_leafc_to_litr1c => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1               
    ccf%m_leafc_to_litr2c => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1               
    ccf%m_leafc_to_litr3c => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1               
    ccf%m_frootc_to_litr1c => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1              
    ccf%m_frootc_to_litr2c => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1              
    ccf%m_frootc_to_litr3c => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1              
    ccf%m_leafc_storage_to_litr1c => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1       
    ccf%m_frootc_storage_to_litr1c => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1      
    ccf%m_livestemc_storage_to_litr1c => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1   
    ccf%m_deadstemc_storage_to_litr1c => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1   
    ccf%m_livecrootc_storage_to_litr1c => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1  
    ccf%m_deadcrootc_storage_to_litr1c => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1  
    ccf%m_leafc_xfer_to_litr1c => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1          
    ccf%m_frootc_xfer_to_litr1c => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1         
    ccf%m_livestemc_xfer_to_litr1c => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1      
    ccf%m_deadstemc_xfer_to_litr1c => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1      
    ccf%m_livecrootc_xfer_to_litr1c => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1     
    ccf%m_deadcrootc_xfer_to_litr1c => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1     
    ccf%m_livestemc_to_cwdc => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1             
    ccf%m_deadstemc_to_cwdc => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1             
    ccf%m_livecrootc_to_cwdc => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1            
    ccf%m_deadcrootc_to_cwdc => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1            
    ccf%m_gresp_storage_to_litr1c => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1       
    ccf%m_gresp_xfer_to_litr1c => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1          
    ccf%m_deadstemc_to_cwdc_fire => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1        
    ccf%m_deadcrootc_to_cwdc_fire => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1   
    ccf%hrv_leafc_to_litr1c => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1                          
    ccf%hrv_leafc_to_litr2c => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1                          
    ccf%hrv_leafc_to_litr3c => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1                          
    ccf%hrv_frootc_to_litr1c => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1                        
    ccf%hrv_frootc_to_litr2c => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1                        
    ccf%hrv_frootc_to_litr3c => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1                        
    ccf%hrv_livestemc_to_cwdc => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1                      
    ccf%hrv_deadstemc_to_prod10c => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1                
    ccf%hrv_deadstemc_to_prod100c => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1              
    ccf%hrv_livecrootc_to_cwdc => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1                    
    ccf%hrv_deadcrootc_to_cwdc => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1                    
    ccf%hrv_leafc_storage_to_litr1c => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1          
    ccf%hrv_frootc_storage_to_litr1c => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1        
    ccf%hrv_livestemc_storage_to_litr1c => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1  
    ccf%hrv_deadstemc_storage_to_litr1c => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1  
    ccf%hrv_livecrootc_storage_to_litr1c => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1
    ccf%hrv_deadcrootc_storage_to_litr1c => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1
    ccf%hrv_gresp_storage_to_litr1c => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1          
    ccf%hrv_leafc_xfer_to_litr1c => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1                
    ccf%hrv_frootc_xfer_to_litr1c => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1              
    ccf%hrv_livestemc_xfer_to_litr1c => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1        
    ccf%hrv_deadstemc_xfer_to_litr1c => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1        
    ccf%hrv_livecrootc_xfer_to_litr1c => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1      
    ccf%hrv_deadcrootc_xfer_to_litr1c => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1      
    ccf%hrv_gresp_xfer_to_litr1c => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1                
    ccf%m_litr1c_to_fire => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1                
    ccf%m_litr2c_to_fire => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1                
    ccf%m_litr3c_to_fire => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1                
    ccf%m_cwdc_to_fire => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1                  
    i3 = i3 + 1
    ! done multi-pointing
    
    if ( crop_prog )then
       ccf%grainc_to_litr1c => oneD_ptr(mn,R3(:,:,i3)); i3 = i3 + 1   
       ccf%grainc_to_litr2c => oneD_ptr(mn,R3(:,:,i3)); i3 = i3 + 1  
       ccf%grainc_to_litr3c => oneD_ptr(mn,R3(:,:,i3)); i3 = i3 + 1  
       ccf%livestemc_to_litr1c => oneD_ptr(mn,R3(:,:,i3));  i3 = i3 + 1
       ccf%livestemc_to_litr2c => oneD_ptr(mn,R3(:,:,i3));  i3 = i3 + 1 
       ccf%livestemc_to_litr3c => oneD_ptr(mn,R3(:,:,i3));  i3 = i3 + 1 
    end if    
    ccf%leafc_to_litr1c => oneD_ptr(mn,R3(:,:,i3));i3=i3+1 
    ccf%leafc_to_litr2c => oneD_ptr(mn,R3(:,:,i3));i3=i3+1 
    ccf%leafc_to_litr3c => oneD_ptr(mn,R3(:,:,i3));i3=i3+1 
    ccf%frootc_to_litr1c => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
    ccf%frootc_to_litr2c => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
    ccf%frootc_to_litr3c => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
    ccf%cwdc_to_litr2c => oneD_ptr(mn,R3(:,:,i3));i3=i3+1  
    ccf%cwdc_to_litr3c => oneD_ptr(mn,R3(:,:,i3));i3=i3+1  
    ccf%litr1_hr => oneD_ptr(mn,R3(:,:,i3));i3=i3+1        
    ccf%litr1c_to_soil1c => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
    ccf%litr2_hr => oneD_ptr(mn,R3(:,:,i3));i3=i3+1        
    ccf%litr2c_to_soil2c => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
    ccf%litr3_hr => oneD_ptr(mn,R3(:,:,i3));i3=i3+1        
    ccf%litr3c_to_soil3c => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
    ccf%soil1_hr => oneD_ptr(mn,R3(:,:,i3));i3=i3+1        
    ccf%soil1c_to_soil2c => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
    ccf%soil2_hr => oneD_ptr(mn,R3(:,:,i3));i3=i3+1        
    ccf%soil2c_to_soil3c => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
    ccf%soil3_hr => oneD_ptr(mn,R3(:,:,i3));i3=i3+1        
    ccf%soil3c_to_soil4c => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
    ccf%soil4_hr => oneD_ptr(mn,R3(:,:,i3));i3=i3+1        
#if (defined CN)
    ccf%dwt_seedc_to_leaf=> oneD_ptr(mn,R3(:,:,i3));i3=i3+1       
    ccf%dwt_seedc_to_deadstem=> oneD_ptr(mn,R3(:,:,i3));i3=i3+1   
    ccf%dwt_conv_cflux=> oneD_ptr(mn,R3(:,:,i3));i3=i3+1          
    ccf%dwt_prod10c_gain=> oneD_ptr(mn,R3(:,:,i3));i3=i3+1        
    ccf%dwt_prod100c_gain=> oneD_ptr(mn,R3(:,:,i3));i3=i3+1       
    ccf%dwt_frootc_to_litr1c=> oneD_ptr(mn,R3(:,:,i3));i3=i3+1    
    ccf%dwt_frootc_to_litr2c=> oneD_ptr(mn,R3(:,:,i3));i3=i3+1    
    ccf%dwt_frootc_to_litr3c=> oneD_ptr(mn,R3(:,:,i3));i3=i3+1    
    ccf%dwt_livecrootc_to_cwdc=> oneD_ptr(mn,R3(:,:,i3));i3=i3+1  
    ccf%dwt_deadcrootc_to_cwdc=> oneD_ptr(mn,R3(:,:,i3));i3=i3+1  
    ccf%dwt_closs=> oneD_ptr(mn,R3(:,:,i3));i3=i3+1               
    ccf%landuseflux=> oneD_ptr(mn,R3(:,:,i3));i3=i3+1             
    ccf%landuptake=> oneD_ptr(mn,R3(:,:,i3));i3=i3+1              
    ccf%prod10c_loss=> oneD_ptr(mn,R3(:,:,i3));i3=i3+1            
    ccf%prod100c_loss=> oneD_ptr(mn,R3(:,:,i3));i3=i3+1           
    ccf%product_closs=> oneD_ptr(mn,R3(:,:,i3));i3=i3+1           
#endif
    ccf%lithr => oneD_ptr(mn,R3(:,:,i3));i3=i3+1           
    ccf%somhr => oneD_ptr(mn,R3(:,:,i3));i3=i3+1           
    ccf%hr => oneD_ptr(mn,R3(:,:,i3));i3=i3+1              
    ccf%sr => oneD_ptr(mn,R3(:,:,i3));i3=i3+1              
    ccf%er => oneD_ptr(mn,R3(:,:,i3));i3=i3+1              
    ccf%litfire => oneD_ptr(mn,R3(:,:,i3));i3=i3+1         
    ccf%somfire => oneD_ptr(mn,R3(:,:,i3));i3=i3+1         
    ccf%totfire => oneD_ptr(mn,R3(:,:,i3));i3=i3+1         
    ccf%nep => oneD_ptr(mn,R3(:,:,i3));i3=i3+1             
    ccf%nbp => oneD_ptr(mn,R3(:,:,i3));i3=i3+1             
    ccf%nee => oneD_ptr(mn,R3(:,:,i3));i3=i3+1             
    ccf%col_cinputs => oneD_ptr(mn,R3(:,:,i3));i3=i3+1     
    ccf%col_coutputs => oneD_ptr(mn,R3(:,:,i3));i3=i3+1    
    ccf%col_fire_closs => oneD_ptr(mn,R3(:,:,i3));i3=i3+1  
#if (defined CN)
    ccf%cwdc_hr => oneD_ptr(mn,R3(:,:,i3));i3=i3+1 
    ccf%cwdc_loss => oneD_ptr(mn,R3(:,:,i3));i3=i3+1 
    ccf%litterc_loss => oneD_ptr(mn,R3(:,:,i3));i3=i3+1 
#endif

  end subroutine init_column_cflux_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_column_nflux_type
!
! !INTERFACE:
  subroutine init_column_nflux_type(beg, end, cnf)
!
! !DESCRIPTION:
! Initialize column nitrogen flux variables
!
! !USES:
    use surfrdMod , only : crop_prog
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (column_nflux_type), intent(inout):: cnf
!
! !REVISION HISTORY:
! Created by Peter Thornton
!
!EOP
!------------------------------------------------------------------------
    cnf%ndep_to_sminn => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
    cnf%nfix_to_sminn => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
    
    cnf%m_leafn_to_litr1n => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1
    cnf%m_leafn_to_litr2n => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1
    cnf%m_leafn_to_litr3n => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1
    cnf%m_frootn_to_litr1n => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1
    cnf%m_frootn_to_litr2n => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1
    cnf%m_frootn_to_litr3n => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1
    cnf%m_leafn_storage_to_litr1n => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1
    cnf%m_frootn_storage_to_litr1n => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1
    cnf%m_livestemn_storage_to_litr1n => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1
    cnf%m_deadstemn_storage_to_litr1n => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1
    cnf%m_livecrootn_storage_to_litr1n => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1
    cnf%m_deadcrootn_storage_to_litr1n => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1
    cnf%m_leafn_xfer_to_litr1n => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1
    cnf%m_frootn_xfer_to_litr1n => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1
    cnf%m_livestemn_xfer_to_litr1n => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1
    cnf%m_deadstemn_xfer_to_litr1n => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1
    cnf%m_livecrootn_xfer_to_litr1n => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1
    cnf%m_deadcrootn_xfer_to_litr1n => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1
    cnf%m_livestemn_to_cwdn => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1
    cnf%m_deadstemn_to_cwdn => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1
    cnf%m_livecrootn_to_cwdn => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1
    cnf%m_deadcrootn_to_cwdn => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1
    cnf%m_retransn_to_litr1n => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1
    cnf%hrv_leafn_to_litr1n => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1             
    cnf%hrv_leafn_to_litr2n => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1             
    cnf%hrv_leafn_to_litr3n => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1             
    cnf%hrv_frootn_to_litr1n => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1            
    cnf%hrv_frootn_to_litr2n => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1            
    cnf%hrv_frootn_to_litr3n => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1            
    cnf%hrv_livestemn_to_cwdn => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1           
    cnf%hrv_deadstemn_to_prod10n => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1        
    cnf%hrv_deadstemn_to_prod100n => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1       
    cnf%hrv_livecrootn_to_cwdn => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1          
    cnf%hrv_deadcrootn_to_cwdn => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1          
    cnf%hrv_retransn_to_litr1n => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1          
    cnf%hrv_leafn_storage_to_litr1n => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1     
    cnf%hrv_frootn_storage_to_litr1n => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1    
    cnf%hrv_livestemn_storage_to_litr1n => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1 
    cnf%hrv_deadstemn_storage_to_litr1n => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1 
    cnf%hrv_livecrootn_storage_to_litr1n => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1
    cnf%hrv_deadcrootn_storage_to_litr1n => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1
    cnf%hrv_leafn_xfer_to_litr1n => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1        
    cnf%hrv_frootn_xfer_to_litr1n => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1       
    cnf%hrv_livestemn_xfer_to_litr1n => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1    
    cnf%hrv_deadstemn_xfer_to_litr1n => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1    
    cnf%hrv_livecrootn_xfer_to_litr1n => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1   
    cnf%hrv_deadcrootn_xfer_to_litr1n => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1   
    cnf%m_deadstemn_to_cwdn_fire => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1
    cnf%m_deadcrootn_to_cwdn_fire => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1
    cnf%m_litr1n_to_fire => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1
    cnf%m_litr2n_to_fire => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1
    cnf%m_litr3n_to_fire => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1
    cnf%m_cwdn_to_fire => oneD_ptr(mn,R3(:,:,i3)); i3=i3+1
    
    i3 = i3 + 1
    if ( crop_prog )then
       cnf%grainn_to_litr1n => oneD_ptr(mn,R3(:,:,i3));    i3 = i3 + 1
       cnf%grainn_to_litr2n => oneD_ptr(mn,R3(:,:,i3));    i3 = i3 + 1
       cnf%grainn_to_litr3n => oneD_ptr(mn,R3(:,:,i3));    i3 = i3 + 1
       cnf%livestemn_to_litr1n => oneD_ptr(mn,R3(:,:,i3));    i3 = i3 + 1
       cnf%livestemn_to_litr2n => oneD_ptr(mn,R3(:,:,i3));    i3 = i3 + 1
       cnf%livestemn_to_litr3n => oneD_ptr(mn,R3(:,:,i3));    i3 = i3 + 1
    end if
    ! done multi-pointing
    
    cnf%leafn_to_litr1n => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
    cnf%leafn_to_litr2n => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
    cnf%leafn_to_litr3n => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
    cnf%frootn_to_litr1n => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
    cnf%frootn_to_litr2n => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
    cnf%frootn_to_litr3n => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
    cnf%cwdn_to_litr2n => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
    cnf%cwdn_to_litr3n => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
    cnf%litr1n_to_soil1n => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
    cnf%sminn_to_soil1n_l1 => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
    cnf%litr2n_to_soil2n => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
    cnf%sminn_to_soil2n_l2 => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
    cnf%litr3n_to_soil3n => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
    cnf%sminn_to_soil3n_l3 => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
    cnf%soil1n_to_soil2n => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
    cnf%sminn_to_soil2n_s1 => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
    cnf%soil2n_to_soil3n => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
    cnf%sminn_to_soil3n_s2 => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
    cnf%soil3n_to_soil4n => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
    cnf%sminn_to_soil4n_s3 => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
    cnf%soil4n_to_sminn => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
    cnf%sminn_to_denit_l1s1 => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
    cnf%sminn_to_denit_l2s2 => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
    cnf%sminn_to_denit_l3s3 => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
    cnf%sminn_to_denit_s1s2 => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
    cnf%sminn_to_denit_s2s3 => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
    cnf%sminn_to_denit_s3s4 => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
    cnf%sminn_to_denit_s4 => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
    cnf%sminn_to_denit_excess => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
    cnf%sminn_leached => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
    cnf%dwt_seedn_to_leaf => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
    cnf%dwt_seedn_to_deadstem => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
    cnf%dwt_conv_nflux => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
    cnf%dwt_prod10n_gain => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
    cnf%dwt_prod100n_gain => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
    cnf%dwt_frootn_to_litr1n => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
    cnf%dwt_frootn_to_litr2n => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
    cnf%dwt_frootn_to_litr3n => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
    cnf%dwt_livecrootn_to_cwdn => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
    cnf%dwt_deadcrootn_to_cwdn => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
    cnf%dwt_nloss => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
    cnf%prod10n_loss => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
    cnf%prod100n_loss => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
    cnf%product_nloss => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
    cnf%potential_immob => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
    cnf%actual_immob => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
    cnf%sminn_to_plant => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
    cnf%supplement_to_sminn => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
    cnf%gross_nmin => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
    cnf%net_nmin => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
    cnf%denit => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
    cnf%col_ninputs => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
    cnf%col_noutputs => oneD_ptr(mn,R3(:,:,i3));i3=i3+1
    cnf%col_fire_nloss => oneD_ptr(mn,R3(:,:,i3));i3=i3+1

  end subroutine init_column_nflux_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_landunit_pstate_type
!
! !INTERFACE:
  subroutine init_landunit_pstate_type(beg, end, lps)
!
! !DESCRIPTION:
! Initialize landunit physical state variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (landunit_pstate_type), intent(inout):: lps
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(lps%t_building(beg:end))
    allocate(lps%t_building_max(beg:end))
    allocate(lps%t_building_min(beg:end))
    allocate(lps%tk_wall(beg:end,nlevurb))
    allocate(lps%tk_roof(beg:end,nlevurb))
    allocate(lps%tk_improad(beg:end,nlevgrnd))
    allocate(lps%cv_wall(beg:end,nlevurb))
    allocate(lps%cv_roof(beg:end,nlevurb))
    allocate(lps%cv_improad(beg:end,nlevgrnd))
    allocate(lps%thick_wall(beg:end))
    allocate(lps%thick_roof(beg:end))
    allocate(lps%nlev_improad(beg:end))
    allocate(lps%vf_sr(beg:end))
    allocate(lps%vf_wr(beg:end))
    allocate(lps%vf_sw(beg:end))
    allocate(lps%vf_rw(beg:end))
    allocate(lps%vf_ww(beg:end))
    allocate(lps%taf(beg:end))
    allocate(lps%qaf(beg:end))
    allocate(lps%sabs_roof_dir(beg:end,1:numrad))
    allocate(lps%sabs_roof_dif(beg:end,1:numrad))
    allocate(lps%sabs_sunwall_dir(beg:end,1:numrad))
    allocate(lps%sabs_sunwall_dif(beg:end,1:numrad))
    allocate(lps%sabs_shadewall_dir(beg:end,1:numrad))
    allocate(lps%sabs_shadewall_dif(beg:end,1:numrad))
    allocate(lps%sabs_improad_dir(beg:end,1:numrad))
    allocate(lps%sabs_improad_dif(beg:end,1:numrad))
    allocate(lps%sabs_perroad_dir(beg:end,1:numrad))
    allocate(lps%sabs_perroad_dif(beg:end,1:numrad))

    lps%t_building(beg:end) = nan
    lps%t_building_max(beg:end) = nan
    lps%t_building_min(beg:end) = nan
    lps%tk_wall(beg:end,1:nlevurb) = nan
    lps%tk_roof(beg:end,1:nlevurb) = nan
    lps%tk_improad(beg:end,1:nlevgrnd) = nan
    lps%cv_wall(beg:end,1:nlevurb) = nan
    lps%cv_roof(beg:end,1:nlevurb) = nan
    lps%cv_improad(beg:end,1:nlevgrnd) = nan
    lps%cv_improad(beg:end,1:5) = nan
    lps%thick_wall(beg:end) = nan
    lps%thick_roof(beg:end) = nan
    lps%nlev_improad(beg:end) = bigint
    lps%vf_sr(beg:end) = nan
    lps%vf_wr(beg:end) = nan
    lps%vf_sw(beg:end) = nan
    lps%vf_rw(beg:end) = nan
    lps%vf_ww(beg:end) = nan
    lps%taf(beg:end) = nan
    lps%qaf(beg:end) = nan
    lps%sabs_roof_dir(beg:end,1:numrad) = nan
    lps%sabs_roof_dif(beg:end,1:numrad) = nan
    lps%sabs_sunwall_dir(beg:end,1:numrad) = nan
    lps%sabs_sunwall_dif(beg:end,1:numrad) = nan
    lps%sabs_shadewall_dir(beg:end,1:numrad) = nan
    lps%sabs_shadewall_dif(beg:end,1:numrad) = nan
    lps%sabs_improad_dir(beg:end,1:numrad) = nan
    lps%sabs_improad_dif(beg:end,1:numrad) = nan
    lps%sabs_perroad_dir(beg:end,1:numrad) = nan
    lps%sabs_perroad_dif(beg:end,1:numrad) = nan

  end subroutine init_landunit_pstate_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_landunit_eflux_type
!
! !INTERFACE:
  subroutine init_landunit_eflux_type(beg, end, lef)
!
! !DESCRIPTION: 
! Initialize landunit energy flux variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end 
    type (landunit_eflux_type), intent(inout):: lef 
!
! !REVISION HISTORY:
! Created by Keith Oleson
!
!EOP
!------------------------------------------------------------------------

    allocate(lef%eflx_traffic(beg:end))
    allocate(lef%eflx_traffic_factor(beg:end))
    allocate(lef%eflx_wasteheat(beg:end))
    allocate(lef%eflx_heat_from_ac(beg:end))

    lef%eflx_traffic(beg:end) = nan
    lef%eflx_traffic_factor(beg:end) = nan
    lef%eflx_wasteheat(beg:end) = nan
    lef%eflx_heat_from_ac(beg:end) = nan

  end subroutine init_landunit_eflux_type

#if (defined CNDV)
!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_gridcell_dgvstate_type
!
! !INTERFACE:
  subroutine init_gridcell_dgvstate_type(beg, end, gps)
!
! !DESCRIPTION:
! Initialize gridcell DGVM variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (gridcell_dgvstate_type), intent(inout):: gps
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(gps%agdd20(beg:end))
    allocate(gps%tmomin20(beg:end))
    allocate(gps%t10min(beg:end))
    gps%agdd20(beg:end) = nan
    gps%tmomin20(beg:end) = nan
    gps%t10min(beg:end) = nan

  end subroutine init_gridcell_dgvstate_type
#endif

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_gridcell_pstate_type
!
! !INTERFACE:
  subroutine init_gridcell_pstate_type(beg, end, gps)
!
! !DESCRIPTION:
! Initialize gridcell physical state variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (gridcell_pstate_type), intent(inout):: gps
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------
    
    
    !allocate(gps%bcphiwet2t(beg:end,1:2))
    !allocate(gps%bcphidry2t(beg:end,1:2))
    !allocate(gps%bcphodry2t(beg:end,1:2))
    !allocate(gps%ocphiwet2t(beg:end,1:2))
    !allocate(gps%ocphidry2t(beg:end,1:2))
    !allocate(gps%ocphodry2t(beg:end,1:2))
    !allocate(gps%dstx01wd2t(beg:end,1:2))
    !allocate(gps%dstx01dd2t(beg:end,1:2))
    !allocate(gps%dstx02wd2t(beg:end,1:2))
    !allocate(gps%dstx02dd2t(beg:end,1:2))
    !allocate(gps%dstx03wd2t(beg:end,1:2))
    !allocate(gps%dstx03dd2t(beg:end,1:2))
    !allocate(gps%dstx04wd2t(beg:end,1:2))
    !allocate(gps%dstx04dd2t(beg:end,1:2))
    
    !gps%bcphiwet2t(beg:end,1:2) = nan
    !gps%bcphidry2t(beg:end,1:2) = nan
    !gps%bcphodry2t(beg:end,1:2) = nan
    !gps%ocphiwet2t(beg:end,1:2) = nan
    !gps%ocphidry2t(beg:end,1:2) = nan
    !gps%ocphodry2t(beg:end,1:2) = nan
    !gps%dstx01wd2t(beg:end,1:2) = nan
    !gps%dstx01dd2t(beg:end,1:2) = nan
    !gps%dstx02wd2t(beg:end,1:2) = nan
    !gps%dstx02dd2t(beg:end,1:2) = nan
    !gps%dstx03wd2t(beg:end,1:2) = nan
    !gps%dstx03dd2t(beg:end,1:2) = nan
    !gps%dstx04wd2t(beg:end,1:2) = nan
    !gps%dstx04dd2t(beg:end,1:2) = nan

  end subroutine init_gridcell_pstate_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_gridcell_efstate_type
!
! !INTERFACE:
  subroutine init_gridcell_efstate_type(beg, end, gve)
!
! !DESCRIPTION:
! Initialize gridcell isoprene emission factor variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (gridcell_efstate_type), intent(inout) :: gve
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein (heald)
!
!EOP
!------------------------------------------------------------------------

    allocate(gve%efisop(6,beg:end))
    gve%efisop(:,beg:end) = nan

  end subroutine init_gridcell_efstate_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_gridcell_wflux_type
!
! !INTERFACE:
  subroutine init_gridcell_wflux_type(beg, end, gwf)
!
! !DESCRIPTION:
! Initialize gridcell water flux variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (gridcell_wflux_type), intent(inout):: gwf
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(gwf%qflx_runoffg(beg:end))
    allocate(gwf%qflx_snwcp_iceg(beg:end))
    allocate(gwf%qflx_liq_dynbal(beg:end))
    allocate(gwf%qflx_ice_dynbal(beg:end))

    gwf%qflx_runoffg(beg:end) = nan
    gwf%qflx_snwcp_iceg(beg:end) = nan
    gwf%qflx_liq_dynbal(beg:end) = nan
    gwf%qflx_ice_dynbal(beg:end) = nan

  end subroutine init_gridcell_wflux_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_gridcell_eflux_type
!
! !INTERFACE:
  subroutine init_gridcell_eflux_type(beg, end, gef)
!
! !DESCRIPTION:
! Initialize gridcell energy flux variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (gridcell_eflux_type), intent(inout):: gef
!
! !REVISION HISTORY:
! Created by David Lawrence
!
!EOP
!------------------------------------------------------------------------
    allocate(gef%eflx_sh_totg(beg:end))
    allocate(gef%eflx_dynbal(beg:end))

    gef%eflx_sh_totg(beg:end) = nan
    gef%eflx_dynbal(beg:end) = nan

  end subroutine init_gridcell_eflux_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_gridcell_wstate_type
!
! !INTERFACE:
  subroutine init_gridcell_wstate_type(beg, end, gws)
!
! !DESCRIPTION:
! Initialize gridcell water state variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (gridcell_wstate_type), intent(inout):: gws
!
! !REVISION HISTORY:
! Created by David Lawrence
!
!EOP
!------------------------------------------------------------------------
    allocate(gws%gc_liq1(beg:end))
    allocate(gws%gc_liq2(beg:end))
    allocate(gws%gc_ice1(beg:end))     
    allocate(gws%gc_ice2(beg:end))    

    gws%gc_liq1(beg:end) = nan
    gws%gc_liq2(beg:end) = nan
    gws%gc_ice1(beg:end) = nan     
    gws%gc_ice2(beg:end) = nan    

  end subroutine init_gridcell_wstate_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_gridcell_estate_type
!
! !INTERFACE:
  subroutine init_gridcell_estate_type(beg, end, ges)
!
! !DESCRIPTION:
! Initialize gridcell energy state variables     
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (gridcell_estate_type), intent(inout):: ges
!
! !REVISION HISTORY:
! Created by David Lawrence
!
!EOP
!------------------------------------------------------------------------
    allocate(ges%gc_heat1(beg:end))     
    allocate(ges%gc_heat2(beg:end))    

    ges%gc_heat1(beg:end) = nan     
    ges%gc_heat2(beg:end) = nan    

  end subroutine init_gridcell_estate_type

end module clmtypeInitMod
