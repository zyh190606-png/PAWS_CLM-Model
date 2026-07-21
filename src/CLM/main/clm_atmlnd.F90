module clm_atmlnd

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: clm_atmlnd
!
! !DESCRIPTION:
! Handle atm2lnd, lnd2atm mapping/downscaling/upscaling/data
!
! !USES:
  use clm_varpar  , only : numrad, ndst   !ndst = number of dust bins.
                                          !only used # ifdef DUST
  use clm_varcon  , only : rair, grav, cpair, hfus, tfrz
  use clm_varctl  , only : iulog
  use decompMod   , only : get_proc_bounds, get_proc_bounds_atm
  use shr_kind_mod, only : r8 => shr_kind_r8
  use nanMod      , only : nan
  use spmdMod     , only : masterproc
  use abortutils  , only : endrun
  !use seq_drydep_mod, only : n_drydep, drydep_method, DD_XLND !CP
  use clm_varpar  , only : nvoc
!
! !PUBLIC TYPES:
  implicit none
PRIVATE
!----------------------------------------------------
! atmosphere -> land variables structure
!----------------------------------------------------
  type,PUBLIC:: atm2lnd_type
     real(r8), pointer :: forc_t(:)       !atmospheric temperature (Kelvin)
     real(r8), pointer :: forc_u(:)       !atm wind speed, east direction (m/s)
     real(r8), pointer :: forc_v(:)       !atm wind speed, north direction (m/s)
     real(r8), pointer :: forc_wind(:)    !atmospheric wind speed   
     real(r8), pointer :: forc_q(:)       !atmospheric specific humidity (kg/kg)
     real(r8), pointer :: forc_hgt(:)     !atmospheric reference height (m)
     real(r8), pointer :: forc_hgt_u(:)   !obs height of wind [m] (new)
     real(r8), pointer :: forc_hgt_t(:)   !obs height of temperature [m] (new)
     real(r8), pointer :: forc_hgt_q(:)   !obs height of humidity [m] (new)
     real(r8), pointer :: forc_pbot(:)    !atmospheric pressure (Pa)
     real(r8), pointer :: forc_th(:)      !atm potential temperature (Kelvin)
     real(r8), pointer :: forc_vp(:)      !atmospheric vapor pressure (Pa) 
     real(r8), pointer :: forc_rho(:)     !density (kg/m**3)
     real(r8), pointer :: forc_rh(:)      !atmospheric relative humidity (%)
     real(r8), pointer :: forc_psrf(:)    !surface pressure (Pa)
     real(r8), pointer :: forc_pco2(:)    !CO2 partial pressure (Pa)
     real(r8), pointer :: forc_lwrad(:)   !downwrd IR longwave radiation (W/m**2)
     real(r8), pointer :: forc_solad(:,:) !direct beam radiation (numrad) (vis=forc_sols , nir=forc_soll )
     real(r8), pointer :: forc_solai(:,:) !diffuse radiation (numrad) (vis=forc_solsd, nir=forc_solld)
     real(r8), pointer :: forc_solar(:)   !incident solar radiation
     real(r8), pointer :: forc_rain(:)    !rain rate [mm/s]
     real(r8), pointer :: forc_snow(:)    !snow rate [mm/s]
     real(r8), pointer :: forc_ndep(:)    !nitrogen deposition rate (gN/m2/s)
     real(r8), pointer :: rainf(:)        !ALMA rain+snow [mm/s]
#ifdef C13
     real(r8), pointer :: forc_pc13o2(:)  !C13O2 partial pressure (Pa)
#endif
     real(r8), pointer :: forc_po2(:)     !O2 partial pressure (Pa)
     real(r8), pointer :: forc_aer(:,:)   ! aerosol deposition array
  end type atm2lnd_type

!----------------------------------------------------
! land -> atmosphere variables structure
!----------------------------------------------------
  type, PUBLIC:: lnd2atm_type
     real(r8), pointer :: t_rad(:)        !radiative temperature (Kelvin)
     real(r8), pointer :: t_ref2m(:)      !2m surface air temperature (Kelvin)
     real(r8), pointer :: q_ref2m(:)      !2m surface specific humidity (kg/kg)
     real(r8), pointer :: h2osno(:)       !snow water (mm H2O)
     real(r8), pointer :: albd(:,:)       !(numrad) surface albedo (direct)
     real(r8), pointer :: albi(:,:)       !(numrad) surface albedo (diffuse)
     real(r8), pointer :: taux(:)         !wind stress: e-w (kg/m/s**2)
     real(r8), pointer :: tauy(:)         !wind stress: n-s (kg/m/s**2)
     real(r8), pointer :: eflx_lh_tot(:)  !total latent HF (W/m**2)  [+ to atm]
     real(r8), pointer :: eflx_sh_tot(:)  !total sensible HF (W/m**2) [+ to atm]
     real(r8), pointer :: eflx_lwrad_out(:) !IR (longwave) radiation (W/m**2)
     real(r8), pointer :: qflx_evap_tot(:)!qflx_evap(_soi + _veg) + qflx_tran_veg
     real(r8), pointer :: fsa(:)          !solar rad absorbed (total) (W/m**2)
     real(r8), pointer :: nee(:)          !net CO2 flux (kg CO2/m**2/s) [+ to atm]
     real(r8), pointer :: ram1(:)         !aerodynamical resistance (s/m)
     real(r8), pointer :: fv(:)           !friction velocity (m/s) (for dust model)
     real(r8), pointer :: flxdst(:,:)       !dust flux (size bins)
     real(r8), pointer :: ddvel(:,:)        !dry deposition velocities
     real(r8), pointer :: flxvoc(:,:)       ! VOC flux (size bins)
  end type lnd2atm_type
  
  type(atm2lnd_type),public,target :: atm_a2l      ! a2l fields on atm grid
  type(lnd2atm_type),public,target :: atm_l2a      ! l2a fields on atm grid

  type(atm2lnd_type),public,target :: clm_a2l      ! a2l fields on clm grid
  type(lnd2atm_type),public,target :: clm_l2a      ! l2a fields on clm grid

  real(r8), pointer, public :: adiag_arain(:)
  real(r8), pointer, public :: adiag_asnow(:)
  real(r8), pointer, public :: adiag_aflux(:)
  real(r8), pointer, public :: adiag_lflux(:)

! !PUBLIC MEMBER FUNCTIONS:
  !public :: init_adiag_type
  !public :: init_atm2lnd_type !CP: moved to clmtypeInitMod, so that when modifying vdata, CLM need not be compiled again
  public :: init_lnd2atm_type
  !public :: clm_mapa2l !CP cut the linkages here
  !public :: clm_map2gcell
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein and Tony Craig, 2006-01-10
!
! !PRIVATE MEMBER FUNCTIONS:

!EOP
!----------------------------------------------------

contains

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_adiag_type
!
! !INTERFACE:

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_atm2lnd_type
!
! !INTERFACE:

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_lnd2atm_type
!
! !INTERFACE:
  subroutine init_lnd2atm_type(beg, end, l2a)
!
! !DESCRIPTION:
! Initialize land variables required by the atmosphere
!
! !ARGUMENTS:
  implicit none
  integer, intent(in) :: beg, end
  type (lnd2atm_type), intent(inout):: l2a

end subroutine init_lnd2atm_type
!CP cut the endless linkages here
end module clm_atmlnd

