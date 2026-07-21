!===============================================================================
! CVS $Id: shr_const_mod.F90,v 1.2.6.1 2004/01/07 20:50:01 mvertens Exp $
! CVS $Source: /fs/cgd/csm/models/CVS.REPOS/shared/csm_share/shr/shr_const_mod.F90,v $
! CVS $Name: clm3_deva_21 $
!===============================================================================

MODULE shr_const_mod

   use shr_kind_mod, only: SHR_KIND_R8, SHR_KIND_IN

   !----------------------------------------------------------------------------
   ! physical constants (all data public) !CP commented meaning and unit
   ! http://www.cesm.ucar.edu/models/ccsm2.0.1/shr_const.html
   !----------------------------------------------------------------------------
   integer(SHR_KIND_IN ), parameter :: 	R8 = SHR_KIND_R8
real(R8), parameter :: 	SHR_CONST_PI = 3.14159265358979323846_R8
real(R8), parameter :: 	SHR_CONST_CDAY = 86400.0_R8
real(R8), parameter :: 	SHR_CONST_SDAY = 86164.0_R8                       ! sidereal day,s
real(R8), parameter :: 	SHR_CONST_OMEGA = 2.0_R8*SHR_CONST_PI/SHR_CONST_SDAY ! angular frequency of earth rotation, radians/s
real(R8), parameter :: 	SHR_CONST_REARTH = 6.37122e6_R8                   ! m
real(R8), parameter :: 	SHR_CONST_G = 9.80616_R8                          ! m/s^2
real(R8), parameter :: 	SHR_CONST_STEBOL = 5.67e-8_R8                     ! stephen-boltzmann constant W/m^2/K^4
real(R8), parameter :: 	SHR_CONST_BOLTZ = 1.3806504e-23_R8                ! boltzman constant J/K
real(R8), parameter :: 	SHR_CONST_AVOGAD = 6.02214179e26_R8               ! 1/kmol, if 1/mol, it should be 6.02214179e23_R8
real(R8), parameter :: 	SHR_CONST_RGAS = SHR_CONST_AVOGAD*SHR_CONST_BOLTZ ! Gas constant J/K/kmol, divide by 1000 to get J/K/mol
real(R8), parameter :: 	SHR_CONST_MWDAIR = 28.966_R8                      ! molar weight of dry air, g/mol
real(R8), parameter :: 	SHR_CONST_MWWV = 18.016_R8                        ! molar weight of water vapor, g/mol
real(R8), parameter :: 	SHR_CONST_RDAIR = SHR_CONST_RGAS/SHR_CONST_MWDAIR ! R/, J/K/kg -- dry air
real(R8), parameter :: 	SHR_CONST_RWV = SHR_CONST_RGAS/SHR_CONST_MWWV     ! J/K/kg -- water vapor
real(R8), parameter :: 	SHR_CONST_ZVIR = (SHR_CONST_RWV/SHR_CONST_RDAIR)-1.0_R8
real(R8), parameter :: 	SHR_CONST_KARMAN = 0.4_R8                         ! von Karman constant
real(R8), parameter :: 	SHR_CONST_PSTD = 101325.0_R8                      ! STD atmospheric pressure, Pa
real(R8), parameter :: 	SHR_CONST_PDB = 0.0112372_R8                      ! Pee Dee Belemnite, 13C/12C
real(R8), parameter :: 	SHR_CONST_TKTRIP = 273.16_R8                      ! Kelvin
real(R8), parameter :: 	SHR_CONST_TKFRZ = 273.15_R8                       ! Kelvin
real(R8), parameter :: 	SHR_CONST_TKFRZSW = SHR_CONST_TKFRZ - 1.8_R8
real(R8), parameter :: 	SHR_CONST_RHODAIR = SHR_CONST_PSTD/(SHR_CONST_RDAIR*SHR_CONST_TKFRZ) ! dry air density, kg/m^3
real(R8), parameter :: 	SHR_CONST_RHOFW = 1.000e3_R8                      ! density of fresh water,  kg/m^3
real(R8), parameter :: 	SHR_CONST_RHOSW = 1.026e3_R8                      ! density of sea water,  kg/m^3
real(R8), parameter :: 	SHR_CONST_RHOICE = 0.917e3_R8                     ! ice density kg/m3
real(R8), parameter :: 	SHR_CONST_CPDAIR = 1.00464e3_R8                   ! specific heat of dry air, J/kg/K
real(R8), parameter :: 	SHR_CONST_CPWV = 1.810e3_R8                       ! specific heat of water vap ~ J/kg/K
real(R8), parameter :: 	SHR_CONST_CPVIR = (SHR_CONST_CPWV/SHR_CONST_CPDAIR)-1.0_R8
real(R8), parameter :: 	SHR_CONST_CPFW = 4.188e3_R8                       ! specific heat of fresh h2o, J/kg/K 
real(R8), parameter :: 	SHR_CONST_CPSW = 3.996e3_R8                       ! specific heat of sea h2o, J/kg/K 
real(R8), parameter :: 	SHR_CONST_CPICE = 2.11727e3_R8                    ! specific heat of fresh ice, J/kg/K 
real(R8), parameter :: 	SHR_CONST_LATICE = 3.337e5_R8                     ! latent heat of fusion, J/kg
real(R8), parameter :: 	SHR_CONST_LATVAP = 2.501e6_R8                     ! latent heat of evaporation  J/kg
real(R8), parameter :: 	SHR_CONST_LATSUB = SHR_CONST_LATICE + SHR_CONST_LATVAP ! latent heat of sublimation  J/kg
real(R8), parameter :: 	SHR_CONST_OCN_REF_SAL = 34.7_R8                   
real(R8), parameter :: 	SHR_CONST_ICE_REF_SAL = 4.0_R8
real(R8), parameter :: 	SHR_CONST_SPVAL = 1.0e30_R8

END MODULE shr_const_mod
