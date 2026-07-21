#include <misc.h>
#include <preproc.h>

! This function is modified by Chaopeng Shen to link in PAWS ==============
! CHAOPENG modifications are signaled by !CP. Any line without !CP has not been modified by CP
! Mainly, this is the driver program for the CN Ecosystem dynamics in CLM4.0 (or Biome-BGC)
! Mortality rate will be implemented ****; psn will be provided by RUE concept modeling; vegetation module will be run on a daily basis
! data pointers in the structures will be linked to a portion of a big memory storage
! functions to be modified include 
! Initialization subroutine (now turned a linking subroutine)
! driver program (to remove some functions)
! gap mortality update routine (modify CNCStateUpdate2,CNNStateUpdate2)
! removed fire
! Data/Classification/Parameters interfacing
! End CP Comments region ===================================================

! Questions:
! litterfall rate, turn over rate and respiration rate rescaling to time step???????
! also, snowcover need to be updated

module CNEcosystemDynMod
#ifdef CN

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: CNEcosystemDynMod
!
! !DESCRIPTION:
! Ecosystem dynamics: phenology, vegetation
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
  use clm_varctl  , only: fpftdyn !CP
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: CNEcosystemDynInit   ! Ecosystem dynamics initialization
  public :: CNEcosystemDyn   ! Ecosystem dynamics: phenology, vegetation
!
! !REVISION HISTORY:
! Created by Peter Thornton
! 19 May 2009: PET - modified to include call to harvest routine
!
!
! !PRIVATE MEMBER FUNCTIONS:
!
! !PRIVATE TYPES:
!EOP
!-----------------------------------------------------------------------

contains


!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CNEcosystemDynInit
!
! !INTERFACE:
  subroutine CNEcosystemDynInit(lbc, ubc, lbp, ubp )
!
! !DESCRIPTION:
! Initialzation of the CN Ecosystem dynamics.
!
! !USES:
    use CNAllocationMod, only : CNAllocationInit
    use CNPhenologyMod , only : CNPhenologyInit
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: lbc, ubc        ! column bounds
    integer, intent(in) :: lbp, ubp        ! pft bounds
!
! !CALLED FROM:
!
! !REVISION HISTORY:
! 04/05/11, Erik Kluzek creation
!
! !LOCAL VARIABLES:
!EOP
!-----------------------------------------------------------------------
     call CNAllocationInit ( lbc, ubc, lbp, ubp )
     call CNPhenologyInit  ( lbp, ubp )

  end subroutine CNEcosystemDynInit

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CNEcosystemDyn
!
! !INTERFACE:
  subroutine CNEcosystemDyn(lbc, ubc, lbp, ubp, num_soilc, filter_soilc, &
                     num_soilp, filter_soilp, num_pcropp, filter_pcropp, doalb)
!
! !DESCRIPTION:
! The core CN code is executed here. Calculates fluxes for maintenance
! respiration, decomposition, allocation, phenology, and growth respiration.
! These routines happen on the radiation time step so that canopy structure
! stays synchronized with albedo calculations.
!
! !USES:
    use clmtype
    !use spmdMod              , only: masterproc !CP
    use CNSetValueMod        , only: CNZeroFluxes
    use CNNDynamicsMod       , only: CNNDeposition,CNNFixation, CNNLeaching
    use CNMRespMod           , only: CNMResp
    use CNDecompMod          , only: CNDecompAlloc
    use CNPhenologyMod       , only: CNPhenology
    use CNGRespMod           , only: CNGResp
    use CNCStateUpdate1Mod   , only: CStateUpdate1,CStateUpdate0
    use CNNStateUpdate1Mod   , only: NStateUpdate1
    use CNGapMortalityMod    , only: CNGapMortality
    use CNCStateUpdate2Mod   , only: CStateUpdate2, CStateUpdate2h
    use CNNStateUpdate2Mod   , only: NStateUpdate2, NStateUpdate2h
    use CNFireMod            , only: CNFireArea, CNFireFluxes !CP
    use CNCStateUpdate3Mod   , only: CStateUpdate3
    use CNNStateUpdate3Mod   , only: NStateUpdate3
    use CNBalanceCheckMod    , only: CBalanceCheck, NBalanceCheck
    use CNPrecisionControlMod, only: CNPrecisionControl
    use CNVegStructUpdateMod , only: CNVegStructUpdate
    use CNAnnualUpdateMod    , only: CNAnnualUpdate
    use CNSummaryMod         , only: CSummary, NSummary !CP need to modify
#if (defined C13)
    use CNC13StateUpdate1Mod , only: C13StateUpdate1,C13StateUpdate0
    use CNC13StateUpdate2Mod , only: C13StateUpdate2, C13StateUpdate2h
    use CNC13StateUpdate3Mod , only: C13StateUpdate3
    use CNC13FluxMod         , only: C13Flux1, C13Flux2, C13Flux2h, C13Flux3
    use C13SummaryMod        , only: C13Summary
#endif
    !use pftdynMod               , only: CNHarvest !CP need to change
    use CNHarvestMod         , only: CNHarvest
    use CNWoodProductsMod    , only: CNWoodProducts
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: lbc, ubc        ! column bounds
    integer, intent(in) :: lbp, ubp        ! pft bounds
    integer, intent(in) :: num_soilc       ! number of soil columns in filter
    integer, intent(in) :: filter_soilc(ubc-lbc+1) ! filter for soil columns
    integer, intent(in) :: num_soilp       ! number of soil pfts in filter
    integer, intent(in) :: filter_soilp(ubp-lbp+1) ! filter for soil pfts
    integer, intent(in) :: num_pcropp      ! number of prog. crop pfts in filter
    integer, intent(in) :: filter_pcropp(:)! filter for prognostic crop pfts
    logical, intent(in) :: doalb           ! true = surface albedo calculation time step
    type(model_type),pointer::clm
    type(column_type),pointer::c
    type(pft_type),pointer::p
!
! !CALLED FROM:
!
! !REVISION HISTORY:
! 10/22/03, Peter Thornton: created from EcosystemDyn during migration to
!                           new vector code.
! 11/3/03, Peter Thornton: removed update of elai, esai, frac_veg_nosno_alb.
!     These are now done in CNVegStructUpdate(), which is called
!     prior to SurfaceAlbedo().
! 11/13/03, Peter Thornton: switched from nolake to soil filtering.
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in arguments
!
! local pointers to implicit out arguments
!
! !OTHER LOCAL VARIABLES:
!
!EOP
!-----------------------------------------------------------------------

 !   if (doalb) then
       clm => clm3
       c   => clm%g%l%c
       p   => c%p
       ! Call the main CN routines
       call CNZeroFluxes(num_soilc, filter_soilc, num_soilp, filter_soilp) !CP removed m_*** and hrv_*** fluxes

       call CNNDeposition(lbc, ubc)

       call CNNFixation(num_soilc,filter_soilc)

       call CNMResp(lbc, ubc, num_soilc, filter_soilc, num_soilp, filter_soilp)

       call CNDecompAlloc(lbp, ubp, lbc, ubc, num_soilc, filter_soilc, &
                          num_soilp, filter_soilp, num_pcropp)

       ! CNphenology needs to be called after CNdecompAlloc, because it
       ! depends on current time-step fluxes to new growth on the last
       ! litterfall timestep in deciduous systems

       call CNPhenology(num_soilc, filter_soilc, num_soilp, filter_soilp, &
                        num_pcropp, filter_pcropp, doalb)

       call CNGResp(num_soilp, filter_soilp)
       
       call CStateUpdate0(num_soilp, filter_soilp)

#if (defined C13)
       call C13StateUpdate0(num_soilp, filter_soilp)

       call C13Flux1(num_soilc, filter_soilc, num_soilp, filter_soilp)
#endif

       call CStateUpdate1(num_soilc, filter_soilc, num_soilp, filter_soilp)

#if (defined C13)
       call C13StateUpdate1(num_soilc, filter_soilc, num_soilp, filter_soilp)
#endif
       
       call NStateUpdate1(num_soilc, filter_soilc, num_soilp, filter_soilp)

       call CNGapMortality(num_soilc, filter_soilc, num_soilp, filter_soilp)

#if (defined C13)
       call C13Flux2(num_soilc, filter_soilc, num_soilp, filter_soilp)
#endif

       call CStateUpdate2(num_soilc, filter_soilc, num_soilp, filter_soilp)

#if (defined C13)
       call C13StateUpdate2(num_soilc, filter_soilc, num_soilp, filter_soilp)
#endif

       call NStateUpdate2(num_soilc, filter_soilc, num_soilp, filter_soilp)
       
       !if (fpftdyn /= ' ') then
          call CNHarvest(num_soilc, filter_soilc, num_soilp, filter_soilp)
       !end if 

#if (defined C13)
       call C13Flux2h(num_soilc, filter_soilc, num_soilp, filter_soilp)
#endif

       call CStateUpdate2h(num_soilc, filter_soilc, num_soilp, filter_soilp)

#if (defined C13)
       call C13StateUpdate2h(num_soilc, filter_soilc, num_soilp, filter_soilp)
#endif

       call NStateUpdate2h(num_soilc, filter_soilc, num_soilp, filter_soilp)
       
       call CNWoodProducts(num_soilc, filter_soilc)
       
       call CNFireArea(num_soilc, filter_soilc) 

       call CNFireFluxes(num_soilc, filter_soilc, num_soilp, filter_soilp)

       call CNNLeaching(lbc, ubc, num_soilc, filter_soilc)

#if (defined C13)
       call C13Flux3(num_soilc, filter_soilc, num_soilp, filter_soilp)
#endif

       call CStateUpdate3(num_soilc, filter_soilc, num_soilp, filter_soilp)

#if (defined C13)
       call C13StateUpdate3(num_soilc, filter_soilc, num_soilp, filter_soilp)
#endif

       call NStateUpdate3(num_soilc, filter_soilc, num_soilp, filter_soilp) !LEACHING IS INSIDE

       call CNPrecisionControl(num_soilc, filter_soilc, num_soilp, filter_soilp)

    if (doalb) then   
       call CNVegStructUpdate(num_soilp, filter_soilp)
    end if

!       call CNAnnualUpdate(num_soilc, filter_soilc, num_soilp, filter_soilp)
       
       call CSummary(num_soilc, filter_soilc, num_soilp, filter_soilp)
       
#if (defined C13)
       call C13Summary(num_soilc, filter_soilc, num_soilp, filter_soilp)
#endif
       
       call NSummary(num_soilc, filter_soilc, num_soilp, filter_soilp)

!    end if  !end of if-doalb block

  end subroutine CNEcosystemDyn
#endif
!-----------------------------------------------------------------------
end  module CNEcosystemDynMod
