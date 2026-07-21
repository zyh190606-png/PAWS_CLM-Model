module CNHarvestMod

  use spmdMod
  use clmtype
  use decompMod   , only : get_proc_bounds
  use clm_varsur  , only : pctspec
  use clm_varpar  , only : max_pft_per_col
  use clm_varctl  , only : iulog
  use shr_sys_mod , only : shr_sys_flush
  use shr_kind_mod, only : r8 => shr_kind_r8
  use abortutils  , only : endrun
  
  public :: CNHarvest
  public :: CNHarvestPftToColumn
  
Contains
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CNHarvest
!
! !INTERFACE:
subroutine CNHarvest (num_soilc, filter_soilc, num_soilp, filter_soilp)
!
! !DESCRIPTION:
! Harvest mortality routine for coupled carbon-nitrogen code (CN)
!
! !USES:
   use clmtype
   use pftvarcon       , only : noveg, nbrdlf_evr_shrub, pprodharv10
   use clm_varcon      , only : secspday
   use clm_time_manager, only : get_days_per_year
!
! !ARGUMENTS:
   implicit none
   integer, intent(in) :: num_soilc       ! number of soil columns in filter
   integer, intent(in) :: filter_soilc(:) ! column filter for soil points
   integer, intent(in) :: num_soilp       ! number of soil pfts in filter
   integer, intent(in) :: filter_soilp(:) ! pft filter for soil points
!
! !CALLED FROM:
! subroutine CNEcosystemDyn
!
! !REVISION HISTORY:
! 3/29/04: Created by Peter Thornton
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in arrays
   integer , pointer :: pgridcell(:)   ! pft-level index into gridcell-level quantities
   integer , pointer :: ivt(:)         ! pft vegetation type

   real(r8), pointer :: leafc(:)              ! (gC/m2) leaf C
   real(r8), pointer :: frootc(:)             ! (gC/m2) fine root C
   real(r8), pointer :: livestemc(:)          ! (gC/m2) live stem C
   real(r8), pointer :: deadstemc(:)          ! (gC/m2) dead stem C
   real(r8), pointer :: livecrootc(:)         ! (gC/m2) live coarse root C
   real(r8), pointer :: deadcrootc(:)         ! (gC/m2) dead coarse root C
   real(r8), pointer :: xsmrpool(:)           ! (gC/m2) abstract C pool to meet excess MR demand
   real(r8), pointer :: leafc_storage(:)      ! (gC/m2) leaf C storage
   real(r8), pointer :: frootc_storage(:)     ! (gC/m2) fine root C storage
   real(r8), pointer :: livestemc_storage(:)  ! (gC/m2) live stem C storage
   real(r8), pointer :: deadstemc_storage(:)  ! (gC/m2) dead stem C storage
   real(r8), pointer :: livecrootc_storage(:) ! (gC/m2) live coarse root C storage
   real(r8), pointer :: deadcrootc_storage(:) ! (gC/m2) dead coarse root C storage
   real(r8), pointer :: gresp_storage(:)      ! (gC/m2) growth respiration storage
   real(r8), pointer :: leafc_xfer(:)         ! (gC/m2) leaf C transfer
   real(r8), pointer :: frootc_xfer(:)        ! (gC/m2) fine root C transfer
   real(r8), pointer :: livestemc_xfer(:)     ! (gC/m2) live stem C transfer
   real(r8), pointer :: deadstemc_xfer(:)     ! (gC/m2) dead stem C transfer
   real(r8), pointer :: livecrootc_xfer(:)    ! (gC/m2) live coarse root C transfer
   real(r8), pointer :: deadcrootc_xfer(:)    ! (gC/m2) dead coarse root C transfer
   real(r8), pointer :: gresp_xfer(:)         ! (gC/m2) growth respiration transfer
   real(r8), pointer :: leafn(:)              ! (gN/m2) leaf N
   real(r8), pointer :: frootn(:)             ! (gN/m2) fine root N
   real(r8), pointer :: livestemn(:)          ! (gN/m2) live stem N
   real(r8), pointer :: deadstemn(:)          ! (gN/m2) dead stem N
   real(r8), pointer :: livecrootn(:)         ! (gN/m2) live coarse root N
   real(r8), pointer :: deadcrootn(:)         ! (gN/m2) dead coarse root N
   real(r8), pointer :: retransn(:)           ! (gN/m2) plant pool of retranslocated N
   real(r8), pointer :: leafn_storage(:)      ! (gN/m2) leaf N storage
   real(r8), pointer :: frootn_storage(:)     ! (gN/m2) fine root N storage
   real(r8), pointer :: livestemn_storage(:)  ! (gN/m2) live stem N storage
   real(r8), pointer :: deadstemn_storage(:)  ! (gN/m2) dead stem N storage
   real(r8), pointer :: livecrootn_storage(:) ! (gN/m2) live coarse root N storage
   real(r8), pointer :: deadcrootn_storage(:) ! (gN/m2) dead coarse root N storage
   real(r8), pointer :: leafn_xfer(:)         ! (gN/m2) leaf N transfer
   real(r8), pointer :: frootn_xfer(:)        ! (gN/m2) fine root N transfer
   real(r8), pointer :: livestemn_xfer(:)     ! (gN/m2) live stem N transfer
   real(r8), pointer :: deadstemn_xfer(:)     ! (gN/m2) dead stem N transfer
   real(r8), pointer :: livecrootn_xfer(:)    ! (gN/m2) live coarse root N transfer
   real(r8), pointer :: deadcrootn_xfer(:)    ! (gN/m2) dead coarse root N transfer
!
! local pointers to implicit in/out arrays
!
! local pointers to implicit out arrays
   real(r8), pointer :: hrv_leafc_to_litter(:)
   real(r8), pointer :: hrv_frootc_to_litter(:)
   real(r8), pointer :: hrv_livestemc_to_litter(:)
   real(r8), pointer :: hrv_deadstemc_to_prod10c(:)
   real(r8), pointer :: hrv_deadstemc_to_prod100c(:)
   real(r8), pointer :: hrv_livecrootc_to_litter(:)
   real(r8), pointer :: hrv_deadcrootc_to_litter(:)
   real(r8), pointer :: hrv_xsmrpool_to_atm(:)
   real(r8), pointer :: hrv_leafc_storage_to_litter(:)
   real(r8), pointer :: hrv_frootc_storage_to_litter(:)
   real(r8), pointer :: hrv_livestemc_storage_to_litter(:)
   real(r8), pointer :: hrv_deadstemc_storage_to_litter(:)
   real(r8), pointer :: hrv_livecrootc_storage_to_litter(:)
   real(r8), pointer :: hrv_deadcrootc_storage_to_litter(:)
   real(r8), pointer :: hrv_gresp_storage_to_litter(:)
   real(r8), pointer :: hrv_leafc_xfer_to_litter(:)
   real(r8), pointer :: hrv_frootc_xfer_to_litter(:)
   real(r8), pointer :: hrv_livestemc_xfer_to_litter(:)
   real(r8), pointer :: hrv_deadstemc_xfer_to_litter(:)
   real(r8), pointer :: hrv_livecrootc_xfer_to_litter(:)
   real(r8), pointer :: hrv_deadcrootc_xfer_to_litter(:)
   real(r8), pointer :: hrv_gresp_xfer_to_litter(:)
   real(r8), pointer :: hrv_leafn_to_litter(:)
   real(r8), pointer :: hrv_frootn_to_litter(:)
   real(r8), pointer :: hrv_livestemn_to_litter(:)
   real(r8), pointer :: hrv_deadstemn_to_prod10n(:)
   real(r8), pointer :: hrv_deadstemn_to_prod100n(:)
   real(r8), pointer :: hrv_livecrootn_to_litter(:)
   real(r8), pointer :: hrv_deadcrootn_to_litter(:)
   real(r8), pointer :: hrv_retransn_to_litter(:)
   real(r8), pointer :: hrv_leafn_storage_to_litter(:)
   real(r8), pointer :: hrv_frootn_storage_to_litter(:)
   real(r8), pointer :: hrv_livestemn_storage_to_litter(:)
   real(r8), pointer :: hrv_deadstemn_storage_to_litter(:)
   real(r8), pointer :: hrv_livecrootn_storage_to_litter(:)
   real(r8), pointer :: hrv_deadcrootn_storage_to_litter(:)
   real(r8), pointer :: hrv_leafn_xfer_to_litter(:)
   real(r8), pointer :: hrv_frootn_xfer_to_litter(:)
   real(r8), pointer :: hrv_livestemn_xfer_to_litter(:)
   real(r8), pointer :: hrv_deadstemn_xfer_to_litter(:)
   real(r8), pointer :: hrv_livecrootn_xfer_to_litter(:)
   real(r8), pointer :: hrv_deadcrootn_xfer_to_litter(:)
!
! !OTHER LOCAL VARIABLES:
   integer :: p                         ! pft index
   integer :: g                         ! gridcell index
   integer :: fp                        ! pft filter index
   real(r8):: am                        ! rate for fractional harvest mortality (1/yr)
   real(r8):: m                         ! rate for fractional harvest mortality (1/s)
   real(r8):: days_per_year             ! days per year
   logical ::do_harvest = .true.
!EOP
!-----------------------------------------------------------------------

   ! assign local pointers to pft-level arrays
   pgridcell                      => clm3%g%l%c%p%gridcell
   
   ivt                            => clm3%g%l%c%p%itype
   leafc                          => clm3%g%l%c%p%pcs%leafc
   frootc                         => clm3%g%l%c%p%pcs%frootc
   livestemc                      => clm3%g%l%c%p%pcs%livestemc
   deadstemc                      => clm3%g%l%c%p%pcs%deadstemc
   livecrootc                     => clm3%g%l%c%p%pcs%livecrootc
   deadcrootc                     => clm3%g%l%c%p%pcs%deadcrootc
   xsmrpool                       => clm3%g%l%c%p%pcs%xsmrpool
   leafc_storage                  => clm3%g%l%c%p%pcs%leafc_storage
   frootc_storage                 => clm3%g%l%c%p%pcs%frootc_storage
   livestemc_storage              => clm3%g%l%c%p%pcs%livestemc_storage
   deadstemc_storage              => clm3%g%l%c%p%pcs%deadstemc_storage
   livecrootc_storage             => clm3%g%l%c%p%pcs%livecrootc_storage
   deadcrootc_storage             => clm3%g%l%c%p%pcs%deadcrootc_storage
   gresp_storage                  => clm3%g%l%c%p%pcs%gresp_storage
   leafc_xfer                     => clm3%g%l%c%p%pcs%leafc_xfer
   frootc_xfer                    => clm3%g%l%c%p%pcs%frootc_xfer
   livestemc_xfer                 => clm3%g%l%c%p%pcs%livestemc_xfer
   deadstemc_xfer                 => clm3%g%l%c%p%pcs%deadstemc_xfer
   livecrootc_xfer                => clm3%g%l%c%p%pcs%livecrootc_xfer
   deadcrootc_xfer                => clm3%g%l%c%p%pcs%deadcrootc_xfer
   gresp_xfer                     => clm3%g%l%c%p%pcs%gresp_xfer
   leafn                          => clm3%g%l%c%p%pns%leafn
   frootn                         => clm3%g%l%c%p%pns%frootn
   livestemn                      => clm3%g%l%c%p%pns%livestemn
   deadstemn                      => clm3%g%l%c%p%pns%deadstemn
   livecrootn                     => clm3%g%l%c%p%pns%livecrootn
   deadcrootn                     => clm3%g%l%c%p%pns%deadcrootn
   retransn                       => clm3%g%l%c%p%pns%retransn
   leafn_storage                  => clm3%g%l%c%p%pns%leafn_storage
   frootn_storage                 => clm3%g%l%c%p%pns%frootn_storage
   livestemn_storage              => clm3%g%l%c%p%pns%livestemn_storage
   deadstemn_storage              => clm3%g%l%c%p%pns%deadstemn_storage
   livecrootn_storage             => clm3%g%l%c%p%pns%livecrootn_storage
   deadcrootn_storage             => clm3%g%l%c%p%pns%deadcrootn_storage
   leafn_xfer                     => clm3%g%l%c%p%pns%leafn_xfer
   frootn_xfer                    => clm3%g%l%c%p%pns%frootn_xfer
   livestemn_xfer                 => clm3%g%l%c%p%pns%livestemn_xfer
   deadstemn_xfer                 => clm3%g%l%c%p%pns%deadstemn_xfer
   livecrootn_xfer                => clm3%g%l%c%p%pns%livecrootn_xfer
   deadcrootn_xfer                => clm3%g%l%c%p%pns%deadcrootn_xfer
   hrv_leafc_to_litter              => clm3%g%l%c%p%pcf%hrv_leafc_to_litter
   hrv_frootc_to_litter             => clm3%g%l%c%p%pcf%hrv_frootc_to_litter
   hrv_livestemc_to_litter          => clm3%g%l%c%p%pcf%hrv_livestemc_to_litter
   hrv_deadstemc_to_prod10c         => clm3%g%l%c%p%pcf%hrv_deadstemc_to_prod10c
   hrv_deadstemc_to_prod100c        => clm3%g%l%c%p%pcf%hrv_deadstemc_to_prod100c
   hrv_livecrootc_to_litter         => clm3%g%l%c%p%pcf%hrv_livecrootc_to_litter
   hrv_deadcrootc_to_litter         => clm3%g%l%c%p%pcf%hrv_deadcrootc_to_litter
   hrv_xsmrpool_to_atm              => clm3%g%l%c%p%pcf%hrv_xsmrpool_to_atm
   hrv_leafc_storage_to_litter      => clm3%g%l%c%p%pcf%hrv_leafc_storage_to_litter
   hrv_frootc_storage_to_litter     => clm3%g%l%c%p%pcf%hrv_frootc_storage_to_litter
   hrv_livestemc_storage_to_litter  => clm3%g%l%c%p%pcf%hrv_livestemc_storage_to_litter
   hrv_deadstemc_storage_to_litter  => clm3%g%l%c%p%pcf%hrv_deadstemc_storage_to_litter
   hrv_livecrootc_storage_to_litter => clm3%g%l%c%p%pcf%hrv_livecrootc_storage_to_litter
   hrv_deadcrootc_storage_to_litter => clm3%g%l%c%p%pcf%hrv_deadcrootc_storage_to_litter
   hrv_gresp_storage_to_litter      => clm3%g%l%c%p%pcf%hrv_gresp_storage_to_litter
   hrv_leafc_xfer_to_litter         => clm3%g%l%c%p%pcf%hrv_leafc_xfer_to_litter
   hrv_frootc_xfer_to_litter        => clm3%g%l%c%p%pcf%hrv_frootc_xfer_to_litter
   hrv_livestemc_xfer_to_litter     => clm3%g%l%c%p%pcf%hrv_livestemc_xfer_to_litter
   hrv_deadstemc_xfer_to_litter     => clm3%g%l%c%p%pcf%hrv_deadstemc_xfer_to_litter
   hrv_livecrootc_xfer_to_litter    => clm3%g%l%c%p%pcf%hrv_livecrootc_xfer_to_litter
   hrv_deadcrootc_xfer_to_litter    => clm3%g%l%c%p%pcf%hrv_deadcrootc_xfer_to_litter
   hrv_gresp_xfer_to_litter         => clm3%g%l%c%p%pcf%hrv_gresp_xfer_to_litter
   hrv_leafn_to_litter              => clm3%g%l%c%p%pnf%hrv_leafn_to_litter
   hrv_frootn_to_litter             => clm3%g%l%c%p%pnf%hrv_frootn_to_litter
   hrv_livestemn_to_litter          => clm3%g%l%c%p%pnf%hrv_livestemn_to_litter
   hrv_deadstemn_to_prod10n         => clm3%g%l%c%p%pnf%hrv_deadstemn_to_prod10n
   hrv_deadstemn_to_prod100n        => clm3%g%l%c%p%pnf%hrv_deadstemn_to_prod100n
   hrv_livecrootn_to_litter         => clm3%g%l%c%p%pnf%hrv_livecrootn_to_litter
   hrv_deadcrootn_to_litter         => clm3%g%l%c%p%pnf%hrv_deadcrootn_to_litter
   hrv_retransn_to_litter           => clm3%g%l%c%p%pnf%hrv_retransn_to_litter
   hrv_leafn_storage_to_litter      => clm3%g%l%c%p%pnf%hrv_leafn_storage_to_litter
   hrv_frootn_storage_to_litter     => clm3%g%l%c%p%pnf%hrv_frootn_storage_to_litter
   hrv_livestemn_storage_to_litter  => clm3%g%l%c%p%pnf%hrv_livestemn_storage_to_litter
   hrv_deadstemn_storage_to_litter  => clm3%g%l%c%p%pnf%hrv_deadstemn_storage_to_litter
   hrv_livecrootn_storage_to_litter => clm3%g%l%c%p%pnf%hrv_livecrootn_storage_to_litter
   hrv_deadcrootn_storage_to_litter => clm3%g%l%c%p%pnf%hrv_deadcrootn_storage_to_litter
   hrv_leafn_xfer_to_litter         => clm3%g%l%c%p%pnf%hrv_leafn_xfer_to_litter
   hrv_frootn_xfer_to_litter        => clm3%g%l%c%p%pnf%hrv_frootn_xfer_to_litter
   hrv_livestemn_xfer_to_litter     => clm3%g%l%c%p%pnf%hrv_livestemn_xfer_to_litter
   hrv_deadstemn_xfer_to_litter     => clm3%g%l%c%p%pnf%hrv_deadstemn_xfer_to_litter
   hrv_livecrootn_xfer_to_litter    => clm3%g%l%c%p%pnf%hrv_livecrootn_xfer_to_litter
   hrv_deadcrootn_xfer_to_litter    => clm3%g%l%c%p%pnf%hrv_deadcrootn_xfer_to_litter


   days_per_year = get_days_per_year()

   ! pft loop
   do fp = 1,num_soilp
      p = filter_soilp(fp)
      g = pgridcell(p)
      
      ! If this is a tree pft, then
      ! get the annual harvest "mortality" rate (am) from harvest array
      ! and convert to rate per second
      if (ivt(p) > noveg .and. ivt(p) < nbrdlf_evr_shrub) then

         if (do_harvest) then
            !am = harvest(g)  !CP
            am = 0.01D0 !CP
            m  = am/(days_per_year * secspday)
         else
            m = 0._r8
         end if   

         ! pft-level harvest carbon fluxes
         ! displayed pools
         hrv_leafc_to_litter(p)               = leafc(p)               * m
         hrv_frootc_to_litter(p)              = frootc(p)              * m
         hrv_livestemc_to_litter(p)           = livestemc(p)           * m
         hrv_deadstemc_to_prod10c(p)          = deadstemc(p)           * m * &
                                                pprodharv10(ivt(p))
         hrv_deadstemc_to_prod100c(p)         = deadstemc(p)           * m * &
                                                (1.0_r8 - pprodharv10(ivt(p)))
         hrv_livecrootc_to_litter(p)          = livecrootc(p)          * m
         hrv_deadcrootc_to_litter(p)          = deadcrootc(p)          * m
         hrv_xsmrpool_to_atm(p)               = xsmrpool(p)            * m

         ! storage pools
         hrv_leafc_storage_to_litter(p)       = leafc_storage(p)       * m
         hrv_frootc_storage_to_litter(p)      = frootc_storage(p)      * m
         hrv_livestemc_storage_to_litter(p)   = livestemc_storage(p)   * m
         hrv_deadstemc_storage_to_litter(p)   = deadstemc_storage(p)   * m
         hrv_livecrootc_storage_to_litter(p)  = livecrootc_storage(p)  * m
         hrv_deadcrootc_storage_to_litter(p)  = deadcrootc_storage(p)  * m
         hrv_gresp_storage_to_litter(p)       = gresp_storage(p)       * m

         ! transfer pools
         hrv_leafc_xfer_to_litter(p)          = leafc_xfer(p)          * m
         hrv_frootc_xfer_to_litter(p)         = frootc_xfer(p)         * m
         hrv_livestemc_xfer_to_litter(p)      = livestemc_xfer(p)      * m
         hrv_deadstemc_xfer_to_litter(p)      = deadstemc_xfer(p)      * m
         hrv_livecrootc_xfer_to_litter(p)     = livecrootc_xfer(p)     * m
         hrv_deadcrootc_xfer_to_litter(p)     = deadcrootc_xfer(p)     * m
         hrv_gresp_xfer_to_litter(p)          = gresp_xfer(p)          * m

         ! pft-level harvest mortality nitrogen fluxes
         ! displayed pools
         hrv_leafn_to_litter(p)               = leafn(p)               * m
         hrv_frootn_to_litter(p)              = frootn(p)              * m
         hrv_livestemn_to_litter(p)           = livestemn(p)           * m
         hrv_deadstemn_to_prod10n(p)          = deadstemn(p)           * m * &
                                                pprodharv10(ivt(p))
         hrv_deadstemn_to_prod100n(p)         = deadstemn(p)           * m * &
                                                (1.0_r8 - pprodharv10(ivt(p)))
         hrv_livecrootn_to_litter(p)          = livecrootn(p)          * m
         hrv_deadcrootn_to_litter(p)          = deadcrootn(p)          * m
         hrv_retransn_to_litter(p)            = retransn(p)            * m

         ! storage pools
         hrv_leafn_storage_to_litter(p)       = leafn_storage(p)       * m
         hrv_frootn_storage_to_litter(p)      = frootn_storage(p)      * m
         hrv_livestemn_storage_to_litter(p)   = livestemn_storage(p)   * m
         hrv_deadstemn_storage_to_litter(p)   = deadstemn_storage(p)   * m
         hrv_livecrootn_storage_to_litter(p)  = livecrootn_storage(p)  * m
         hrv_deadcrootn_storage_to_litter(p)  = deadcrootn_storage(p)  * m

         ! transfer pools
         hrv_leafn_xfer_to_litter(p)          = leafn_xfer(p)          * m
         hrv_frootn_xfer_to_litter(p)         = frootn_xfer(p)         * m
         hrv_livestemn_xfer_to_litter(p)      = livestemn_xfer(p)      * m
         hrv_deadstemn_xfer_to_litter(p)      = deadstemn_xfer(p)      * m
         hrv_livecrootn_xfer_to_litter(p)     = livecrootn_xfer(p)     * m
         hrv_deadcrootn_xfer_to_litter(p)     = deadcrootn_xfer(p)     * m
         
      end if  ! end tree block

   end do ! end of pft loop

   ! gather all pft-level litterfall fluxes from harvest to the column
   ! for litter C and N inputs

   call CNHarvestPftToColumn(num_soilc, filter_soilc)

end subroutine CNHarvest
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CNHarvestPftToColumn
!
! !INTERFACE:
subroutine CNHarvestPftToColumn (num_soilc, filter_soilc)
!
! !DESCRIPTION:
! called at the end of CNHarvest to gather all pft-level harvest litterfall fluxes
! to the column level and assign them to the three litter pools
!
! !USES:
  use clmtype
  use clm_varpar, only : max_pft_per_col, maxpatch_pft
!
! !ARGUMENTS:
  implicit none
  integer, intent(in) :: num_soilc       ! number of soil columns in filter
  integer, intent(in) :: filter_soilc(:) ! soil column filter
!
! !CALLED FROM:
! subroutine CNphenology
!
! !REVISION HISTORY:
! 9/8/03: Created by Peter Thornton
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in scalars
   integer , pointer :: ivt(:)      ! pft vegetation type
   real(r8), pointer :: wtcol(:)    ! pft weight relative to column (0-1)
   real(r8), pointer :: pwtgcell(:) ! weight of pft relative to corresponding gridcell
   real(r8), pointer :: lf_flab(:)  ! leaf litter labile fraction
   real(r8), pointer :: lf_fcel(:)  ! leaf litter cellulose fraction
   real(r8), pointer :: lf_flig(:)  ! leaf litter lignin fraction
   real(r8), pointer :: fr_flab(:)  ! fine root litter labile fraction
   real(r8), pointer :: fr_fcel(:)  ! fine root litter cellulose fraction
   real(r8), pointer :: fr_flig(:)  ! fine root litter lignin fraction
   integer , pointer :: npfts(:)    ! number of pfts for each column
   integer , pointer :: pfti(:)     ! beginning pft index for each column
   real(r8), pointer :: hrv_leafc_to_litter(:)
   real(r8), pointer :: hrv_frootc_to_litter(:)
   real(r8), pointer :: hrv_livestemc_to_litter(:)
   real(r8), pointer :: phrv_deadstemc_to_prod10c(:)
   real(r8), pointer :: phrv_deadstemc_to_prod100c(:)
   real(r8), pointer :: hrv_livecrootc_to_litter(:)
   real(r8), pointer :: hrv_deadcrootc_to_litter(:)
   real(r8), pointer :: hrv_leafc_storage_to_litter(:)
   real(r8), pointer :: hrv_frootc_storage_to_litter(:)
   real(r8), pointer :: hrv_livestemc_storage_to_litter(:)
   real(r8), pointer :: hrv_deadstemc_storage_to_litter(:)
   real(r8), pointer :: hrv_livecrootc_storage_to_litter(:)
   real(r8), pointer :: hrv_deadcrootc_storage_to_litter(:)
   real(r8), pointer :: hrv_gresp_storage_to_litter(:)
   real(r8), pointer :: hrv_leafc_xfer_to_litter(:)
   real(r8), pointer :: hrv_frootc_xfer_to_litter(:)
   real(r8), pointer :: hrv_livestemc_xfer_to_litter(:)
   real(r8), pointer :: hrv_deadstemc_xfer_to_litter(:)
   real(r8), pointer :: hrv_livecrootc_xfer_to_litter(:)
   real(r8), pointer :: hrv_deadcrootc_xfer_to_litter(:)
   real(r8), pointer :: hrv_gresp_xfer_to_litter(:)
   real(r8), pointer :: hrv_leafn_to_litter(:)
   real(r8), pointer :: hrv_frootn_to_litter(:)
   real(r8), pointer :: hrv_livestemn_to_litter(:)
   real(r8), pointer :: phrv_deadstemn_to_prod10n(:)
   real(r8), pointer :: phrv_deadstemn_to_prod100n(:)
   real(r8), pointer :: hrv_livecrootn_to_litter(:)
   real(r8), pointer :: hrv_deadcrootn_to_litter(:)
   real(r8), pointer :: hrv_retransn_to_litter(:)
   real(r8), pointer :: hrv_leafn_storage_to_litter(:)
   real(r8), pointer :: hrv_frootn_storage_to_litter(:)
   real(r8), pointer :: hrv_livestemn_storage_to_litter(:)
   real(r8), pointer :: hrv_deadstemn_storage_to_litter(:)
   real(r8), pointer :: hrv_livecrootn_storage_to_litter(:)
   real(r8), pointer :: hrv_deadcrootn_storage_to_litter(:)
   real(r8), pointer :: hrv_leafn_xfer_to_litter(:)
   real(r8), pointer :: hrv_frootn_xfer_to_litter(:)
   real(r8), pointer :: hrv_livestemn_xfer_to_litter(:)
   real(r8), pointer :: hrv_deadstemn_xfer_to_litter(:)
   real(r8), pointer :: hrv_livecrootn_xfer_to_litter(:)
   real(r8), pointer :: hrv_deadcrootn_xfer_to_litter(:)
!
! local pointers to implicit in/out arrays
   real(r8), pointer :: hrv_leafc_to_litr1c(:)
   real(r8), pointer :: hrv_leafc_to_litr2c(:)
   real(r8), pointer :: hrv_leafc_to_litr3c(:)
   real(r8), pointer :: hrv_frootc_to_litr1c(:)
   real(r8), pointer :: hrv_frootc_to_litr2c(:)
   real(r8), pointer :: hrv_frootc_to_litr3c(:)
   real(r8), pointer :: hrv_livestemc_to_cwdc(:)
   real(r8), pointer :: chrv_deadstemc_to_prod10c(:)
   real(r8), pointer :: chrv_deadstemc_to_prod100c(:)
   real(r8), pointer :: hrv_livecrootc_to_cwdc(:)
   real(r8), pointer :: hrv_deadcrootc_to_cwdc(:)
   real(r8), pointer :: hrv_leafc_storage_to_litr1c(:)
   real(r8), pointer :: hrv_frootc_storage_to_litr1c(:)
   real(r8), pointer :: hrv_livestemc_storage_to_litr1c(:)
   real(r8), pointer :: hrv_deadstemc_storage_to_litr1c(:)
   real(r8), pointer :: hrv_livecrootc_storage_to_litr1c(:)
   real(r8), pointer :: hrv_deadcrootc_storage_to_litr1c(:)
   real(r8), pointer :: hrv_gresp_storage_to_litr1c(:)
   real(r8), pointer :: hrv_leafc_xfer_to_litr1c(:)
   real(r8), pointer :: hrv_frootc_xfer_to_litr1c(:)
   real(r8), pointer :: hrv_livestemc_xfer_to_litr1c(:)
   real(r8), pointer :: hrv_deadstemc_xfer_to_litr1c(:)
   real(r8), pointer :: hrv_livecrootc_xfer_to_litr1c(:)
   real(r8), pointer :: hrv_deadcrootc_xfer_to_litr1c(:)
   real(r8), pointer :: hrv_gresp_xfer_to_litr1c(:)
   real(r8), pointer :: hrv_leafn_to_litr1n(:)
   real(r8), pointer :: hrv_leafn_to_litr2n(:)
   real(r8), pointer :: hrv_leafn_to_litr3n(:)
   real(r8), pointer :: hrv_frootn_to_litr1n(:)
   real(r8), pointer :: hrv_frootn_to_litr2n(:)
   real(r8), pointer :: hrv_frootn_to_litr3n(:)
   real(r8), pointer :: hrv_livestemn_to_cwdn(:)
   real(r8), pointer :: chrv_deadstemn_to_prod10n(:)
   real(r8), pointer :: chrv_deadstemn_to_prod100n(:)
   real(r8), pointer :: hrv_livecrootn_to_cwdn(:)
   real(r8), pointer :: hrv_deadcrootn_to_cwdn(:)
   real(r8), pointer :: hrv_retransn_to_litr1n(:)
   real(r8), pointer :: hrv_leafn_storage_to_litr1n(:)
   real(r8), pointer :: hrv_frootn_storage_to_litr1n(:)
   real(r8), pointer :: hrv_livestemn_storage_to_litr1n(:)
   real(r8), pointer :: hrv_deadstemn_storage_to_litr1n(:)
   real(r8), pointer :: hrv_livecrootn_storage_to_litr1n(:)
   real(r8), pointer :: hrv_deadcrootn_storage_to_litr1n(:)
   real(r8), pointer :: hrv_leafn_xfer_to_litr1n(:)
   real(r8), pointer :: hrv_frootn_xfer_to_litr1n(:)
   real(r8), pointer :: hrv_livestemn_xfer_to_litr1n(:)
   real(r8), pointer :: hrv_deadstemn_xfer_to_litr1n(:)
   real(r8), pointer :: hrv_livecrootn_xfer_to_litr1n(:)
   real(r8), pointer :: hrv_deadcrootn_xfer_to_litr1n(:)
!
! local pointers to implicit out arrays
!
!
! !OTHER LOCAL VARIABLES:
   integer :: fc,c,pi,p               ! indices
!EOP
!-----------------------------------------------------------------------

   ! assign local pointers
   lf_flab                        => pftcon%lf_flab
   lf_fcel                        => pftcon%lf_fcel
   lf_flig                        => pftcon%lf_flig
   fr_flab                        => pftcon%fr_flab
   fr_fcel                        => pftcon%fr_fcel
   fr_flig                        => pftcon%fr_flig

   ! assign local pointers to column-level arrays
   npfts                          => clm3%g%l%c%npfts
   pfti                           => clm3%g%l%c%pfti
   hrv_leafc_to_litr1c              => clm3%g%l%c%ccf%hrv_leafc_to_litr1c
   hrv_leafc_to_litr2c              => clm3%g%l%c%ccf%hrv_leafc_to_litr2c
   hrv_leafc_to_litr3c              => clm3%g%l%c%ccf%hrv_leafc_to_litr3c
   hrv_frootc_to_litr1c             => clm3%g%l%c%ccf%hrv_frootc_to_litr1c
   hrv_frootc_to_litr2c             => clm3%g%l%c%ccf%hrv_frootc_to_litr2c
   hrv_frootc_to_litr3c             => clm3%g%l%c%ccf%hrv_frootc_to_litr3c
   hrv_livestemc_to_cwdc            => clm3%g%l%c%ccf%hrv_livestemc_to_cwdc
   chrv_deadstemc_to_prod10c        => clm3%g%l%c%ccf%hrv_deadstemc_to_prod10c
   chrv_deadstemc_to_prod100c       => clm3%g%l%c%ccf%hrv_deadstemc_to_prod100c
   hrv_livecrootc_to_cwdc           => clm3%g%l%c%ccf%hrv_livecrootc_to_cwdc
   hrv_deadcrootc_to_cwdc           => clm3%g%l%c%ccf%hrv_deadcrootc_to_cwdc
   hrv_leafc_storage_to_litr1c      => clm3%g%l%c%ccf%hrv_leafc_storage_to_litr1c
   hrv_frootc_storage_to_litr1c     => clm3%g%l%c%ccf%hrv_frootc_storage_to_litr1c
   hrv_livestemc_storage_to_litr1c  => clm3%g%l%c%ccf%hrv_livestemc_storage_to_litr1c
   hrv_deadstemc_storage_to_litr1c  => clm3%g%l%c%ccf%hrv_deadstemc_storage_to_litr1c
   hrv_livecrootc_storage_to_litr1c => clm3%g%l%c%ccf%hrv_livecrootc_storage_to_litr1c
   hrv_deadcrootc_storage_to_litr1c => clm3%g%l%c%ccf%hrv_deadcrootc_storage_to_litr1c
   hrv_gresp_storage_to_litr1c      => clm3%g%l%c%ccf%hrv_gresp_storage_to_litr1c
   hrv_leafc_xfer_to_litr1c         => clm3%g%l%c%ccf%hrv_leafc_xfer_to_litr1c
   hrv_frootc_xfer_to_litr1c        => clm3%g%l%c%ccf%hrv_frootc_xfer_to_litr1c
   hrv_livestemc_xfer_to_litr1c     => clm3%g%l%c%ccf%hrv_livestemc_xfer_to_litr1c
   hrv_deadstemc_xfer_to_litr1c     => clm3%g%l%c%ccf%hrv_deadstemc_xfer_to_litr1c
   hrv_livecrootc_xfer_to_litr1c    => clm3%g%l%c%ccf%hrv_livecrootc_xfer_to_litr1c
   hrv_deadcrootc_xfer_to_litr1c    => clm3%g%l%c%ccf%hrv_deadcrootc_xfer_to_litr1c
   hrv_gresp_xfer_to_litr1c         => clm3%g%l%c%ccf%hrv_gresp_xfer_to_litr1c
   hrv_leafn_to_litr1n              => clm3%g%l%c%cnf%hrv_leafn_to_litr1n
   hrv_leafn_to_litr2n              => clm3%g%l%c%cnf%hrv_leafn_to_litr2n
   hrv_leafn_to_litr3n              => clm3%g%l%c%cnf%hrv_leafn_to_litr3n
   hrv_frootn_to_litr1n             => clm3%g%l%c%cnf%hrv_frootn_to_litr1n
   hrv_frootn_to_litr2n             => clm3%g%l%c%cnf%hrv_frootn_to_litr2n
   hrv_frootn_to_litr3n             => clm3%g%l%c%cnf%hrv_frootn_to_litr3n
   hrv_livestemn_to_cwdn            => clm3%g%l%c%cnf%hrv_livestemn_to_cwdn
   chrv_deadstemn_to_prod10n        => clm3%g%l%c%cnf%hrv_deadstemn_to_prod10n
   chrv_deadstemn_to_prod100n       => clm3%g%l%c%cnf%hrv_deadstemn_to_prod100n
   hrv_livecrootn_to_cwdn           => clm3%g%l%c%cnf%hrv_livecrootn_to_cwdn
   hrv_deadcrootn_to_cwdn           => clm3%g%l%c%cnf%hrv_deadcrootn_to_cwdn
   hrv_retransn_to_litr1n           => clm3%g%l%c%cnf%hrv_retransn_to_litr1n
   hrv_leafn_storage_to_litr1n      => clm3%g%l%c%cnf%hrv_leafn_storage_to_litr1n
   hrv_frootn_storage_to_litr1n     => clm3%g%l%c%cnf%hrv_frootn_storage_to_litr1n
   hrv_livestemn_storage_to_litr1n  => clm3%g%l%c%cnf%hrv_livestemn_storage_to_litr1n
   hrv_deadstemn_storage_to_litr1n  => clm3%g%l%c%cnf%hrv_deadstemn_storage_to_litr1n
   hrv_livecrootn_storage_to_litr1n => clm3%g%l%c%cnf%hrv_livecrootn_storage_to_litr1n
   hrv_deadcrootn_storage_to_litr1n => clm3%g%l%c%cnf%hrv_deadcrootn_storage_to_litr1n
   hrv_leafn_xfer_to_litr1n         => clm3%g%l%c%cnf%hrv_leafn_xfer_to_litr1n
   hrv_frootn_xfer_to_litr1n        => clm3%g%l%c%cnf%hrv_frootn_xfer_to_litr1n
   hrv_livestemn_xfer_to_litr1n     => clm3%g%l%c%cnf%hrv_livestemn_xfer_to_litr1n
   hrv_deadstemn_xfer_to_litr1n     => clm3%g%l%c%cnf%hrv_deadstemn_xfer_to_litr1n
   hrv_livecrootn_xfer_to_litr1n    => clm3%g%l%c%cnf%hrv_livecrootn_xfer_to_litr1n
   hrv_deadcrootn_xfer_to_litr1n    => clm3%g%l%c%cnf%hrv_deadcrootn_xfer_to_litr1n

   ! assign local pointers to pft-level arrays
   ivt                            => clm3%g%l%c%p%itype
   wtcol                          => clm3%g%l%c%p%wtcol
   pwtgcell                       => clm3%g%l%c%p%wtgcell  
   hrv_leafc_to_litter              => clm3%g%l%c%p%pcf%hrv_leafc_to_litter
   hrv_frootc_to_litter             => clm3%g%l%c%p%pcf%hrv_frootc_to_litter
   hrv_livestemc_to_litter          => clm3%g%l%c%p%pcf%hrv_livestemc_to_litter
   phrv_deadstemc_to_prod10c        => clm3%g%l%c%p%pcf%hrv_deadstemc_to_prod10c
   phrv_deadstemc_to_prod100c       => clm3%g%l%c%p%pcf%hrv_deadstemc_to_prod100c
   hrv_livecrootc_to_litter         => clm3%g%l%c%p%pcf%hrv_livecrootc_to_litter
   hrv_deadcrootc_to_litter         => clm3%g%l%c%p%pcf%hrv_deadcrootc_to_litter
   hrv_leafc_storage_to_litter      => clm3%g%l%c%p%pcf%hrv_leafc_storage_to_litter
   hrv_frootc_storage_to_litter     => clm3%g%l%c%p%pcf%hrv_frootc_storage_to_litter
   hrv_livestemc_storage_to_litter  => clm3%g%l%c%p%pcf%hrv_livestemc_storage_to_litter
   hrv_deadstemc_storage_to_litter  => clm3%g%l%c%p%pcf%hrv_deadstemc_storage_to_litter
   hrv_livecrootc_storage_to_litter => clm3%g%l%c%p%pcf%hrv_livecrootc_storage_to_litter
   hrv_deadcrootc_storage_to_litter => clm3%g%l%c%p%pcf%hrv_deadcrootc_storage_to_litter
   hrv_gresp_storage_to_litter      => clm3%g%l%c%p%pcf%hrv_gresp_storage_to_litter
   hrv_leafc_xfer_to_litter         => clm3%g%l%c%p%pcf%hrv_leafc_xfer_to_litter
   hrv_frootc_xfer_to_litter        => clm3%g%l%c%p%pcf%hrv_frootc_xfer_to_litter
   hrv_livestemc_xfer_to_litter     => clm3%g%l%c%p%pcf%hrv_livestemc_xfer_to_litter
   hrv_deadstemc_xfer_to_litter     => clm3%g%l%c%p%pcf%hrv_deadstemc_xfer_to_litter
   hrv_livecrootc_xfer_to_litter    => clm3%g%l%c%p%pcf%hrv_livecrootc_xfer_to_litter
   hrv_deadcrootc_xfer_to_litter    => clm3%g%l%c%p%pcf%hrv_deadcrootc_xfer_to_litter
   hrv_gresp_xfer_to_litter         => clm3%g%l%c%p%pcf%hrv_gresp_xfer_to_litter
   hrv_leafn_to_litter              => clm3%g%l%c%p%pnf%hrv_leafn_to_litter
   hrv_frootn_to_litter             => clm3%g%l%c%p%pnf%hrv_frootn_to_litter
   hrv_livestemn_to_litter          => clm3%g%l%c%p%pnf%hrv_livestemn_to_litter
   phrv_deadstemn_to_prod10n        => clm3%g%l%c%p%pnf%hrv_deadstemn_to_prod10n
   phrv_deadstemn_to_prod100n       => clm3%g%l%c%p%pnf%hrv_deadstemn_to_prod100n
   hrv_livecrootn_to_litter         => clm3%g%l%c%p%pnf%hrv_livecrootn_to_litter
   hrv_deadcrootn_to_litter         => clm3%g%l%c%p%pnf%hrv_deadcrootn_to_litter
   hrv_retransn_to_litter           => clm3%g%l%c%p%pnf%hrv_retransn_to_litter
   hrv_leafn_storage_to_litter      => clm3%g%l%c%p%pnf%hrv_leafn_storage_to_litter
   hrv_frootn_storage_to_litter     => clm3%g%l%c%p%pnf%hrv_frootn_storage_to_litter
   hrv_livestemn_storage_to_litter  => clm3%g%l%c%p%pnf%hrv_livestemn_storage_to_litter
   hrv_deadstemn_storage_to_litter  => clm3%g%l%c%p%pnf%hrv_deadstemn_storage_to_litter
   hrv_livecrootn_storage_to_litter => clm3%g%l%c%p%pnf%hrv_livecrootn_storage_to_litter
   hrv_deadcrootn_storage_to_litter => clm3%g%l%c%p%pnf%hrv_deadcrootn_storage_to_litter
   hrv_leafn_xfer_to_litter         => clm3%g%l%c%p%pnf%hrv_leafn_xfer_to_litter
   hrv_frootn_xfer_to_litter        => clm3%g%l%c%p%pnf%hrv_frootn_xfer_to_litter
   hrv_livestemn_xfer_to_litter     => clm3%g%l%c%p%pnf%hrv_livestemn_xfer_to_litter
   hrv_deadstemn_xfer_to_litter     => clm3%g%l%c%p%pnf%hrv_deadstemn_xfer_to_litter
   hrv_livecrootn_xfer_to_litter    => clm3%g%l%c%p%pnf%hrv_livecrootn_xfer_to_litter
   hrv_deadcrootn_xfer_to_litter    => clm3%g%l%c%p%pnf%hrv_deadcrootn_xfer_to_litter

   do pi = 1,maxpatch_pft
      do fc = 1,num_soilc
         c = filter_soilc(fc)

         if (pi <=  npfts(c)) then
            p = pfti(c) + (pi - 1)*clm3%ngridcells !CP from PAWS p for a column is not consecutive

            if (pwtgcell(p)>0._r8) then

               ! leaf harvest mortality carbon fluxes
               hrv_leafc_to_litr1c(c) = hrv_leafc_to_litr1c(c) + &
                  hrv_leafc_to_litter(p) * lf_flab(ivt(p)) * wtcol(p)
               hrv_leafc_to_litr2c(c) = hrv_leafc_to_litr2c(c) + &
                  hrv_leafc_to_litter(p) * lf_fcel(ivt(p)) * wtcol(p)
               hrv_leafc_to_litr3c(c) = hrv_leafc_to_litr3c(c) + &
                  hrv_leafc_to_litter(p) * lf_flig(ivt(p)) * wtcol(p)

               ! fine root harvest mortality carbon fluxes
               hrv_frootc_to_litr1c(c) = hrv_frootc_to_litr1c(c) + &
                  hrv_frootc_to_litter(p) * fr_flab(ivt(p)) * wtcol(p)
               hrv_frootc_to_litr2c(c) = hrv_frootc_to_litr2c(c) + &
                  hrv_frootc_to_litter(p) * fr_fcel(ivt(p)) * wtcol(p)
               hrv_frootc_to_litr3c(c) = hrv_frootc_to_litr3c(c) + &
                  hrv_frootc_to_litter(p) * fr_flig(ivt(p)) * wtcol(p)

               ! wood harvest mortality carbon fluxes
               hrv_livestemc_to_cwdc(c)  = hrv_livestemc_to_cwdc(c)  + &
                  hrv_livestemc_to_litter(p)  * wtcol(p)
               chrv_deadstemc_to_prod10c(c)  = chrv_deadstemc_to_prod10c(c)  + &
                  phrv_deadstemc_to_prod10c(p)  * wtcol(p)
               chrv_deadstemc_to_prod100c(c)  = chrv_deadstemc_to_prod100c(c)  + &
                  phrv_deadstemc_to_prod100c(p)  * wtcol(p)
               hrv_livecrootc_to_cwdc(c) = hrv_livecrootc_to_cwdc(c) + &
                  hrv_livecrootc_to_litter(p) * wtcol(p)
               hrv_deadcrootc_to_cwdc(c) = hrv_deadcrootc_to_cwdc(c) + &
                  hrv_deadcrootc_to_litter(p) * wtcol(p)

               ! storage harvest mortality carbon fluxes
               hrv_leafc_storage_to_litr1c(c)      = hrv_leafc_storage_to_litr1c(c)      + &
                  hrv_leafc_storage_to_litter(p)      * wtcol(p)
               hrv_frootc_storage_to_litr1c(c)     = hrv_frootc_storage_to_litr1c(c)     + &
                  hrv_frootc_storage_to_litter(p)     * wtcol(p)
               hrv_livestemc_storage_to_litr1c(c)  = hrv_livestemc_storage_to_litr1c(c)  + &
                  hrv_livestemc_storage_to_litter(p)  * wtcol(p)
               hrv_deadstemc_storage_to_litr1c(c)  = hrv_deadstemc_storage_to_litr1c(c)  + &
                  hrv_deadstemc_storage_to_litter(p)  * wtcol(p)
               hrv_livecrootc_storage_to_litr1c(c) = hrv_livecrootc_storage_to_litr1c(c) + &
                  hrv_livecrootc_storage_to_litter(p) * wtcol(p)
               hrv_deadcrootc_storage_to_litr1c(c) = hrv_deadcrootc_storage_to_litr1c(c) + &
                  hrv_deadcrootc_storage_to_litter(p) * wtcol(p)
               hrv_gresp_storage_to_litr1c(c)      = hrv_gresp_storage_to_litr1c(c)      + &
                  hrv_gresp_storage_to_litter(p)      * wtcol(p)

               ! transfer harvest mortality carbon fluxes
               hrv_leafc_xfer_to_litr1c(c)      = hrv_leafc_xfer_to_litr1c(c)      + &
                  hrv_leafc_xfer_to_litter(p)      * wtcol(p)
               hrv_frootc_xfer_to_litr1c(c)     = hrv_frootc_xfer_to_litr1c(c)     + &
                  hrv_frootc_xfer_to_litter(p)     * wtcol(p)
               hrv_livestemc_xfer_to_litr1c(c)  = hrv_livestemc_xfer_to_litr1c(c)  + &
                  hrv_livestemc_xfer_to_litter(p)  * wtcol(p)
               hrv_deadstemc_xfer_to_litr1c(c)  = hrv_deadstemc_xfer_to_litr1c(c)  + &
                  hrv_deadstemc_xfer_to_litter(p)  * wtcol(p)
               hrv_livecrootc_xfer_to_litr1c(c) = hrv_livecrootc_xfer_to_litr1c(c) + &
                  hrv_livecrootc_xfer_to_litter(p) * wtcol(p)
               hrv_deadcrootc_xfer_to_litr1c(c) = hrv_deadcrootc_xfer_to_litr1c(c) + &
                  hrv_deadcrootc_xfer_to_litter(p) * wtcol(p)
               hrv_gresp_xfer_to_litr1c(c)      = hrv_gresp_xfer_to_litr1c(c)      + &
                  hrv_gresp_xfer_to_litter(p)      * wtcol(p)

               ! leaf harvest mortality nitrogen fluxes
               hrv_leafn_to_litr1n(c) = hrv_leafn_to_litr1n(c) + &
                  hrv_leafn_to_litter(p) * lf_flab(ivt(p)) * wtcol(p)
               hrv_leafn_to_litr2n(c) = hrv_leafn_to_litr2n(c) + &
                  hrv_leafn_to_litter(p) * lf_fcel(ivt(p)) * wtcol(p)
               hrv_leafn_to_litr3n(c) = hrv_leafn_to_litr3n(c) + &
                  hrv_leafn_to_litter(p) * lf_flig(ivt(p)) * wtcol(p)

               ! fine root litter nitrogen fluxes
               hrv_frootn_to_litr1n(c) = hrv_frootn_to_litr1n(c) + &
                  hrv_frootn_to_litter(p) * fr_flab(ivt(p)) * wtcol(p)
               hrv_frootn_to_litr2n(c) = hrv_frootn_to_litr2n(c) + &
                  hrv_frootn_to_litter(p) * fr_fcel(ivt(p)) * wtcol(p)
               hrv_frootn_to_litr3n(c) = hrv_frootn_to_litr3n(c) + &
                  hrv_frootn_to_litter(p) * fr_flig(ivt(p)) * wtcol(p)

               ! wood harvest mortality nitrogen fluxes
               hrv_livestemn_to_cwdn(c)  = hrv_livestemn_to_cwdn(c)  + &
                  hrv_livestemn_to_litter(p)  * wtcol(p)
               chrv_deadstemn_to_prod10n(c)  = chrv_deadstemn_to_prod10n(c)  + &
                  phrv_deadstemn_to_prod10n(p)  * wtcol(p)
               chrv_deadstemn_to_prod100n(c)  = chrv_deadstemn_to_prod100n(c)  + &
                  phrv_deadstemn_to_prod100n(p)  * wtcol(p)
               hrv_livecrootn_to_cwdn(c) = hrv_livecrootn_to_cwdn(c) + &
                  hrv_livecrootn_to_litter(p) * wtcol(p)
               hrv_deadcrootn_to_cwdn(c) = hrv_deadcrootn_to_cwdn(c) + &
                  hrv_deadcrootn_to_litter(p) * wtcol(p)

               ! retranslocated N pool harvest mortality fluxes
               hrv_retransn_to_litr1n(c) = hrv_retransn_to_litr1n(c) + &
                  hrv_retransn_to_litter(p) * wtcol(p)

               ! storage harvest mortality nitrogen fluxes
               hrv_leafn_storage_to_litr1n(c)      = hrv_leafn_storage_to_litr1n(c)      + &
                  hrv_leafn_storage_to_litter(p)      * wtcol(p)
               hrv_frootn_storage_to_litr1n(c)     = hrv_frootn_storage_to_litr1n(c)     + &
                  hrv_frootn_storage_to_litter(p)     * wtcol(p)
               hrv_livestemn_storage_to_litr1n(c)  = hrv_livestemn_storage_to_litr1n(c)  + &
                  hrv_livestemn_storage_to_litter(p)  * wtcol(p)
               hrv_deadstemn_storage_to_litr1n(c)  = hrv_deadstemn_storage_to_litr1n(c)  + &
                  hrv_deadstemn_storage_to_litter(p)  * wtcol(p)
               hrv_livecrootn_storage_to_litr1n(c) = hrv_livecrootn_storage_to_litr1n(c) + &
                  hrv_livecrootn_storage_to_litter(p) * wtcol(p)
               hrv_deadcrootn_storage_to_litr1n(c) = hrv_deadcrootn_storage_to_litr1n(c) + &
                  hrv_deadcrootn_storage_to_litter(p) * wtcol(p)

               ! transfer harvest mortality nitrogen fluxes
               hrv_leafn_xfer_to_litr1n(c)      = hrv_leafn_xfer_to_litr1n(c)      + &
                  hrv_leafn_xfer_to_litter(p)      * wtcol(p)
               hrv_frootn_xfer_to_litr1n(c)     = hrv_frootn_xfer_to_litr1n(c)     + &
                  hrv_frootn_xfer_to_litter(p)     * wtcol(p)
               hrv_livestemn_xfer_to_litr1n(c)  = hrv_livestemn_xfer_to_litr1n(c)  + &
                  hrv_livestemn_xfer_to_litter(p)  * wtcol(p)
               hrv_deadstemn_xfer_to_litr1n(c)  = hrv_deadstemn_xfer_to_litr1n(c)  + &
                  hrv_deadstemn_xfer_to_litter(p)  * wtcol(p)
               hrv_livecrootn_xfer_to_litr1n(c) = hrv_livecrootn_xfer_to_litr1n(c) + &
                  hrv_livecrootn_xfer_to_litter(p) * wtcol(p)
               hrv_deadcrootn_xfer_to_litr1n(c) = hrv_deadcrootn_xfer_to_litr1n(c) + &
                  hrv_deadcrootn_xfer_to_litter(p) * wtcol(p)

            end if
         end if

      end do

   end do

end subroutine CNHarvestPftToColumn
!-----------------------------------------------------------------------

end module