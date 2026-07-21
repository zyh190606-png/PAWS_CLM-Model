#include <misc.h>
#include <preproc.h>
module pftvarcon

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: pftvarcon
!
! !DESCRIPTION:
! Module containing vegetation constants and method to
! read and initialize vegetation (PFT) constants.
!
! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use abortutils  , only : endrun
  use clm_varpar  , only : mxpft, numpft, numrad, ivis, inir
  use clm_varctl  , only : iulog
!
! !PUBLIC TYPES:
  implicit none
  save
!
! Vegetation type constants
!
  character(len=40) pftname(0:mxpft) !PFT description

  integer :: noveg                  !value for not vegetated 
  integer :: ndllf_evr_tmp_tree     !value for Needleleaf evergreen temperate tree
  integer :: ndllf_evr_brl_tree     !value for Needleleaf evergreen boreal tree
  integer :: ndllf_dcd_brl_tree     !value for Needleleaf deciduous boreal tree
  integer :: nbrdlf_evr_trp_tree    !value for Broadleaf evergreen tropical tree
  integer :: nbrdlf_evr_tmp_tree    !value for Broadleaf evergreen temperate tree
  integer :: nbrdlf_dcd_trp_tree    !value for Broadleaf deciduous tropical tree
  integer :: nbrdlf_dcd_tmp_tree    !value for Broadleaf deciduous temperate tree
  integer :: nbrdlf_dcd_brl_tree    !value for Broadleaf deciduous boreal tree
  integer :: ntree                  !value for last type of tree
  integer :: nbrdlf_evr_shrub       !value for Broadleaf evergreen shrub
  integer :: nbrdlf_dcd_tmp_shrub   !value for Broadleaf deciduous temperate shrub
  integer :: nbrdlf_dcd_brl_shrub   !value for Broadleaf deciduous boreal shrub
  integer :: nc3_arctic_grass       !value for C3 arctic grass
  integer :: nc3_nonarctic_grass    !value for C3 non-arctic grass
  integer :: nc4_grass              !value for C4 grass
  integer :: npcropmin              !value for first crop
  integer :: ncorn                  !value for corn
  integer :: nscereal               !value for spring temperate cereal
  integer :: nwcereal               !value for winter temperate cereal
  integer :: nsoybean               !value for soybean
  integer :: npcropmax              !value for last prognostic crop in list
  integer :: nc3crop                !value for generic crop
  integer :: nirrig                 !value for irrigated generic crop

  real(r8):: dleaf(0:mxpft)       !characteristic leaf dimension (m)
  real(r8):: c3psn(0:mxpft)       !photosynthetic pathway: 0. = c4, 1. = c3
  real(r8):: mp(0:mxpft)          !slope of conductance-to-photosynthesis relationship
  real(r8):: qe25(0:mxpft)        !quantum efficiency at 25C (umol CO2 / umol photon)
  real(r8):: xl(0:mxpft)          !leaf/stem orientation index
  real(r8):: rhol(0:mxpft,numrad) !leaf reflectance: 1=vis, 2=nir
  real(r8):: rhos(0:mxpft,numrad) !stem reflectance: 1=vis, 2=nir
  real(r8):: taul(0:mxpft,numrad) !leaf transmittance: 1=vis, 2=nir
  real(r8):: taus(0:mxpft,numrad) !stem transmittance: 1=vis, 2=nir
  real(r8):: z0mr(0:mxpft)        !ratio of momentum roughness length to canopy top height (-)
  real(r8):: displar(0:mxpft)     !ratio of displacement height to canopy top height (-)
  real(r8):: roota_par(0:mxpft)   !CLM rooting distribution parameter [1/m]
  real(r8):: rootb_par(0:mxpft)   !CLM rooting distribution parameter [1/m]
  real(r8):: crop(0:mxpft)        ! crop pft: 0. = not crop, 1. = crop pft
  real(r8):: irrigated(0:mxpft)   ! irrigated pft: 0. = not, 1. = irrigated
  real(r8):: smpso(0:mxpft)       !soil water potential at full stomatal opening (mm)
  real(r8):: smpsc(0:mxpft)       !soil water potential at full stomatal closure (mm)
  real(r8):: fnitr(0:mxpft)       !foliage nitrogen limitation factor (-)
  ! begin new pft parameters for CN code
  real(r8):: slatop(0:mxpft)      !SLA at top of canopy [m^2/gC]
  real(r8):: dsladlai(0:mxpft)    !dSLA/dLAI [m^2/gC]
  real(r8):: leafcn(0:mxpft)      !leaf C:N [gC/gN]
  real(r8):: flnr(0:mxpft)        !fraction of leaf N in Rubisco [no units]
  real(r8):: woody(0:mxpft)       !woody lifeform flag (0 or 1)
  real(r8):: lflitcn(0:mxpft)      !leaf litter C:N (gC/gN)
  real(r8):: frootcn(0:mxpft)      !fine root C:N (gC/gN)
  real(r8):: livewdcn(0:mxpft)     !live wood (phloem and ray parenchyma) C:N (gC/gN)
  real(r8):: deadwdcn(0:mxpft)     !dead wood (xylem and heartwood) C:N (gC/gN)
  real(r8):: grperc(0:mxpft)       !growth respiration parameter
  real(r8):: grpnow(0:mxpft)       !growth respiration parameter

! for crop
  real(r8):: graincn(0:mxpft)      !grain C:N (gC/gN)
  real(r8):: mxtmp(0:mxpft)        !parameter used in accFlds
  real(r8):: baset(0:mxpft)        !parameter used in accFlds
  real(r8):: declfact(0:mxpft)     !parameter used in CNAllocation
  real(r8):: bfact(0:mxpft)        !parameter used in CNAllocation
  real(r8):: aleaff(0:mxpft)       !parameter used in CNAllocation
  real(r8):: arootf(0:mxpft)       !parameter used in CNAllocation
  real(r8):: astemf(0:mxpft)       !parameter used in CNAllocation
  real(r8):: arooti(0:mxpft)       !parameter used in CNAllocation
  real(r8):: fleafi(0:mxpft)       !parameter used in CNAllocation
  real(r8):: allconsl(0:mxpft)     !parameter used in CNAllocation
  real(r8):: allconss(0:mxpft)     !parameter used in CNAllocation
  real(r8):: ztopmx(0:mxpft)       !parameter used in CNVegStructUpdate
  real(r8):: laimx(0:mxpft)        !parameter used in CNVegStructUpdate
  real(r8):: gddmin(0:mxpft)       !parameter used in CNPhenology
  real(r8):: hybgdd(0:mxpft)       !parameter used in CNPhenology
  real(r8):: lfemerg(0:mxpft)      !parameter used in CNPhenology
  real(r8):: grnfill(0:mxpft)      !parameter used in CNPhenology
  integer :: mxmat(0:mxpft)        !parameter used in CNPhenology
  integer :: mnNHplantdate(0:mxpft)!minimum planting date for NorthHemisphere (YYYYMMDD)
  integer :: mxNHplantdate(0:mxpft)!maximum planting date for NorthHemisphere (YYYYMMDD)
  integer :: mnSHplantdate(0:mxpft)!minimum planting date for SouthHemisphere (YYYYMMDD)
  integer :: mxSHplantdate(0:mxpft)!maximum planting date for SouthHemisphere (YYYYMMDD)
  real(r8):: planttemp(0:mxpft)    !planting temperature used in CNPhenology (K)
  real(r8):: minplanttemp(0:mxpft) !mininum planting temperature used in CNPhenology (K)
  real(r8):: froot_leaf(0:mxpft)   !allocation parameter: new fine root C per new leaf C (gC/gC) 
  real(r8):: stem_leaf(0:mxpft)    !allocation parameter: new stem c per new leaf C (gC/gC)
  real(r8):: croot_stem(0:mxpft)   !allocation parameter: new coarse root C per new stem C (gC/gC)
  real(r8):: flivewd(0:mxpft)      !allocation parameter: fraction of new wood that is live (phloem and ray parenchyma) (no units)
  real(r8):: fcur(0:mxpft)         !allocation parameter: fraction of allocation that goes to currently displayed growth, remainder to storage
  real(r8):: fcurdv(0:mxpft)       !alternate fcur for use with cndv
  real(r8):: lf_flab(0:mxpft)      !leaf litter labile fraction
  real(r8):: lf_fcel(0:mxpft)      !leaf litter cellulose fraction
  real(r8):: lf_flig(0:mxpft)      !leaf litter lignin fraction
  real(r8):: fr_flab(0:mxpft)      !fine root litter labile fraction
  real(r8):: fr_fcel(0:mxpft)      !fine root litter cellulose fraction
  real(r8):: fr_flig(0:mxpft)      !fine root litter lignin fraction
  real(r8):: leaf_long(0:mxpft)    !leaf longevity (yrs)
  real(r8):: evergreen(0:mxpft)    !binary flag for evergreen leaf habit (0 or 1)
  real(r8):: stress_decid(0:mxpft) !binary flag for stress-deciduous leaf habit (0 or 1)
  real(r8):: season_decid(0:mxpft) !binary flag for seasonal-deciduous leaf habit (0 or 1)
  real(r8):: pconv(0:mxpft)        !proportion of deadstem to conversion flux
  real(r8):: pprod10(0:mxpft)      !proportion of deadstem to 10-yr product pool
  real(r8):: pprod100(0:mxpft)     !proportion of deadstem to 100-yr product pool
  real(r8):: pprodharv10(0:mxpft)  !harvest mortality proportion of deadstem to 10-yr pool

  ! new pft parameters for CN-fire code
  real(r8):: resist(0:mxpft)       !resistance to fire (no units)

  ! pft parameters for CNDV code
  ! from LPJ subroutine pftparameters
  real(r8) pftpar20(0:mxpft)       !tree maximum crown area (m2)
  real(r8) pftpar28(0:mxpft)       !min coldest monthly mean temperature
  real(r8) pftpar29(0:mxpft)       !max coldest monthly mean temperature
  real(r8) pftpar30(0:mxpft)       !min growing degree days (>= 5 deg C)
  real(r8) pftpar31(0:mxpft)       !upper limit of temperature of the warmest month (twmax)
  real(r8), parameter :: reinickerp = 1.6_r8 !parameter in allometric equation
  real(r8), parameter :: dwood  = 2.5e5_r8   !cn wood density (gC/m3); lpj:2.0e5
  real(r8), parameter :: allom1 = 100.0_r8   !parameters in
  real(r8), parameter :: allom2 =  40.0_r8   !...allometric
  real(r8), parameter :: allom3 =   0.5_r8   !...equations
  real(r8), parameter :: allom1s = 250.0_r8  !modified for shrubs by
  real(r8), parameter :: allom2s =   8.0_r8  !X.D.Z
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: pftconrd ! Read and initialize vegetation (PFT) constants
  PRIVATE :: r8
!
! !REVISION HISTORY:
! Created by Sam Levis (put into module form by Mariana Vertenstein)
! 10/21/03, Peter Thornton: Added new variables for CN code
! 06/24/09, Erik Kluzek: Add indices for all pft types, and add expected_pftnames array and comparision
! 09/17/10, David Lawrence: Modified code to read in netCDF pft physiology file
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: pftconrd
!
! !INTERFACE:
  subroutine pftconrd
!
! !DESCRIPTION:
! Read and initialize vegetation (PFT) constants
!
! !USES:
    use fileutils , only : getfil
!    use ncdio_pio , only : ncd_io, ncd_pio_closefile, ncd_pio_openfile, file_desc_t, &
!                           ncd_inqdid, ncd_inqdlen
    use clm_varctl, only : fpftcon
    use clm_varcon, only : tfrz
    use spmdMod   , only : masterproc, mpicom, MPI_REAL8,MPI_INTEGER !MPI_CHARACTER, 
    use nanMod    , only : nan
    use matTable
    
!
! !ARGUMENTS:
    implicit none
!
! !CALLED FROM:
! routine initialize in module initializeMod
!
! !REVISION HISTORY:
! Created by Gordon Bonan
!
!
! !LOCAL VARIABLES:
!EOP
    character(len=256) :: locfn ! local file name
    integer :: i,n              ! loop indices
    integer :: ier              ! error code
    !type(file_desc_t) :: ncid   ! pio netCDF file id
    integer :: dimid            ! netCDF dimension id
    integer :: npft             ! number of pfts on pft-physiology file
    logical :: readv            ! read variable in or not
    character(len=32) :: subname = 'pftconrd'              ! subroutine name
    !
    ! Expected PFT names: The names expected on the fpftcon file and the order they are expected to be in.
    ! NOTE: similar types are assumed to be together, first trees (ending with broadleaf_deciduous_boreal_tree
    !       then shrubs, ending with broadleaf_deciduous_boreal_shrub, then grasses starting with c3_arctic_grass
    !       and finally crops, ending with soybean
    ! DO NOT CHANGE THE ORDER -- WITHOUT MODIFYING OTHER PARTS OF THE CODE WHERE THE ORDER MATTERS!
    !
    character(len=40), parameter :: expected_pftnames(0:mxpft) = (/ &
                 'not_vegetated                      '  &
               , 'needleleaf_evergreen_temperate_tree'  &
               , 'needleleaf_evergreen_boreal_tree   '  &
               , 'needleleaf_deciduous_boreal_tree   '  &
               , 'broadleaf_evergreen_tropical_tree  '  &
               , 'broadleaf_evergreen_temperate_tree '  &
               , 'broadleaf_deciduous_tropical_tree  '  &
               , 'broadleaf_deciduous_temperate_tree '  &
               , 'broadleaf_deciduous_boreal_tree    '  &
               , 'broadleaf_evergreen_shrub          '  &
               , 'broadleaf_deciduous_temperate_shrub'  &
               , 'broadleaf_deciduous_boreal_shrub   '  &
               , 'c3_arctic_grass                    '  &
               , 'c3_non-arctic_grass                '  &
               , 'c4_grass                           '  &
               , 'c3_crop                            '  &
               , 'c3_irrigated                       '  &
               , 'corn                               '  &
               , 'spring_temperate_cereal            '  &
               , 'winter_temperate_cereal            '  &
               , 'soybean                            '  &
    /)
!-----------------------------------------------------------------------

    ! Set specific vegetation type values
#if (defined PAWS)
    real*8,pointer    :: ptr1(:),ptr2(:,:),ptr3(:,:,:)
    integer,pointer    :: ptr1L(:),ptr2L(:,:),ptr3L(:,:,:)
     
    gd_ptr => gd_base .G. fpftcon
!    fieldname =(/'pftname','z0mr','displar','dleaf','c3psn','mp','qe25','rholvis', &
!   &  'rholnir','rhosvis','rhosnir','taulvis','taulnir','tausvis','tausnir','xl',&
!   & 'roota_par','rootb_par','slatop','dsladlai','leafcn','flnr','smpso','smpsc',&
!   & 'fnitr','woody','lflitcn','frootcn','livewdcn','deadwdcn','grperc','grpnow',&
!   & 'froot_leaf','stem_leaf','croot_stem','flivewd','fcur','fcurdv','lf_flab',  &
!   & 'lf_fcel','lf_flig','fr_flab','fr_fcel','fr_flig','leaf_long','evergreen',  &
!   & 'stress_decid','season_decid','resist','pftpar20','pftpar28','pftpar29',    &
!   & 'pftpar30','pftpar31','pconv','pprod10','pprodharv10','pprod100','graincn', &
!   & 'mxtmp','baset','declfact','bfact', 'aleaff','arootf','astemf', 'arooti',   &
!   & 'fleafi', 'allconsl', 'allconss', 'crop', 'irrigated', 'ztopmx', 'laimx',   &
!   & 'gddmin', 'hybgdd', 'lfemerg', 'grnfill', 'mxmat', 'planting_temp',         &
!   &  'min_planting_temp', 'min_NH_planting_date', 'min_SH_planting_date',       &
!   &  'max_NH_planting_date', 'max_SH_planting_date' /)
    
    pftname = (/                              &
   & 'not_vegetated                      ',       &
   & 'needleleaf_evergreen_temperate_tree',   &
   & 'needleleaf_evergreen_boreal_tree   ',      &
   & 'needleleaf_deciduous_boreal_tree   ',      &
   & 'broadleaf_evergreen_tropical_tree  ',     &
   & 'broadleaf_evergreen_temperate_tree ',    &
   & 'broadleaf_deciduous_tropical_tree  ',     &
   & 'broadleaf_deciduous_temperate_tree ',    &
   & 'broadleaf_deciduous_boreal_tree    ',       &
   & 'broadleaf_evergreen_shrub          ',             &
   & 'broadleaf_deciduous_temperate_shrub',   &
   & 'broadleaf_deciduous_boreal_shrub   ',      &
   & 'c3_arctic_grass                    ',                       &
   & 'c3_non-arctic_grass                ',                   &
   & 'c4_grass                           ',                              &
   & 'c3_crop                            ',                               &
   & 'c3_irrigated                       ',                          &
   & 'corn                               ',                                  &
   & 'spring_temperate_cereal            ',               &
   & 'winter_temperate_cereal            ',               &
   & 'soybean                            '/)
    
    call getPtr(gd_ptr,'z0mr',ptr1); z0mr = ptr1;
    call getPtr(gd_ptr,'displar',ptr1); displar = ptr1;
    call getPtr(gd_ptr,'dleaf',ptr1); dleaf = ptr1;
    call getPtr(gd_ptr,'c3psn',ptr1); c3psn = ptr1;
    call getPtr(gd_ptr,'mp',ptr1); mp = ptr1;
    call getPtr(gd_ptr,'qe25',ptr1); qe25 = ptr1;
    call getPtr(gd_ptr,'rholvis',ptr1); rhol(:,ivis) = ptr1;
    call getPtr(gd_ptr,'rholnir',ptr1); rhol(:,inir) = ptr1;
    call getPtr(gd_ptr,'rhosvis',ptr1); rhos(:,ivis) = ptr1;
    call getPtr(gd_ptr,'rhosnir',ptr1); rhos(:,inir) = ptr1;
    call getPtr(gd_ptr,'taulvis',ptr1); taul(:,ivis) = ptr1;
    call getPtr(gd_ptr,'taulnir',ptr1); taul(:,inir) = ptr1;
    call getPtr(gd_ptr,'tausvis',ptr1); taus(:,ivis) = ptr1;
    call getPtr(gd_ptr,'tausnir',ptr1); taus(:,inir) = ptr1;
    call getPtr(gd_ptr,'xl',ptr1); xl = ptr1;
    call getPtr(gd_ptr,'roota_par',ptr1); roota_par = ptr1;
    call getPtr(gd_ptr,'rootb_par',ptr1); rootb_par = ptr1;
    call getPtr(gd_ptr,'slatop',ptr1); slatop = ptr1;
    call getPtr(gd_ptr,'dsladlai',ptr1); dsladlai = ptr1;
    call getPtr(gd_ptr,'leafcn',ptr1); leafcn = ptr1;
    call getPtr(gd_ptr,'flnr',ptr1); flnr = ptr1;
    call getPtr(gd_ptr,'smpso',ptr1); smpso = ptr1;
    call getPtr(gd_ptr,'smpsc',ptr1); smpsc = ptr1;
    call getPtr(gd_ptr,'fnitr',ptr1); fnitr = ptr1;
    call getPtr(gd_ptr,'woody',ptr1); woody = ptr1;
    call getPtr(gd_ptr,'lflitcn',ptr1); lflitcn = ptr1;
    call getPtr(gd_ptr,'frootcn',ptr1); frootcn = ptr1;
    call getPtr(gd_ptr,'livewdcn',ptr1); livewdcn = ptr1;
    call getPtr(gd_ptr,'deadwdcn',ptr1); deadwdcn = ptr1;
    call getPtr(gd_ptr,'grperc',ptr1); grperc = ptr1;
    call getPtr(gd_ptr,'grpnow',ptr1); grpnow = ptr1;
    call getPtr(gd_ptr,'froot_leaf',ptr1); froot_leaf = ptr1;
    call getPtr(gd_ptr,'stem_leaf',ptr1); stem_leaf = ptr1;
    call getPtr(gd_ptr,'croot_stem',ptr1); croot_stem = ptr1;
    call getPtr(gd_ptr,'flivewd',ptr1); flivewd = ptr1;
    call getPtr(gd_ptr,'fcur',ptr1); fcur = ptr1;
    call getPtr(gd_ptr,'fcurdv',ptr1); fcurdv = ptr1;
    call getPtr(gd_ptr,'lf_flab',ptr1); lf_flab = ptr1;
    call getPtr(gd_ptr,'lf_fcel',ptr1); lf_fcel = ptr1;
    call getPtr(gd_ptr,'lf_flig',ptr1); lf_flig = ptr1;
    call getPtr(gd_ptr,'fr_flab',ptr1); fr_flab = ptr1;
    call getPtr(gd_ptr,'fr_fcel',ptr1); fr_fcel = ptr1;
    call getPtr(gd_ptr,'fr_flig',ptr1); fr_flig = ptr1;
    call getPtr(gd_ptr,'leaf_long',ptr1); leaf_long = ptr1;
    call getPtr(gd_ptr,'evergreen',ptr1); evergreen = ptr1;
    call getPtr(gd_ptr,'stress_decid',ptr1); stress_decid = ptr1;
    call getPtr(gd_ptr,'season_decid',ptr1); season_decid = ptr1;
    call getPtr(gd_ptr,'resist',ptr1); resist = ptr1;
    call getPtr(gd_ptr,'pftpar20',ptr1); pftpar20 = ptr1;
    call getPtr(gd_ptr,'pftpar28',ptr1); pftpar28 = ptr1;
    call getPtr(gd_ptr,'pftpar29',ptr1); pftpar29 = ptr1;
    call getPtr(gd_ptr,'pftpar30',ptr1); pftpar30 = ptr1;
    call getPtr(gd_ptr,'pftpar31',ptr1); pftpar31 = ptr1;
    call getPtr(gd_ptr,'pconv',ptr1); pconv = ptr1;
    call getPtr(gd_ptr,'pprod10',ptr1); pprod10 = ptr1;
    call getPtr(gd_ptr,'pprodharv10',ptr1); pprodharv10 = ptr1;
    call getPtr(gd_ptr,'pprod100',ptr1); pprod100 = ptr1;
    call getPtr(gd_ptr,'graincn',ptr1); graincn = ptr1;
    call getPtr(gd_ptr,'mxtmp',ptr1); mxtmp = ptr1;
    call getPtr(gd_ptr,'baset',ptr1); baset = ptr1;
    call getPtr(gd_ptr,'declfact',ptr1); declfact = ptr1;
    call getPtr(gd_ptr,'bfact',ptr1); bfact = ptr1;
    call getPtr(gd_ptr,'aleaff',ptr1); aleaff = ptr1;
    call getPtr(gd_ptr,'arootf',ptr1); arootf = ptr1;
    call getPtr(gd_ptr,'astemf',ptr1); astemf = ptr1;
    call getPtr(gd_ptr,'arooti',ptr1); arooti = ptr1;
    call getPtr(gd_ptr,'fleafi',ptr1); fleafi = ptr1;
    call getPtr(gd_ptr,'allconsl',ptr1); allconsl = ptr1;
    call getPtr(gd_ptr,'allconss',ptr1); allconss = ptr1;
    call getPtr(gd_ptr,'crop',ptr1); crop = ptr1;
    call getPtr(gd_ptr,'irrigated',ptr1); irrigated = ptr1;
    call getPtr(gd_ptr,'ztopmx',ptr1); ztopmx = ptr1;
    call getPtr(gd_ptr,'laimx',ptr1); laimx = ptr1;
    call getPtr(gd_ptr,'gddmin',ptr1); gddmin = ptr1;
    call getPtr(gd_ptr,'hybgdd',ptr1); hybgdd = ptr1;
    call getPtr(gd_ptr,'lfemerg',ptr1); lfemerg = ptr1;
    call getPtr(gd_ptr,'grnfill',ptr1); grnfill = ptr1;
    call getPtr(gd_ptr,'mxmat',ptr1L); mxmat = ptr1L;
    call getPtr(gd_ptr,'planting_temp',ptr1); planttemp = ptr1;
    call getPtr(gd_ptr,'min_planting_temp',ptr1); minplanttemp = ptr1;
    call getPtr(gd_ptr,'min_NH_planting_date',ptr1L); mnNHplantdate = ptr1L;
    call getPtr(gd_ptr,'min_SH_planting_date',ptr1L); mnSHplantdate = ptr1L;
    call getPtr(gd_ptr,'max_NH_planting_date',ptr1L); mxNHplantdate = ptr1L;
    call getPtr(gd_ptr,'max_SH_planting_date',ptr1L); mxSHplantdate = ptr1L;
    if (laimx(mxpft) .ne. 6D0 .or. dsladlai(1).ne. 0.001250D0) then
       call endrun('alignment problem with reading pftcon.mat')
    endif
    
#else   
    if (masterproc) then
       write(iulog,*) 'Attempting to read PFT physiological data .....'
    end if
    call getfil (fpftcon, locfn, 0)
    call ncd_pio_openfile (ncid, trim(locfn), 0)
    call ncd_inqdid(ncid,'pft',dimid)
    call ncd_inqdlen(ncid,dimid,npft)

    call ncd_io('pftname',pftname, 'read', ncid, readvar=readv, posNOTonfile=.true.) 
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('z0mr',z0mr, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('displar',displar, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('dleaf',dleaf, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('c3psn',c3psn, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('mp',mp, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('qe25',qe25, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('rholvis',rhol(:,ivis), 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('rholnir',rhol(:,inir), 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('rhosvis',rhos(:,ivis), 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('rhosnir', rhos(:,inir), 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('taulvis',taul(:,ivis), 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('taulnir',taul(:,inir), 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('tausvis',taus(:,ivis), 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('tausnir',taus(:,inir), 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('xl',xl, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('roota_par',roota_par, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('rootb_par',rootb_par, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('slatop',slatop, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('dsladlai',dsladlai, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('leafcn',leafcn, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('flnr',flnr, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('smpso',smpso, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('smpsc',smpsc, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('fnitr',fnitr, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('woody',woody, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('lflitcn',lflitcn, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('frootcn',frootcn, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('livewdcn',livewdcn, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('deadwdcn',deadwdcn, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('grperc',grperc, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('grpnow',grpnow, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('froot_leaf',froot_leaf, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('stem_leaf',stem_leaf, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('croot_stem',croot_stem, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('flivewd',flivewd, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('fcur',fcur, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('fcurdv',fcurdv, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('lf_flab',lf_flab, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('lf_fcel',lf_fcel, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('lf_flig',lf_flig, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('fr_flab',fr_flab, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('fr_fcel',fr_fcel, 'read', ncid, readvar=readv, posNOTonfile=.true.)    
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('fr_flig',fr_flig, 'read', ncid, readvar=readv, posNOTonfile=.true.)    
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('leaf_long',leaf_long, 'read', ncid, readvar=readv, posNOTonfile=.true.)    
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('evergreen',evergreen, 'read', ncid, readvar=readv, posNOTonfile=.true.)    
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('stress_decid',stress_decid, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('season_decid',season_decid, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('resist',resist, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('pftpar20',pftpar20, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('pftpar28',pftpar28, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('pftpar29',pftpar29, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('pftpar30',pftpar30, 'read', ncid, readvar=readv, posNOTonfile=.true.)  
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('pftpar31',pftpar31, 'read', ncid, readvar=readv, posNOTonfile=.true.)  
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('pconv',pconv, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('pprod10',pprod10, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('pprodharv10',pprodharv10, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('pprod100',pprod100, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('graincn',graincn, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('mxtmp',mxtmp, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('baset',baset, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('declfact',declfact, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('bfact',bfact, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('aleaff',aleaff, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('arootf',arootf, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('astemf',astemf, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('arooti',arooti, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('fleafi',fleafi, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('allconsl',allconsl, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('allconss',allconss, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('crop',crop, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('irrigated',irrigated, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('ztopmx',ztopmx, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('laimx',laimx, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('gddmin',gddmin, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('hybgdd',hybgdd, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('lfemerg',lfemerg, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('grnfill',grnfill, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('mxmat',mxmat, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('planting_temp',planttemp, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('min_planting_temp',minplanttemp, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('min_NH_planting_date',mnNHplantdate, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('min_SH_planting_date',mnSHplantdate, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('max_NH_planting_date',mxNHplantdate, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('max_SH_planting_date',mxSHplantdate, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_pio_closefile(ncid)
#endif
    do i = 0, mxpft
       if ( trim(adjustl(pftname(i))) /= trim(expected_pftnames(i)) )then
          write(iulog,*)'pftconrd: pftname is NOT what is expected, name = ', &
                        trim(pftname(i)), ', expected name = ', trim(expected_pftnames(i))
          call endrun( 'pftconrd: bad name for pft on fpftcon dataset' )
       end if
       if ( trim(pftname(i)) == 'not_vegetated'                       ) noveg               = i
       if ( trim(pftname(i)) == 'needleleaf_evergreen_temperate_tree' ) ndllf_evr_tmp_tree  = i
       if ( trim(pftname(i)) == 'needleleaf_evergreen_boreal_tree'    ) ndllf_evr_brl_tree  = i
       if ( trim(pftname(i)) == 'needleleaf_deciduous_boreal_tree'    ) ndllf_dcd_brl_tree  = i
       if ( trim(pftname(i)) == 'broadleaf_evergreen_tropical_tree'   ) nbrdlf_evr_trp_tree  = i
       if ( trim(pftname(i)) == 'broadleaf_evergreen_temperate_tree'  ) nbrdlf_evr_tmp_tree  = i
       if ( trim(pftname(i)) == 'broadleaf_deciduous_tropical_tree'   ) nbrdlf_dcd_trp_tree  = i
       if ( trim(pftname(i)) == 'broadleaf_deciduous_temperate_tree'  ) nbrdlf_dcd_tmp_tree  = i
       if ( trim(pftname(i)) == 'broadleaf_deciduous_boreal_tree'     ) nbrdlf_dcd_brl_tree  = i
       if ( trim(pftname(i)) == 'broadleaf_evergreen_shrub'           ) nbrdlf_evr_shrub     = i
       if ( trim(pftname(i)) == 'broadleaf_deciduous_temperate_shrub' ) nbrdlf_dcd_tmp_shrub = i
       if ( trim(pftname(i)) == 'broadleaf_deciduous_boreal_shrub'    ) nbrdlf_dcd_brl_shrub = i
       if ( trim(pftname(i)) == 'c3_arctic_grass'                     ) nc3_arctic_grass     = i
       if ( trim(pftname(i)) == 'c3_non-arctic_grass'                 ) nc3_nonarctic_grass  = i
       if ( trim(pftname(i)) == 'c4_grass'                            ) nc4_grass            = i
       if ( trim(pftname(i)) == 'c3_crop'                             ) nc3crop              = i
       if ( trim(pftname(i)) == 'c3_irrigated'                        ) nirrig               = i
       if ( trim(pftname(i)) == 'corn'                                ) ncorn                = i
       if ( trim(pftname(i)) == 'spring_temperate_cereal'             ) nscereal             = i
       if ( trim(pftname(i)) == 'winter_temperate_cereal'             ) nwcereal             = i
       if ( trim(pftname(i)) == 'soybean'                             ) nsoybean             = i
    end do

    ntree                = nbrdlf_dcd_brl_tree  ! value for last type of tree
    npcropmin            = ncorn                ! first prognostic crop
    npcropmax            = nsoybean             ! last prognostic crop in list

#if (defined CNDV)
    fcur(:) = fcurdv(:)
#endif

    !
    ! Do some error checking
    !
    if ( npcropmax /= mxpft )then
       call endrun( trim(subname)//' ERROR: npcropmax is NOT the last value' )
    end if
    do i = 0, mxpft
       if (     (irrigated(i) == 1.0_r8) .and. i == nirrig )then
          ! correct
       else if ( irrigated(i) == 0.0_r8 )then
          ! correct
       else
          call endrun( trim(subname)//' ERROR: irrigated has wrong values' )
       end if
       if (      crop(i) == 1.0_r8 .and. (i >= nc3crop .and. i <= npcropmax) )then
          ! correct
       else if ( crop(i) == 0.0_r8 )then
          ! correct
       else
          call endrun( trim(subname)//' ERROR: crop has wrong values' )
       end if
       if ( (i /= noveg) .and. (i < npcropmin) .and. &
            abs(pconv(i)+pprod10(i)+pprod100(i) - 1.0_r8) > 1.e-7_r8 )then
          call endrun( trim(subname)//' ERROR: pconv+pprod10+pprod100 do NOT sum to one.' )
       end if
       if ( pprodharv10(i) > 1.0_r8 .or. pprodharv10(i) < 0.0_r8 )then
          call endrun( trim(subname)//' ERROR: pprodharv10 outside of range.' )
       end if
    end do

    if (masterproc) then
       write(iulog,*) 'Successfully read PFT physiological data'
       write(iulog,*)
    end if

  end subroutine pftconrd

end module pftvarcon
