#include <misc.h>
#include <preproc.h>
#define BUILDPIO
#undef  BUILDPIO
#undef  SWITCH_DIMS
#define SWITCH_DIMS

module ncdio
! originally contains parallel IO routines. all that is removed by CP. CP is dark-minded.
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: ncdioMod
!
! !DESCRIPTION:
! Generic interfaces to write fields to netcdf files
!
! !USES:
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use spmdMod        , only : masterproc, mpicom, MPI_REAL8, MPI_INTEGER, &
                              MPI_LOGICAL
  !use clmtype        , only : gratm, grlnd, nameg, namel, namec, namep, allrof !CP
  use clm_varcon     , only : spval,ispval
  use shr_sys_mod    , only : shr_sys_flush
  use abortutils     , only : endrun
  use clm_varctl     , only : single_column
  use clm_varctl     , only : iulog
  use clm_varctl     , only : ncd_lowmem2d, ncd_pio_def
  use clm_varctl     , only : ncd_pio_UseRearranger, ncd_pio_useBoxRearr
  use clm_varctl     , only : ncd_pio_SerialCDF, ncd_pio_IODOF_rootonly
  use clm_varctl     , only : ncd_pio_DebugLevel, ncd_pio_num_iotasks
  !use clm_mct_mod
  !use spmdGathScatMod
  !use decompMod      , only : get_clmlevel_gsize,get_clmlevel_dsize !CP
  !use perf_mod       , only : t_startf, t_stopf !CP
#if (defined BUILDPIO)
  use piolib_mod     ! _EXTERNAL
  use pio_types      ! _EXTERNAL
  use pio_kinds      , only : pio_offset
#endif
!
! !PUBLIC TYPES:
  implicit none
  include 'netcdf.inc'
  save

!
! !PUBLIC MEMBER FUNCTIONS:
!
  public :: check_ret   ! checks return status of netcdf calls

  integer,parameter,public :: ncd_int       = nf_int
  integer,parameter,public :: ncd_float     = nf_float
  integer,parameter,public :: ncd_double    = nf_double
  integer,parameter,public :: ncd_char      = nf_char
  integer,parameter,public :: ncd_global    = nf_global
  integer,parameter,public :: ncd_write     = nf_write
  integer,parameter,public :: ncd_nowrite   = nf_nowrite
  integer,parameter,public :: ncd_clobber   = nf_clobber
  integer,parameter,public :: ncd_noclobber = nf_noclobber
  integer,parameter,public :: ncd_share     = nf_share
  integer,parameter,public :: ncd_fill      = nf_fill
  integer,parameter,public :: ncd_nofill    = nf_nofill
  integer,parameter,public :: ncd_unlimited = nf_unlimited

!
! !REVISION HISTORY:
!
!
! !PRIVATE MEMBER FUNCTIONS:
!
  
!EOP
!-----------------------------------------------------------------------

contains

!--------------------------------------------------------------------
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: check_ret
!
! !INTERFACE:
  subroutine check_ret(ret, cstring)
!
! !DESCRIPTION:
! Check return status from netcdf call
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: ret
    character(len=*) :: cstring
!
! !REVISION HISTORY:
!
!EOP
    character(len=*),parameter :: subname='check_ret' ! subroutine name
!-----------------------------------------------------------------------

    if (ret /= NF_NOERR) then
       write(iulog,*)'netcdf error from ',trim(subname),':',trim(cstring),':',trim(NF_STRERROR(ret))
       call endrun()
    end if

  end subroutine check_ret

end module ncdio
