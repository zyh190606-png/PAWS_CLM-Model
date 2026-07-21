module paws_varctl

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: paws_varctl
!
! !DESCRIPTION:
! Module containing run control variables
! The purpose of this module is to hold public data for:
! control information (e.g. DBGFLG, MethodsChoosingInfo)
! file sharing info (e.g. input/output filenames, file open IDs)
! other information
!
! !USES:
!
!
! !PUBLIC MEMBER FUNCTIONS:
  implicit none

  private
  save
!
! Filenames:
!
  integer, PARAMETER :: clen = 500
  character(len=clen), public :: drvfile     = ' '                            ! driver control file
  character(len=clen), public :: parfile     = ' '                            ! parameter change file
  character(len=clen), public :: matinfile   = ' '                            ! input mat file
  character(len=clen), public :: matmonth    = ' '                            ! monthly save mat file
  character(len=clen), public :: budgeterf   = ' '                            ! budgeter output files
  character(len=clen), public :: recf        = ' '                            ! Rec file
  character(len=clen), public :: resultf     = ' '                            ! result file (objective function) for optimization
!
! Run control variables:
!
  logical,             public :: dbgflg                                       ! whether this is a debug run
  logical,             public :: writeout    = .false.                        ! some writing flag
  logical,             public :: writeoutL   = .false.                        ! long term monitoring writing flag
  character(len=clen), public :: prj         = ' '                            ! project name
  integer,             public :: nyrs        = 0                              ! internal variable to store which year of sim this is
  integer,             public :: ndays       = 0                              ! internal variable to store which day of sim this is
  integer,             public :: nmonths     = 0                              ! internal variable to store which month of sim this is
  character*100,       public :: mode        = ''                             ! ='normal' for full sim; ='column' for column only run
  integer,             public :: nrepeat     = 1                              ! how many repeats to do simulation run
  integer,             public :: nMatOutputFreq = 1                           ! output frequency. =0 no month.mat output, =1, month.mat output
  
end module paws_varctl
