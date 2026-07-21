module abortutils

!-----------------------------------------------------------------------
!BOP
! !MODULE: abortutils
!
! !DESCRIPTION:
! Abort the model for abnormal termination
!
! !USES:
  use clm_varctl, only : iulog
!
! !REVISION HISTORY:
! Author: CCM Core group
!
!EOP
!-----------------------------------------------------------------------

   private
   save

   public :: endrun

CONTAINS

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: endrun
!
! !INTERFACE:
   subroutine endrun(msg)
!
! !DESCRIPTION:
! Abort the model for abnormal termination
!
! Author: CCM Core group
!
!-----------------------------------------------------------------------
   !use spmdMod,     only: mpicom
   use shr_sys_mod, only: shr_sys_flush
!
! !ARGUMENTS:
   implicit none
   character(len=*), intent(in), optional :: msg    ! string to be printed
   CHARACTER*1000::msg2
!
! !REVISION HISTORY:
! Author: CCM Core group
!
!EOP
!-----------------------------------------------------------------------
call ch_assert(msg)
!CP modified dummy module
! IF (PRESENT(msg)) then
!   msg2 = msg
! ELSE
!   msg2 = 'endrun hit without msg'
! ENDIF
! WRITE(*,*) msg2
! PAUSE
!STOP ! go out and see where it is!
end subroutine endrun

end module abortutils
