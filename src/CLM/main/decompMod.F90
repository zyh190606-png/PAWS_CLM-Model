! These functions have been modified significantly to compile in a non-MPI mode
! all pretty dummy right now. for HPC, should get back the original MPI procedures

MODULE decompMod
PUBLIC:: get_proc_bounds,get_proc_bounds_atm
CONTAINS
subroutine get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)
USE memoryChunk
IMPLICIT NONE
     integer, optional, intent(out) :: begp, endp  ! proc beginning and ending
                                                   ! pft indices
     integer, optional, intent(out) :: begc, endc  ! proc beginning and ending
                                                   ! column indices
     integer, optional, intent(out) :: begl, endl  ! proc beginning and ending
                                                   ! landunit indices
     integer, optional, intent(out) :: begg, endg  ! proc beginning and ending
                                                   ! gridcell indices


     if (present(begp)) then
        begp = 1
     endif
     if (present(endp)) then
        endp = nn
     endif
     if (present(begc)) then
        begc = 1
     endif
     if (present(endc)) then
        endc = mn
     endif
     if (present(begl)) then
        begl = 1
     endif
     if (present(endl)) then
        endl = mn
     endif
     if (present(begg)) then
        begg = 1
     endif
     if (present(endg)) then
        endg = mn
     endif
!begg=1 
!endg=mn
!begl=1
!endl=mn
!begc=1
!endc=mn
!begp=1
!endp=nn

end subroutine

subroutine get_proc_bounds_atm(begg, endg, begl, endl, begc, endc, begp, endp)
USE memoryChunk
IMPLICIT NONE
integer begg, endg, begl, endl, begc, endc, begp, endp


begg=1 
endg=mn
begl=1
endl=mn
begc=1
endc=mn
begp=1
endp=nn

end subroutine

function gsMap_lnd_gdc2glo()
use nanMOD, only: NaN
IMPLICIT NONE
real*8 gsMap_lnd_gdc2glo
gsMap_lnd_gdc2glo = NAN
end function

subroutine get_proc_global(numg, numl, numc, nump)
USE memoryChunk
IMPLICIT NONE
integer numg, numl, numc, nump
numg = mn
numl = mn
numc = mn
nump = nn

end subroutine


function get_proc_clumps()
USE memoryChunk
IMPLICIT NONE
integer get_proc_clumps
!CP dummy function
get_proc_clumps = 1

end function

subroutine get_clump_bounds(nc, begg, endg, begl, endl, begc, endc, begp, endp)
USE memoryChunk
IMPLICIT NONE
integer begg, endg, begl, endl, begc, endc, begp, endp
integer, intent(in) :: nc
!CP dummy subroutine
!just use the whole thing
begg = 1
endg = mn
begl = 1
endl = mn
begc = 1
endc = mn
begp = 1
endp = nn

end subroutine
END MODULE
