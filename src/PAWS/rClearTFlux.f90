
!
      subroutine rClearTFlux(rid) 
      ! (E) Clear Transition Fluxes
      ! Output:
      ! r.Riv.Qoc,Qgc,evap,trib
      use Vdata
      implicit none
      integer :: rid
      
      !r(rid).Riv.prcp(:) = 0;
      r(rid)%Riv%Qoc(:) = 0.0D0
      r(rid)%Riv%Qgc(:) = 0.0D0
      !r(rid).Riv.Ql(:) = 0;
      r(rid)%Riv%evap(:) = 0.0D0
      r(rid)%Riv%trib(:) = 0.0D0
      end subroutine rClearTFlux
