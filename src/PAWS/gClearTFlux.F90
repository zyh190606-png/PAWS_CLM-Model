!

      subroutine gClearTFlux()
      !% (E) This is to clean up transitional fluxes, applied after each land cell step
      use Vdata
      implicit none
      integer :: i, kk, iGW
      integer :: rid
      
      do i = 1, size(g%Riv%ID)
        rid = g%Riv%ID(i)
        r(rid)%Riv%fG = 0.0D0    !% Clear Accumulative Transition Fluxes from river to groundwater
      enddo
      
      kk = 1
      do iGW = 2, size(g%GW)
        if (abs(g%GW(iGW)%tt - w%t) < 1D-4) then
          g%GW(iGW-1)%Qperc = 0.0D0   !% Clean the accumulative percolation flux
          if (kk==1) then
            w%Rec%Mem(5) = w%Rec%Mem(5)-w%Rec%Mem(6)    !% This is to get mass balance in budgeter right
            w%Rec%Mem(1) = w%Rec%Mem(1)-w%Rec%Mem(6)
            kk = 0
          endif
        endif
      enddo
      
      !% ##### THINK ABOUT MINIMIZING OPERATIONS
      !% Some of these are unnecessary
      g%VDZ%Dperc = 0.0D0
      g%VDZ%DF = 0.0D0
      g%Veg%Gprcp = 0.0D0
      g%Veg%EvapG = 0.0D0
      g%Veg%PET = 0.0D0
      g%Veg%SNOM = 0.0D0
      g%OVN%Rf = 0.0D0
      g%OVN%Qoc = 0.0D0
      g%OVN%Ex = 0.0D0
      w%Stats%tGprcp = 0.0D0
      !%for i=1:length(g(gid).GW)-1,g(gid).GW(1).Qperc(:)=0; end      
      
      end subroutine gClearTFlux
