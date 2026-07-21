!

      subroutine wSyncTime(T)
      !% Synchronize every component of a watershed (and all its subwatersheds)
      !% to be at time T
      use Vdata
      implicit none
      integer :: T, i, k, c
      real*8 :: gid
      
      w%t = T
      !for i = 1:length(w(wid).child)
      !    c = w(wid).child(i);
      !    wSyncTime(c,T);
      !end
      !do i = 1, size(w%g)   !% there should only be 1
      do i = 1, int(w%g)
          gid = w%g
          g%t = T + 0.0D0
          g%OVN%t = T + 0.0D0
          g%VDZ%t = T + 0.0D0
          do k = 1, size(g%GW)
            g%GW(k)%tt = T + 0.0D0
          enddo
          g%Wea%t = T + 0.0D0
          !g%ET%t = T + 0.0D0
          g%Veg%t = T + 0.0D0
      enddo
      do i = 1, size(w%CR)
          c = w%CR(i)
          r(c)%Riv%t = T + 0.0D0
      enddo
      
      end subroutine wSyncTime

      subroutine solversSetDT(dT)
      !% Synchronize every component of a watershed (and all its subwatersheds)
      !% to be at time T
      use Vdata
      implicit none
      integer :: T, i, k, c
      real*8 :: gid, dt
      
      w%dt = dT
      do i = 1, int(w%g)
          gid = w%g
          g%t = dT
          g%OVN%dt = dT
          g%VDZ%dtP = (/dT/3D0, dT, dT /)
          do k = 1, size(g%GW)
            g%GW(k)%dt = dT + 0.0D0
          enddo
          !g%Wea%dt = T + 0.0D0
          !g%ET%t = T + 0.0D0
          !g%Veg%dt = T + 0.0D0
      enddo
      do i = 1, size(w%CR)
          c = w%CR(i)
          r(c)%Riv%dt = dT + 0.0D0
      enddo
      
      end subroutine solversSetDT
