!

      subroutine rMinDtSet(RL,RDT,dtNew,dtOld)
      !% (E) Check for maximum velocity and COU2
      !% The rational is that extreme cases happen when flow is very high (This-
      !% may be wrong) and thus velocity is high. 
      !% The most correct way is to identify if exchange term is too big, but-
      !% there is no good way to check for that prior to computation of river flow
      !% Two Criteria: (1). Umax (2) Qmax
      !% Input:
      !% RL: list of river IDs;-
      !% RDT: remaining dt to march in this river network loop step

      !global r w
      use Vdata
      implicit none
      integer :: RL(size(w%CR))
      real*8 :: RDT, dtNew, dtOld
      real*8 :: COU2, Qmax, COUMAX, mtimes
      integer i, n, rid
      
      COU2 = w%Tolerances%COU
      Qmax = w%Tolerances%Qmax   !% Tolerances
      COUMAX = 0.0D0
      do i = 1, size(RL)
        rid = RL(i)
        COUMAX= max(COUMAX,r(rid)%Riv%umax*r(rid)%Riv%dt*86400.0D0/r(rid)%DM%d(1))
      enddo

      mtimes = COUMAX / COU2
      !% Q = r(RL(end)).Riv.Qx(end-1);
      !% if Q>Qmax
      !%     mtimes = max(2,mtimes);
      !% %     w(1).Tolerances.Qmax = w(1).Tolerances.Qmax;
      !% end
      if (size(RL)>0) then
         dtOld = r(RL(size(RL)))%Riv%dt
         dtNew = dtOld
      else
         dtOld  = g%OVN%dt
         dtNew  = g%OVN%dt
      endif

      if (mtimes>1) then
        mtimes = max(mtimes*1.2D0, 2.0D0)
        dtNew = min(dtNew / mtimes, RDT)
        do i = 1, size(RL)
          r(RL(i))%Riv%dt = dtNew   !%r(RL(i)).Riv.dt/mtimes;
        enddo
      else
        if (dtNew > RDT) then
          dtNew = min(dtNew, RDT)
          do i = 1, size(RL)
            r(RL(i))%Riv%dt = dtNew
          enddo
        endif
      endif
      
      end subroutine rMinDtSet
