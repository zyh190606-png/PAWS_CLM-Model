!

      subroutine F_oc_dw(j)
      !% (E) Diffusive Wave Approximation to Overland-Channel Exchange-
      !% Used an explicit DW formulation, not self-similar. Maybe not good for-
      !% big dt. Also compute groundwater level for the channel cell

      !% #######################################################################
      !% May need to think of a new method: what does this flow in have to do with
      !% h depth, as a flow rate? its nonsense. The real factor is the distance to
      !% stream. we can compute a contribution fraction (similar to variable
      !% source area)
      !global g r w
      use Vdata
      implicit none
      integer,intent(in)::j
      integer :: i,k,nr,nc,n,jj,kf
      integer :: rid, k1, k2
      real*8 :: DTE, GRAV, dx, dt_l, dt_r, Ar, hhh, m_err, t1, ff
      integer,dimension(:),pointer::idx=>NULL(),Rep=>NULL(),Cidx=>NULL()
      real*8,dimension(:),POINTER::Len=>NULL(),ZBank=>NULL(),Zch=>NULL()
      real*8,dimension(:),POINTER::h0=>NULL(),Qoc1D=>NULL()
      real*8,dimension(:),pointer::gOVNh1D=>NULL(),gTopoE1D=>NULL(),gOVNMann1D=>NULL(),gGWh1D=>NULL()
      real*8,dimension(:),allocatable:: Mann
      real*8,dimension(size(g%Riv%Lidx(j)%ptr)) :: ho, Zo
      real*8,dimension(size(g%Riv%Stages(j)%ptr)) :: Moc, EstE
      real*8,dimension(size(g%Riv%Len(j)%ptr)) :: Ab
      real*8,dimension(size(g%Riv%TM(j)%ptr,1)) :: Qr
      real*8,dimension(size(r(g%Riv%ID(j))%Riv%h)) :: hg, tempAc

      real*8,dimension(:,:),pointer::M=>NULL()
      
      GRAV = 9.81D0
      rid = g%Riv%ID(j)
      Len => g%Riv%Len(j)%ptr
      ![nr,nc]=size(Len);
      nr = size(Len,1)
      nc = 1    !nc = size(Len,2)
      idx => g%Riv%Lidx(j)%ptr  !% idx is 1d (column) array
      ZBank => g%Riv%ZBank(j)%ptr
      Zch => g%Riv%Stages(j)%ptr    !% Zch should also be 1d (column) array
      !ho = max(g(gid).OVN.h(idx),0);
      gOVNh1D => OneD_ptr(size(g%OVN%h),g%OVN%h)
      Qoc1D   => OneD_ptr(size(g%OVN%Qoc),g%OVN%Qoc)
      ho = gOVNh1D(idx)
      where (ho<0.0D0)
        ho = 0.0D0
      end where
      gTopoE1D => OneD_ptr(size(g%Topo%E),g%Topo%E)
      Zo = ho + gTopoE1D(idx)

      !% debug
      if (debugg) then
      hhh = sum(g%OVN%h)/size(g%OVN%h)+sum(r(rid)%Riv%Ac(2:r(rid)%DM%msize(1)-1))*r(rid)%DM%d(1)/w%Stats%TA;
      endif

      !% These may be wrong:
      !if numel(g(gid).OVN.Mann)==1, Mann=repmat(g(gid).OVN.Mann,1,length(Zch));
      if (size(g%OVN%Mann)==1) then
        if(allocated(Mann)) deallocate(Mann)
        allocate(Mann(size(g%OVN%Mann,1)*size(g%OVN%Mann,2)*size(Zch)))
        n = 0
        do k = 1, size(Zch)
          do jj = 1, size(g%OVN%Mann,2)
            do i = 1, size(g%OVN%Mann,1)
              n = n + 1
              Mann(n) = g%OVN%Mann(i,jj)
            enddo
          enddo
        enddo
      else
        if(allocated(Mann)) deallocate(Mann)
        allocate(Mann(size(idx)))
        gOVNMann1D => OneD_ptr(size(g%OVN%Mann),g%OVN%Mann)
        Mann = gOVNMann1D(idx)
      endif
      !%S0=g(gid).Topo.S0;
      dx = g%DM%d(2)
      !% Here, S0 stands for the mean slope in this cell. This can be
      !% computed from integrating elevation and divide by area (obtain
      !% height) and then divide by the length to the channel. Assuming
      !% the lowest point is where the channel is.
      !% F_oc(n,ZBank,Zch,Zo,ho,Ab,L,dx,dt,Mann,M)
      dt_l = g%OVN%dt
      dt_r = r(rid)%Riv%dt
      DTE = max(dt_l,dt_r)
      !Moc = zeros(size(Zch));
      Moc = 0.0D0
      Ab = g%Riv%W(j)%ptr * Len
      !F_oc(ZBank,Zch,Zo,ho,Ab,Len,dx,dt_l*86400,Mann,Moc)
      call F_oc(size(ZBank),ZBank,Zch,Zo,ho,Ab,Len,dx,dt_l*86400,Mann,Moc)
      Ar = dx**2


      !% debug
      !%g.OVN.Qoc(idx) = g.OVN.Qoc(idx)+Moc/Ar;
      !% debug over

      !% THIS WILL NOT GET BACK Moc AS SOME idx ARE REPETITIVE!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!
      !Rep => g%Riv%RepLidx(j)%ptr
      !if (Rep(1)/=0) then
      !    do i = 1, size(Rep)
      !        k1 = Rep(i)
      !        k2 = idx(k1)
      !        !g%OVN%h(k2) = g%OVN%h(k2) - Moc(k1) / Ar
      !        gOVNh1D(k2) = gOVNh1D(k2) - Moc(k1) / Ar
      !        !%g.OVN.Qoc(k2)=g.OVN.Qoc(k2)+Moc(k1)/Ar   !% debug
      !        Qoc1D(k2) = Qoc1D(k2)+Moc(k1)/Ar;
      !    enddo
      !endif
      !!!!!!!!!!!!!!!!!!!!
      !%if any(g.OVN.h<-1e-1)
      !%    4;
      !%end
      !%h1 = g.OVN.h(:);
      !% Add the mass to river lateral inflow
      M => g%Riv%TM(j)%ptr
      Qr = matmul(M, Moc)    !%Ab = r(rid).Riv.Ab(Cidx);
      !%DDT = max(DTE,dt_r); % Possible Upscaling of flux.
      Cidx => g%Riv%Cidx(j)%ptr

      !%r(rid).Riv.Qoc(Cidx)=r(rid).Riv.Qoc(Cidx)+Qr'./(r(rid).Riv.dt*86400); % ASSUME Cidx SHALL NOT REPEAT
      !% #################### BE VERY CAREFUL ABOUT THIS FLUX TERM TIME STEP HERE
      !% IT SHOULD BE Qoc*(THE TIME THIS FLUX WILL BE APPLIED)=MASS REMOVED FROM LAND
      !%}
      !% Previously converted to a source term, it is can cause instability
      !% (together with tributaries. now directly take the water away
!      if (any(r(rid)%Riv%Ac<0D0)) then
!        call display('negative Ac before Qr')
!      endif
      tempAc = 0D0; kf = 0
      tempAc(Cidx) = r(rid)%Riv%Ac(Cidx) + Qr / r(rid)%Riv%DX
      do i=1,size(Cidx)
         if (tempAc(Cidx(i))<0D0) then
             kf = kf + 1
             ! This really happens very rarely
             ff = -r(rid)%Riv%Ac(Cidx(i))/(tempAc(Cidx(i))-r(rid)%Riv%Ac(Cidx(i)))
             do jj=1,size(M,2)
                 if (M(i,jj) .ne. 0D0) then
                     Moc(jj) = Moc(jj)*(ff*0.97D0)
                 endif
             enddo
         endif
      enddo
      if (kf > 0) Qr = matmul(M, Moc)
      tempAc(Cidx) = r(rid)%Riv%Ac(Cidx) + Qr / r(rid)%Riv%DX
      r(rid)%Riv%Ac = tempAc
      r(rid)%Riv%vrv(10,Cidx) = r(rid)%Riv%vrv(10,Cidx) + Qr ! m3
      
      
      !% Deduct the mass from OVN
      if(.not. associated(g%Riv%Qoc(j)%ptr)) allocate(g%Riv%Qoc(j)%ptr(size(Moc)))
      g%Riv%Qoc(j)%ptr = Moc / DTE   !% Save as a rates #####################
      !% positive as flow from land to riv
      EstE = Zch + Moc / Ab    !% Estimated Stage (due to overland alone)
      !h0 = g.OVN.h(:);
      h0 => OneD_Ptr(size(g%OVN%h), g%OVN%h)
      !g%OVN%h(idx) = g%OVN%h(idx) - Moc / Ar
      ! Above is the vectorized statement, but when there are repeated indices in idx, only one of the repeated
      !values will be assigned, leading to MBE. Therefore we do the RepLidx fix below. However, different compilers may
      !have different interpretations of what to do under this scenario. In addition, this performance gain is not
      !significant in Fortran. Therefore we do this loop, and remove the RepLidx fix
      ! In GW_source??
      do i = 1,size(Moc)
      !gOVNh1D(idx) = gOVNh1D(idx) - Moc / Ar
      !Qoc1D(idx) = Qoc1D(idx)+Moc/Ar;
         t1=Moc(i) / Ar
         gOVNh1D(idx(i) ) = gOVNh1D(idx(i) ) - t1
         Qoc1D(idx(i) ) = Qoc1D(idx(i) )+t1;
      enddo
      
!      if (any(r(rid)%Riv%Ac<0D0)) then
!        call display('negative Ac after Qr')
!      endif
      r(rid)%Riv%h(Cidx) = r(rid)%Riv%Ac(Cidx) / r(rid)%Riv%w(Cidx)
      w%Rec%acc%Qoc(j) = w%Rec%acc%Qoc(j) + sum(Qr)

      r(rid)%Riv%EstE(Cidx) = matmul(M, EstE)
      !%test = mean(h0-h1)-sum(r(rid).Riv.Qoc)*(g(gid).OVN.dt*86400)/w.Stats.TA;
      !% test2= mean(h0)-mean(h1)-sum(r(rid).Riv.Qoc)*(g(gid).OVN.dt*86400)/w.Stats.TA;
      !% if abs(test)>1e-12 || abs(test2)>1e-8
      !%     555
      !% end
      !% r(rid).DM.d convert the unit into m3/s/m. But notice dx may be set to
      !% variable
      !hg = zeros(size(r(rid).Riv.h));
      hg = 0.0D0
      gGWh1D => OneD_ptr(size(g%GW(1)%h),g%GW(1)%h)
      hg(Cidx) = matmul(gGWh1D(idx), g%Riv%TN(j)%ptr)    !% groundwater elevation
      hg(1) = hg(2)
      hg(size(hg)) = hg(size(hg)-1)
      r(rid)%Riv%hg = hg

      !% debug
      if (debugg) then
      m_err = sum(g%OVN%h)/size(g%OVN%h)+sum(r(rid)%Riv%Ac(2:r(rid)%DM%msize(1)-1))*r(rid)%DM%d(1)/w%Stats%TA -  hhh;
      if (abs(m_err)>1D-14) then
          call display('F_oc_dw:: MBE')
      endif
      !% also prcp
      !%{
      !if length(g(gid).OVN.prcp)>1
      !    r(rid).Riv.prcp(Cidx)= r(rid).Riv.prcp(Cidx)+(M*g(gid).OVN.prcp(idx))'*(dt_l/DDT/86400);
      !else
      !    r(rid).Riv.prcp(Cidx)= r(rid).Riv.prcp(Cidx)+g(gid).OVN.prcp*(dt_l/DDT/86400);
      !end
      !%} 
      endif
      !%commented out on May 01
      
      if(allocated(Mann)) deallocate(Mann)
      
      end subroutine F_oc_dw
