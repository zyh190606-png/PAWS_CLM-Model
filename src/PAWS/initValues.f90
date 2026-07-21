
!
      subroutine initValues()
      !  % (E) These are some values initiated after parameter change has been applied
      !  % This is to avoid values being modified by changepar
      !  % Applying multiplier parameters set in adhoc_set (and modified in
      !  % changepar)
      !  % VDZ initial values: heads, theta, WP, FC, they are determined from here
      !  % because equilibrium profile THE is a function of K
      !  % This equilibrium is an initial guess of soil water content to quickly
      !  % bring the model into reality
      !  % Assume hydrostatic equilibrium and evaluate THE and h
      !  % This assumption may not be valid in arid regions, we need to make
      !  % modifications if necessary, and this is only an initial guess
      !  % also modifying River E and K previously
      use Vdata
      use gd_Mod
      implicit none
      type(VDZ_type) VVDZ
      type(sta_type),pointer:: day
      real*8,dimension(maxval(w%CR))::rs2

      real*8,dimension(size(g%GW(1)%h,1),size(g%GW(1)%h,2))::WTE
      real*8,dimension(size(g%GW(1)%EB,1),size(g%GW(1)%EB,2))::EB, TempH
      real*8,dimension(:,:,:),pointer::h=>NULL(), RR=>NULL(), dp=>NULL()
      !real*8,dimension(size(g%VDZ%h,1),size(g%VDZ%h,2),size(g%VDZ%h,3))::h
      real*8,dimension(:,:,:),pointer::S=>NULL()
      real*8,dimension(size(g%GW(1)%h,1),size(g%GW(1)%h,2),10)::NLC
      !real*8,dimension(size(g%VDZ%h,1),size(g%VDZ%h,2),size(g%VDZ%h,3))::RR,dP
      !real*8,dimension(size(g%VDZ%REP(:,:,:),1),size(g%VDZ%REP(:,:,:),2),size(g%VDZ%REP(:,:,:),3))::RR
      real*8,dimension(size(g%GW(1)%h,1),size(g%GW(1)%h,2))::AWW
      real*8,dimension(size(g%VDZ%E(1,1,:)))::E
      real*8,dimension(size(E))::f
      real*8,dimension(size(g%GW(1)%K,1),size(g%GW(1)%K,2))::KK
      real*8,dimension(size(g%GW(1)%botK))::botK
      real*8,dimension(size(g%GW(1)%EB,1),size(g%GW(1)%EB,2))::DD
      real*8,dimension(size(DD))::DD1D
      real*8,dimension(:),pointer::control=>null()
      !real*8,dimension(size(g%OVN%dP,1),size(g%OVN%dP,2),size(g%OVN%dP,3))::dP

      integer :: nr,nc,nz,i,j,k,m,n,icase, ii, rid
      integer :: ns,ny,nt
      integer,dimension(:),allocatable :: un
      real*8 :: minST,minS,mp,vp, tt = 0D0, ttm = 0D0
      logical :: initializeTHE = .true.
 
 
      !write(777,*) 'initvalues here'
      VVDZ = g%VDZ
      WTE = g%GW(1)%h
      EB = g%GW(1)%EB
      h => VVDZ%h
      allocate(S(size(h,1),size(h,2),size(h,3)))
      
      ![nr,nc,nz]=size(h);
      nr = size(h,1)
      nc = size(h,2)
      nz = size(h,3)
      !h = repmat(WTE,[1,1,nz])-VDZ.E;
      
      if (size(w%control)>1) then
        if (w%control(2)>0) then
          initializeTHE = .false.
        endif
      endif

      if (initializeTHE) then
          do i = 1, nz
            h(:,:,i) = WTE - VVDZ%E(:,:,i)
          end do
      endif

      !un=find(WTE<EB);
      n = 0
      do j = 1, size(WTE,2)
        do i = 1, size(WTE,1)
          if(WTE(i,j)<EB(i,j)) then
            n = n + 1
          end if
        end do
      end do
      allocate(un(n))
      n = 0
      m = 0
      do j = 1, size(WTE,2)
        do i = 1, size(WTE,1)
          n = n + 1
          if(WTE(i,j)<EB(i,j)) then
            m = m + 1
            un(m) = n
          end if
        end do
      end do
      h(:,:,1) = 0
      ![h,THE,S,K]=vanGenuchten(h,VDZ.THES,VDZ.THER,VDZ.ALPHA,VDZ.LAMBDA,VDZ.N,VDZ.KS,0);
      !g%VDZ%THE = THE
      !g%VDZ%K = KKK
      !g%VDZ%h = h
      call vanGenuchten(h,VVDZ%THES,VVDZ%THER,VVDZ%ALPHA,VVDZ%LAMBDA,VVDZ%N,VVDZ%KS,int(0),g%VDZ%THE,S,g%VDZ%K,g%VDZ%h)

      ![h,g(gid).VDZ.WP] = vanGenuchten(-15*10.2,VDZ.THES,VDZ.THER,VDZ.ALPHA,VDZ.LAMBDA,VDZ.N,VDZ.KS,0);
      call vanGenuchten(-15.0D0*10.2D0,VVDZ%THES,VVDZ%THER,VVDZ%ALPHA,VVDZ%LAMBDA,VVDZ%N,VVDZ%KS,int(0),g%VDZ%WP)
      ![h,g(gid).VDZ.FC] = vanGenuchten(-0.33*10.2,VDZ.THES,VDZ.THER,VDZ.ALPHA,VDZ.LAMBDA,VDZ.N,VDZ.KS,0);
      call vanGenuchten(-0.33D0*10.2D0,VVDZ%THES,VVDZ%THER,VVDZ%ALPHA,VVDZ%LAMBDA,VVDZ%N,VVDZ%KS,int(0),g%VDZ%FC)
      !g%VDZ%WP = max(cat(4,g.VDZ.WP,g.VDZ.THER),[],4);
      do k = 1, size(g%VDZ%WP,3)
        do j = 1, size(g%VDZ%WP,2)
          do i = 1, size(g%VDZ%WP,1)
            if(g%VDZ%WP(i,j,k)<g%VDZ%THER(i,j,k)) then
              g%VDZ%WP(i,j,k)=g%VDZ%THER(i,j,k)
            end if
          end do
        end do
      end do

      !% use the dominant soil property to specific dynamic storage parameters
      !NLC = zeros(nr,nc,10); % nonlinear coefficient for the drainable porosity
      NLC = 0.0D0
      g%VEG%hI = 0D0

      !if isfield(g.VDZ,'REP')
      if(associated(g%VDZ%REP)) then
        RR => g%VDZ%REP(:,:,:)
        !% ![1ES,2THES-THER,3ALPHA,4nn,5:-1/nn,6:-(1+1/nn),7:Kl,8DZL,9HL,10HB]
        NLC(:,:,1)=g%Topo%E
        NLC(:,:,2)=(RR(:,:,2)-RR(:,:,3))*g%VDZ%VParm(2)
        !%REP(i,j,:)=[k,THES THER ALPHA N Lambda KS Ko];
        NLC(:,:,3)=RR(:,:,4)
        NLC(:,:,4)=RR(:,:,5)
        NLC(:,:,5)=-1.0D0/RR(:,:,5)
        NLC(:,:,6)=-(1.0D0+1.0D0/RR(:,:,5))
        NLC(:,:,7)=g%GW(1)%botK
        NLC(:,:,8)=g%GW(1)%E-g%GW(2)%E
        NLC(:,:,9)=0.0D0    ! % NEED TO BE UPDATED
        NLC(:,:,10)=g%GW(1)%EB
        g%GW(1)%NLC = NLC
        !g.GW(1).STO = zeros(nr,nc);
        g%GW(1)%STO = 0.0D0
      end if

      !% groundwater unconfined aquifer specific yield
      !%
      !AWW = zeros(nr,nc);
      AWW = 0.0D0
      do j = 1, nc
        do i = 1, nr
          k = g%VDZ%NZC(i,j)-1_8
          AWW(i,j)=sum((g%VDZ%THES(i,j,2:k)-g%VDZ%FC(i,j,2:k))*g%VDZ%DZ(i,j,2:k))/sum(g%VDZ%DZ(i,j,2:k))
        end do
      end do
      !g%GW(1)%ST = AWW*g%VDZ%VParm(2)
      minST = 5D-4
      if (maxval(g%GW(1)%ST)>100) then
         call display( 'water body BC detected, max ST:'//num2str(maxval(g%GW(1)%ST)) )
      endif
      
      WHERE (g%GW(1)%ST < 100)
         g%GW(1)%ST = minST ! Now uniformly set ST = 5e-4, because this doesn't matters, it is a transient conduit
         ! explanation: if ST>100, this means it is a water body boundary condition (close to fixed head)
      ENDWHERE
      g%Veg%VDB%ho = g%Veg%VDB%ho * g%VDZ%VParm(5)
      
      call subsurfaceExtrapolate()
      !%}
      !% soil conductivity decay with depth: not necessarily useful
      !if (size(g%VDZ%VParm)>=7 .and. g%VDZ%VParm(7)/=0) then
      !  vp = g%VDZ%VParm(7)
      !  mp = g%VDZ%VParm(8)
      !  do j = 1, nc
      !    do i = 1, nr
      !          E = g%VDZ%E(i,j,1)-g%VDZ%E(i,j,:)
      !          f = dexp(-E*vp)
      !          !f(f<mp)=mp; f(f>10)=10; % should give them a limit
      !          where(f<mp)
      !            f=mp
      !          end where
      !          where(f>10)
      !            f=10
      !          end where
      !          g%VDZ%KS(i,j,:)=g%VDZ%KSSAVE(i,j,:)*f
      !          
      !    end do
      !  end do
      !end if
      !% we shall do the head derivative calculation at the bottom of the aquitard
      !%K = sqrt(g.GW(1).K.*g.GW(2).K);
      !KK = (g%GW(1)%K * g%VDZ%KS(:,:,10) * g%GW(2)%K)**(1.0D0/3.0D0)
      !g%GW(1)%botK = KK * g%VDZ%VParm(9)
      ! changed 2013/Mar/04: now use harmonic mean
      KK =  (1D0/g%GW(1)%K + 1D0/g%VDZ%KS(:,:,10) + 1D0/g%GW(2)%K)
      WHERE (g%GW(1)%K .eq. 0D0 .OR. g%GW(2)%K .eq. 0D0 .OR. g%VDZ%KS(:,:,10) .eq. 0D0)
        g%GW(1)%botK = 0D0;
      ELSEWHERE
        g%GW(1)%botK = g%VDZ%VParm(9)/ KK
      ENDWHERE
      if (maxval(g%GW(1)%botK) > 0.0005) then
        g%GW(2)%dt = g%GW(1)%dt
        call display( 'Decreased second layer GW dt to the same as first layer')
      end if

      botK = oneD_Ptr(size(g%GW(1)%botK),g%GW(1)%botK)
      do i =1, size(g%GW(1)%Cond)
        g%GW(1)%Cond(i)%K = botK(g%GW(1)%Cond(i)%idx)       !% very small K for the aquitard
      enddo
      DD = 0.5D0*(g%GW(1)%EB+g%GW(1)%ES) - g%GW(2)%E
      DD1D = oneD_Ptr(size(DD),DD)
      do i = 1, size(g%GW(1)%Cond)
        g%GW(1)%Cond(i)%D = DD1D(g%GW(1)%Cond(i)%idx)     !% thickness to lower layer
      enddo
      g%GW(2)%ST = g%GW(2)%ST * g%GW(2)%D           !% Specific storage * thickness = storage coefficient!!!

      !% slope rescaleing
      !%g.Topo.S0 = g.Topo.S0*(g.Topo.S0p/meanALL(g.Topo.S0(3:end-2,3:end-2)));
      !for i=w.CR, rs2(i)=-(r(i).Riv.E(end)-r(i).Riv.E(1))/r(i).Riv.DDist(end); end
      do i = 1, size(w%CR)
        ii = w%CR(i)
        rs2(ii)=-(r(ii)%Riv%E(size(r(ii)%Riv%E))-r(ii)%Riv%E(1))/r(ii)%Riv%DDist(size(r(ii)%Riv%DDist))
      end do
      minS = maxval(rs2)        !% minimum microslope for outflow: maximum river slope (fractal theory??)
      !g.Topo.S0(g.Topo.S0<minS)=minS;
      where (g%Topo%S0 < minS)
        g%Topo%S0 = minS
      end where
      
      ttm = 0D0
      if (ttm .eq. 0) then ! this is just to create a debugging hook. normally we should always do this
        do j = 1, nc
          do i = 1, nr
                tt = 0D0;  ! flow depth
                ttm = 0D0  ! MANNING
                ! Hard-wiring here for the moment
                ! Currently, hard-wiring Manning's cofficient from TR-55
                ! no matter for what unit, N is always the same (1/m^(1/3)), just add 1.486 on top of it if using US Customary system
                do k = 1,size(g%Veg%RPT,3)
                  if (g%Veg%RPT(I,J,K) .eq. g%Veg%LUidx(1)) then 
                    ! is urban
                    ! BECAUSE URBAN is handled separately, we should not average its N here
                    !tt  = tt  + g.Veg.Frac(I,J,K) * g.VDZ.VParm(10)
                    !ttm = ttm + g.Veg.Frac(I,J,K) * 0.011D0
                  elseif (g%Veg%RPT(I,J,K) > g%Veg%LUidx(1)) then 
                    ! water or bare
                    tt  = tt  + g%Veg%Frac(I,J,K) * g%VDZ%VParm(10)
                    ttm = ttm + g%Veg%Frac(I,J,K) * 0.025D0 
                  elseif (g%Veg%RPT(I,J,K) <= 8) then
                    ! trees
                    tt  = tt  + g%Veg%Frac(I,J,K) * g%VDZ%VParm(10)
                    ttm = ttm + g%Veg%Frac(I,J,K) * 0.7D0  ! between light under brush and dense brush
                  else
                    ! grass and bush
                    tt  = tt  + g%Veg%Frac(I,J,K) * g%VDZ%VParm(10)
                    ttm = ttm + g%Veg%Frac(I,J,K) * 0.24D0 ! dense grass
                  endif
                enddo
                IF (g%Veg%Imp(I,J)<1D0) then
                  g%OVN%Mann(I,J) = ttm/(1D0 - g%Veg%Imp(I,J))
                ELSE
                  g%OVN%Mann(I,J) = 0.011D0
                ENDIF 
                
                if ((isnan(g%OVN%Mann(I,J)) .or. g%OVN%Mann(I,J) < 1e-12) .and. (.not. w%DM%Mask(I,J)) ) then
                  ! already outside the domain. just give it a value
                  g%OVN%Mann(I,J) = 0.5D0
                endif
          enddo
        enddo
      endif
                
      g%OVN%SN = (86400.0D0*g%Topo%S0**0.5D0)/(g%OVN%Mann*g%OVN%dist)
      g%OVN%SNSAVE = g%OVN%SN + 0.0D0
      
      ! a new change: ffrac is reserved exclusively for urban, and therefore should be using urban Mann
      KK = (86400.0D0*g%Topo%S0**0.5D0)/(0.011D0*g%OVN%dist)*(g%VDZ%VParm(10)**(2.0D0/3.0D0))
      !if(.not. associated(g%OVN%ffrac)) allocate(g%OVN%ffrac(size(KK,1),size(KK,2)))
      !forall(i=1:size(KK,1),j=1:size(KK,2)) g%OVN%ffrac(i,j) = dexp(-KK(i,j)*g%VDZ%dtP(3))
      g%OVN%ffrac = dexp(-KK*g%VDZ%dtP(3))        !% fraction that stayed (not runoff)

      !% lowland depression storage characterstics
      !% ! [1.hback, 2.hod, 3.Kdg, 4.dfrac, 5.gwST, 6.mr, 7.g1, 8.g2, 9.g3, 10.EB, 11.E]
      !% ! g1 = K*dt/mr, g2=g1*Ad/(Ac*S), g3=g1/(1+g2)
      dP => g%OVN%dP
      dP(:,:,1) = g%OVN%hback
      dP(:,:,2) = 0.002D0     !%g.OVN.ho;
      !%dP(:,:,2)=g.OVN.dP(:,:,2); % dP(2) is updated in dailyVegUpdate
      !%dP(:,:,4); original dP(:,:,4) is dfrac
      KK = dsqrt(g%GW(1)%K*g%VDZ%GA%KBar(:,:,10))*(1.0D0-g%Veg%Imp)      !% impervious
      KK = KK/g%OVN%FParm(2)
      dP(:,:,3) = KK
      dP(:,:,5) = g%GW(1)%ST
      dP(:,:,11) = max(dP(:,:,11), g%GW(1)%EB+0.4D0)
      dP(:,:,11) = max(dP(:,:,11), min(g%GW(1)%h-0.5D0,g%Topo%E))
      if (any(dP(:,:,11)<g%GW(1)%EB)) then
        call display( 'negative mr may be generated')
      end if
      dP(:,:,6) = min(2.0D0,(dP(:,:,11)-g%GW(1)%EB)/2.0D0)
      dP(:,:,7) = KK*g%VDZ%dtP(3)/dP(:,:,6)
      dP(:,:,8) = dP(:,:,7)*g%OVN%dP(:,:,4)/g%GW(1)%ST
      icase = 1
      if (icase==1) then
        dP(:,:,9) = dP(:,:,7)/(1.0D0+dP(:,:,8))
      else
        dP(:,:,9) = 0.0D0
      end if
      dP(:,:,10) = g%GW(1)%EB
      !%dP(:,:,11)=g.GW(1).Ed;
      !g%OVN%dP = dP

      call rivBedK()
      !%g.OVN.SN(:)=0; % if to use F_contribute
      do i = 1, size(w%CR)
      ![1. hg_accu (m); 2. Qgc_n (m/day); 3. gamma (-); 4. beta (m/s) ; 5 .ST; 6. zeta=beta*dt/mr/(1+K*dt/mr)] -- see notes
        rid = w%CR(i)
        r(rid)%Riv%E = r(rid)%Riv%E+g%OVN%FParm(1)      !% River bed elevation greatly affect baseflow -5+1=-4
        r(rid)%Riv%K = r(rid)%Riv%K/g%OVN%FParm(3)
        r(rid)%Riv%vrv(1,:)=0D0    ! H^*
        r(rid)%Riv%vrv(2,:)=0D0    ! q_{gc}^n
        r(rid)%Riv%vrv(3,:)=r(rid)%Riv%w*r(rid)%DM%d(1)/(g%DM%d(1)*g%DM%d(2)) ! gamma, storage factor (<< 1)
        r(rid)%Riv%vrv(4,:)=r(rid)%Riv%vrv(3,:)*r(rid)%Riv%K/86400D0 ! [m/s] used inside riv_cv
        r(rid)%Riv%vrv(5,:)=minST    ! storage
        r(rid)%Riv%vrv(6,:)=0D0 ! zeta
        g%Riv%ZBank(i)%ptr = g%Riv%ZBank(i)%ptr + g%OVN%FParm(1)
        !%g.Riv.ZBank{i}=g.Riv.ZBank{i}+1; % no need to do this two
        !%anymore!!!!!!!!!!!!!!
      end do
      
      


      end subroutine initValues

!========================================================================================
      subroutine rivBedK()
      use Vdata
      implicit none
      integer :: gid, n, i, nrr, ncc, rid
      integer, dimension(:), allocatable :: idx, Cidx
      real*8, dimension(:), allocatable :: KK1, KK2, KK, EH
      real*8, dimension(size(g%GW(1)%K))::K1D
      real*8, dimension(size(g%GW(1)%h))::h1D
      real*8, dimension(size(g%VDZ%GA%KBar(:,:,10)))::KV
      real*8,dimension(:),pointer:: mask1D


      !global g r
      !% Groundwater Interaction with river object if applicable
      !if isfield(g(gid),'Riv') && ~isempty(g(gid).Riv) && ~isempty(g(gid).Riv.ID)
      if (associated(g%Riv%ID)) then
          n = size(g%Riv%ID)
          K1D = oneD_ptr(size(g%GW(1)%K),g%GW(1)%K)
          KV = oneD_Ptr(size(g%VDZ%GA%KBar(:,:,10)),g%VDZ%GA%KBar(:,:,10))
          h1D = oneD_Ptr(size(g%GW(1)%h),g%GW(1)%h)
          mask1D => oneD_Ptr(size(w%DM%mask2),w%DM%Map)
          do i = 1, n
              !idx= g(gid).Riv.Lidx{i}; Cidx = g(gid).Riv.Cidx{i};
              if(allocated(idx)) deallocate(idx)
              allocate(idx(size(g%Riv%Lidx(i)%ptr)))
              idx = g%Riv%Lidx(i)%ptr
              if(allocated(Cidx)) deallocate(Cidx)
              allocate(Cidx(size(g%Riv%Cidx(i)%ptr)))
              Cidx = g%Riv%Cidx(i)%ptr
              rid = g%Riv%ID(i)
              if(allocated(KK1)) deallocate(KK1)
              allocate(KK1(size(idx)))
              KK1 = K1D(idx)
              ![nrr,ncc]=size(KK1);
              if(allocated(KK2)) deallocate(KK2)
              allocate(KK2(size(idx)))
              KK2 = KV(idx)
              if(allocated(KK)) deallocate(KK)
              allocate(KK(size(idx)))
              where (mask1D(idx) .eq. 0) KK1 = 0D0
              KK = dsqrt(KK1*KK2)       !% average of vertical and lateral
              !if nrr>1, KK=KK'; end
              !r(rid).Riv.K(Cidx) = KK*g(gid).Riv.TN{i};
              r(rid)%Riv%K(Cidx) = matmul(KK, g%Riv%TN(i)%ptr)
              !EH = (g(gid).GW(1).h(idx))'*g(gid).Riv.TN{i};
              if(allocated(EH)) deallocate(EH)
              allocate(EH(size(g%Riv%TN(i)%ptr,2)))
              EH = matmul(h1D(idx), g%Riv%TN(i)%ptr)
              r(rid)%Riv%mr(Cidx) = max(r(rid)%Riv%E(Cidx) - EH, 2.0D0);
              !r(rid).Riv.mr(r(rid).Riv.mr==0)=2;
              where(r(rid)%Riv%mr==0) r(rid)%Riv%mr=2.0D0
          enddo
      endif

      if(allocated(idx)) deallocate(idx)
      if(allocated(Cidx)) deallocate(Cidx)
      if(allocated(KK1)) deallocate(KK1)
      if(allocated(KK2)) deallocate(KK2)
      if(allocated(KK)) deallocate(KK)
      if(allocated(EH)) deallocate(EH)
      
    end subroutine rivBedK

!========================================================================================
    
      subroutine subsurfaceExtrapolate()
      use Vdata
      use matTable
      implicit none
      integer :: i,j,k,nr,nc,nz,ioe
      integer :: kStrat, vStrat, nStrat
      real*8 ::  kParam, vParam, nParam
      real*8 t1, t2, t3,f
      real*8,dimension(size(g%VDZ%h,3)):: r1, r2, dl, du, fact
      real*8,dimension(:),pointer::z
      real*8,dimension(:,:,:),pointer:: h3
      real*8,dimension(:,:),pointer:: lastZ
      TYPE(gd), pointer:: gdVDZ
      character*100, dump
      
      h3 => g%VDZ%h
      nr = size(h3,1)
      nc = size(h3,2)
      nz = size(h3,3)
      ! process extrapolate options.
      ! currently, let it be a light-weight program.
      if (size(g%VDZ%VParm)>=18 .AND. g%VDZ%VParm(18)>0) then
        kStrat = 1; kParam = g%VDZ%VParm(18)>0
      else      
        open(unit=17,file='subsurfaceExtrapolate.txt',ACTION='READ',status='unknown',iostat=ioe)
        if (ioe /= 0) then
          kStrat = 0; kParam = 0
        else
          read(17,*) dump, kStrat, kParam
        endif
      close(17)
      endif
      
      gdVDZ => gdG .G. 'VDZ'
      if (.not. isempty(gdVDZ,'lastZ')) then
          call getPtr(gdVDZ, 'lastZ', lastZ)
      else
          allocate(lastZ(size(h3,1),size(h3,2)) )
          lastZ = 5D0*0.3048D0
      endif    
      ! apply modifications
      do i=1,nr
          do j=1,nc
              if (w%DM%mask(i,j)) then
                  z => g%VDZ%E(i,j,:)
                  
                  if (kStrat .eq. 1) then  ! e-folding
                      dl = g%Topo%E(i,j) - lastZ(i,j) - g%VDZ%EB(i,j,:)
                      du = dl - g%VDZ%DZ(i,j,:)
                      f = kParam/(1D0+150D0*g%Topo%S0(i,j))  ! Ying 2013 Science Supporting Materials kParam==100
                      fact = f*(exp(-du/f)-exp(-dl/f))/(dl-du);
                      do k=1,nz
                          if (du(k)>=0D0) then
                              t1 = f*(exp(-du(k)/f)-exp(-dl(k)/f))/(dl(k)-du(k)); ! average K in layer
                              g%VDZ%Ko(i,j,k)=t1*g%VDZ%Ko(i,j,k)
                              g%VDZ%KS(i,j,k)=t1*g%VDZ%KS(i,j,k)
                          elseif (dl(k)>0D0) then
                              t1 = f*(exp(-0D0/f)-exp(-dl(k)/f))/dl(k);
                              g%VDZ%Ko(i,j,k)=t1*g%VDZ%Ko(i,j,k)
                              g%VDZ%KS(i,j,k)=t1*g%VDZ%KS(i,j,k)
                          endif
                      enddo
                  endif
              endif
          enddo
      enddo
      
      
      end subroutine subsurfaceExtrapolate