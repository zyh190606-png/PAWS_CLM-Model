!

      subroutine r_MacCormack(rid)
      !% (E) RKFV scheme for 1D diffusive wave equation
      !% For historical reasons this is named MacCormack, in fact this function-
      !% is not MacCormack scheme. Future revision should modify its name-
      !% for river flow, the upper and lower boundary conditions are stored in-
      !% r.UPST,UPSV,DST and DSV. They are wrapped in rLateral
      !% Input:  r(rid).Riv states
      !% Output: updated r(rid).Riv
      !global r
      use Vdata
      use flow
      implicit none
      integer rid
      REAL*8,DIMENSION(size(r(rid)%Riv%Qx))::Qx0    ! unused output in original MatLab code
      REAL*8,DIMENSION(:),POINTER:: hGacc ! hG is temporarily accrued until it is used by groundwater flow
      REAL*8 factor
      Integer,save:: status = 0
      integer i, j, nr, id, m, k, n1,n2
      integer :: nv = 9, opt = 0
      type localP
        real*8, dimension(:), pointer :: p  => null()
      end type
      type localRivs
        type(localP), dimension(:), pointer :: r  => null()  ! each element is for a river
      end type
      type(localRivs),dimension(:),pointer,save:: tempVars ! each element is a variable that is to be saved
      type(localRivs),pointer:: tempV ! each element is a variable
      real*8,DIMENSION(:),pointer :: LoadFact,Rbmem,Rcmem
      real*8,DIMENSION(:),pointer :: Rb,Rc,Hg,St,Phi,Fact
      integer,dimension(:),pointer:: Cidx
      Real*8 Ar, T

      type(rRiv_type), pointer :: Riv      
      !% Think of a better way for time step unit coordination
      r(rid)%Riv%dt = r(rid)%Riv%dt * 86400.0D0
      !%MacCormack(r(rid).Riv)
      !%{
      !fG_old = sum(r(rid).Riv.fG(:));
      !A_old  = sum(r(rid).Riv.Ac(:))*r(rid).Riv.DX;
      !%}
      !riv_cv(r(rid).Riv);
      Riv => r(rid)%Riv
         
      nr = size(w%CR)

      DO I=1,nr
        IF (w%CR(I) .eq. rid) then
          k = I
          EXIT
        endif
      ENDDO
      Ar = g%DM%d(1)*g%DM%d(2)
      
      !write(99,*) 'here'
!      if (status .eq. 0) then
!         allocate(tempVars(nv) )
!         DO J=1,nv
!           tempV => tempVars(J)
!           allocate(tempV%r(nr) )
!             DO I=1,nr
!               id = w%CR(I)
!               m =  size(r(id)%Riv%h)
!               allocate(tempV%r(I)%p(m) )
!               tempV%r(I)%p = 0D0
!             ENDDO
!         ENDDO
!
!         DO J=1,nr
!            id = w%CR(J)
!            Cidx => g%Riv%Cidx(j)%ptr            
!            
!            LoadFact => tempVars(1)%r(J)%p ! LoadingFactor is how much (fraction) land cell area the river is interacting with (due to discretization)
!            St => tempVars(2)%r(J)%p ! LoadingFactor is how much (fraction) land cell area the river is interacting with (due to discretization)
!            Fact => tempVars(3)%r(J)%p ! LoadingFactor is how much (fraction) land cell area the river is interacting with (due to discretization)
!            LoadFact = 1D0
!            DO I=1,size(Cidx)
!              LoadFact(Cidx(I) ) = sum(g%Riv%TM(j)%ptr(I,:) )
!            ENDDO
!      
!            call raster2Riv(g%VDZ%m, g%GW(1)%ST, J, St);
!            St(1)=St(2); St(size(St))=St(size(St)-1)
!            
!            ! if w varies in time then this needs to vary in time, too
!            Fact = (r(id)%Riv%w * r(id)%Riv%dx)/(St*LoadFact*Ar)
!         ENDDO
!         status = 1
!      endif

!      LoadFact => tempVars(1)%r(k)%p
!      St => tempVars(2)%r(k)%p
!      Fact => tempVars(3)%r(k)%p
!      Rbmem => tempVars(4)%r(k)%p
!      Rcmem => tempVars(5)%r(k)%p
!      Hg => tempVars(6)%r(k)%p
!      Rb => tempVars(7)%r(k)%p
!      Rc => tempVars(8)%r(k)%p
!      phi => tempVars(9)%r(k)%p
!      
!      n1=37; n2=38
!      IF (abs(Riv%t - w%t)<1D-8) THEN
!        call raster2Riv(g%VDZ%m, g%VDZ%DPerc, k, Rc); 
!        Rc = Rc/w%dt  ! Current loading of cell recharge
!        phi = Rc - Rcmem
!        phi = (Rc - Rcmem - Rbmem)*(Riv%dt/86400D0)/St
!        call raster2Riv(g%VDZ%m, g%GW(1)%h, k, Hg); 
!        Rb = 0D0
!        write(444,*) 'Calculating phi', phi(n1:n2)
!      ENDIF  

      T = Riv%t - w%t
      !Riv%hg =  Hg - (Riv%fG/(Ar * LoadFact)-Rbmem*T + (Rc - Rcmem)*T)/St ! Predicted current Hg level
      !Riv%hg =  Hg - (Riv%fG/(Ar * LoadFact))/ST ! Predicted current Hg level
!      Riv%hg(1)=Riv%hg(2); Riv%hg(size(Hg)) = Riv%hg(size(hg)-1)
!
!      write(444,*) 'Transformed Hg before', Hg(n1:n2)
!      write(444,*) 'Predicted Hg before', Riv%Hg(n1:n2)
!      write(444,*) 'River stage before', Riv%h(n1:n2)+Riv%E(n1:n2)
!      write(444,*) 'Rc before', Rc(n1:n2)
!      write(444,*) 'Rcmem before', Rcmem(n1:n2)
!      write(444,*) 'Rbmem before', Rbmem(n1:n2)

      IF (abs(Riv%t - w%t)<1D-8) THEN
        Riv%vrv(1,:) = Riv%hg  ! temporary - anticipated variable. it will be changed from inside riv_cv
      ENDIF
      if (ANY(Riv%hg<Riv%E-100)) THEN
        call display( num2str(rid)//' '// num2str(w%t)// ' river GW suspiciously low')
      endif
      !write(99,*) 'here2'
      
      call riv_cv(size(Riv%h),Riv%h,Riv%hx,Riv%Ac,Riv%u,Riv%W,Riv%WF,Riv%R,Riv%E,Riv%Qx,   &
     &   Riv%dx,Riv%Mann,Riv%S,Riv%SM,Riv%S0,Riv%SE,Riv%dt,1,Riv%Code,Riv%CType,  &
     &   Riv%rtr(1),Riv%UPST,Riv%UPSV,Riv%DST,Riv%DSV,Riv%STATE,Riv%umax,Riv%hg,Riv%mr,Riv%K,Riv%fG,  &
     &   Riv%t,Riv%vrv)   

      r(rid)%Riv%Eta = r(rid)%Riv%E + r(rid)%Riv%h
      r(rid)%Riv%dt = r(rid)%Riv%dt / 86400.0D0
      r(rid)%Riv%t = r(rid)%Riv%t + r(rid)%Riv%dt

      IF (abs(Riv%t - w%t - w%dt)<1D-8) THEN
         Riv%vrv(2,:) = 0.5D0*(Riv%vrv(2,:)+Riv%fG/(g%DM%d(1)*g%DM%d(2)*w%dt)) ! Inside riv_cv, Riv%vrv(2,:)*dt should yield (m)
      ENDIF
      
!      write(444,*) 'fG term', Riv%fG(n1:n2)/(Ar * LoadFact(n1:n2))/T
!      T = Riv%dt
!      Riv%hg =  Hg - (Riv%fG/(Ar * LoadFact)-Rbmem*T + (Rc - Rcmem)*T)/St ! Predicted Hg level
!      write(444,*) 'Predicted Hg after', Riv%Hg(n1:n2)
!      write(444,*) 'delta mass after', St(n1:n2)*(Riv%Hg(n1:n2) - Hg(n1:n2))
      
!      IF (abs(Riv%t - w%t - w%dt)<1D-8) THEN
!         Rcmem = Rc
!         Rbmem = Rb/w%dt
!      ENDIF

      !% Compute Q to downstream
      r(rid)%Net%DSQ = r(rid)%Net%DSQ + r(rid)%Riv%Qx(size(r(rid)%Riv%Qx)-1)
      r(rid)%Riv%EstE = r(rid)%Riv%Eta
      Qx0 = r(rid)%Riv%Qx
      
      end subroutine r_MacCormack

      subroutine raster2Riv(m, dat, j, output, option)
      use vData
      integer m(2)
      integer j, rid, opt
      real*8 dat(m(1),m(2))
      real*8, dimension(size(r(g%Riv%ID(j) )%Riv%h)) :: output
      integer, optional :: option
      real*8, dimension(:), pointer :: dat1d

      opt = 0
      !if (present(option)) opt = option

      rid = g%Riv%ID(j)
      if (opt .eq. 0) then
        dat1D => OneD_ptr(size(dat),dat)
        output(g%Riv%Cidx(j)%ptr) = matmul(dat1d(g%Riv%Lidx(j)%ptr), g%Riv%TN(j)%ptr)
      endif

      end subroutine