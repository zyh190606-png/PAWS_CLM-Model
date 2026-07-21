!

      SUBROUTINE F_oc(n,ZBank,Zch,Zo,ho,Ab,L,dx,dt,Mann,M)
      IMPLICIT NONE
      INTEGER n
      REAL*8 dx,dt,Ar
      REAL*8, DIMENSION(n)::ZBank,Zch,Ab,L,Mann,M,ho,Zo
      ! Zo is eta on land, ho is depth
      INTEGER I,J,K,MAXIT
      REAL*8,DIMENSION(n)::hc,C
      REAL*8 d0,PREC,dhdx,sgn,As,Bs,Q,Ms,ds,T,Ass,GAM,d ! As new land elevation, Bs new riv elevation
      REAL*8 A,B,Bss,Bm,U,MM,dh,hN,Mmx,hh, signs(2)
      logical flipflop
      
      PREC = 1e-12
      Ar = dx**2.D0
      GAM = 0.2D0
      ! For places where landwater flows to river, explicit, elsewhere, implicit
      hc = Zch - ZBank
      C = L*2.D0*dt/Mann
      MAXIT = 15
      WHERE (hc<PREC) hc = 0.D0
      DO I = 1,n
        !IF (I .eq. 60) THEN
        !print*,3
        !endif
      
        IF (Zo(I)>=Zch(I)) THEN  ! From land to River
          d = ho(I)
          IF (d < PREC .OR. Zo(I)<ZBank(I)) THEN
            M(I) = 0.D0
          ELSE
            ! Explicitly Compute Flux
            dhdx= (Zo(I) - max(Zch(I),ZBank(I)))/(dx/2.D0)
            !U = (d**(2.D0/3.D0))*(dhdx**0.5D0)/Mann(I)
            !M(I) = U*2*L(I)*dt*d
            ! predictor
            M(I)= (d**(5.D0/3.D0))*(dhdx**0.5D0)*C(I)
            ! 
            hh = d - M(I)/Ar
            if (hh>0) THEN
              MM = (hh**(5.D0/3.D0))*(dhdx**0.5D0)*C(I)
            else
              MM = 0D0
            endif
            M(I)=0.5D0*(M(I)+MM)
            ! maximum acceptable rise of water level in stream
            Mmx=(Zo(I)-Zch(I))*Ab(I)/(1D0+Ab(I)/Ar) ! (Zo-Q/Ar=Zs, Q=(Zs-Zch)*Ab, equilibrium stage)
            
            M(I)=min(M(I),Ar*ho(I),Mmx)
            
          ENDIF
        ELSE  ! From River to Land
          d = hc(I)
          IF (ZBank(I)>Zch(I) .OR. d<PREC) THEN ! No flow
            M(I) = 0.D0
          ELSE
            ! River Flow to Land, have to implicitly solve (iteratively)
            A = Zo(I); As = A+(Zch(I)-Zo(I))*Ab(I)/(Ar+Ab(I));
            MM = max(ZBank(I),Zo(I))
            B= Zch(I); Bs = 0.5D0*(MM+Zch(I)); 
            signs=0  ! adaptive picard: function is not smooth. curve is wierd
            flipflop=.false.
            DO J = 1,MAXIT
              Ms = (B - Bs)*Ab(I)
              As = Ms/Ar + A
              ds = max(Bs - ZBank(I),1D-5)
              !!T = 1/(C(I)*(ds**(5.D0/3.D0)))
              !!Bss = As + (Ms*T*(dx/2.D0))**2D0
              T = C(I)*(ds**(5.D0/3.D0))/Ab(I)
              Bss = B - T*sign(1D0,Bs - As)*(abs(Bs - As)/(dx/2D0))**0.5D0

              signs(2:3)=signs(1:2)
              IF (abs(Bss-Bs)<1D-4) EXIT
              IF (Bss < Bs) THEN
                GAM = 0.5D0
                signs(1) = 1
              ELSE
                GAM = 0.2D0
                signs(1) = -1
              ENDIF
              if (signs(1)*signs(2)<0 .AND. signs(2)*signs(3)<0) flipflop=.true.
              if (flipflop) then
                GAM = 0.08D0
              endif

              BM = max(MM, Bss)
              Bs = Bs + GAM*(BM-Bs)
              !IF (BM<MM) THEN
              ! Bs = MM+0.0001D0
              !ELSE
              ! Bs = BM
              !ENDIF
            ENDDO

            IF (J .eq. MAXIT+1) THEN
              Bs = BM ! debug stop
            ELSE
              Bs = BM ! debug stop
            ENDIF

            M(I) = -abs((As- A)*Ar)
            
            Mmx = -(Zch(I)-max(Zo(I),ZBank(I)))*Ar/(1D0+Ar/Ab(I))
            M(I) = max(M(I),Mmx) ! most less negative
            ! HOW COME U DO NOT CONVERGE?
            
            !M(I) = (B - BS)*Ab(I)
            !U = M(I)/(2*L(I)*ds*dt)
            !dhdx= (Bs-As)/(dx/2.D0)
            !U = (ds**(2.D0/3.D0))*(dhdx**0.5D0)/Mann(I)
            !M(I)=U*2*L(I)*dt*ds
            !As = Ms/Ar + A
          ENDIF
        ENDIF
      ENDDO
      END SUBROUTINE
      
      
