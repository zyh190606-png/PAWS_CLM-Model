#include <preproc.h> 
    ! Flow library for the PAWS watershed model
    ! Michigan State UniV% CEE
    ! Chaopeng Shen
    ! PhD, Environmental Engineering

      MODULE FLOW
      !USE VData
      IMPLICIT NONE
      PRIVATE
      PUBLIC :: GGA,vdz1c,Aquifer2D,applyCond,surface
      PUBLIC :: aquifer2dNL,runoff,ovn_cv,lowland,riv_cv
      !PUBLIC :: Aquifer2D
      !SAVE
      CONTAINS

      SUBROUTINE Aquifer2D(ny,nx,nG,h,E,D,K,Wm,S,MAPP,BBC,NBC,                       &
     &             dd,dt,hN,DR)
      ! Compute Lateral Movement of the Unconfined Aquifer
      ! PDE Unit: [L]/[T]
      ! The difference from Confined Aquifer is that the aquifer thickness is changing
      ! Also the Specific Yield (differential water capacity) is changing
      ! ATYPE=0, confined aquifer, =1, unconfined aquifer
      ! h is hydraulic head (L), E is bottom elevation(L), D is aquifer thickness (for ATYPE=0,D=h-E)
      ! K is conductivity(L/T), W is the source term (L/T)
      ! S is specific water capacity (specific yield for unconfined or storativity for confined)
      ! Mapp is a map of active cells (or percentage active)
      ! BBC is boundary conditions, NBC is number of BBCs.
      
      USE VData
      IMPLICIT NONE
      INTEGER ny,nx,NBC,MTD,nG(2),nCond ! NBC is the size arrays for BC specifications
      INTEGER ATYPE 
      REAL*8,DIMENSION(ny,nx)::MAPP,WTD,K,WTDNEW,D,hN,S
      REAL*8,DIMENSION(ny,nx),Target::E,HH,DR,h,Wm ! E for a 2D aquifer is its bottom
      !REAL*8,DIMENSION(nz,nt) :: THES,FC ! Database items. For computing Specific Yield
      TYPE(BC_type):: BBC(NBC)
      REAL*8 DD(2),DX, DY,dt       ! Assuming DX, DY is the same (uniform grid)
      
      REAL*8,DIMENSION(ny,nx),target:: T,Thalfx,Thalfy                       &
     &           ,Tempx,Tempy,Dperc2
      INTEGER I,J,nGx,nGy,MD(5)
      REAL*8,DIMENSION(ny-2*nG(1),nx-2*nG(2)),target::                       &
     &           A1b,A1,A2,A3,A4,A5,Rm,TEMP
      REAL*8,DIMENSION((ny-2*nG(1))*(nx-2*nG(2)))::RHS
      REAL*8,DIMENSION((ny-2*nG(1))*(nx-2*nG(2)),5)::M
      REAL*8,DIMENSION(:,:),Pointer:: hNew,hE,hW,hNN,hS,DRR
      INTEGER N1,NM,M1,MM,NNY,NNX,NN,n,iter,itmax,np,IND(2)
      REAL*8 err,meanW,meanH0,meanH,meanE,H0(ny,nx),DP2,PK(2)
      REAL*8,DIMENSION(:),POINTER::pp1,pp2,pp3,pp4,pp5 ! print 
      REAL*8,DIMENSION(:),ALLOCATABLE::pk1,pk2,pk3,pk4,pk5
      REAL*8,DIMENSION(3,3):: TEMPI1,TEMPI2
      !meanW = sum(W*Mapp)/size(h); meanH0=sum((H-E)*S*Mapp)/size(h) ! print *, debugging
      H0 = H
      
      nGy = int(nG(1)); nGx = int(nG(2))
      DY  = DD(1); DX = DD(2)
      N1=nGy+1;NM=ny-nGy;M1=nGx+1;MM=nx-nGx
      NNY = ny-2*nG(1); NNX=nx-2*nG(2); NN = NNY*NNX; n=NN
      Dperc2 = 0D0
      ! This solver itself is indiscriminate about the method used for the modeling. 
      ! The actually method used is determined from discretization and Boundary Condition
      WHERE (Mapp<=0) K=0.0D0
      
      !IF (ATYPE .eq. 0) THEN
      !ELSEIF (ATYPE .eq. 1) THEN
      !  D = h - E  ! Linearize Transmissivity Term
      !ENDIF
      HH = h ! Make a copy
      ! Apply Inactive Cell Condition, if any more than indicated by MAPP
      CALL applyCond(0,BBC,NBC,n,A1,A2,A3,A4,A5,M,RHS,K,HH)
      T = D * K
      Tempx = T + CSHIFT(T,dim=2,shift=1)
      Thalfx = (T*CSHIFT(T,dim=2,shift=1)/Tempx) *(2.D0/DX)
      !hS => T(8:10,15:17) ! print debug
      !hS => Tempx(8:10,15:17)
      
      Tempy = T + CSHIFT(T,dim=1,shift=1)
      Thalfy = (T*CSHIFT(T,dim=1,shift=1)/Tempy) *(2.D0/DY)
      WHERE (Tempx<=0) Thalfx = 0.0D0
      WHERE (Tempy<=0) Thalfy = 0.0D0
      
      ! Calculate Coefficients
      ! Here A1,A2,A3,A4,A5 meaning have been updated and it's different from Matlab code GW_sol_fd
      ! A2 is relation to its SN, A3 to its NN, A4 to its WN, A5 to its EN
      ! This is consistent with Ghost Cell numbering
      Tempx = Thalfx/(DX)
      Tempy = Thalfy/(DY)
      A2 = -Tempy(N1-1:NM-1,M1:MM)
      A3 = -Tempy(N1:NM,M1:MM)
      A4 = -Tempx(N1:NM,M1-1:MM-1)
      A5 = -Tempx(N1:NM,M1:MM)
      A1b= -(A2+A3+A4+A5)
      Temp = S(N1:NM,M1:MM)/dt
      A1 = Temp + A1b
      Rm  = Temp*HH(N1:NM,M1:MM) + Wm(N1:NM,M1:MM)
      ! Implicit Source Terms, leakance to lower aquifer, etc
      !hS => A1(7:9,14:16) ! print debug
      !hS => R(7:9,14:16)
      !hS => A2(7:9,14:16)
      !hS => A3(7:9,14:16)
      !hS => A4(7:9,14:16)
      !hS => A5(7:9,14:16)
      !hS => A1b(7:9,14:16)
      ! Exfiltration if possible
      ! Apply Dirichlet Boundary Condition
      ! Flux BC should be converted to a source term and a no flow BC
      call applyCond(1,BBC,NBC,n,A1,A2,A3,A4,A5,M,Rm,K,HH)
      ! Apply Lower/Upper Boundary conditions
      call applyCond(2,BBC,NBC,n,A1,A2,A3,A4,A5,M,Rm,K,HH)
      ! Assembly Final Matrices, Save to rank 1 format
      !NNY = ny - 2*nGy; NNX = nx - 2*nGx
      !DO J = nGx+1,nx-nGx
      !  IS=(J-1)*NNY+1; IE = J*NNY; 
      !  M(IS:IE,1)=A4(:,J)
      !  M(IS:IE,2)=A2(:,J)
      !  M(IS:IE,3)=A1(:,J)
      !  M(IS:IE,4)=A3(:,J)
      !  M(IS:IE,5)=A5(:,J)
      !  RHS(IS:IE)= R(:,J)
      !ENDDO
      ! M=[A4 A2 A1 A3 A5]
      CALL elem(NN,A4,M(:,1))
      CALL elem(NN,A2,M(:,2))
      CALL elem(NN,A1,M(:,3))
      CALL elem(NN,A3,M(:,4))
      CALL elem(NN,A5,M(:,5))
      CALL elem(NN,Rm,RHS)
      
      ! Solve equation
      MD = (/-NNY,-1,0,1,NNY/)
      itmax = 100
      hNew => HH(N1:NM,M1:MM) ! Use hNew as a pointer to reduce copying data
      !DO I = N1,NM
      !DO J = M1,MM
      ! if (ISNAN(W(I,J))) THEN
      ! ENDIF
      !ENDDO
      !ENDDO
      CALL pcg(n,5,M,MD,RHS,hNew,4,1d-8,itmax,iter,err) !hNew is hydraulic head
      !hNN => HH(8:10,15:17) ! print debug
      !TEMPI1 =  hNN - hE
      !TEMPI2 = hW*dt/S(8:10,15:17) ! print debug
      !hN(N1:NM,M1:MM) = reshape(hNew,(/NNY,NNX/))-E(N1:NM,M1:MM)
      hN = HH
      ! Compute Lateral Fluxes
      hE => HH(N1:NM,M1+1:MM+1)
      hW => HH(N1:NM,M1-1:MM-1)
      hNN=> HH(N1+1:NM+1,M1:MM)
      hS => HH(N1-1:NM-1,M1:MM)
      hNew=> HH(N1:NM,M1:MM)
      DRR=> DR(N1:NM,M1:MM)
      DRR= -(A4*(hW-hNew)+A5*(hE-hNew)+A2*(hS-hNew)+A3*(hNN-hNew))
      ! recharge to lower
      np = size(h)
      !np = size(BBC(1)%Cidx)    ! Note: in the new Fortran version of the
      ! model, BBC is actually g%GW(x)%Cond(x)
      
      !ALLOCATE(pk1(np),pk2(np),pk3(np),pk4(np),pk5(np))
      !pp1 => oneD_Ptr(nn,hNew); pk1 = pp1(BBC(1)%Cidx)
      !pp2 => oneD_Ptr(ny*nx,h0); pk2 = pp2(BBC(1)%idx)
      !pp3 => oneD_Ptr(ny*nx,HH); pk3 = pp3(BBC(1)%idx)
      !pp3 => oneD_Ptr(ny*nx,Dperc2); pk3 = pp3(BBC(1)%idx)
      
      !pk3 = BBC(1)%K*(pk1 - BBC(1)%h)/BBC(1)%D
      !DP2 = sum(pk3)/size(h)
      
      !DEALLOCATE(pk1,pk2,pk3,pk4,pk5)
      
      !meanH = sum((HH-E)*S*Mapp)/size(HH)
      !meanE = meanH-(meanH0 + meanW*dt - DP2*dt)
      !if (abs(meanE)>1e-7) THEN
      !print*,3
      !endif
      !TEMP = (hNew - h(N1:NM,M1:MM))*S(N1:NM,M1:MM)/dt
      ! Percolation to lower layer is computed outside
      END SUBROUTINE Aquifer2D

!====================================================================================

      SUBROUTINE Aquifer2DNL(ny,nx,nG,h,E,D,K,Wm,S,MAPP,BBC,NBC,                       &
     &             dd,dt,hN,DR,FR,NLC,STO,nC)
      ! Non-Linear version of the GW code. 
      ! FR is to multiply (THES-THER), currently not used, applied outside
      ! The nonlinearity comes from using Hilberts' drainable porosity as the ST
      ! PDE Unit: [L]/[T]
      ! Added three input:
      ! HB: bottom of aquifer
      ! NLC, nonlinearity coefficient as: 
      ![1ES,2THES-THER,3ALPHA,4nn,5:-1/nn,6:-(1+1/nn)] 
      ! STO, current storage, put zero for initialization
      ! nc: the number of coefficient
      
      
      USE VData
      IMPLICIT NONE
      INTEGER ny,nx,NBC,MTD,nG(2),nCond,nC ! NBC is the size arrays for BC specifications
      INTEGER ATYPE 
      REAL*8,DIMENSION(ny,nx)::MAPP,WTD,K,WTDNEW,D,hN,S
      REAL*8,DIMENSION(ny,nx),Target::E,HH,DR,h,Wm ! E for a 2D aquifer is its bottom
      REAL*8,DIMENSION(ny,nx),Target::STO ! E for a 2D aquifer is its bottom
      REAL*8,Target:: NLC(ny,nx,nC)
      !REAL*8,DIMENSION(nz,nt) :: THES,FC ! Database items. For computing Specific Yield
      TYPE(BC_type):: BBC(NBC)
      REAL*8 DD(2),DX, DY,dt,FR       ! Assuming DX, DY is the same (uniform grid)
      
      REAL*8,DIMENSION(ny,nx),target:: T,Thalfx,Thalfy                       &
     &           ,Tempx,Tempy,Dperc2,hP,DRT,STT
      INTEGER I,J,nGx,nGy,MD(5)
      REAL*8,DIMENSION(ny-2*nG(1),nx-2*nG(2)),target::                       &
     &           A1b,A1,A2,A3,A4,A5,Rm,TEMP
      REAL*8,DIMENSION((ny-2*nG(1))*(nx-2*nG(2)))::RHS
      REAL*8,DIMENSION((ny-2*nG(1))*(nx-2*nG(2)),5)::M
      REAL*8,DIMENSION(:,:),Pointer:: hNew,hE,hW,hNN,hS,DRR
      INTEGER N1,NM,M1,MM,NNY,NNX,NN,n,iter,itmax,np,IND(2),IIT
      INTEGER I1,J1
      REAL*8 err,meanW,meanH0,meanH,meanE,H0(ny,nx),DP2,PK(2),errP,DL
      REAL*8,DIMENSION(:),POINTER::pp1,pp2,pp3,pp4,pp5 ! print 
      REAL*8,DIMENSION(:),ALLOCATABLE::pk1,pk2,pk3,pk4,pk5
      REAL*8,DIMENSION(3,3):: TEMPI1,TEMPI2
      REAL*8,DIMENSION(:,:),POINTER::ES,RP
      REAL*8,DIMENSION(:),POINTER::NL
      REAL*8 HD,HHI
      meanW = sum(Wm*Mapp)/size(h); meanH0=sum((H-E)*S*Mapp)/size(h) ! print *, debugging
      
      nGy = int(nG(1)); nGx = int(nG(2))
      DY  = DD(1); DX = DD(2)
      N1=nGy+1;NM=ny-nGy;M1=nGx+1;MM=nx-nGx
      NNY = ny-2*nG(1); NNX=nx-2*nG(2); NN = NNY*NNX; n=NN
      Dperc2 = 0D0
      ! This solver itself is indiscriminate about the method used for the modeling. 
      ! The actually method used is determined from discretization and Boundary Condition
      WHERE (Mapp<=0) K=0.0D0
      
      ES => NLC(:,:,1)
      ! Calculate S(h) if needed
      DO I=1,ny
      DO J=1,nx
        IF (Mapp(i,j)>0) THEN
          ! if STO is initialized, none of it should be <=0
          ![1ES,2THES-THER,3ALPHA,4nn,5:-1/nn,6:-(1+1/nn),7:Kl,8DZL,9HL,10HB] 
          IF (STO(I,J)<=0D0) THEN
           NL => NLC(I,J,:)
           STO(I,J)=(NL(1)-h(I,J))*NL(2)*                                         &
     &       (1+abs(NL(3)*(NL(1)-h(I,J)))**NL(4))**NL(5)+                         &
     &       NL(2)*(h(I,J)-NL(10))
          ELSE
           EXIT
          ENDIF
        ENDIF
      ENDDO
      ENDDO
      
      
      T = D * K
      Tempx = T + CSHIFT(T,dim=2,shift=1)
      Thalfx = (T*CSHIFT(T,dim=2,shift=1)/Tempx) *(2.D0/DX)
      !hS => T(8:10,15:17) ! print debug
      !hS => Tempx(8:10,15:17)
      
      Tempy = T + CSHIFT(T,dim=1,shift=1)
      Thalfy = (T*CSHIFT(T,dim=1,shift=1)/Tempy) *(2.D0/DY)
      WHERE (Tempx<=0) Thalfx = 0.0D0
      WHERE (Tempy<=0) Thalfy = 0.0D0
      ERRP = 1D3; IIT = 0
      DRT = 0D0; TEMP = 0D0
      Rm = 0D0; A1=0D0;STT=0D0

      ! Calculate Coefficients
      ! Here A1,A2,A3,A4,A5 meaning have been updated and it's different from Matlab code GW_sol_fd
      ! A2 is relation to its SN, A3 to its NN, A4 to its WN, A5 to its EN
      ! This is consistent with Ghost Cell numbering
      Tempx = Thalfx/(DX)
      Tempy = Thalfy/(DY)
      A2 = -Tempy(N1-1:NM-1,M1:MM)
      A3 = -Tempy(N1:NM,M1:MM)
      A4 = -Tempx(N1:NM,M1-1:MM-1)
      A5 = -Tempx(N1:NM,M1:MM)
      A1b= -(A2+A3+A4+A5)
      
      hP = h; ! hp is iteration storage
      
      !hh = 0D0 ! hh is newton update
      ![1ES,2THES-THER,3ALPHA,4nn,5:-1/nn,6:-(1+1/nn),7:Kl,8DZL,9HL,10HB] 
      DO WHILE (ERRP > 5D-5 .AND. IIT<=5)
      IF (IIT .eq. 20) THEN
       DO I=1,ny
       DO J=1,nx
        IF (Mapp(i,j)>0) THEN
        NL => NLC(I,J,:)
        TEMP(I,J)=NL(2)*(1-(1+abs(NL(3)*(HP(I,J)-NL(1)))**NL(4))             &
     &           **NL(6))
        A1(I,J) = Temp(I,J)/dt + A1b(I,J)
        Rm(I,J)  = -DR(I,J)-Wm(I,J)
        ENDIF
       ENDDO
       ENDDO
      ELSE
      ! remember R,A1,A2,A3,A4,A5 are one cell less
       DO I=N1,NM
       DO J=M1,MM
        IF (Mapp(i,j)>0) THEN
        I1 = I-nGy; J1=J-nGx
        NL => NLC(I,J,:)
        HD = max(NL(1)-HP(I,J),0D0); 
        DL = NL(7)/NL(8)
        HHI=min(max(HP(I,J)-NL(10),0D0),NL(1)-NL(10))
        TEMP(I,J)=NL(2)*                                                         &
     &       (1-(1+abs(NL(3)*HD)**NL(4))**NL(6))
        TEMP(I,J)=max(5D-4,TEMP(I,J))
        STT(I,J)=HD*NL(2)*                                                       &
     &       (1+abs(NL(3)*HD)**NL(4))**NL(5)+                                     &
     &       NL(2)*HHI
        A1(I1,J1) = Temp(I,J)/dt + A1b(I1,J1) + DL
        DRT(I,J)= -(A4(I1,J1)*(hP(I,J-1)-HP(I,J))+A5(I1,J1)*(hP(I,J+1)-             &
     &  HP(I,J))+A2(I1,J1)*(hp(I-1,J)-hp(I,J))+A3(I1,J1)*(hp(I+1,J)                 &
     &      -hp(I,J)))
        Rm(I1,J1)=-((STT(I,J)-STO(I,J))/dt+NL(7)*(hp(I,J)-NL(9))/NL(8)               &
     &      -Wm(I,J)-DRT(I,J))
        ENDIF
       ENDDO
       ENDDO
      ENDIF
      
      ! in this scheme BC is implemented in the explicit term
      !call applyCond(1,BBC,NBC,n,A1,A2,A3,A4,A5,M,R,K,HP)
      
      ! M=[A4 A2 A1 A3 A5]
      CALL elem(NN,A4,M(:,1))
      CALL elem(NN,A2,M(:,2))
      CALL elem(NN,A1,M(:,3))
      CALL elem(NN,A3,M(:,4))
      CALL elem(NN,A5,M(:,5))
      CALL elem(NN,Rm,RHS)
      
      !n=MD; for i=1:NN, for j=1:5, jj=i+n(j); if jj>=1 & jj<=60, M(i,jj)=fort(i,j); end; b(i)=fort(i,end); end; end
      ! Solve equation
      MD = (/-NNY,-1,0,1,NNY/)
      itmax = 100
      hNew => HH(N1:NM,M1:MM) ! Use hNew as a pointer to reduce copying data
      !DO I = N1,NM
      !DO J = M1,MM
      ! if (ISNAN(W(I,J))) THEN
      !  print*, 3
      ! ENDIF
      !ENDDO
      !ENDDO
      !RHS =0D0
      hNew = 0D0 ! Initial guess to the residual is always 0
      CALL pcg(n,5,M,MD,RHS,hNew,4,1d-7,itmax,iter,err)
      !if (iter >=100) then
      !do i=1,NN
      !write(122,"(6g16.9)") M(i,1:5),RHS(i)
      !enddo
      !print *,3
      !endif
      
      ERRP=MAXVAL(ABS(hNew))
      hP(N1:NM,M1:MM) = hP(N1:NM,M1:MM) + hNew
      IIT = IIT + 1
      ENDDO
      
      S(N1:NM,M1:MM) = TEMP(N1:NM,M1:MM)
      DR = DRT
      STO = STT
      ! recharge to lower
      np = size(h)
      np = size(BBC(1)%Cidx)
      hN = hP
      !ALLOCATE(pk1(np),pk2(np),pk3(np),pk4(np),pk5(np))
      !pp1 => oneD_Ptr(nn,hNew); pk1 = pp1(BBC(1)%Cidx)
      !pp2 => oneD_Ptr(ny*nx,h0); pk2 = pp2(BBC(1)%idx)
      !pp3 => oneD_Ptr(ny*nx,HH); pk3 = pp3(BBC(1)%idx)
      !pp3 => oneD_Ptr(ny*nx,Dperc2); pk3 = pp3(BBC(1)%idx)
      
      !pk3 = BBC(1)%K*(pk1 - BBC(1)%h)/BBC(1)%D
      !DP2 = sum(pk3)/size(h)
      
      !DEALLOCATE(pk1,pk2,pk3,pk4,pk5)
      
      !meanH = sum((HH-E)*S*Mapp)/size(HH)
      !meanE = meanH-(meanH0 + meanW*dt - DP2*dt)
      !if (abs(meanE)>1e-7) THEN
      !print*,3
      !endif
      !TEMP = (hNew - h(N1:NM,M1:MM))*S(N1:NM,M1:MM)/dt
      ! Percolation to lower layer is computed outside
      END SUBROUTINE Aquifer2DNL

      FUNCTION STNL(H,NL)
      ![1ES,2THES-THER,3ALPHA,4nn,5:-1/nn,6:-(1+1/nn),7,HB]
      REAL*8 STNL,H,NL(*)
      STNL=(NL(1)-h)*NL(2)*                                                        &
     &       (1+abs(NL(3)*(NL(1)-h))**NL(4))**NL(5)+                               &
     &       NL(2)*(h-NL(7))
      END FUNCTION STNL

      FUNCTION STNLD(H,NL)
      ![1ES,2THES-THER,3ALPHA,4nn,5:-1/nn,6:-(1+1/nn),7,HB]
      REAL*8 STNLD,H,NL(*)
      STNLD=NL(2)*(1D0-(1D0+abs(NL(3)*(H-NL(1)))**NL(4))**NL(6))
      END FUNCTION STNLD

      SUBROUTINE applyCond(TYP,BBC,NBC,n,A1,A2,A3,A4,A5,M,RHS,K,H)
      ! Apply boundary condition for the groundwater model
      ! CP Shen
      USE VDATA
      integer TYP,n,NBC
      TYPE(BC_type) :: BBC(NBC)
      REAL*8,DIMENSION(n)::RHS,A1,A2,A3,A4,A5,K,H,TEMP
      REAL*8 M(n,5)
      integer I,NI
      integer,DIMENSION(:),POINTER :: idxy,Cidx,Nidx,Sidx,Widx,Eidx
      REAL*8, POINTER:: P
      ! M=[A4 A2 A1 A3 A5]

      DO I=1,NBC
       IF (BBC(I)%JType .eq. TYP) THEN
         SELECT CASE (BBC(I)%JType)
          CASE (0)
           idxy=>BBC(I)%idx
           K(idxy) = 0.0D0
          CASE (1)
          ! Indices are row sequential indices
           idxy=>BBC(I)%idx ! In the original matrix
           H(idxy) = BBC(I)%h
           ! These indices are all for the Shrinked Matrix, which is for final matrix building:
           ! These indices give the location of BC Cells' neighbors
           Cidx=>BBC(I)%Cidx
           Nidx=>BBC(I)%Nidx; Sidx=>BBC(I)%Sidx
           Widx=>BBC(I)%Widx; Eidx=>BBC(I)%Eidx
           select case (BBC(I)%bnd)
           ! 0-inside, 1-S,2-N,3-W,4-E
            case (0) ! Internal Cells, for modeling of lakes or streams
             RHS(Sidx)=RHS(Sidx)-A3(Sidx)*BBC(I)%h;M(Sidx,4)=0.D0 ! SN
             RHS(Nidx)=RHS(Nidx)-A2(Nidx)*BBC(I)%h;M(Nidx,2)=0.D0 ! NN
             RHS(Widx)=RHS(Widx)-A5(Widx)*BBC(I)%h;M(Widx,5)=0.D0 ! WN
             RHS(Eidx)=RHS(Eidx)-A4(Eidx)*BBC(I)%h;M(Eidx,1)=0.D0 ! EN
             M(Cidx,:)=0.0D0; M(Cidx,3)=1.0D0
             RHS(Cidx)=BBC(I)%h
            ! Boundary Cells:
            !! CURRENTLY SOME PROBLEMS
            !! IF GHOST CELL IS MARKED 0 in MAPP, their K will be set 0
            !! THEN the A values are all wrong!
            case (1) ! This is a South Boundary, Correct its North Neighbor (To NN's S)
             RHS(Nidx)=RHS(Nidx)-A2(Nidx)*BBC(I)%h
            case (2) ! This is a Nouth Boundary, Correct its South Neighbor (To SN's N)
             RHS(Sidx)=RHS(Sidx)-A3(Sidx)*BBC(I)%h
            case (3) ! This is a West Boundary, Correct its East Neighbor (To EN's W)
             RHS(Eidx)=RHS(Eidx)-A4(Eidx)*BBC(I)%h
            case (4) ! This is a East Boundary, Correct its West Neighbor (To WN's E)
             RHS(Widx)=RHS(Widx)-A5(Widx)*BBC(I)%h
           end select
          CASE (2) ! Upper/Lower Layers Implicit Condition
          ! CONSIDER ARRAY SIZES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !TEMP = 0D0
          NI = size(BBC(I)%K)
             idxy=>BBC(I)%Cidx
             TEMP(1:NI) = BBC(I)%K/BBC(I)%D
             RHS(idxy)=RHS(idxy)+TEMP(1:NI)*BBC(I)%h
             A1(idxy)=A1(idxy)+TEMP(1:NI)
         END SELECT
       ENDIF
      ENDDO
      END SUBROUTINE applyCond

      SUBROUTINE vdz1c(NZ,NA,E,DF,Dperc,h,THE,K,CI,SRC,ddt,dtP,                &
     &   MAXIT,DRAIN,WTL,THER,THES,Lambda,alpha,n,KS,DZ,DDZ,BBCI,                  &
     &   Rm,STATE,GWL,CC,PRatio,ho,SN,Rff,m_err)
      ! Vadose Zone Model using Celia Picard iteration scheme
      ! CP Shen. CEE, MSU
      ! PDE Unit: T^-1
      ! Boundary Condition is implemented with a ghost cell on the top or in the bottom
      ! For Equation BC the top layer is included in the computation.
      ! Input:
      ! ddt: an initial guess of time step (recorded from last time step)
      ! dtP: [minDT maxDT macroDT]
      ! BBCI: a derived type containing boundary conditions
      ! BBCI(2)%BV: boundary condition values. explained inside
      ! SN = (S0^0.5)/(ManningN*ldist). If turning off runoff is desired, set this to 0.
      ! R,CC,DRAIN not used in current version
      ! now, use R to record steps solved by the botdrain and saturation status
      ! ho: initial groundwater abstraction for runoff
      ! NZ: number of cells, NA: number of active cells
      ! WTL: the first layer whose top is below water table. This is an initial guess, may be modified in this subroutine
      ! other variables self-explanatory, see variable reference
      ! Output:
      ! modified h, THE, K (unimportant), WTL, GWL, Rff, DF, Dperc

       USE VData
       IMPLICIT NONE
       integer NZ, NA,DRAINL,WTL ! Matrix Size, Active Cells
       REAL*8 :: DF,Dperc,DM,ddt,dtP(3),DRAIN,DRAINZ,DRoot,DUnsat,SN ! Infiltration and percolation in the time step
       REAL*8, DIMENSION(NZ) :: h, K, THE, CI, E,Rm(*),THESAVE ! E is elevation of cell centers
       REAL*8, DIMENSION(NZ) :: THER,THES,KS,DZ,DDZ,Lambda,alpha,n,SRC
       REAL*8  err, h_err, t, GWL, GWL0, RR, CC, PRatio, ho, Rf, DR,Rff ! Rff is passed in, Rf is local
       REAL*8, DIMENSION(NZ) :: Khalf,temp,a,b,c,d,THEP,CP,HP,m ! K(i=1/2) and DZ(i)
       integer :: ITER, NBC, MAXIT, STATE, N1,NM, JType1, JType2,NMAX
       REAL*8,DIMENSION(nz) :: HTEMP,THETEMP, temp1, temp2, tn1,tn2,tn3          &
     &                         ,tn4,tn5,DTHES,f,THEOLD
       REAL*8 :: PREC , temp3(NZ-2), temp4, tempDF, DTHE, THE0,x,xL,THE2
       TYPE(BC_type) :: BBC(2), BBCI(2) ! Boundary Condition!
       !ddt(1):current dt, ddt(2):dtmin, ddt(3):dtmax, ddt(4):MacroTimeStep(before going out of program)
       REAL*8 :: dt, RDT, dtmin, dtmax, SDT,tf,S(NZ),WTD,WTDL,SS(NZ)
       REAL*8 :: dtPREC,HH,hh1,hh2,hh3,t1,t2,t3,q1,FULL,FULL1,km,EUC
       REAL*8 :: ERRP,RS1,RS2,SS1,m0,qbot,qtop,qm,DFMX,Imx,P(6)
       integer I,J,IMAX,ITOTAL,MAXIT2,NM1 !,IWTL,NM2,WTL0,W1,W2,EN1
       integer ROUTE(10),sigBC,BC_OLD,NAA,WTL0
       REAL*8 ::m_err,m_err2,THE1,DFF,h0(nz),SF,Dpercc,tt,Qin,GST
       REAL*8 :: BC1,BC2,BC3,BC4,THEB,DZB,gm(10),sumS,ZT,KT ! g is temporary storage variable
       REAL*8,POINTER:: BV(:)
       REAL*8 Dperc2,hh4
       LOGICAL:: SATH,BC_NC,DEBUG,CVG
       PARAMETER(DEBUG = .TRUE.) ! CHANGE THIS TO ACTIVATE MASS CONSERVATION DEBUGGING

       ! f is the fraction of the Dupuit flow


       !if (h(2)<-30 .OR. h(2)>0 .OR. SRC(1)>1) THEN
       !print *,3
       !endif

       !IF (SRC(2)<0) THEN
       !PRINT*,3
       !endif

       ERRP = 5.D-5
       BBC = BBCI ! Do not modify input BBC
       t = 0.0D0
       dt = ddt;dtmin=dtP(1); dtmax=dtP(2);tf = dtP(3);
       dtPREC=0.01*dtmin
       RDT = tf-t;
       MAXIT = 10
       MAXIT2= 11
       HTEMP = h
       SDT = dt
       JType1 = BBCI(1)%JType
       JType2 = BBCI(2)%JType
       PREC = 1.0D-12
       N1=1; NM=NA-1; ! Will be modified later
       NM1=NA-1
       DFF = 0D0
       Rm(1:2)=0D0;
       !if (JType1 .eq. 5) Rf = 0D0
       !MMASS = .FALSE.
       ROUTE = 0
       ZT = 0D0
       t1 = 0D0
       NMAX = min(WTL,na)
       !THE_SAVE=THE
       !IWTL = 0;
       !hh1 = 0.D0; EN1 = 1;
       SATH = .FALSE.
       CVG = .TRUE.
       Rf =0D0
       BC_OLD = 0
       IF (JType2 .eq. 3) THEN
        h(NA) = BBC(2)%Bdata - BBC(2)%BV(6) !! assume EQ Lower BC!!!!!!!!!!!!!!!!!!!!!
       ELSE
        h(NA) = BBC(2)%BData - BBC(2)%BData2
       ENDIF

       ! ###### print mask this debugging statement
       hh2 = h(1)
       SS1 = SRC(1)
       h0 = h;  !
       ! !
       IF (DEBUG) THEN
       NAA = NA
       THE0 = sum(THE(2:NA-1)*DZ(2:NA-1))  !
       WTL0 = WTL
       GWL0 = GWL
       THESAVE = THE
       ENDIF
       THEB = THE(NA)
       ! Make things easier
       ! 00000000000000000000000000000000000000000000000000000000000000000000000000000000
       IF (SRC(1)<0D0 .AND. h(1)>0D0) THEN
        HH = min(h(1),abs(SRC(1)*dtP(3)))
        h(1)=h(1)-HH
        SRC(1)=SRC(1)+HH/dtP(3)
       ENDIF

       IF (h(1)<0D0 .AND. SRC(1)>0D0) THEN
        tt = min(-h(1),SRC(1)*dtP(3))
        h(1)=h(1)+tt
        SRC(1)=SRC(1)-tt/dtP(3)
       ENDIF

       IF (h(1)<0D0 .and. JType1 .eq. 5) THEN
        hh3 = h(1) ! This will be added back to the final solution
        h(1) = 0D0
       ELSE
        hh3 = 0D0
       ENDIF

       ! IF Fully Saturated, GWL should be where the lower BC table is
       ! disable this as lower bottom is fixed, recharge is considered as passing the interface
       ! between na-1 and na
       !IF (BBC(2)%JType==1) THEN
       ! does it matter where this head is put?
        ! NM2 = insertLoc(NZ,E,BBC(2)%BData,1,1)
        ! NA = max(NM2+1,NA)
        ! if (NA>NZ) NA = NZ
        ! WTL = NA
       !ENDIF

       ! ########11111111111111111111111111111111111111111111111111111111111111111111111111
       ! time loop, input h, SRC, mass balance must be conserved
       tn1 = n/(n-1D0)
       tn2 = 1D0/tn1
       tn4 = (2D0*n-1D0)/n
       DTHES = THES-THER
       S = SRC
       DO WHILE (t < tf-1e-12)
       RDT = tf-t;
       dt = min(dt,RDT,dtmax)
       dt = max(dt,min(tf-t,dtmin))
       !THEOLD = THE ! save for computing recharge
       !S(1:NA) = SRC(1:NA);  !temporary source term storage for the fractional time steps
       S(1) = SRC(1)  ! only the first one may get changed
       IF (.NOT. CVG) h(1) = hh1 ! h(1) may have been changed during the time step, if not convergence, then get it back
       hh1 = h(1)

       IF (h(1)<0D0 .AND. S(1)>0D0) THEN
        tt = min(-h(1),S(1)*dt)
        h(1)=h(1)+tt
        S(1)=S(1)-tt/dt
        ! If this modification wasn't successful,i.e, not convergent, and time step was reduced and sent back here
        ! then it will get back hh1 as the statement above
       ENDIF
       hh4 = h(1) ! will be used for BBC(1) upper BC

       sumS = sum(S(2:NA-1)*DZ(2:NA-1))
       ! ###################BBBBBBBBBBBBBBBBBBBBBBBBBBBB

       IF (abs(STATE)<1.D-5) THEN
         ! SHOULD START WITH H!!!!
         !SS = (THE - THER)/(THES - THER)
         !m = 1 - 1/n
         !h = (SS**(-1/m)-1)**(1/n)/(-abs(alpha))
         ! BELOW IS WITH h initialization
         temp1 = abs(alpha*h)**n
         temp2 = 1+temp1
         tn3 = (THE-THER)/(DTHES)
         THE = THER + (DTHES)/(temp2**(tn2))
         CI = -DTHES/(temp2 **tn4)                                          &
     &    * temp1 *(n-1D0)/h
         !CI = -(THES-THER)/(temp2 **((2*n-1)/n))                                   &
         !&    * temp1 *(n-1)/h
       ENDIF
       K = KS*(tn3)**Lambda*                                                &
     &    (1-(1-(tn3)**tn1)**tn2)**2.0D0
       !K = KS*((THE-THER)/(THES-THER))**Lambda*                                   &
       !&    (1-(1-((THE-THER)/(THES-THER))**tn1)**tn2)**2.0D0
       WHERE (h>-1d-5)
         K=KS
         CI = 0.0D0
         THE = THES
       ENDWHERE
       !IF (minval(K)<=0D0) THEN
       !print*,4
       !endif

       THEP = THE
       CP = CI

       ! Upper and Lower Boundary Cell for Dirichlet and no flow
       ! K(1) is set as KS(2) because if interaction with ground is necessary
       ! assuming ground surface is saturated
       K(1) = KS(2)*PRatio  ! K(1) ADJUSTMENT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !K(NA)=KS(NA)
       K(NA)=K(NA-1)
       err = 1.0D3
       h_err = 1.0D3
       !Khalf = K * CSHIFT(K,1)/(DZ*CSHIFT(K,1)+ CSHIFT(DZ,1)*K)*2.0D0 ! K_(i+1/2), valid from i=1 to N-1. Last one is invalid
       Khalf = (K + CSHIFT(K,1))/2.D0
       ! see what happens with the geometric mean
       !Khalf = sqrt(K * CSHIFT(K,1))
       !Khalf(1) = sqrt(K(1) * K(2))

       ! Do we need the area of the cell? UNIT??
       temp = Khalf/DDZ  ! valid 1:N-1. If DDZ(1)=0, valid 2:N-1 (in this case, solve from 2:N-1)
       ! a and c will not change from iteration to iteration:
       !a = - CSHIFT(temp,-1)/DZ   ! relation with i-1, valid 2:N
       !c = - temp/DZ ! relation with i+1, valid 1:N-1
       a(2:NM1) = - temp(1:NM1-1)/DZ(2:NM1)
       c(2:NM1) = - temp(2:NM1)/DZ(2:NM1)
       ITER = 0
       HP = h;
       IF (BBC(2)%JType .eq. 3) THEN
       ! BV is an array providing data that is needed for Equation lower BC
       ! [1 ST,2 dzl (=DZ(uc)/2+D(aquitard)),3 DR, 4, Kl, 5, hl(pressure head), 6 E(center) uc aquifer]
       ! From now on (09/21/09), use center of unconfined aquifer for this cell's position
       ! Bdata= HB
         BV => BBC(2)%BV
         NM = NA ! WATER TABLE CELL
         !h(NM) = BBC(2)%Bdata - BV(6) ! h(NM) -- pressure head --- this is a bug
         ! this statement here causes model to ignore the raise of water table
         ! from previous fractional small time steps and thus over-estimate recharge

         DZB = E(NM-1)-BV(6) ! DDZ(NM-1)
         BC1 =  -Khalf(NM-1)/DZB  ! a(nm)
         BC2 = -BC1+BV(4)/BV(2) ! b(nm) less the storage term
         !BC3 = BV(1)*h(NM)/dt+BV(4)*BV(5)/BV(2)                                    &
         !&             +BV(3)+(Khalf(NM-1)-BV(4))  ! d(nm)
         BC3 =                  BV(4)*BV(5)/BV(2)                                    &
     &            + (Khalf(NM-1)-BV(4))  ! d(nm)-Dupuit flow term less the storage term
         if (E(1)-BBC(2)%Bdata<0.2D0) then ! be careful
           GST = (1D-4+BV(1))/2D0 ! Groundwater storage controlled
         else
           GST = BV(1)
         endif
       ELSEIF (BBC(2)%JType .eq. 4) THEN
         BV => BBC(2)%BV
         ! BV meanings are different from JType=3, pertaining to the confined aquifer
         ! [1 ST,2 DZ(bottom cell),3 DR, 4, Kl, 5, hl(pressure head), 6 E(center) c aquifer, 7 DR2,8 EB uc]
         NM = NA
         ! first determine location of bottom cell
         ! BData: GWL of last time step
         if (E(1)-BBC(2)%Bdata<0.2D0) then ! be careful
           GST = (1D-4+BV(1))/2D0 ! Groundwater storage controlled
         else
           GST = BV(1)
         endif
         EUC = BBC(2)%BData2 ! E center uc
         DZB = EUC - BV(6) ! to c aq.
         BC3 = E(NM-1)-EUC ! to bottom cell center
         BC2 = BV(4)/DZB
         BC1 = GST/dt
         P(3)= -Khalf(NM-1)/BC3 ! a(NM)
         IF (BC2 > PREC) THEN
         P(1)= BC1/BC2+1D0
         P(2)= BC1*P(5)+BV(4)+BV(7) ! may change as hl change
         P(4)= -P(3)+BC2-BC2/P(1)
         P(5)= BV(5) !! hl, may change !!!
         ELSE
         P(1)= 1D0  ! p1
         P(2)= 0D0  ! p2
         P(4)= -P(3)  ! b(NM)
         ENDIF
         ! GST?????
       ELSE
         DZB = DDZ(NM-1)
       ENDIF

       ! #######################
       ! ###### BBBBBBBBBBBBBBBBBBCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        ! Top Boundary Conditions
        !IF (JType1 .eq. 5) CALL topBound(nz,na-1,dz,ddz,h,HP,BBC,                         &
        !&        Khalf,S,sigBC,DFMX,qbot,qsur,dt,THE,THES,KS,hh3,DZB,sumS) ! modify HP instead of h
       !JType = 1; Bdata = 0D0
       IF (JType1 .eq. 5) THEN
       RS1 = S(1)+h(1)/dt
       Imx = Khalf(1)*((max(h(1),0D0)-h(2))/(DZ(2)/2D0)+1D0)
       qtop = RS1
       DFMX = sum((THES(2:NM1)-THE(2:NM1))*DZ(2:NM1))
       IF (DFMX<1d-5) THEN ! already saturated
       ! we must be very careful selecting BC2 when it is saturated, it may explode!
          ZT = E(1)-E(NA-1)+DZB
          KT = ZT/(sum(DDZ(1:NM1-1)/Khalf(1:NM1-1))+DZB/Khalf(NM1))
          gm(1) = GST;
          gm(2) = KT/ZT;
          tt = (gm(1)*h(NM1)+KT+BV(3))/(gm(1)+gm(2))
          !qbot = KT*(-tt/ZT+1)/20D0  ! must consider that at this time storage is very small!
          qbot = KT*(-tt/ZT+1)/50D0
          Qin = (sumS+qtop-qbot)*dt
          IF (Qin>0 .or. qtop-BV(3)>0D0) THEN
            !IF (swHead .eq. 1) THEN
              BBC(1)%JType = 3 ! ponding
            !ELSE
            !  BBC(1)%JType = 1
            !  BBC(1)%Bdata = Qin !????????????????????????
            !ENDIF
            ROUTE(9)=ROUTE(9)+1
          ELSE
            BBC(1)%JType = 2 ! this cannot be exfiltration because SATH and Qin<0, this can be evaporation
            BBC(1)%Bdata = qtop
            h(1) = 0D0; hp(1)=0D0
            ROUTE(10)=ROUTE(10)+1
          ENDIF
       ELSE
          qbot = Khalf(NM1)*((h(NM1)-h(NM1+1))/DZB+1D0)
          Qin = (sumS+qtop-qbot)*dt
          IF (Qin>DFMX) THEN ! will be saturated
            !IF (swHead .eq. 1) THEN
            if (qbot<0D0 .and. -qbot>qtop*10D0) THEN
              ! this 'saturation is dominated by groundwater discharge
              ! need to verify that GWL is indeed higher than surface
              if (BV(6)+h(NA)>E(1)) THEN
                BBC(1)%JType = 3
              else
                BBC(1)%JType = 2
                BBC(1)%BData = qtop
                HP(1) = 0D0
              endif
            else
              BBC(1)%JType = 3 ! ponding
            endif
            !ELSE
            !  BBC(1)%JType = 1
            !  BBC(1)%Bdata = Qin - DFMX
            !ENDIF
          ELSE
            IF (qtop<=0D0) THEN
              ! assuming this controlled and should not destroy the solution
              BBC(1)%JTYpe = 2
              IF (t .eq. 0D0 .and. qtop>-1D-2) THEN
                 BBC(1)%Bdata = qtop
              ELSE
                 BBC(1)%Bdata = 0D0
              ENDIF
            ELSE
              !q1 = Khalf(2)*((h(2)-h(3))/(DDZ(2)/2D0)+1D0) ! why divide by 2? this is a bug
              !q1 = Khalf(2)*((h(2)-h(3))/DDZ(2)+1D0)  ! worked okay, but underestimate h(2) ending value
              !q1 = Khalf(2)*((0D0-h(3))/DDZ(2)+1D0)
              km = Khalf(2)
              DO I = 2,nm
                km = min(km,khalf(I))
                IF (THE(I)<THES(I)-PREC) EXIT
              ENDDO
              !FULL=(qtop - max(0D0,Khalf(2)))*dt-(THES(2)-THE(2))*DZ(2)
              !FULL1= full/dtp(3)/khalf(2) ! index from experience
              FULL=(qtop - km)*dt-(THES(2)-THE(2))*DZ(2)
              FULL1= full/dtp(3)/km
              IF (FULL1>5D0 .or. (qtop>Imx .AND. qtop>KS(2))) THEN
              !IF ((qtop>Imx .AND. qtop>KS(2))) THEN
                BBC(1)%JType = 3 ! 3 should work here
                ROUTE(7)=ROUTE(7)+1
              ELSE  ! I or 2????
                b(2) = CP(2)/dt-c(2)
                d(2) = CP(2)*H(2)/dt+(qtop-Khalf(2))/DZ(2)+S(2)
                rr = d(2)/b(2)
                if (d(2)>0 .and. rr>2000) then
                 ! may lead to convergence issue, now we first infiltrate that water
                 THE(2)=THE(2)+qtop*dt/dz(2)
                 t2 = (THE(2)-THER(2))/(THES(2)-THER(2))
                 t3 = -N(2)/(N(2)-1D0)
                 h(2)=((t2**t3-1D0)**(1D0/n(2)))/(-abs(alpha(2)))
                 BBC(1)%Bdata = 0D0
                 DF = DF + BBC(1)%BData*dt
                 ROUTE(5)=ROUTE(5)+10
                else
                 BBC(1)%Bdata = qtop
                endif
                BBC(1)%JType = 2
                ROUTE(8)=ROUTE(8)+1
              ENDIF
            ENDIF
          ENDIF
       ENDIF
       IF (BBC(1)%JType .eq. 2 .and. h(1)>0D0) THEN
         S(1)=0D0
         BBC(1)%BData2 = h(1)
         HP(1) = 0D0; h(1) = 0D0
       ELSE
         BBC(1)%BData2 = 0D0
       ENDIF
       ENDIF ! IF (JType3 .eq. 5) THEN
       ! ###### BBBBBBBBBBBBBBBBBBCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       ! now I've taken the approach to pre-decide BC, and correct for too much infiltration
       ! ###### 222222222222222222222222222222222222222222222222222222222222222222222222
       ! iteration loop, input HP,S

       ! ################# QQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ
       ! ########### fast exit for saturation excess
       ! reference Chaopeng's thesis Chapter of modeling
       ! for saturation excess, we exit the program via a faster way
       !IF ((BBC(1)%JType .eq. 3) .AND. (DFMX<1d-7)) THEN
       IF (.FALSE.) THEN ! closed for test
         if (ZT .eq. 0D0) THEN
          ZT = E(1)-E(NA-1)+DZB
          KT = ZT/(sum(DDZ(1:NA-2)/Khalf(1:NA-2))+DZB/Khalf(NA-1))
         ENDIF
       !g(1) = BV(1)/dt; g(2) = KT/ZT;  ! computed inside topbound
       gm(3)=BV(4)/BV(2);
       gm(4)=gm(2)/(gm(1)+gm(2)+gm(3))
       gm(5) = KT + BV(3) - BV(4)  ! not doing Wnz here
       gm(6) = 1D0/(1/dt+gm(2)-gm(2)*gm(4))    ! p1
       gm(8) = gm(1)*h(NM)+gm(3)*BV(5)
       gm(7) = h(1)/dt+S(1) + sumS -KT+gm(4)*(gm(8)+gm(5))  ! p2
       !CALL satExcess(g(6),g(7),h(1),g(9),SN,ho) !  Picard
       gm(9) = rtnobd(satExcessFunc,h(1),1D-7,SN,ho,gm(6),gm(7),0D0,0D0)
       !SUBROUTINE satExcessFunc(h1N,f,df,SN,ho,p1,p2,d5,d6)
       !rtnobd(funcd,x0,xacc,INP1,INP2,INP3,INP4,INP5,INP6)
       IF (gm(9)>0) THEN ! successful
       ITER = 1
       h(1) = gm(9)
       h(NM) = (gm(8)+gm(2)*h(1)+gm(5))/(gm(1)+gm(2)+gm(3))
       SF = KT*((h(1)-h(NM))/ZT+1D0)
       tt = SF*dt
       DF = DF + tt - sumS*dt ! must take out the sumS portion.
       Dperc = Dperc + tt
       t = t + dt
       Rf = Rf + SN*(max(h(1)-ho,0D0))**(5D0/3D0)*dt
       route(4)=route(4)+1
       dt = dtP(3)
       h(NM-1) = DZB*SF/Khalf(NM-1)-DZB+h(NM)
       Rm(2) = 1D0 ! must be saturated
       DO I=NM-2,2,-1
         ! update the heads
         h(I) = DDZ(I)*SF/Khalf(I)-DDZ(I)+h(I+1)
       ENDDO
       GOTO 100 ! swear not, this is the single place in the model I've used GOTO. This is to avoid too many embedding--CHAOPENG
       ELSE
        ! saturated column would drain dry
        BBC(1)%JType = 2
        BBC(1)%BData = qtop
        route(4)=route(4)+100
       ENDIF
       ENDIF
       ! END ############# QQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ
       ! ########## end of faster exit
       
       ! Common iterations
       DO WHILE ((.NOT. BC_NC) .OR. (err >=ERRP .OR. h_err>ERRP) .AND.                                    &
     &   (ITER < MAXIT2) .AND. (ITER<MAXIT .OR. abs(dt-dtmin)<=dtPREC))
         ITER = ITER+1
         WHERE (HP>-1e-5)
            THEP = THES
            CP = 0.0D0
         ENDWHERE
         !b = -(a+c) + CP/dt
         b(2:NM1) = -(a(2:NM1)+c(2:NM1)) + CP(2:NM1)/dt
         d = CP*HP/dt + S + (CSHIFT(Khalf,-1)-Khalf)/DZ-(THEP-THE)/dt
         !d(2:NM1) = CP(2:NM1)*HP(2:NM1)/dt + S(2:NM1) +                          &
         !&    (Khalf(1:NM1-1)-Khalf(2:NM1))/DZ(2:NM1)-                              &
         !&    (THEP(2:NM1)-THE(2:NM1))/dt


         ! then h = hp when convergence has been reached, this way we do not modify
         ! mass prematurely

         ! BC(1) is for upper side. BC(2) is for lower side
         ! ################
         ! IMPLEMENT BOUNDARY CONDITIONS
         SELECT CASE(BBC(1)%JType)
         CASE(1) ! Head BC
           N1 = 2
           !HTEMP(1) = BBC(1)%BData
           HH = BBC(1)%BData - E(1)
           HTEMP(1) = HH
           !b(1) = 1 NOT NECESSARY
           !d(1) = BBC(1)%BData
           d(N1) = d(N1) - a(N1)*HH  ! RHS Throw-Over
         CASE(2) ! Flux BC
           IF (JType1 .eq. 5) THEN
             N1 = 2
           ELSE
             N1 = FLOOR(BBC(1)%BData2)+1 ! can be cells other than 2, to work with GGA
           ENDIF
           b(N1) = CP(N1)/dt-c(N1)
           d(N1) = CP(N1)*HP(N1)/dt-(THEP(N1)-THE(N1))/dt +                          &
     &     (BBC(1)%BData-Khalf(N1))/DZ(N1) + S(N1)  ! q sign(direction)!
           HTEMP(1) = HP(1)
         CASE(3) ! Equation BC
           ! This means that the top is actually the surface ground
           ! An equation of source and infiltration will be solved for this layer
           N1 = 1
           IF (hp(1)>ho) THEN
              !qm = SN*((hp(1)-ho)**(2D0/3D0)) ! partial linearization
              !qm = SN*((hp(1)-ho)**(5D0/3D0)) a bug
              qm = SN*(max(h(1)-ho,0D0)**(2D0/3D0))
           ELSE
              qm = 0D0
           ENDIF
           ! SN = 86400*SurfaceSlope^(1/2)/manningN/ldist
           a(1) = 0.D0
           c(1) = - Khalf(1)/(DZ(2)/2.D0)
           b(1) = - c(1) + 1.D0/dt + qm   ! partial linearization
           !d(1) = - Khalf(1)+ S(1)+hh1/dt +qm*ho  ! partial linearization
           d(1) = - Khalf(1)+ S(1)+hh4/dt +qm*ho  ! partial linearization
           !b(1) = -c(1) + 1.D0/dt
           !d(1) = - Khalf(1)+ S(1)+hh1/dt -qm
           a(2) = - Khalf(1)/((DZ(2)*DZ(2))/2.D0) ! used to be a(2) = - K(2)/(DZ(2)/2), a bug!!
           b(2) = -(a(2)+c(2)) + CP(2)/dt !! Unnecessary????
           !d(2) = CP(2)*HP(2)/dt + S(2) + (CSHIFT(Khalf,-1)-Khalf)/DZ-(THEP-THE)/dt
         END SELECT

         ! Lower BC, mostly likely JType=3
         SELECT CASE(BBC(2)%JType)
         CASE(0)
         CASE(1) ! Head BC
           !HTEMP(NZ) = BBC(2)%BData
           NM = NA-1
           !b(NZ) = 1  ! NOT NECESSARY
           !d(NZ) = BBC(2)%BData
           d(NM) = d(NM) - c(NM)*BBC(2)%BData
         CASE(2) ! Flux BC
           NM = NA-1
           b(NM) = CP(NM)/dt-a(NM)
           d(NM) = CP(NM)*HP(NM)/dt-(THEP(NM)-THE(NM))/dt +                    &
     &     (Khalf(NM-1)-BBC(2)%BData)/DZ(NM) + S(NM)
         CASE(3) ! Equation BC--consider lower
           NM = NA
           ! h(NM)=BC1  ! h(NM) may change through  the fractional steps
           ! [1 ST,2 zl (=DZ(uc)/2+D(aquitard)),3 DR, 4, Kl, 5, hl(pressure head), 6 E(center) uc aquifer, 7,GW.DZ]

           ! decide the water table location
           CALL waterTable(nm,h,THE,THES,DZ,0.95D0,WTL,x,BV(7))
           ! now compute the Dupuit flow fraction
           if (WTL >= nm) THEN
             xL = BV(7)+x*DZ(WTL-1)
           else
             xL = E(WTL)-BV(6)+(DZ(WTL)+BV(7))/2D0+x*DZ(WTL-1)
           endif
           NMAX = max(NMAX,WTL)
           DR = BV(3)/xL
           a(NM)=BC1
           !tt = 1/(1D0 + BV(3)/Khalf(NM-1));
           !IF (abs(BV(3))>Khalf(NM-1)) THEN
             tt = 1D0
           !ENDIF
           !tt = min(max(tt,0.3D0),2D0)
           !b(NM)=BC2 + BV(1)*tt/dt
           b(NM)=BC2 + GST/dt
           d(WTL:NM-1) = d(WTL:NM-1) + DR
           d(NM)=BC3 + BV(3)*(BV(7)/xL) + GST*h(NM)*tt/dt ! unit is different for this equation

           IF (WTL>2) d(WTL-1) = d(WTL-1) + x*DR
           GWL = E(WTL)+DZ(WTL)/2D0+x*DZ(WTL-1)

           ! correct coefficients for cell (NM-1)
           ! since for this BC we used uc aquifer center as the cell center, the corresponding coefficients of NM-1 must also be modified
           ! later include a changing C term for NM!
           c(NM-1)=BC1/DZ(NM-1)
           b(NM-1)=-(c(NM-1)+a(NM-1))+CP(NM-1)/dt
         CASE(4) ! bedrock bottom BC
         ! [1 ST,2 DZ(bottom cell),3 DR, 4, Kl, 5, hl(pressure head), 6 E(center) c aquifer, 7 DR2,8 EB uc]
           NM = NA
           CALL waterTable(nm,h,THE,THES,DZ,0.95D0,WTL,x,BV(2))
           if (WTL >= nm) THEN
             xL = BV(2)+x*DZ(WTL-1)
           else
             xL = E(WTL)-BBC(2)%Bdata2+(DZ(WTL)+BV(2))/2D0+x*DZ(WTL-1)
           endif
           ! P(5) is current hl
           NMAX = max(NMAX,WTL)
           DR = BV(3)/xL ! in fact ql
           a(NM)=P(3)
           b(NM)=P(4)
           d(NM)= P(2)/P(1)+BV(2)/xL*BV(3)+(Khalf(NM-1)-BV(4))
           d(WTL:NM-1) = d(WTL:NM-1) + DR
           IF (WTL>2) d(WTL-1) = d(WTL-1) + x*DR
           GWL = E(WTL)+DZ(WTL)/2D0+x*DZ(WTL-1)
           c(NM-1)=P(3)/DZ(NM-1)
           b(NM-1)=-(c(NM-1)+a(NM-1))+CP(NM-1)/dt
         END SELECT

         ! Solve EQ
         CALL THOMAS(NZ,N1,NM,a,b,c,d,HTEMP)

         ! Error Convergence
         temp1 = abs(alpha*HTEMP)**n
         temp2 = 1+temp1
         THETEMP = THER + DTHES/(temp2**tn2)
         WHERE (HTEMP > -1e-5) THETEMP = THES
         !temp3=THETEMP(N1:NM) - THEP(N1:NM)
         err = sum(abs(THETEMP - THEP))
        !temp = HTEMP - HP
         h_err = sum(abs(HTEMP(N1:NM) - HP(N1:NM)))
         HP = HTEMP
         THEP = THETEMP
         CP = - DTHES/(temp2 **tn4) * temp1 *(n-1)/HP
         !CP = -(THES-THER)/(temp2 **((2.D0*n-1.D0)/n))                              &
         !&      * temp1 *(n-1)/HP
         IF (BBC(1)%JType .eq. BC_OLD) THEN
            BC_NC = .TRUE.
         ELSE
            BC_NC = .FALSE.
         ENDIF
         BC_OLD = BBC(1)%JType
       ! NOTES: HERE NOTICE HP or HTEMP IS USED!!!!

         t1 = (h(1)-h(2))/DDZ(1)+1.0D0 ! wrong for passed in JType1=1
         t2 = K(1)
         if (t1 > 0) then
           K(1) = PRatio*KS(2)
         else
           K(1) = K(2)
         endif

         IF (abs(K(1)-t2)>1d-10) THEN
           BC_NC = .FALSE.
           Khalf(1)=0.5D0*(K(1)+K(2))
           !Khalf(1) = sqrt(K(1) * K(2))
         ENDIF

       ENDDO
       ! ##############222222222222222222222222222222222222222222222222222222222222

       IF (ITER >= MAXIT .AND. dt>dtmin) THEN ! Adjust Time Step***********
         dt = max(dtmin,dt/3.0D0)
         SDT = dt
         !IF (BBC(1)%JType .eq. 2 .AND. MMASS) h(1)=BBC(1)%BData2
         ROUTE(1)=ROUTE(1)+1
         CVG = .FALSE.
         !IF (MMass) h(1:NA) = h0(1:NA) ! return to old states 'coz may have been changed by DFMX step
       ELSEIF (ITER < MAXIT) THEN
         !IF (ITER .eq. 50) THEN
         !print *,3
         !ENDIF
         CVG = .TRUE.
       !BV(4)*((h(NA)-BV(5))/BV(2)+1)*dt
         t = t + dt
         !DM = DM + sum((THEP(2:NM)-THE(2:NM))*DZ(2:NM))  ! temporary, check mass balance
         !THE_SAVE=THE; H_SAVE=h
         !DTHE = DTHE + sum((THEP(2:NM)-THE(2:NM))*DZ(2:NM))
         THE = THEP
         h = HP
         CI = CP
         state = 0
         SELECT CASE(BBC(1)%JType)
         CASE(2)
           ! GGA BBC should be differently treated
           IF (JType1 .eq. 5) THEN
             !h(1) = 0D0 ! If marched a step using flux BC, this ponding depth should be 0
             DF = DF + BBC(1)%BData*dt
             BBC(1)%BData= 0.0D0 ! Free surface type BC, will drain in one time step
           ENDIF
         CASE(1)
           HH = BBC(1)%BData - E(1)
           DF = DF + Khalf(1)*t1*dt  ! Don't forget gravity ! Don't forget the BC
         CASE(3)
           DF = DF + Khalf(1)*((h(1)-h(2))/DDZ(1)+1.0D0)*dt
           IF (h(1)<ho) THEN
             !Rf = 0D0
           ELSE
             Rf = Rf + (h(1)-ho)*qm*dt ! partial linearization
             !Rf = qm*dt
           ENDIF
         END SELECT
         SELECT CASE(BBC(2)%JType)
         CASE(4)
         ! [1 ST,2 zl (=DZ(uc)/2+D(aquitard)),3 DR, 4, Kl, 5, hl(pressure head), 6 E(center) uc aquifer, 7,GW.DZ]
         ! [1 ST,2 DZ(bottom cell),3 DR, 4, Kl, 5, hl(pressure head), 6 E(center) c aquifer, 7 DR2,8 EB uc]
         IF (BV(4) > PREC) THEN
         BC4 = P(2)/(BC2*P(1))+h(NM)/P(1)
         Dperc2 = BV(4)*((h(NA)-BC4)/DZB+1D0)*dt
         ELSE
         Dperc2 = 0D0
         ENDIF
         Dperc = Dperc + Dperc2 - BV(3)*dt
         !Dperc=Dperc+Khalf(NA-1)*                                                & ! This is wrong for mass balance
         !&  ((h(NA-1)-h(NA))/BC3+1.0D0)*dt !Dperc2+(BV(3)-DR*BV(2))*dt
         CASE(3)
         !IF (BBC(2)%JType .eq. 3) THEN
         !m_err2 = hh2 - DF + SRC(1)*t - h(1)
         !Dperc=Dperc+Khalf(NA-1)*                                                &
         !&  ((h(NA-1)-h(NA))/DZB+1.0D0)*dt ! NA-1 is to accomodate BBC=3
         Dperc=Dperc+Khalf(NA-1)*                                                &
     &  ((h(NA-1)-h(NA))/DZB+1.0D0)*dt-(BV(3)-DR*BV(7))*dt ! NA-1 is to accomodate BBC=3
         !   Dperc=DF-sum((THE(2:NMAX)-THEOLD(2:NMAX)-                            &
         !&           S(2:NMAX)*dt)*DZ(2:NMAX))+Dperc
         CASE(2)
         !ELSEIF (BBC(2)%JType .eq. 2) THEN
           Dperc = Dperc + BBC(2)%BData*dt
           !print *, Dperc,BBC(2)%BData,dt,qbot,qsur
         CASE DEFAULT
         !ELSE
           Dperc=Dperc+Khalf(NA-1)*                                                &
     &  ((h(NA-1)-BBC(2)%BData)/DDZ(NA-1)+1.0D0)*dt
         END SELECT
         ! Adjust Time Step
         IF (abs(dt-RDT)<1.d-10) THEN
         ELSE
           IF (ITER<=3) THEN
             dt = dt * 1.25D0
           ELSEIF (ITER>5) THEN ! use to be 4 which works
             dt = dt/1.25D0
           ENDIF
           SDT = dt
         ENDIF
         ROUTE(2)=ROUTE(2)+1
       ELSE
         CVG = .TRUE.
         IF (JType1 .eq. 5) THEN
         CALL BOTDRAIN(nz,NA,h,THE,THES,THER,KS,KHALF,ALPHA,N,                   &
     &            LAMBDA,S,DZ,DDZ,BBC,DFF,DPercc,dt,E,DZB,BV)
         DF = DF + DFF
         ! may consider runoff
         Dperc = Dperc + Dpercc
         Rm(1) = Rm(1) + dt
         endif
         IF (h(1)>ho) THEN
          DFF= (h(1)-ho)*(1-exp(-SN*max(0D0,h(1)**(2/3D0))*dt))
          Rf = Rf + DFF
          h(1)=h(1)-DFF
         ENDIF
         ROUTE(3)=ROUTE(3)+1
         STATE = 0
         t = t + dt

       ENDIF !  IF (ITER >= MAXIT) THEN************************
       ddt = SDT
       !MMASS = .FALSE.
 100   CONTINUE
       ENDDO !  TIME LOOP

       Rff = Rff + Rf
       h(1) = h(1)+hh3
       THE(NA) = THEB
       if (h(1)<-1D-11 .AND. JType1 .eq. 5) THEN
         CALL h1negative(nz,na-1,h,THE,THER,THES,DZ,DF,ALPHA,LAMBDA,n)
         ROUTE(6)=ROUTE(6)+1
         K(1) = 0D0 ! debugging info
         K(2) = 1D0
       ENDIF
       CALL waterTable(nm,h,THE,THES,DZ,0.95D0,WTL,x,BV(2))
       GWL = E(WTL)+DZ(WTL)/2D0+x*DZ(WTL-1)
       !K
       !error indices
       !mask this portion out after debugging is complete
       IF (DEBUG) THEN
       if (JType1 .eq. 5) then
         dff = df
       else
         dff = 0D0 ! only for GGA sealed top
         hh2 = BBC(1)%JIndex ! h(1)_n in fact
       endif
       if (abs(BBC(2)%JINDEX)>PREC) THEN
         t2 =0D0
       else
         t2 = h(1)
       endif
       m_err = DTHE+h(1)-Rf-hh1-(DFF-Dperc+sum(SRC(2:NM)*DZ(2:NM)))
       THE2 = sum(THE(2:NA-1)*DZ(2:NA-1))
       DTHE = THE2 - THE0
       t1 = sum(SRC(2:NAA-1)*DZ(2:NAA-1))*dtP(3)
       m_err = DTHE-(DFF-Dperc+t1)
       m_err2 = hh2 - DF + ss1*dtP(3) - t2 - Rf ! IF BBC(1)=2 prescribed this need to add that as well.
       K(5) = m_err ! print debug
       K(6) = m_err2
       IF (abs(m_err)>ERRP*0.1 .OR. abs(m_err2)>1e-9                             &
       !IF (abs(m_err)>ERRP*0.1                                               &
     &       .OR. isnan(m_err)) THEN
       !print *,'m_err VDZ=',m_err
       endif
       ENDIF
       ! print debug
       ! record for output
       !PRatio = m_err
       !K(1) = DF
       !K(2) = Dperc
       !K(3) = t1
       !K(4) = THE0
       !K(6) = SS1
       !K(7) = DF
       IF (ROUTE(1)>0) THEN
         K(1)=0D0
         K(2)=2D0;
       ENDIF
       !IF (BBC(1)%JType .eq. 3) THEN
       !m_err=DTHE+h(1)-hh2-(sum(SRC(2:NM)*DZ(2:NM))+hh1)*dtP(3)+Dperc
       !if (abs(m_err)>1e-5) THEN
       !h(1) = h(1) - m_err
       !print *,'m_err=',m_err
       !endif
       !endif
       !########1111111111111111111111111111111111111111111111111111111111
       ! Compute Recharge Flux to saturated zone
       BBC(1)%BDATA2=0D0
       if (minval(h)>=0D0) THEN
         Rm(2) = 1D0
       else
         Rm(2) = 0D0
       endif
       !R = GWL ! R is not used, here used to store GWL

       END SUBROUTINE vdz1c
       ! #################################################################### END SUBROUTINE
      SUBROUTINE Runoff(RF,frac,ho,hc,dt,hI,h,VParm,nzl,hback,imp,typ)
      !  Computes runoff from impervious region and
       implicit none
       integer nzl
       real*8 RF,fk,ho,h(nzl),dt,hI,hc,hback,imp
       real*8 VParm(*)
       REAL*8 frac,f,k,sep,f2,fI,h0,typ
       ! typ<0 -- only calculate run-on, not going to do runoff from h(1) to OVN.h (do it in vdz1c)
       !parameter(sep = 0.5D0)
       ! function to calculate runoff
       ! hc is h in flow domain, hI is h in impervious ponding, h is vadose zone h
       ! ho is initial abstraction of ponding layer
       ! hback is the barrier height for flow domain to flow back (an offset)
       ! which should be the difference between average elevation and bottom elevation of the cell
       !RETURN
       ! imp is the impervious cover (fraction). higher the imp, higher direct runoff (higher f2)
       !Rf = 0D0
       !k = SN*(VParm(10)**(2D0/3D0))
       !frac = exp(-k*dt)

       IF (frac <= 0D0) RETURN

       if (hI > 0D0) then
         fI = max((1D0 - frac)*hI,0D0)
         f2 = fI*imp ! fraction that goes to channel filaments
         !hc= hc + f2
         h(1)=h(1)+(fI-f2)
         hI = hI - fI
         Rf = Rf + f2
       endif

       if (h(1)> ho .and. h(1)>hc-hback .and. typ>0D0) then
          h0 = max(h(1)-max(hc-hback,0D0),0D0)
          f = max((1D0 - frac)*h(1),0D0)
          f = min(f,(h(1)-hc+hback)/2D0,h(1)-ho)
          !hc= hc + f  ! g.OVN.h uses Rf as a source term, so, not adding it here
          h(1)=h(1)-f
          Rf = Rf + f
       elseif (h(1)<hc-hback .and. hc>hback) then
          ! allow flowing back from flow domain to ponding domain??
          ! comment out if not allowed
          ! run-on
          f = 0.5D0*(hc - hback + h(1))
          h(1) = h(1) + f
          !hc = hc - f
          Rf = Rf - f
       endif

       !R(2) =0D0
       END SUBROUTINE Runoff

       SUBROUTINE Lowland(RF,hd,hff,hI,hG,dt,VParm,dP,imp,m)
       implicit none
       integer nzl,ICASE,ICASE2
       real*8 RF,fk,ho,dt,hI,hc,hg,hd,hff,hback,imp
       real*8,Target:: VParm(10),dP(11)
       REAL*8 frac,f,k,sep,f2,fI,h0,typ, ht,PREC,hG0
       REAL*8 hN,hgN,m1,m2,m,hf0,mratio,g23,m_err
       REAL*8,DIMENSION(:),POINTER :: gg
       integer ROUTE(5)

       ! Lowland is the subroutine for depression storage
       ! treats overland flow as first aggregating to depression storage of the cell
       ! which will then flow out, if exceeding capacity of the depression storage
       ! This is similar to original flow domain, but with a storage and exchange with groundwater
       ! dP contains parameters of the depression storage which are:
       ! [1.hback, 2.hod, 3.Kdg, 4.dfrac, 5.gwST, 6.mr, 7.g1, 8.g2, 9.g3, 10.EB(GW), 11.Ed]
       ! hff is the volumetric depth of flow for the flow domain (subtracting hc by hod)
       ! hod is the storage depth of the depression storage which is related to topography and land use
       ! Kdg is the leakance between local depression storage and groundwater, related to landuse and soil types
       ! dfrac is the fraction of depression storage in the cell
       ! mr is the thickness of leakance bed
       ! g1 = K*dt/mr, g2=g1*Ad/(Ac*S), g3=g1/(1+g2)
       ! ICASE = 0

       ! typ<0 -- only calculate run-on, not going to do runoff from h(1) to OVN.h (do it in vdz1c)
       !parameter(sep = 0.5D0)
       ! function to calculate runoff
       ! hc is h in flow domain, hI is h in impervious ponding, h is vadose zone h
       ! ho is initial abstraction of ponding layer
       ! hG is the mean groundwater elevation in the cell
       ! hback is the barrier height for flow domain to flow back (an offset)
       ! which should be the difference between average elevation and bottom elevation of the cell
       !RETURN
       ! imp is the impervious cover (fraction). higher the imp, higher direct runoff (higher f2)
       !Rf = 0D0
       !k = SN*(VParm(10)**(2D0/3D0))
       !frac = exp(-k*dt)

       PREC = 1D-14
       ! Input is Rf
       ! Calculate hc, hff
       ! If hff is negative, don't make it 0, but keep mass balance
       !h0 = hff
       !ht = hd + Rf   first do groundwater exchange
       ! reason: otherwise this is too fast. vdz, runoff, storage and exchange all happen in the same time step.
       gg => dP(7:9)
       !mratio = 0.5D0

       ICASE = 2
       ROUTE = 0 ! debug

       if (hff<0D0) THEN
         hf0 = hff ! save this negative number
         ht = hd
       else
         hf0 = 0D0
         ht = hd + hff/dP(4) ! add moving water back
       endif
       IF (ht <= PREC .AND. hg<dP(11)) THEN
         NULLIFY(gg)
         RETURN ! no exchange shall happen
       ENDIF

       hG0 = hg
       g23=gg(2)*gg(3)


       IF (hg<dP(11)-dP(6)-1) THEN
         hN = (ht-gg(1)*dP(6))/(1D0+gg(1))
       ELSE
       IF (ICASE .eq. 1) THEN
         ! explicit for gw
         hN = 1D0/(1D0+gg(1))*(ht+gg(1)*hg)
       ELSEIF (ICASE .eq. 2) THEN
         ! implicit for gw
         hN = (ht+gg(3)*hg+(g23-gg(1))*dP(11))/(1+gg(1)-g23)
       ENDIF
       ENDIF
       hgN = (hg+gg(2)*hN+gg(2)*dP(11))/(1D0+gg(2))

       ! limitations
       ! hN>0; hg>=EB; if hgN<E-mr, use other formulation
       if (hgN >= hG) THEN ! recharge
         hN = max(hN,0D0)
         m  = dP(4)*(hN - ht) ! should be negative
         hg = hG - m/dP(5)
         ROUTE(1)=1
       else ! exfiltration
         hgN = max(dP(10)+0.1D0,hgN)
         m  = dP(5)*(hG - hgN)
         hN = ht + m/dP(4)
         hG = hGN
         ROUTE(2)=1
       endif

       m_err = (hN-ht)*dP(4)+(hG-hG0)*DP(5)
       if (abs(m_err)>1D-8) THEN
       print *, 3 ! debug
       endif

       hff = max(hN-dP(2),0D0)
       hd  = hN - hff
       hff = hff*dP(4)+hf0
       NULLIFY(gg)

       ! If want to allow it flow one step first
       ! then do nothing here and add Rf to Ovn.h in ovn_cv
       ! If water is first put into depression storage, then
       ! add it up here

       !R(2) =0D0
       END SUBROUTINE Lowland


       SUBROUTINE waterTable(nz,h,THE,THES,DZ,THRSD,WTL,x,DZB)
       ! recognize the position of the water table, both water table layer (WTL) and position (GWL)
       ! WTL: the first fully saturated layer (water table higher than its upper edge)
       implicit none
       integer nz, WTL
       real*8 h(nz),x,the(nz),thrsd,THES(nz),DZ(nz),x0,DZB
       integer I,J

       ! tried to enable the capillary zone with THE criterion, definitely too large K for the Vauclin case

       ! x is the fraction of the first unsaturated that is saturated from linear interpolation in layer WTL-1

       if (WTL >= nz) WTL = nz-1
       DO WHILE (WTL<nz .AND. (.NOT. isnan(DZ(WTL+1))) ) ! bottom layer is groundwaters
         IF (h(WTL)<0D0) THEN
         !IF (THE(WTL)<THES(WTL)*THRSD) THEN
            WTL = WTL + 1
         ELSE
            EXIT
         ENDIF
       ENDDO
       DO WHILE (WTL>2) ! top layer is ground surface
         IF (h(WTL-1)>0D0) THEN
         !IF (THE(WTL)>THES(WTL)*THRSD) THEN
            WTL = WTL - 1
         ELSE
            EXIT
         ENDIF
       END DO
       if (WTL .eq. 2) THEN
         x = 1.D0
       elseif (WTL .eq. nz) THEN
         !x = max((-h(WTL)/(h(WTL-1)-h(WTL))-0.5D0),0D0)
         x = -h(WTL)*(DZB+DZ(WTL-1))/(h(WTL-1)-h(WTL)) 
         if (x>=DZB) THEN
           x = 0.5D0*(x  - DZB)/DZ(WTL-1)
         else
           x = 0D0
         endif
       else
         !x = -h(WTL)*(DZ(WTL)+DZ(WTL-1))/(2D0*(h(WTL-1)-h(WTL)))
         x = -h(WTL)*(DZ(WTL)+DZ(WTL-1))/(h(WTL-1)-h(WTL)) ! there is a factor of 2 but that can be avoided, faster
         !x = (THE(WTL+1)-THES(WTL)*THRSD)/(THE(WTL+1)-THE(WTL))
         if (x>=DZ(WTL)-1D-15) THEN
           x = max(0.5D0*(x  - DZ(WTL))/DZ(WTL-1),0D0)
         else
           x = 0.5D0*x/DZ(WTL) + 0.5D0
           WTL = WTL + 1
         endif
       endif
       if (isnan(x)) x = 0D0

       END SUBROUTINE waterTable

       SUBROUTINE BOTDRAIN(nz,NA,h,THE,THES,THER,KS,KHALF,ALPHA,N,            &
     &            LAMBDA,S,DZ,DDZ,BBC,DDF,DPerc,dt,E,DZB,BV)
       USE Vdata
       IMPLICIT NONE
       integer nz,na
       REAL*8,DIMENSION(nz):: h,THE,THES,KS,ALPHA,LAMBDA,KHALF,S,DDZ
       REAL*8,DIMENSION(nz):: THER,DZ,E,N
       REAL*8 DDF,Dperc,dt,m,RF,DZB,BV(*),C
       TYPE(BC_type) :: BBC(2)
       ! This is a linearized fail-safe bottom drainage model that RE can fall back on
       ! it is simplified as it only means to help the RE solver (occasionally) in case of convergence problems
       ! It computes flux from bottom up, the flux is the min of
       ! 1.explicit flux 2.saturation of receiving cell and 3.available moisture of providing cell
       ! may be less accurate, but more stable
       REAL*8 f1,f2,f3,KU,Q,FU,FL,THE0,THE1,test,SS,test2,SSS,hh2
       real*8 tt(5)
       integer I,J,NM,RC
       NM = NA-1
       hh2 = h(1)
       !h(NA) = BBC(2)%Bdata - E(NA)
       THE0 = sum(THE(2:NA)*DZ(2:NA))
       DDF = 0D0
       Dperc  = 0D0; RF= 0D0
       C = 0.4D0

       ! now we use an implicit formulation for botdrain for stability
       ! we assume h(NM-1) not changing, use implicit method
       ! [1 ST,2 dzl (=DZ(uc)/2+D(aquitard)),3 DR, 4, Kl, 5, hl(pressure head), 6 E(center) uc aquifer]

       tt(1)=BV(1)/dt; tt(2)=Khalf(NA-1)/DZB;
       tt(3)=tt(1)+tt(2)
       tt(4)=tt(2)*h(NA-1)+BV(3)+tt(1)*h(NA)+Khalf(NA-1)
       h(NA)=tt(4)/tt(3)

       ! Bottom Flux
       I = NM
       f1 = Khalf(NM)*((h(NM)-h(NA))/DZB+1)
       IF (f1>0) THEN
         f2 = 10000
         f3 = C*(THE(NM)-THER(NM))*DZ(NM)/dt
         FU = min(f1,f2,f3)
       ELSE
         f2 = 10000
         f3 = 10000
         FU = max(f1,-f2,-f3)
       ENDIF

       Dperc = FU*dt


       DO I = NM,3,-1
         FL = FU; FU = 0D0
         Q = 0D0
         THE(I)=THE(I)-FL*dt/DZ(I)+S(I)*dt
         IF (THE(I)>THES(I)) THEN
           Q = (THE(I)-THES(I))*DZ(I)/dt
           THE(I) = THES(I)
           FU = -Q ! flux should go up
         ELSE
           !IF (FL > 0D0) THEN
           ! cell already gives water to below, be careful, don't provide up too much
             SS=(THE(I)-THER(I))/(THES(I)-THER(I))
             m = 1D0 - 1D0/n(I);
             h(I)=((SS**(-1.D0/m)-1D0)**(1.D0/n(I)))/(-abs(alpha(I)))
           !ENDIF
           f1 = Khalf(I-1)*((h(I-1)-h(I))/DDZ(I-1)+1) ! when alpha is large this can be erraneous
           !IF (abs(ALPHA(I))>2D0)
           if (f1<0D0) THEN
             ss = 0.5D0*(THE(I)+THE(I-1))
             f1 = max(f1,-(ss - THE(I-1))*DZ(I-1)/dt)
           endif
           if (f1>0D0) THEN ! .or. FL>0D0) THEN
           ! FL>0 condition is to prevent cell from giving water to both directions
            if (F1<0D0) f1 = 10000
             !f2 = (THES(I)-THE(I))*DZ(I-1)/dt
             !f3 = (THE(I-1)-THER(I-1))*DZ(I)*C/dt
             f2 = (THES(I)-THE(I))*DZ(I)/dt
             f3 = (THE(I-1)-THER(I-1))*DZ(I-1)*C/dt
             FU = min(f1,f2,f3)
           else
             !f2 = (THES(I-1)-THE(I-1))*DZ(I-1)/dt
             f2 = 100000D0
             f3 = (THE(I)-THER(I))*DZ(I)*C/dt
             FU = max(f1,-f2,-f3)
           endif
           THE(I) = THE(I)+FU*dt/DZ(I)
         ENDIF
       ENDDO

       ! Top Cell
       FL = FU; I=2;
       THE(2) = THE(2)-FL*dt/DZ(2)+S(I)*dt

       IF (THE(2)>=THES(2)) THEN
         Q = (THE(I)-THES(I))*DZ(I)/dt
         THE(2) = THES(2)
         !h(1) = h(1) + Q*dt + S(1)*dt
         IF (BBC(1)%JType .eq. 2) THEN
           h(1) = h(1) + Q*dt + BBC(1)%Bdata*dt
         ELSE
           h(1) = h(1) + Q*dt + S(1)*dt
         ENDIF
         DDF = -Q*dt
       ELSE
          IF (BBC(1)%JType .eq. 2) THEN ! flux
           Q = 0D0
           THE(2)=THE(2)+(BBC(1)%BData)*dt/DZ(2)
           !h(1)=0D0
           IF (THE(2)>THES(2)) THEN ! impossible to squeeze the water inside the column
             h(1) = h(1)+(THE(2)-THES(2))*DZ(2) + S(1)*dt
             Q = (THE(I)-THES(I))*DZ(I)/dt
             THE(2)=THES(2)
           ENDIF
           DDF = DDF + BBC(1)%BData*dt - Q*dt !!!!!!!!!!PROBLEM WITH MASS BALANCE!!!!
          ELSEIF (BBC(1)%JType .eq. 3) THEN
           IF (FL > 0D0) THEN
           ! cell already gives water to below, be careful, don't provide up too much
             SS=(THE(I)-THER(I))/(THES(I)-THER(I))
             m = 1D0 - 1D0/n(I);
             h(I)=((SS**(-1.D0/m)-1D0)**(1.D0/n(I)))/(-abs(alpha(I)))
           ENDIF
           f1 = Khalf(1)*((h(1)-h(2))/DDZ(2)+1)
           if (f1>0D0) THEN ! .or. FL>0D0) THEN
           ! FL>0 condition is to prevent cell from giving water to both directions
             if (f1<0D0) f1=100000
             f2 = (THES(2)-THE(2))*DZ(2)/dt
             f3 = h(1)/dt + S(1)
             FU = min(f1,f2,f3)
           else
             !f2 = (THES(I-1)-THE(I-1))*DZ(I-1)/dt
             f2 = 100000D0
             f3 = (THE(2)-THER(2))*DZ(2)*0.8D0/dt
             FU = max(f1,-f2,-f3)
           endif
           DDF = FU*dt
           h(1)=h(1)-FU*dt + S(1)*dt
           THE(2)=THE(2)+DDF/DZ(2)
         ENDIF
       ENDIF
       ! Update h using THE
       DO I = 2,NM
         IF (THE(I)<THES(I)) THEN
           SS=(THE(I)-THER(I))/(THES(I)-THER(I))
           m = 1D0 - 1D0/n(I);
           h(I)=((SS**(-1.D0/m)-1D0)**(1.D0/n(I)))/(-abs(alpha(I)))
         ELSEIF (h(I-1)>=0D0) THEN
           h(I)=h(I-1)+DDZ(I-1)
         ELSE
           h(I)=0D0
         ENDIF
       ENDDO
       IF (BBC(1)%Jtype .eq. 2) RF = BBC(1)%Bdata*dt
       THE1 = sum(THE(2:NA)*DZ(2:NA))
       SSS =  sum(S(2:NA)*DZ(2:NA))*dt
       !test = THE1 + h(1)- (THE0 + hh2+S(1)*dt+SSS                        &
       !&       - Dperc + RF)
       !test2 = THE1 -(THE0+SSS+DDF-Dperc)
       !IF (max(abs(test),abs(test2))>1d-7) THEN
       !print*,3
       !ENDIF
       END SUBROUTINE

       SUBROUTINE satExcess(p1,p2,h1,h1N,SN,ho)
       ! picard iteration
       REAL*8 p1,p2,h1,hN,h1N,hNN,SN,ho,ERRP,X
       integer I,IMAX
       REAL*8 Fg,ERR
       ERRP = 1D-10
       h1N = h1
       hNN = hN
       ERR = 1000
       I = 1
       IMAX = 20
       DO WHILE (abs(ERR)>ERRP .and. I<=IMAX)
       IF (h1N>ho) THEN
         Fg = SN*(h1N-ho)**(5D0/3D0)
       ELSE
         Fg = 0D0
       ENDIF
       X = p1*(p2-Fg)
       ERR = X - h1N
       h1N = h1N + 0.5*ERR
       I = I + 1
       ENDDO
       IF (I > IMAX) h1N = -1 ! failed

       END SUBROUTINE satExcess

       SUBROUTINE satExcessFunc(h1N,f,DDF,SN,ho,p1,p2,d5,d6)
       !CALL funcd(rtnobd,f,DDF,INP1,INP2,INP3,INP4,INP5,INP6)
       !Ft = rtsafe(curveF,0.0D0,Fr,PREC, t, Am, Bm, Fex,km,0.0D0)
       implicit none
       real*8 h1N,f,DDF,SN,ho,p1,p2,d5,d6
       real*8 Fg,Dg
       IF (h1N>ho) THEN
         Fg = SN*(h1N-ho)**(5D0/3D0)
         Dg = (h1N-ho)**(2D0/3D0)
       ELSE
         Fg = 0D0
         Dg = 0D0
       ENDIF
       f = p1*(p2-Fg) - h1N
       DDF = -(5D0*SN*p1/3D0)*Dg-1
       !                        5/3
       !f = p1 (p2 - SN (h - ho)   ) - h
       !                        2/3
       !        5 SN p1 (h - ho)
       !DDF = - ------------------- - 1
       !                3
       END SUBROUTINE satExcessFunc

       SUBROUTINE ovn_cv(m,h,U,V,S,d,E,Ex,Ey,Mann,ddt)
       USE Vdata
       IMPLICIT NONE
       integer m(2),nG(2)
       REAL*8 d(2),ddt
       REAL*8,DIMENSION(m(1),m(2)),TARGET::h,E,Ex,Ey,Mann,U,V,S
       REAL*8,DIMENSION(m(1),m(2),0:2),TARGET::dh
       REAL*8,DIMENSION(m(1),m(2)),TARGET::hx,hy,Eta,hnew,h1,hN
       REAL*8,DIMENSION(m(1),m(2)) :: fhuX,fhvY,UU,VV
       LOGICAL,DIMENSION(m(1),m(2)) :: MM
       REAL*8,POINTER,DIMENSION(:,:):: dhh
       REAL*8 t,rdt,dt,ft,PREC,minDt
       REAL*8 temp,sgn,dx,dy,t23,t12
       LOGICAL:: DRY
       integer i,j,k,n
       ! hN is input, hNew is output (after RK3)
       t = 0D0
       rdt = ddt
       ft = ddt
       dt = ddt*86400D0
       minDt = dt/8D0
       PREC = 1D-10
       dy = d(1); dx = d(2)
       DRY = .TRUE.
       t23 = 2D0/3D0
       t12 = 0.5D0

       do while (t<ft-PREC)
        dh = 0D0;
        rdt = ft - t; dt = min(rdt*86400D0,dt)
        ! whole function
        ! #################################
        do n = 0,2
          dhh => dh(:,:,n)
          if (n .eq. 0) then
            hN = h
          elseif (n .eq. 1) then
            hN = h + dt*dh(:,:,0)
          elseif (n .eq. 2) then
            hN = h + 0.25D0*dt*dh(:,:,0) + 0.25D0*dt*dh(:,:,1)
          endif
          ! input to one_step function hN,output dhh
          ! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@
          ! notice memory effect - remnant value unchanged from last operation
          MM = .FALSE.
          DO i = 2,m(1)-1
          DO j = 2,m(2)-1 ! not going across boundary
            IF (hN(I,J)>PREC .OR. S(I,J)>0D0) THEN
              MM(I,J) = .TRUE.
              MM(I-1,J) = .TRUE.;MM(I,J-1) = .TRUE.;
              MM(I+1,J) = .TRUE.;MM(I,J+1) = .TRUE.;
              ! neighbors. the interfaces of neighbors don't matter except the one to this cell
              ! (other interfaces are all zero)
              DRY = .FALSE.
            ENDIF
          ENDDO
          ENDDO
          !MM = .TRUE.; DRY = .FALSE. ! debug. verified that with this statement it is exactly
          ! the same as the matlab code

          IF (DRY) THEN
            ! surface is totally dry
            ! let's get out of here
            RETURN
          ENDIF

          Eta = 0D0
          hx =0D0; hy = 0D0
          fhuX=0D0; fhvY = 0D0
          Eta = hN + E
          !UU = 0D0; VV = 0D0
          CALL interfaceH(m(1),m(2),hN,Eta,Ex,Ey,hx,hy,MM)
          ! assume no boundary condition needs to be applied
          ! use diffusive wave
          DO i = 2,m(1)-1
          DO j = 2,m(2)-1
            IF (MM(I,J)) THEN
              temp = (Eta(i,j)-Eta(i,j+1))/dx
              if (temp>0) then
                sgn = 1D0
              else
                sgn = -1D0
              endif
          UU(i,j)=(sgn/mann(i,j))*(hx(i,j)**t23)*(abs(temp)**t12)
              temp = (Eta(i,j)-Eta(i+1,j))/dy
              if (temp>0) then
                sgn = 1D0
              else
                sgn = -1D0
              endif
          VV(i,j)=(sgn/mann(i,j))*(hy(i,j)**t23)*(abs(temp)**t12)
              fhux(I,J)=UU(i,j)*hx(i,j)
              fhvy(I,J)=VV(i,j)*hy(i,j)
            ENDIF
          ENDDO
          ENDDO
          fhux(:,m(2)-1) = 0D0;
          fhvy(m(1)-1,:) = 0D0;

          DO i = 2,m(1)-1
          DO j = 2,m(2)-1
            IF (MM(I,J)) THEN
         dhh(I,J)=-(fhuX(I,J)-fhuX(I,J-1))/dx-(fhvY(I,J)-fhvY(I-1,J))/dy
              dhh(I,J)=dhh(I,J)+S(I,J)
            ENDIF
          ENDDO
          ENDDO
          ! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        enddo
        ! #################################
        hNew = h +dt/6D0*(dh(:,:,0)+4D0*dh(:,:,2)+dh(:,:,1))

        if (minval(hnew)<-0.3D0 .and. dt>mindt) THEN
         dt = min(rdt*86400D0,max(dt/2d0,mindt))
        else
         t = t + dt/86400D0
         !h1 = h ! print debug
         h = hnew
        endif
       enddo
       t = ft ! avoid runoff error
       NULLIFY(dhh)
       END SUBROUTINE ovn_cv





       SUBROUTINE h1negative(nz,na,h,THE,THER,THES,DZ,DF,                 &
     &             ALPHA,LAMBDA,n)
       ! this subroutine tries to fill up upper cell hh3<0
       integer nz,na
       REAL*8,DIMENSION(nz)::h,THE,THES,DZ,THER,ALPHA,                      &
     &             LAMBDA,n
       integer I,J,K
       REAL*8 TH,PREC,mm,SS,DF,C

       I = 2
       C = 0.25D0 ! how much you can extract
       DO WHILE (I <= na)
         TH = min((THE(I)-THER(I))*DZ(I)*C,-h(1))
         h(1) = h(1)+TH
         DF = DF - TH ! too much infiltration, reduce correspondingly
         THE(I) = THE(I)-TH/DZ(I)
         mm = 1D0-1/(n(I))
         SS = (THE(I)-THER(I))/(THES(I)-THER(I))
         h(I) = ((SS**(-1D0/mm)-1)**(1D0/n(I)))/(-abs(ALPHA(I)))
         IF (h(1)>=0D0) EXIT
         I = I+1
       ENDDO


       END SUBROUTINE h1negative

       !!!!! GENERALIZED GREEN AND AMPT METHOD with source term
       SUBROUTINE GGA(NZ,WFP,DF,THE,THES,THER,KS,h,DZ,P,                         &
     &            S,E,KBar,FM,dt,alpha,N,LAMBDA,FL,Code,THEW)
       ! Solve 1D GGA
       USE Vdata
       USE locfunc
       IMPLICIT NONE
       integer :: NZ
       REAL*8, DIMENSION(NZ) :: THE,THES,THER,KS,DZ,E,KBar,FM,h
       REAL*8, DIMENSION(NZ) :: ALPHA,N,LAMBDA,S
       REAL*8 :: prcp,DF,DPerc,FM1,dt,WFP,P,FL,Code
       REAL*8 ::THEW
       !WFP wetting front position, = n.x means WF is in the n-th layer, x fraction from above
       REAL*8 :: h0 ! TEMPORARY
       integer :: jt,ROUTE(10)
       ! Input.
       ! THE, THES: Current and saturation moisture content
       ! KS: Saturated Conductivity
       ! jt: soil types
       ! WFP :: Wetting Front Location (n-th layer, the distance from layer Upper Boundary to Wetting front)
       ! Output
       ! DF: cumulative infiltration amount in this time step
       ! SWF: Suction.. It should be a function, but may also be an array input
       ! BC1: sum(L_i*k_m/k_i), precomputed
       ! Kbar: harmonic mean of conductivities
       ! FM: sum of available pore space from top down, sum(L_i*(THES-THER)_i)
       REAL*8 :: SW, tr, Fr, LEN, MxF ,BC1, x,kmm,WFP0
       REAL*8 :: t, Am, Bm, Fp, tp, tf, tp1, tp2, km, dTHE ! TEMPORARY
       REAL*8 :: L, IFC, Ft, Fex, IFSlow, temp1,temp2,err
       ! Wetting Front Length, Infiltration Capacity, Cumulative Infiltration
       integer :: IWFP, ICASE, K, IU,I,J,SCODE
       LOGICAL :: T_F
       REAL*8 :: PREC,FMM,t1,t2,t3,ss,x1,THE0,M1,ME,M2,M3,ME2,hh0
       ! Local Variables
       ! IWFP: layer # of wetting front
       ! LEN, Length from layer upper boundary to wetting front, SW: Suction at wetting front
       ! tr,Fr, remaining time and infiltration amount to end current fractional step
       ! MxF : Maximum Infiltration rate (Infiltration Capacity)
       ! Fex = existing wetting front moisture in the cell
       ! BC: coefficient of B. B = sum(L_i/K_i)
       ! (1) Determine If ponding will happen in this time step

       ! Mass balance:
       ! GA's description of THETA will cause mass balance issues if not properly handled
       ! output THE: cell average THE, but THEW is used to save dry portion THE
       ! If code = 1 when passed in, last time step GGA was also used, thus WFP was correct
       ! We use WFP and THE(IWFP) (cell average theta) to find out the THE in the dry portion
       ! of the IWFP cell
       ! There may be inconsistencies, in which case WFP need to be modified

       ! print debug
       !THE0 = sum(THE*DZ); hh0 = h(1)
       !M1 = sum(DZ*S)*dt
       !M2 = DF

       ! Local Variable Initialization
       dTHE =0.0D0;Am=0.0D0;Bm=0.0D0;
       Fp=0.0D0;
       tp=0.0D0;
       ROUTE = 0
       tf=0.0D0;tp1=0.0D0;tp2=0.0D0;
       km=0.0D0;dTHE=0.0D0;
       L=0.0D0;IFC=0.0D0;Ft=0.0D0;Fex=0.0D0;IFSlow=0.0D0;
       !SW=0.0D0;Fr=0.0D0;LEN=0.0D0;MxF=0.0D0; kmm=0D0
       !temp1=0.0D0; temp2=0.0D0
       PREC = 1.0D-12
       K = 0
       !DF = 0.0D0
       DPerc = 0.0D0
       T_F = .TRUE.

       ! Note that first layer is the ponding layer! Soil starts from the second one!
       IF (Code .eq. 0) THEN ! last time step was not GGA, so WFP cannot be used
         DO I = 2,NZ
           IF (THE(I)>=THES(I)-0.001D0) THEN
           ELSE
             EXIT
           ENDIF
         ENDDO
         IF (I<nz .OR. (I .eq. nz .and. THE(I)<=THES(I)-0.001D0)) THEN
           WFP = DBLE(I-1)
         ELSE
           WFP = nz
         ENDIF
         IWFP = max(FLOOR(WFP+PREC)+1,2) ! WFP = 2.2 means it is in the third layer, 0.2 cell from the top
       ELSE
         ! need to inherit THEW, which is the saved dry half THE in the IWFP cell
         IWFP = max(FLOOR(WFP+PREC)+1,2) ! WFP = 2.2 means it is in the third layer, 0.2 cell from the top
         !THE(IWFP) = THES(IWFP)*x+(1-x)*THE(IWFP)
         x = WFP - (IWFP - 1D0)
         t1 = (THE(IWFP)-THES(IWFP)*x)/(1D0-x)
         if (t1>THES(IWFP) .or. t1<THER(IWFP)) THEN
           ! reposition fix -- this should not happen often
           ! place WFP to fit the average THE
           t1 = max(min(THEW,THE(IWFP)/2D0),(THER(IWFP)+THE(IWFP))/2D0)
           x1 = (THE(IWFP)-t1)/(THES(IWFP)-t1)
           THE(IWFP)=t1
           WFP = (IWFP-1D0) + x1
           ROUTE(5)=ROUTE(5)+1
         else
           THE(IWFP)=t1 ! safely transfer
         endif
       ENDIF

       tr = dt
       !h0 = h(1)
       LEN = DZ(IWFP)*max(WFP+1-IWFP,0D0)
       Fex = LEN * (THES(IWFP) - THE(IWFP))
       WFP0 = WFP

       prcp = max(1D-20,P)

       ! consider source term
       ! source term is assumed to be immediately supplied by rainfall and ponding infiltration
       ! so at each step directly modify rainfall amount
       ! (1) remove h0; (2) reduce prcp
       IF (IWFP>=2) THEN
         CALL GGA_Source(2,IWFP,DZ,S,h,WFP,prcp,dt,dt,DF,SCODE)
         IF (SCODE .eq. 1) RETURN
       ENDIF

       DO WHILE (tr > PREC .AND. IWFP <= NZ .AND. SCODE < 1)

         K = K + 1 ! Number of Fractional Steps
         h0= h(1)
         !if (h(1)<-PREC) THEN
         !print*,'debug'
         !endif
         IU = IWFP - 1
         !IF (T_F) THEN ! Need to obtain new parameters --- also need to be done if tp is reached -- KMM is different
           IF (IWFP <=NZ .AND. ABS(LEN-DZ(IWFP))>PREC) THEN
           ! need to provide: Fex, Len, IWFP, tr
             IF (T_F) THEN ! save some effort if not new layer
             dTHE = THES(IWFP) - THE(IWFP)
             do while (dTHE .eq. 0D0 .and. IWFP <=NZ)
               ! wetting front meet perched water table
               ! move the wetting front
               L = L + DZ(IWFP)
               IWFP = IWFP + 1
               dTHE = THES(IWFP) - THE(IWFP)
           IF (abs(S(IWFP))>PREC .and. h(1)>=0D0) THEN
           CALL GGA_Source(IWFP,IWFP,DZ,S,h,WFP,prcp,dt,tr,DF,SCODE)
           ENDIF
             enddo
             IF (IWFP >NZ .OR. (IWFP .eq. nz .and.                                  &
     &          ABS(LEN-DZ(IWFP))<=PREC)) THEN
               ! exit, update variables

               h(IWFP-1)=LEN + DZ(IWFP-1)/2D0
               DO I=IWFP-2,2,-1
                h(I) = h(I+1) + E(I) - E(I+1)
               ENDDO
               THE(2:IWFP-1)=THES(2:IWFP-1)
               WFP = DBLE(IWFP-1)
               RETURN
             ENDIF

             !SW = -h(IWFP)
             SW = swFunc(h(IWFP),ALPHA(IWFP),N(IWFP),LAMBDA(IWFP))
             SW = max(SW,1D-3)
             km = KS(IWFP)
             ENDIF
             Fr = dTHE*DZ(IWFP) - Fex !!!!!!!!!!!!!!!!
             IF (IWFP .EQ. 2) THEN
              L = DMAX1(LEN,1.0E-20)
              Am = dTHE * SW
              Bm = 0.0D0
              Kmm=km
             ELSE
              L = SUMDZ(IU,E,DZ) + LEN ! SUMDZ(I)=sum(DZ(1:I))
              BC1 = SUMDZ(IU,E,DZ)/Kbar(IU) ! BC1(m)=sum(DZ(1:m)/KS(1:m))
              Kmm = L/(BC1+LEN/km) ! current mean K -- at entrance of a new layer, this is not correct
              !Bm = BC1*km*dTHE - FM(IU) ! FM = sum(L_i*(THES-THER)_i)
              Bm = BC1*km*dTHE ! This is used because now Fex contains only water in this cell.
              Am = dTHE*(SW + SUMDZ(IU,E,DZ) - BC1*km)
             ENDIF
           ELSE IF (IWFP.EQ.NZ .AND. ABS(LEN-DZ(IWFP))<PREC) THEN
           !! FULLY SATURATED
           !! will not come in here
             SW = -h(IWFP)
             temp1 = KS(NZ)*(1 + SW/SUMDZ(NZ,E,DZ))
             Ft = DMIN1(h(1)+prcp*tr, temp1*tr)
             DF = DF + Ft
             DPerc = DPerc + Ft
             MxF= h(1) + tr*prcp - Ft
             h(1) = DMAX1(0.0D0,MxF-PREC)
             RETURN
           ELSE
             PRINT*,4
           ENDIF
         !ENDIF
         IFC = KMM*(1+SW/L)
         ! Compute tp : fractional time to change of infiltrating curve
         IF (prcp > IFC-PREC) THEN
          ! ponding from start to end
           tp = 1e20;
           tf = 1/km*(Fr - Am*log((Am+Bm+DZ(IWFP)*dTHE)/(Am+Bm+Fex)))
           IF (tf<=0) tf = 1e20
           ICASE = 1
         ELSEIF (h(1)<1D-8) THEN
           ! just infiltrate as prcp
           ! but it may naturally pond as F increases
           IF (prcp<=km) THEN
             Fp = 1D20
             tp = 1D20
           ELSE
             Fp = Am/(prcp/km-1)-Bm + FM(IU)
             tp = (Fp - Fex - FM(IU))/prcp
           ENDIF
           tf = Fr / prcp
           ICASE = 2
         ELSE
           ! slightly complicated. h0 may go away in the time step
           ! but prcp may also equal Inf before that happens
           ! before pond goes away, infiltrate as IFC, after that would be prcp
           ! First compute when it shall deplete the pond
           ! The equation is nonlinear (curveT). h0+prcp*tp2 = km*tp2+Am*log((A+B+Fex+h0+tp2*prcp)/(A+B+Fex))
           ! It may have 0,1 or 2 roots
           IF (ABS(prcp - km)<1.0E-10) THEN
              tp2 = (EXP(h(1)/Am)*(Am+Bm+Fex)-(Am+Bm+Fex+h(1)))/prcp
           ELSE
              ! Finding upper bound
              IF (prcp < km) THEN ! line sloping down, find x-axis intercept
                temp1 = h(1)/(km-prcp)
              ELSE ! prcp > km. there may be two roots (or no root), using Rolle's Rule to find middle point, where derivative is zero
           !temp1 = (Am*prcp*(Am+Bm+Fex)/(prcp-km)-(Am+Bm+Fex+h(1)))/prcp
                temp1 = -(Am+Bm+Fex+h(1)+Am*prcp/(km-prcp))/prcp
              ENDIF
              IF (temp1 <0) THEN
                ! inflection point is <0, what is >0 is monotonous
                tp2 = rtsafe(curveT,0.0D0,tr,PREC,h0 ,Am,Bm,                   &
     &               Fex,km,prcp)
              else
                ! inflection point is >0, there may be 0,1,2 roots, first try the smaller segment
                tp2 = rtsafe(curveT,0.0D0,temp1,PREC,h0 ,Am,Bm,                   &
     &               Fex,km,prcp)
                if (tp2<0 .or. isnan(tp2) .or. tp2>temp1) then
                tp2 = rtsafe(curveT,temp1,max(temp1*2D0,tr),PREC,h0 ,                   &
     &               Am,Bm,Fex,km,prcp)
                endif
              endif
              ! there might be root >tr, but we don't care in that case

              if (isnan(tp2)) then
               write(777,'(11F17.13)') Am,Bm,Fex,km,t1,h0,DZ(IWFP),       &
     &               dthe,Fr,prcp,tp2
              tp2 = rtnewt(curveT,temp1,max(temp1*2D0,tr),PREC,h0 ,                      &
     &               Am,Bm,Fex,km,prcp,Err)
              endif ! to use with plotNewton
              if (tp2<0 .or. err>1D-5) then
                tp2 = 1e20
              endif
           ENDIF
           tp = tp2
           !Fp = h0 + tp * prcp ! Infiltration Amount at change-status time
           tf = 1/km*(Fr - Am*log((Am+Bm+DZ(IWFP)*dTHE)/(Am+Bm+Fex)))

           IF (tf<=0) tf = 1e20
           IF (tp<=0) tf = 1e20
           ICASE = 3

         ENDIF

         t = DMIN1(tr,tp,tf)  ! Whichever time comes sooner
         T_F   = (DABS(t - tf)<PREC) ! Filled up cell
         ! Compute Infiltration and new wetting front position
         IF (T_F) THEN
           Ft = Fr
           ROUTE(1)=ROUTE(1)+1
         ELSEIF (t .EQ. tp2) THEN  ! Pond gone
           Ft = h(1) + tp2 * prcp
           ROUTE(2)=ROUTE(2)+1
         ELSEIF (ICASE .EQ. 2) THEN ! f = prcp
           Ft = t * prcp
           ROUTE(3)=ROUTE(3)+1
         ELSE  ! f = IFC
           Ft = rtnewt(curveF,0.0D0,min(IFC*t,Fr),PREC, t, Am, Bm,           &
     &             Fex,km,0.0D0,Err)
           !Ft = rtsafe(curveF,0.0D0,Fr,PREC, t, Am, Bm, Fex,km,0.0D0)
           ROUTE(4)=ROUTE(4)+1
           IF (FT<=0 .or. Ft>h(1) + t*prcp .or. abs(err)>1D-5) then
           Ft = rtsafe(curveF,0.0D0,IFC*t,PREC, t, Am, Bm, Fex,km,1D0)
           ENDIF
           IF (FT<=0 .or. Ft>h(1) + t*prcp .or. isnan(Ft)) then
           !write(999,'(14F17.13)') Am,Bm,Fex,km,t1,h0,DZ(IWFP),             &
           !&               dthe,Fr,prcp,tp2,ft,IFC,t
           Ft = rtsafe(curveF,0.0D0,IFC*t,PREC, t, Am, Bm, Fex,km,1D0)
           Ft = h(1) + t*prcp
           ENDIF
         ENDIF

         t1 = h(1) + t*prcp - Ft

         DF = DF + Ft
         tr = tr - t
         Fr = Fr - Ft
         Fex= DZ(IWFP)*dTHE - Fr
         LEN= LEN + Ft/dTHE
         L = L + Ft/dTHE ! This is for fractional time increment of IFC


         h(1)= t1

         if (t1<0 .or. Fr<0D0 .or. isnan(t1)) then
         ! print debug
          !write(999,'(13F16.13)') Am,Bm,Fex,km,t1,h0,DZ(IWFP),             &
          !&               dthe,Fr,prcp,tp2,ft,IFC
         endif

         IF (T_F .AND. IWFP <NZ) THEN ! Filled up. Proceed to next cell
           IWFP =  IWFP + 1
           LEN = 0.0D0
           Fex = 0.0D0
           ! deal with source term again
           IF (abs(S(IWFP))>PREC .and. h(1)>=0D0) THEN
           CALL GGA_Source(IWFP,IWFP,DZ,S,h,WFP,prcp,dt,tr,DF,SCODE)
           ENDIF
         ENDIF
       ENDDO

       ! Update state variables (h, THE, WPF)
       IF (IWFP>2) THEN
         h(IWFP-1)=LEN + DZ(IWFP-1)/2D0
         DO I=IWFP-2,2,-1
            h(I) = h(I+1) + E(I) - E(I+1)
         ENDDO
         THE(2:IWFP-1)=THES(2:IWFP-1)
       ENDIF


       IF (IWFP>1) THEN
         x = LEN/DZ(IWFP)
         THEW = THE(IWFP)
         THE(IWFP) = THES(IWFP)*x+(1D0-x)*THE(IWFP)
         ! does this transformation ONLY when switching back to RE
         ! remember h(IWFP) needs to be updated outside
         ! update h(IWFP)
         SS = (THE(IWFP)-THER(IWFP))/(THES(IWFP)-THER(IWFP))
         t1 = -N(IWFP)/(N(IWFP)-1D0)
         h(IWFP)=((SS**t1-1D0)**(1D0/n(IWFP)))/(-abs(alpha(IWFP)))
       ENDIF



       WFP = DBLE(IWFP-1) + x
       FL = (WFP - max(IWFP-1D0,WFP0))*(THES(IWFP)-THE(IWFP))/dt  ! Flux into the WFP cell, used for mass balance

       ! print debug
       !M3 = sum(THE*DZ)
       !ME = M3 - (THE0 + DF - M2 + M1)
       !ME2= M3+h(1)-(THE0+hh0+M1+P*dt)
       !if (abs(ME)>1D-8 .or. abs(ME2)>1D-8) THEN
       !Print *,4
       !endif
       !!THE0 = sum(THE*DZ)
       !!M1 = sum(DZ*S)*dt
       !!M2 = DF
       END SUBROUTINE GGA

       !rtnobd(funcd,x0,xacc,INP1,INP2,INP3,INP4,INP5,INP6)
       !CALL funcd(rtnobd,f,df,INP1,INP2,INP3,INP4,INP5,INP6)
       SUBROUTINE curveT(tr,f,df,h0,Am,Bm,Fex,km,prcp)
       REAL*8 tr,f,df,h0,Am,Bm,Fex,km,prcp,temp
       temp=(Am+Bm+Fex+h0+tr*prcp)/(Am+Bm+Fex)
       f = km*tr-h0-tr*prcp+Am*log(temp)
       df= km - prcp + Am*prcp/temp
       END SUBROUTINE curveT

       SUBROUTINE GGA_Source(N1,NM,DZ,S,h,WFP,prcp,dt,tr,DF,SCODE)
       IMPLICIT NONE
       integer N1,NM,SCODE
       REAL*8,DIMENSION(NM)::S,h,DZ
       REAL*8 :: WFP,prcp,dt,DF,tr
       REAL*8 t1,t2,t3,ss
       ! dt is the full time step
       ! tr is remaining
       SCODE = 0
       !ss = SUM(DZ(2:IWFP)*S(2:IWFP)) * tr
       ss = SUM(DZ(N1:NM)*S(N1:NM)) * dt
       t1 = h(1) + ss
       if (t1 < 0) THEN
         t3 = ss
         ss = ss + h(1)
         DF = DF + h(1)
         h(1) = 0D0 ! pond gone
         t2 = prcp*tr + ss
         if (t2 < 0D0) then
           ! both rain and h(1) are consumed by source term, normally should not happen since it got in here
           WFP = 1D0
           ss = ss + prcp*tr ! rain used up
           DF = DF + prcp*tr
           S(N1:NM) = S(N1:NM)*(ss/t3)
           SCODE = 1  ! Source too much to handle
         else
           prcp = t2/tr
           DF = DF - ss
           S(N1:NM)=0D0
         endif
       else
         h(1) = t1
         S(N1:NM)=0D0
         DF = DF - ss
       endif
       END SUBROUTINE GGA_Source

       Function swFunc(h,alpha,n,l)
       ! function to compute average wetting front suction
       ! see Jia 1998, Neuman 1976 for details
       ! for computation simplicity, obtain the value from closest n,l,h values
       ! in computing the matrix, alpha was assumed 1 (or,say, the integration variable was alpha*h)
       USE Vdata
       IMPLICIT NONE
       REAL*8 swFunc,hh
       REAL*8 h,alpha,n,l,sd
       integer I,J,K,ii

       i = round((n - g%VDZ%GA%n(1))/g%VDZ%GA%dn)+1;
       i = min(max(i,1),g%VDZ%GA%m(1));
       j = round((l - g%VDZ%GA%l(1))/g%VDZ%GA%dl)+1;
       j = min(max(j,1),g%VDZ%GA%m(2));

       hh = h*alpha
       if (hh<g%VDZ%GA%hspan(1,1)) then
        k = 1;
       elseif (hh > g%VDZ%GA%hspan(g%VDZ%GA%nspan,3)) THEN
        k = g%VDZ%GA%m(3);
       else
        do ii = 1,g%VDZ%GA%nspan
            if (ii .eq. 1) then
                sd = 1;
            else
                sd = g%VDZ%GA%hspan(ii-1,4);
            endif
            if (hh>= g%VDZ%GA%hspan(ii,1) .and. hh<g%VDZ%GA%hspan(ii,3)) then
                k = round((hh - g%VDZ%GA%hspan(ii,1))/g%VDZ%GA%hspan(ii,2))+sd;
                EXIT
            endif
        enddo
       endif
       swFunc = (g%VDZ%GA%SW(i,j,g%VDZ%GA%m(3))-g%VDZ%GA%SW(i,j,k))/alpha

       ! Treat K(h) as a two parameter function K(a,h)
       ! SW(x)=F(a,0)-F(a,h) in which F(a,x)=int^x_{-inf}K(h)dh = 1/a*int^{x*a}_{-inf}K(u)du, u=a*h
       ! F(a,x)=1/a*F(1,a*h). F(1,x) is what was tabulated

       END FUNCTION swFunc

       FUNCTION round(d)
       integer round
       real*8 d
       round = floor(d+0.5D0)

       END FUNCTION round


       SUBROUTINE curveF(Ft,f,df,tr,Am,Bm,Fex,km,dump)
       REAL*8 Ft,f,df,tr,Am,Bm,Fex,km,dump
       ! rtsafe calls this function, 2nd and 3rd arguments are output
       ! f is function value, df is derivative
       ! Ft is variable
       ! This equation being minimized is
       f = km*tr - Ft + Am*log((Am+Bm+Fex+Ft)/(Am+Bm+Fex))
       df= -1 + Am/(Am+Bm+Ft+Fex)
       END SUBROUTINE curveF


      SUBROUTINE dailySoilT(nz,na,Temp,KT,DZ,DDZ,CT,BBC,dt,a,bbase,c)
      ! solve a simple heat equation to update soil temperature
      USE Vdata
      IMPLICIT NONE
      integer nz ! last cell nz is a boundary cell (i.e. not solved for)
      integer NM,N1,NA
      REAL*8,DIMENSION(nz):: Temp,KT,CT,DZ,DDZ,a,b,c,bbase
      REAL*8 dt
      TYPE(BC_type)::BBC(2)

      REAL*8,DIMENSION(nz)::Khalf,Ttemp,HTemp,d
      REAL*8 HH

      N1 = 2
      NM = NA - 1
      HTemp = 0D0
      d = 0D0

      IF (a(nz-1) .eq. 0D0) THEN  ! only compute when a(nz-1) = 0, this is the signal
        Khalf = KT * CSHIFT(KT,1)/(DZ*CSHIFT(KT,1)+ CSHIFT(DZ,1)*KT) ! harmonic mean, may try other means
        Ttemp = Khalf/DDZ  ! valid 1:N-1. If DDZ(1)=0, valid 2:N-1 (in this case, solve from 2:N-1)
        a = - CSHIFT(Ttemp,-1)/DZ   ! relation with i-1, valid 2:N
        c = - Ttemp/DZ ! relation with i+1, valid 1:N-1
        bbase = -(a+c)
      ENDIF
      b = bbase + CT/dt
      d = (CT*Temp)/dt

      ! apply BC
      ! top BC is always the top daily temperature
      HH = BBC(1)%BData
      d(N1) = d(N1) - a(N1)*HH
      ! bottom BC can be dirichlet, or neumann (same derivative as above or no flux)
      SELECT CASE(BBC(2)%JType)
         CASE(0) ! no heat flow (semi)
           NM = NA-1
           TEMP(NA) = TEMP(NA-1)
           d(NM) = d(NM) - c(NM)*TEMP(NA)
         CASE(1) ! Dirichlet BC
           NM = NA-1
           d(NM) = d(NM) - c(NM)*BBC(2)%BData
         CASE(2) ! no heat flow (true)
           NM = NA-1
           Khalf(NM) = 0D0
           c(NM) = 0D0
           b(NM) = -(a(NM)+c(NM)) + CT(NM)/dt
           d(NM) = (CT(NM)*Temp(NM))/dt
         CASE(4) ! Free Drainage (same derivative as inner cells)
           NM = NM -1
           Temp(na)=max(0.0D0,2*TEMP(na-1)-TEMP(na-2)) ! may need to adjust for cell size
           d(NM) = d(NM) - c(NM)*BBC(2)%BData
      END SELECT

      ! solve for equation
      CALL THOMAS(NZ,N1,NM,a,b,c,d,TEMP)

      END SUBROUTINE


      SUBROUTINE SURFACE(m,n,nz,nt,ns,Veg,V,H,dt,OVH,OVF)
      ! This function takes care of ET, CANOPY, DEPRESSION STORAGE
      ! Rule of Thumb: Should not use pointer for output unless it is very clear given the context
      ! C H A O P E N G    S H E N
      ! Env. Engr.
      ! Code Unit: m

      USE Vdata
      IMPLICIT NONE
      integer m,n,nt,ns,nz,H ! H is hour of day
      integer,DIMENSION(:,:),POINTER:: nRPT
      integer,dimension(:,:),pointer:: SMAP
      integer,DIMENSION(:,:,:),POINTER::RPT
      REAL*8 dt
      REAL*8,DIMENSION(:,:),POINTER::ET,Cpc,hh
      REAL*8,DIMENSION(:,:,:),POINTER::FRAC,LAI,HI,ROOT,KC
      REAL*8,DIMENSION(:,:),POINTER::MAPP,Imp
      REAL*8,DIMENSION(:,:),POINTER::prcp,rad,wnd,hmd,ret,rets,tempe
      REAL*8,DIMENSION(:,:,:),POINTER::EB
      REAL*8,DIMENSION(:),POINTER::EE,THE,THES,THER,VS,FC,WP,hhh
      REAL*8,DIMENSION(m,n)::OVH,OVF
      TYPE(Vege_type)::Veg
      !TYPE(STA)::DAY(2)
      TYPE(VDZ_type)::V
      TYPE(THERMAL_type),POINTER:: th
      integer I,J,K,L,P,RDN,EDN,ETN,II,RMX
      REAL*8 RPET,PREC,NewS,RPrcp,FPET,EvapF
      REAL*8 PT(nt),RDF(nz),RFrac(nt)
      REAL*8,DIMENSION(nz):: PET,TP,EVAP,EDF,F
      REAL*8 THRSD, SS, RATIO, TEMP, temp2,temp3,temp4,temp5,temp6
      REAL*8,DIMENSION(:),POINTER:: TEMPP, tempp2,tempp3,tempp4
      integer NR,NCODE(3)
      REAL*8::EU=0D0,EL=0D0,RSnow=0D0,RS=0D0,RadS,Sa,TT,albedo,SNOMP
      REAL*8 ETS,ST,ETA,test,CPET,VA,GAM,OM,test2,OCOV,Newsnow,tt2

      REAL*8,DIMENSION(7) :: inpt
      REAL*8,DIMENSION(:),POINTER:: sitev,statev,outv
      REAL*8 cump, cume ,cummr, input(7), fa,errmbal,w1
      integer iflag(5),nstepday ! consider change later


      NR = 2; NCODE = 0; nstepday = 24;
      fa = 1000D0/24D0
      NCODE(1) =1; ! impervious treatment
      NCODE(2) =1; ! snowmelt

      th => g%VDZ%th
      RPT => Veg%RPT
      KC => Veg%Kc
      !CS => Veg%CS
      !CSEvap => Veg%CSEvap
      Cpc => Veg%Cpc
      !DELTA => Veg%Delta
      LAI => Veg%LAI
      Imp => Veg%Imp
      PREC = 1e-9
      !RDFN => Veg%RDFN
      FRAC => Veg%Frac
      !BSoil => Veg%BSoil
      hh => V%h(:,:,1)
      !ETC => Veg%ETC
      SMap => Veg%SMAP

      prcp => g%Wea%Day%prcp
      ret => g%Wea%Day%ret
      tempe => g%Wea%Day%temp
      hmd => g%Wea%Day%hmd
      rad => g%Wea%Day%rad
      GAM = V%VParm(1)


      iflag = (/0, 0, 3, 1, 4/) ! iflag(1)=0 means inpt(5)=qsi, inpt(6)=qli
      cump = 0D0; cume = 0D0; cummr = 0D0

      EB => V%EB
      MAPP=> V%MAPP
      ETS = 0D0
      !OPEN(UNIT=3,file='rec.dat')
      PET = 0.D0
      Sa = Veg%SNOCoef(2) ! normally set to 2.6, Qian and Gustafson
      DO I = 1,m
      DO J = 1,n
       IF (MAPP(I,J)>PREC) THEN
        ! ##11111111111111111111111111111111111111111111111111111111111111111111111111
         L = SMap(I,J)
         PET = 0.D0; Evap = 0.D0; RS = 0.D0; RadS=0.D0; Rprcp=0.D0;
         albedo = 0.23D0; TP = 0D0; NEWSNOW = 0D0; OCOV = 1D0
         !temp = Veg%SWE(I,J) ! print
         !ETA = Veg%hI(I,J)+Veg%CS(I,J)+V%h(I,J,1)+OVH(I,J)
         ! (0) Deal with impervious layer and water
      IF ((I .eq. 19) .AND.(J .eq. 42)) THEN !.AND. Veg%t>=731103.499999D0) THEN
         !print *,3
      endif
      hhh => V%h(I,J,:)
         ! (0) Decide if snow or not
         RPET=ret(H,L)*dt
         Veg%PET(I,J)=RPET

         IF (NCODE(2) .eq. 0) THEN
         ! PRECIPITATION and SNOWMELT, Brubaker's mixed indices method 55555555555555555555555555555555555555555555555555555555555555

         IF (prcp(H,L)>0) THEN
           IF (g%Wea%Day%SNOW(L)) THEN ! prcp is snow
             NEWSNOW = g%Wea%Day%prcp(H,L)*dt
             Veg%SNOW(I,J)=NEWSNOW
             Veg%SWE(I,J)=Veg%SWE(I,J)+Veg%SNOW(I,J)
             Rprcp = 0.D0
           ELSE
             Rprcp = prcp(H,L)*dt
           ENDIF
         ENDIF

         ! Calculate Snowmelt if applicable
         IF (Veg%SWE(I,J)>0) THEN
           Veg%Salb(I,J) = max(Veg%Salb(I,J) - 0.05D0*dt,0.35D0) ! albedo change due to snow aging, assumed 0.05 per day--may change later
           IF (g%Wea%Day%SNOW(L)) Veg%Salb(I,J) = 0.8D0 ! New snow albedo=0.8
           ! Update Albedo
           IF (Veg%SWE(I,J)>5d-4) THEN
             albedo = Veg%Salb(I,J)
           ELSE
             albedo=Veg%Salb(I,J)*Veg%SNOCOV(I,J)+                              &
     &             (1D0-Veg%SNOCOV(I,J))*Veg%albedo(I,J)
           ENDIF

           RadS = max(rad(H,L)*(1-Veg%Salb(I,J))-g%Wea%Day%Hb(H,L),0D0) ! Radiation on snow


           ! now modify SNOMP, the equation should give daily, now
           ! (1), times the SNOMP by dt
           ! (2), if temp<0, no snowmelt
           IF (g%Wea%Day%Temp(H,L)<0D0) THEN
             SNOMP = 0D0
           ELSE
             SNOMP = (RadS/334.D0+2D-4*max(g%Wea%Day%Temp(H,L),0.D0))*dt
             SNOMP = min(SNOMP,Veg%SWE(I,J)) !Brubaker (1996, Hydrological processes), See Dingman
                             !However this is on an hourly level. need to verify its effect
            ! Is the unit of RadS correct? From the input, it is MJ/m2/day, 334 is latent heat of
            ! fusion, kJ/kg, so MJ/ton--> MJ/m3, so correct

           ENDIF
           ! 334 is latent heat of melting*rho water, Unit m
           Veg%SNOM(I,J) = SNOMP*Veg%SNOCOV(I,J)
           Veg%SWE(I,J) = Veg%SWE(I,J)-Veg%SNOM(I,J)
           RS = Veg%SWE(I,J)/Veg%SNOCoef(1) !Veg%SNOCoef(1) is SNOCOV100
           ! new snow cover area:
           OCOV = (1D0-Veg%SNOCOV(I,J))
           RPET = RPET*OCOV ! assume this portion of PET is used up by snowmelting
           Rprcp = Rprcp + Veg%SNOM(I,J) ! Melted Snow becomes available precipitation (on the ground)
           Veg%SNOCOV(I,J)=1D0-(exp(-Sa*RS)-RS*exp(-Sa))
           ! Reference:Noah LSM,see Qian and Gustafson 2009, Journal of geophysical research Vol 114
         ENDIF
         !!!!!!!!!!5555555555555555555555555555555555555555555555555555555555555555555555555555555555555

         ELSE ! UEB code
         ! decide if it is necessary to invoke the snow module
         ! temporarily remove snowflag (to compare with UEB original
         IF ((tempe(H,L)<th%SParm(1) .and. prcp(H,L)>0D0)                                 &
     &       .or. Veg%SWE(I,J)>0D0) THEN !Only enter snow model when: there is new snow, or there are remaining snow on ground
         ! rad unit is MJ/m2/day
         ! ime unit need to be in hour
         input=(/tempe(H,L),prcp(H,L)/24D0,g%Wea%Day%wnd(H,L),hmd(H,L),             &
     &      rad(H,L)*fa,g%Wea%Day%hli(H,L)*fa,g%Wea%Day%cosz(H,L)/) ! incoming longwave is zero?, notice rad unit!
         outv => th%Soutv(I,J,:)
         ! to provide to snow model,radiation need to be in kJ/m^2/h
         ! from MJ/m2/day, that's 1000/24
         statev => th%Sstatev(I,J,:)
         sitev => th%Ssitev(I,J,:)
         tempp => th%tss(I,J,:)
         tempp2 => th%tsa(I,J,:)

         ! snow cover will need to be fixed
         !w1 = th%Sstatev(I,J,2)
         CALL SnowLSub(DT*24,1,input,th%Ssitev(I,J,:),th%Sstatev(I,J,:)         &
     &  ,th%SParm,iflag,nstepday, cump,cume,cummr                                &
     &  ,outv,th%tss(I,J,:),th%tsa(I,J,:),th%nd,th%dfc)

         !write(303,"(24g10.2)") V%t,outv(1:23)
         !errmbal=w1+cump-cummr-cume-th%Sstatev(I,J,2)
         !if (abs(cummr-outv(22)+outv(23))>1d-8 .or.                             &
         !&        abs(errmbal)>1D-5) then
         !endif


         RPET = RPET*(1D0-th%Sstatev(I,J,7)) ! assume this portion of PET is used up by snowmelting
         Rprcp = outv(22)
         ! testing phase, remove all extra water:
         !Rprcp = 0D0
         Veg%SWE(I,J) = th%Sstatev(I,J,2)
         Veg%Salb(I,J) = outv(3)
         Veg%SNOCOV(I,J) = outv(16)
         Veg%SNOM(I,J) = outv(23) ! temporarily use this to store sublimation!!!!!!!
         th%state = 1D0
         ELSE
         ! only rain
         Rprcp = prcp(H,L)*dt

         ENDIF

         ENDIF

         ! test snow phase
         !IF (NCODE(3) > 1) THEN

         ! now treat impervious surfaces
         ! it has its own storage: hI
         ! rule of applying ET: each class working on its own portion,
         ! do not use others, except for shared Canopy Storage (for now)
         VA = 1D0 - Veg%Imp(I,J);tt2=0D0
         IF ((NCODE(1) .eq. 1) .and. Veg%Imp(I,J)>0D0) THEN
         ! impervious treatment(1) : disregard
         tt2 = RPrcp*Veg%Imp(I,J)
         Veg%hI(I,J) = Veg%hI(I,J) + tt2
         !RPET = RPET*(1D0-Veg%Imp(I,J)) ! take this portion out of total PET, will only be applied on impervious fraction
         CPET = RPET*Veg%Imp(I,J)
         Evap(1)= max(min(CPET,Veg%hI(I,J)),0D0)
         Veg%hI(I,J) = Veg%hI(I,J) - Evap(1)
         RPET = RPET - CPET

         Rprcp = Rprcp - tt2
         ! now both RPET and Rprcp are for vegetated portion
         ELSE
         CPET = 0D0
         ENDIF

         ! now also evaporate from flow domain
         if (OVH(I,J)>0D0) then
           ! RPET in the flow domain preferrably evaporate the flow water
           ! if extra, will come back join rpet
           FPET = RPET*OVF(I,J)
           EvapF= max(min(OVH(I,J),FPET),0D0)
           OVH(I,J)=OVH(I,J)-EvapF
           RPET = RPET - FPET
           Evap(1)=Evap(1)+EvapF
         else
           FPET = 0D0
         endif


          ! Using station prcp!!!!!
         IF (g%Wea%Day%SNOW(L)) THEN
           Veg%NewS(I,J) = 0D0
         ELSE
           Veg%NewS(I,J) = min(Cpc(I,J)-Veg%CS(I,J),Rprcp)
         ENDIF
         Rprcp = Rprcp - Veg%NewS(I,J) ! subtract canopy storage from remaining prcp
         Veg%CS(I,J)= Veg%CS(I,J)+Veg%NewS(I,J)
         !temp = Veg%CS(I,J)
         Veg%EvapCS(I,J)= min(Veg%CS(I,J),RPET)
         !temp = Veg%EvapCS(I,J)
         Veg%CS(I,J)= Veg%CS(I,J)-Veg%EvapCS(I,J)
         !temp = Veg%CS(I,J)
         RPET=RPET-Veg%EvapCS(I,J)
         !temp = RPET(I,J)


         ! Use Soil ET demand to evaporate surface water
         !SS = max(0D0,(Rprcp(I,J)-RPET(I,J)*Veg%BSoil(I,J))/dt) ! Source term for overland flow
         !Evap(1)=min(RPET(I,J),(hh(I,J)+Rprcp*dt)*Veg%BSoil(I,J))
         ! Previously, this section is behind PT calculation, I think that was a mistake
         ETS=min(RPET,(max(hh(I,J)+Rprcp,0.D0)))
         RPET = RPET-ETS
         Evap(1) = Evap(1)+ETS
         !Rprcp = Rprcp - ETS

         ! print debug
         !test2 = ret(H,L)*dt*OCOV - RPET -                                     &
         !&    (Veg%EvapCS(I,J) + Evap(1))


          ! (1) Compute PET for different crops using crop coefficient
         !IF (Rprcp > PREC) THEN
         !### 222222222222222222222222222222222222222222222222222222222222222222222222222
         DO K=1,nt
           !P = RPT(I,J,K)
           PT(K) = RPET*KC(I,J,K)/VA
           !Veg%DELTA(I,J,K)=(1-exp(-0.5D0*LAI(I,J,K))) ! Beer-Lambert Law
         ENDDO
         ! This PT is a whole cover, theoretical value, 'what would be if entirely consisted of this class'
         ! As such, it needs to be rescaled back to 1, (divided by VA)

         RFrac = FRAC(I,J,1:nt)*Veg%Delta(I,J,1:nt) ! This is repetitive computation (already done in updateveg, but not stored)

         !!!! 444444444444444444444444444444444444444444444444444444444444444444444444444444
         ! Compute EDF and RDF
         ! RDN and RDFN
         RDF = 0.D0
         EE => V%EB(I,J,:)
         RMX = 1;

         IF (sum(PT)>PREC) THEN
         DO K=1,nt
            RDN = Veg%RDFN(I,J,K)
            IF (RDN>1) THEN
              RDF(1)=0.0D0
              RDF(2:RDN-1)=V%DZ(I,J,2:RDN-1)
              RDF(RDN)=Veg%ROOT(I,J,K)-(EE(1)-EE(RDN-1)) ! Assume linear RDF
              RDF(2:RDN)=RDF(2:RDN)/sum(RDF(2:RDN)) ! rescale to 1
              RDF(RDN+1:nz)=0.D0
              PET(2:RDN) = PET(2:RDN)+ RDF(2:RDN)*PT(K)*RFrac(K)
              RMX = max(RMX,RDN)
            ENDIF
         ENDDO
         ELSE
         RMX=1;RDN=1;
         ENDIF

         !EDN = min(3,V%NZC(I,J))


         EDF=0.D0
         IF (RPET > PREC) THEN
          DO K=2,V%NZC(I,J)-1
           EL=EE(1)-EE(K); EU=EE(1)-EE(K-1)
           EDF(K)=EL/(EL+exp(2.374D0-7.13D0*EL)/1000)-                                      &
     &            EU/(EU+exp(2.374D0-7.13D0*EU)/1000) ! SWAT P133 with unit converted to m
      IF ((EDF(K)<0.02D0) .or.(EL > 0.2) .or. (K .eq. V%NZC(I,J)-1))                       &
     &    THEN
             EDN = K
             EXIT
            ENDIF
          ENDDO
         ELSE
          EDN = 1
         ENDIF
         EDF(EDN+1:nz)=0.D0
         ! end computing EDF and RDF $4444444444444444444444444444444444444444444444444444444444

         ! Water Stress
         THE => V%THE(I,J,:)
         THER => V%THER(I,J,:)
         THES => V%THES(I,J,:)
         FC => V%FC(I,J,:)
         WP => V%WP(I,J,:)
         F = 0.0D0; !F2=0.0D0
         !F(2:ETN)=1.D0-                                                                   &
         !&    (WP(2:ETN)*THE(2:ETN)-THER(2:ETN))/(THES(2:ETN)-THER(2:ETN)) ! Plant water stress
         ! now treat evaporation and transpiration differently, each on its own territory
         !temp =  (Veg%BSoil(I,J)/VA)*RPET
         temp =  Veg%BSoil(I,J)*RPET  ! still not sure whether we should do the re-scaling
         !Evap(2:EDN)=Veg%BSoil(I,J)*(RPET/VA)*EDF(2:EDN) ! should we rescale VA here?-- i think so
         ! Impervious is excluded here as BSoil got rid of impervious layer
         DO K = 2,EDN
           IF (THE(K)>FC(K)) THEN
             F(K) = 1D0
           ELSEIF (THE(K)<=WP(K)+0.002D0) THEN
             F(K) = 0D0
           ELSE
             F(K) = exp(2.5D0*(THE(K)-FC(K))/(FC(K)-WP(K)))
           ENDIF
           Evap(K) =temp*F(K)*EDF(K)
         ENDDO

         F = 0D0
         !GAM = 0.03D0

         ETN = min(RMX,V%NZC(I,J)-1)
         DO K=2,ETN
           ! use LK00 paper
           THRSD = .05D0*(FC(K)-WP(K))+WP(K) ! safety threshold not to drain water
           IF (THE(K)<THRSD) THEN
            F(K) = 0D0
           ELSE
            ST = (THE(K)-WP(K))/THES(K)
            F(K) = ST**(GAM/(THE(K)-WP(K)))
           ENDIF
           TP(K) = PET(K)*F(K) ! PET is just potential transpiration (evaporation has been done above)
         ENDDO

         ! Output Evap,Tp,ET, and convert to VDZ source term
         ! Accounting
         Veg%TP(I,J)=sum(TP(2:ETN))
         Veg%Evap(I,J)=sum(Evap(1:EDN)) ! All evaporation from soil and ground surface (not from canopy)
         Veg%EvapG(I,J)=Evap(1)
         Veg%ET(I,J)=Veg%Evap(I,J)+Veg%TP(I,J)+Veg%EvapCS(I,J) ! Actual ET
         Veg%Gprcp(I,J)=Rprcp+tt2 ! prcp on the ground, should add prcp to impervious as well
         !temp2 = Veg%ET(I,J)
         !temp3 = Veg%ET(I,J)-prcp(H,L)*dt
         hhh => V%h(I,J,:)
         RPET= RPET - Veg%TP(I,J) - (Veg%Evap(I,J)-Evap(1))
         Veg%RPET(I,J)=RPET
         !temp4 = RPET(I,J)
         !temp5 = Veg%PET(I,J)
         VS => V%SRC(I,J,:) ! Assume only source term in the VDZ are ET extractions
         VS = 0.D0
         VS(1) = (-ETS+ Rprcp)/dt
         ETN = max(ETN,EDN)
         VS(2:ETN) = - ((TP(2:ETN)+Evap(2:ETN))/V%DZ(I,J,2:ETN))/dt
         !IF (VS(1)>1) THEN
         !ENDIF

         !IF (VS(1)>0 .AND. VS(2)<0) THEN
         !WRITE(3,"(8g16.6)") Veg%EvapCS(I,J),Evap(1),temp,Veg%TP(I,J),                     &
         !&         temp2,Veg%PET(I,J),prcp(H,L)*dt
         !ENDIF

         !## Test mass and ET conservation section
         !tt = Veg%EvapCS(I,J)
         ! mass balance checker
         ! print debug
         !temp2 = Veg%SNOCOV(I,J)
         !temp3 = ret(H,L)*dt*(1D0-Veg%SNOCOV(I,J))
         !temp4 = sum(VS(2:nz)*V%DZ(I,j,2:nz))*dt
         !temp = Veg%SWE(I,J)
         !tt2 = Veg%SNOM(I,J)
         !test = Veg%hI(I,J)+Veg%CS(I,J)+OVH(I,J)+VS(1)*dt + hhh(1)-Eta                              &
         !& -(prcp(H,L)*dt-NEWSNOW + Veg%SNOM(I,J) -Evap(1)-Veg%EvapCS(I,J))
         !temp2 = ret(H,L)*dt*OCOV - RPET -   CPET -FPET-                         &
         !&    (Veg%EvapCS(I,J) + ETS - sum(VS(2:nz)*V%DZ(I,j,2:nz))*dt)
         !if (abs(test)>1e-8 .or. abs(temp2)>1e-6) THEN
         !print *,4
         !endif
         !IF (ETN<1 .or. ETN>100) THEN
         !print *,4
         !endif
         !ETS = ETS + Veg%TP(I,J) + temp

      ! ###22222222222222222222222222222222222222222222222222222222222222222222222222
      ! ###111111111111111111111111111111111111111111111111111111111111111111111111111
      ! testing snow add one if
        !ENDIF
        ENDIF
      ENDDO
      ENDDO
      !CLOSE(3)
      ! print debug
      !ETS = ETS/(67D0*50D0)
      !ST = (sum(V%SRC(:,:,2:12)*V%DZ(:,:,2:12)))/(67D0*50D0)
      !test = ST*dt + ETS
      !IF (abs(test)>1e-10) THEN
      !print *, 4
      !endif

      END SUBROUTINE SURFACE
      
!=====================================================================================
      subroutine riv_cv(nx,h,hx,A,u,W,WF,R,E,Qx,dx,Mann,S,SM,S0,SE,ddt                &
     &   ,nGho,Code,CType,rtr,UPST,UPSV,DST,DSV,STATE,umax,hg,mr,KB,fG,               &
     &   tt)
      ! RK3 river flow
      ! Chaopeng Shen
      ! PhD Environmental Engineering
      ! Michigan State Univ.
      IMPLICIT NONE
      integer nx, nGho,EQQ !EQQ: 1--KW, 2--DW, 3--DyW
      REAL*8,DIMENSION(nx):: h,A,u,E,Mann,S,S0,W,Qx,SM,SE,R,wF,hg,mr,KB
      REAL*8 dx,ddt, rtr, code, fG(nx)
      integer CType(nx) ! Boundary Condition Types
      REAL*8 UPSV, DSV(3),UPST,DST
      INTEGER STATE, ISTATE,ICODE
      
      integer T(nx),I,K,KK,NN, MAXNN !SHIFT
      REAL*8 G, tt, tt1, tt2, sff, Neg, RDT, dt,EL,ER, PREC
      REAL*8,DIMENSION(nx):: dAdt0,dAdt1,dAdt2,fG0,fG1,fG2,hN
      REAL*8,DIMENSION(nx):: QQ,Q0,Q1,Q2,U0,U1,U2,hx,ht,rt,hgg,fgg
      REAL*8 UMAX,COU,sgn(nx),temp(nx),temp2(nx)
      REAL*8 TQ(2),SUMQ,SUMQ2,T1,T2,TA(nx),QQgc(nx)
      REAL*8 AA_old,AAA,QQgc_old,ww,test,test2,Qout
      
      G = 9.81D0      
      PREC = 1.0D-11
      sgn = 1.0D0
      Neg = 1e-6
      QQ = 0.0D0
      Qx = 0.0D0
      RDT = DDT
      DT = DDT
      EQQ = int(code)
      K = 1
      dAdt0=0.0D0;dAdt1=0.0D0;dAdt2=0.0D0;
      U0=0.0d0;U1=0.0d0;U2=0.0d0;
      Q0=0.0D0;Q1=0.0D0;Q2=0.0D0;
      hgg = 0.D0; fgg=0.D0
      COU = 0.4D0
      ! Adaptively Choosing dt according to Courant Condition:
      umax=max(maxval(abs(u))*3.D0/2.D0,1.0D-12) ! Kinematic Wave Speed
      dt = COU*dx/umax*0.8D0 ! 0.8 is for the safety of growing velocity
      dt = min(ddt,dt,RDT)
      
      ! print for debug, comment them
      !AA_old = sum(A(2:nx-1))*dx
      !ww = sum(S)*dx*DDT
      !IF (tt>=733408.4166D0 .AND. abs(AA_old-0.250621D0)<1d-3)THEN
      !print *,3
      !endif
      MAXNN = 5
      

      do while (RDT>PREC .AND. K<=MAXNN)
        
      ! Second Order TVD RK method:
        !STATE = 1
        !CALL one_step(A,dAdt0,Q0,u,hx,h,R)
        
        !STATE = 0
        !TEMP = DT*(KB/86400D0)/mr
        
        CALL one_step(A,dAdt0,Q0,u,hx,h,R,fG0,dx,EQQ) ! TESTING WHAT IS WRONG
        STATE = 0
        
        CALL one_step(A+dt*dAdt0,dAdt1,Q1,u1,hx,h,R,fG1,dx,EQQ)
        
        CALL one_step(A+(dt/4)*(dAdt0+dAdt1),dAdt2,Q2,u2,hx,h,R,                     &
     &               fG2,dx,EQQ)
        
        umax = max(maxval(abs(u)),maxval(abs(u1)))*(3D0/2D0)
        if (umax<PREC) umax = PREC
        if (COU*dx/(umax) < dt .AND. K<=MAXNN-1) THEN
           dt = dt/2.D0
           K  = K + 1
        else
        
        A  = A + (dt/6.D0)*(dAdt0+4.D0*dAdt2+dAdt1)
           
           QQ  = QQ + (dt/6.D0)*(Q0 + Q1 + 4.D0*Q2)
           !fG = fG + (dt*dx/6.D0)*(fG0+4.D0*fG2+fG1) ! Now Unit is m3
           RDT= RDT - dt
           dt = min(RDT,dt)
           ! Correct what may become negative depth: Only Allow that much to flow out
           DO I =2,nx
             IF (A(I)<0) THEN
             TQ=0.0D0; SUMQ=0.0D0;SUMQ2=0.0D0;T1=0.0D0;T2=0.0D0
               IF (S(I)>=0) THEN
               ! Recover Negative Depth from outflow fluxes
                 TQ(1)=QQ(I-1); TQ(2)=QQ(I)
                 SUMQ = abs(A(I))*dx
                 IF (TQ(1)<0) T1 = TQ(1)
                 IF (TQ(2)>0) T2 = TQ(2)
                 SUMQ2= T1 + T2
                 T1 = T1*SUMQ/SUMQ2; T2 = T2*SUMQ/SUMQ2
                 A(I) = 0.0D0
                 A(I-1)=A(I-1)-T1/dx; A(I+1)=A(I+1)-T2/dx
                 QQ(I-1)=TQ(1)+T1
                 QQ(I)=TQ(2)-T2
               ELSE
                 T1 = A(I)+(QQ(I)-QQ(I-1))*(dt/dx) ! Recover A before Q is applied (After S is used)
                 IF (T1>0) THEN
                 ! Use all S, the same as S(I)>0
                   TQ(1)=QQ(I-1); TQ(2)=QQ(I)
                   SUMQ = abs(A(I))*dx
                   IF (TQ(1)<0) THEN
                      T1 = -TQ(1)
                   ELSE
                      T1 = 0D0
                   ENDIF
                   IF (TQ(2)>0) THEN
                      T2 = TQ(2)
                   ELSE
                      T2 = 0D0
                   ENDIF
                   !T1 = -DBLE(TQ(1)<0)*TQ(1); T2=DBLE(TQ(2)>0)*TQ(2) ! gfortran complains about this statement
                   SUMQ2= T1 + T2
                   T1 = T1*SUMQ/SUMQ2; T2 = T2*SUMQ/SUMQ2
                   A(I) = 0.0D0
                   A(I-1)=A(I-1)-T1/dx; A(I+1)=A(I+1)-T2/dx
                   QQ(I-1)=TQ(1)+T1
                   QQ(I)=TQ(2)-T2
                 ELSE ! source term is sufficient to reduce to negative depth, in this case, no outflow at all
                 ! we warrant a little negative depth due to source term mass-balance
                 ! notice the signs are different for Q(I) and Q(I-1)
                   TQ(1)=QQ(I-1); TQ(2)=QQ(I)
                   IF (TQ(1)<0) THEN
                     A(I-1)=A(I-1)+TQ(1)/dx;
                     A(I)=A(I)-TQ(1)/dx
                     QQ(I-1)=0.0D0
                   ENDIF
                   IF (TQ(2)>0) THEN
                     A(I+1)=A(I+1)-TQ(2)/dx
                     A(I)=A(I)+TQ(2)/dx
                     QQ(I)=0.0D0
                   ENDIF
                 ENDIF
               ENDIF
             ENDIF
           ENDDO
           ! End Correction Loop
        endif
      
      enddo
      Qx = QQ/DDT
      h = A/w
      WHERE (hx>0)
        u = Qx /(wF * hx)
      ELSEWHERE
        u = 0.D0
      ENDWHERE
      ! Do Groundwater interaction operator: Backward Euler
      !!!!! GEOMETRY
      ! New updated h
      
      ! print for debug, comment them
      !AAA = sum(A(2:nx-1))*dx
      !Qout = Qx(nx-1)*DDT
      !test2 = AAA - (AA_old + ww - Qout)
      ICODE = 1
      IF (ICODE >0) THEN
      WHERE ((hg - (E - mr))>0) ! Reference: Gunduz Thesis
        hgg = hg
      ELSEWHERE
        hgg = (E - mr)
      ENDWHERE
      TEMP = (DDT*(KB/86400D0))/mr
      !TEMP2 = hgg - E
      !TEMP2 = (h+E-hg)
      hN = (TEMP*hg+h+E)/(1D0+TEMP) - E
      !TEMP2 = (hN - h)/(ddt/86400)
      
      WHERE (hN<0.0D0) hN=0.0D0
      QQgc = (hN-h)*w*dx
      fG = fG + QQgc ! Unit :m3
      A = hN * w      
      ENDIF
      ! print for debug, comment them
      !AAA = sum(A(2:nx-1))*dx
      
      !test = AAA - (AA_old+sum(QQgc)+ ww - Qout)
      !if (abs(test)>1e-8 .OR. abs(test2)>1e-8) THEN
      !print *,tt
      !endif
      
      ! end of do while loop
      !STATE = 0
      !CALL one_step(A,dAdt2,Q2,u,hx,h,R,fG1,dx,EQQ) ! output Qx,u,hx,h
      !STATE = 1  ! if not changed before coming back, then STATE=1
      
      
      CONTAINS
! BEGIN CONTAINS ##########################################################
      SUBROUTINE one_step(AA,dAdt,Q,UU,hhx,hh,RR,fGG,dx,EQ)
      ! Input: AA,(UU,wF,hhx). Output: dAdt,Q,UU,hhx,hh,RR
      ! Does not modify AA
      IMPLICIT NONE
      REAL*8,DIMENSION(nx):: AA,dAdt,Q,hN,hgg,TEMPC
      INTEGER OPT,EQ
      REAL*8 hhx(nx),Ex(nx),Ey(nx),SE(nx),Sgn(nx),UU(nx),hh(nx)
      REAL*8 RR(nx),AF(nx),Eta(nx),hy(nx),fGG(nx),AAA,dx
      
      dAdt = 0.0D0
      Ex = 0.0D0; Ey = 0.0D0
      Sgn= 1.0D0
      Eta = 0.0D0
      Q  = 0.0D0
      fGG = 0.D0
      

      !IF (STATE .ne. 1) Then

        ! Resolve Geometry !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Currently only rectangular geometry      
        hhx = 0.0D0; hy = 0.0D0
        hh = 0.0D0; UU = 0.0D0
        hh = AA/w
        Eta = E+hh
        
        CALL BC(nx,hh,AA,u,w,UPST,UPSV,DST,DSV,dx)
        
        CALL interfaceH(1,nx,hh,Eta,Ex,Ey,hhx,hy)
        
        IF (int(UPST) .eq. 0) hhx(1)=0.D0
        
        RR = hhx
        ! Apply BC
        IF (EQ .eq. 2) THEN  ! Diffusive Wave!
          ! Compute U and Flux and dhdt
          SE = (Eta - CSHIFT(Eta,1))/dx
          !SE = 1.0D-3
          !SE = (S0 - (CSHIFT(h,1)-h)/dx)
          
          WHERE (SE<0) Sgn=-1.0D0  
          UU = Sgn*(RR**(2.D0/3.D0))*(abs(SE)**0.5D0)/Mann
        ENDIF
        
      !ENDIF

      ! Channel Type of the interface should just inherit from upstream channel
      ! here also resolve geometry !!!!!!!!!!!!!!!!!!      
      Q  = wF * hhx * UU ! This Q is Q at this moment, this is different from the output Qx, which is a mass-balance Q (Q in a period of time)
      ! Now compute interchange with groundwater, use an implicit update (operator splitting)
      ! for rectangular channel:
      
      IF (int(UPST) .eq. 2) Q(1)=UPSV  ! UPST type 2 is not working inside BC. it will be overwritten outside
      
      dAdt =  (CSHIFT(Q,-1)-Q)/dx + S
      
      dAdt(1) =0.0D0
      dAdt(nx)=0.0D0
      ! As long as channel geometry is longitidinally homogeneous, this relation holds
      END SUBROUTINE one_step
! END CONTAINS ##############################
      end subroutine riv_cv
      

      !CALL BC(nx,hh,AA,u,w,UPST,UPSV,DST,DSV)
      SUBROUTINE BC(nx,h,A,u,w,UPST,UPSV,DST,DSV,dx)
      integer nx
      REAL*8,DIMENSION(nx):: h,A,u,w
      REAL*8 UPSV, DSV(3), UPST, DST
      INTEGER UPS, DS
      REAL*8 sgn, EL,ER,AAA,dx
      ! BC
      
      UPS = int(UPS)
      DS  = int(DST)
      SELECT CASE(UPS)
      CASE(0) ! WALL
       h(1)=0.0D0
       A(1)=0.0D0
       u(1)=0.0D0
       
      CASE(1) ! Dirichlet Head
       h(1)=UPSV
       w(1)=w(2)
       A(1)=h(1)*w(1)
       u(1)=u(2)
      CASE(2) ! Prescribed Flux
       h(1)=h(2)
       A(1)=A(2)
       u(1)=UPSV/A(1)
      CASE(3) ! Neumann--same as internal cell
       h(1)=h(2)
       A(1)=A(2)
       u(1)=u(2)
      END SELECT
      
      SELECT CASE(DS)
      CASE(0) ! WALL
       h(nx)=0.0D0
       A(nx)=0.0D0
       u(nx)=0.0D0
      CASE(1) ! Dirichlet Head
      
       h(nx)=DSV(1)
       w(nx)=w(nx-1)
       
       A(nx)=h(nx)*w(nx)
       
       u(nx)=u(nx-1)
      CASE(2) ! Prescribed Flux
      CASE(3) ! Free Outflow Condition
       h(nx)=max(0.0D0,min(h(nx-1)-0.01,2*h(nx-1)-h(nx-2))) ! outflow shouldn't have head higher than inside
       u(nx)=max(0.0D0,2*u(nx-1)-u(nx-2))
       A(nx)=h(nx)*w(nx)
      CASE(4) ! Junction BC
      END SELECT
      
      END SUBROUTINE BC
      
      Subroutine interfaceH(ny,nx,h,E,Ex,Ey,hx,hy) ! THis version used a limiter similar to the minmod to control overshooting flux from a shallow cell
      implicit none
      integer ny, nx
      real*8, DIMENSION(ny,nx):: h, E, Ex,Ey,hx,hy
      
      !!! THIS ONE DOES NOT LOOK AT Ex or Ey, DISREGARD EDGE EFFECT
      
      integer i,j
      real*8 ER(ny,nx)
      REAL*8 PREC,t1,t2,phi
      
      hx = 0.0D0
      hy = 0.0D0
      ER = 0.0D0
      PREC = 1.0D-10
      ! hx!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF (nx>1) THEN
      DO I = 1,ny
        DO J = 1,nx-1
          if (E(I,J)>=E(I,J+1)) THEN ! Positive Wind
            if (h(I,J)<PREC) THEN
              hx(I,J)=0.0D0
            elseif (h(I,J+1)>PREC) THEN
              if (h(I,J)*4D0<h(I,J+1)) THEN !!!!!!!!!!!!!!!!!!!!!!!!!!CONSIDER MULTIPLY BY 4, like in interfaceHH.F90!!!!!!!!!!!!!!!!!!
                hx(I,J)=h(I,J)
              ELSE
                hx(I,J)=0.5D0*(h(I,J)+h(I,J+1)) ! modification 05/06/2009
              endif
            else ! half-wet
              hx(I,J)=h(I,J)
            endif
            
          else
            if (h(I,J+1)<PREC) THEN
              hx(I,J) = 0.0D0
            elseif (h(I,J)>PREC) THEN
              !hx(I,J) = 0.5D0*(h(I,J)+h(I,J+1))
              
              if (h(I,J+1)*4D0<h(I,J)) THEN
                hx(I,J)=h(I,J+1)
              ELSE
                hx(I,J)=0.5D0*(h(I,J)+h(I,J+1)) ! modification 05/06/2009
              endif
            else ! half-wet
              hx(I,J) = h(I,J+1)
            endif
            
          endif
        ENDDO
      ENDDO
      hx = max(hx,0.0D0)
      ENDIF
      
      ! hy!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF (ny>1) THEN
      DO I = 1,ny-1
        DO J = 1,nx
          if (E(I,J)>=E(I+1,J)) THEN
            if (h(I,J)<PREC) THEN
              hy(I,J)=0.0D0
            elseif (h(I+1,J)>PREC) THEN
              !hy(I,J)=0.5D0*(h(I,J)+h(I+1,J))
              ! remark, the following selection is to impose stability, avoid having little depth while on the other end depth is large
              ! then the outflow would be too large. This can be improved later
              if (h(I,J)<h(I+1,J)) THEN
                hy(I,J)=h(I,J)
              ELSE
                hy(I,J)=0.5D0*(h(I,J)+h(I+1,J)) ! modification 05/06/2009
              endif
            else ! half-wet
              hy(I,J)=h(I,J)
            endif
          else
            if (h(I+1,J)<PREC) THEN
              hy(I,J) = 0.0D0
            elseif (h(I,J)>PREC) THEN
              !hy(I,J) = 0.5D0*(h(I,J)+h(I+1,J))
              
              if (h(I+1,J)<h(I,J)) THEN
                hy(I,J)=h(I+1,J)
              ELSE
                hy(I,J)=0.5D0*(h(I,J)+h(I+1,J)) ! modification 05/06/2009
              endif
            else ! half-wet
              hy(I,J) = h(I+1,J)
            endif
          endif
        ENDDO
      ENDDO
      hy = max(hy,0.0D0)
      ENDIF      
            
      END SUBROUTINE interfaceH
      
      END MODULE FLOW

