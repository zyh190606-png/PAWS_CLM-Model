#include <preproc.h> 
    ! Flow library for the PAWS watershed model
    ! Michigan State UniV% CEE
    ! Chaopeng Shen
    ! PhD, Environmental Engineering

      MODULE FLOW
      !USE VData
      IMPLICIT NONE
      PRIVATE
      PUBLIC :: GGA,vdz1c,Aquifer2D,Aquifer3D,applyCond,surface
      PUBLIC :: aquifer2dNL,runoff,ovn_cv,lowland,riv_cv,fixedHeadGWupdate
      !PUBLIC :: Aquifer2D
      !SAVE
      
      interface vGh
        module procedure vGh_scalar
        module procedure vGh_array
      end interface
      
      interface vGTHE
        module procedure vGTHE_scalar
        module procedure vGTHE_array
      end interface
      
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
      
      USE VData, Only : oneD_ptr, pcg, elem, BC_type, DEBUGG
      USE displayMod, Only: display
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
      REAL*8,DIMENSION(:,:),Pointer:: hNew=>NULL(),hE=>NULL(),hW=>NULL()
      REAL*8,DIMENSION(:,:),Pointer:: hNN=>NULL(),hS=>NULL(),DRR=>NULL()
      INTEGER N1,NM,M1,MM,NNY,NNX,NN,n,iter,itmax,np,IND(2)
      REAL*8 err,meanW,meanH0,meanH,meanE,H0(ny,nx),DP2,PK(2)
      REAL*8,DIMENSION(:),POINTER::pp1=>NULL(),pp2=>NULL(),pp3=>NULL()
      REAL*8,DIMENSION(:),POINTER::pp4=>NULL(),pp5=>NULL() ! print 
      REAL*8,DIMENSION(:),ALLOCATABLE::pk1,pk2,pk3,pk4,pk5
      REAL*8,DIMENSION(3,3):: TEMPI1,TEMPI2
      IF (DEBUGG) THEN
        meanW = sum(Wm*Mapp)/size(h); meanH0=sum((H-E)*S*Mapp)/size(h) ! print *, debugging
      ENDIF
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
      CALL applyCond(0,BBC,NBC,n,5,A1,A2,A3,A4,A5,M,RHS,K,HH)
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
      
      ! Apply Dirichlet Boundary Condition
      ! Flux BC should be converted to a source term and a no flow BC
      call applyCond(1,BBC,NBC,n,5,A1,A2,A3,A4,A5,M,Rm,K,HH)
      ! Apply Lower/Upper Boundary conditions
      call applyCond(2,BBC,NBC,n,5,A1,A2,A3,A4,A5,M,Rm,K,HH)
      
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
      !np = size(h)
      
      IF (DEBUGG) THEN
          if (size(BBC)>0) then
          np = size(BBC(1)%Cidx)    ! Note: in the new Fortran version of the
          ! model, BBC is actually g%GW(x)%Cond(x)
      
          ALLOCATE(pk1(np),pk2(np),pk3(np),pk4(np),pk5(np))
          pp1 => oneD_Ptr(nn,hNew); pk1 = pp1(BBC(1)%Cidx)
          pp2 => oneD_Ptr(ny*nx,h0); pk2 = pp2(BBC(1)%idx)
          pp3 => oneD_Ptr(ny*nx,HH); pk3 = pp3(BBC(1)%idx)
          pp3 => oneD_Ptr(ny*nx,Dperc2); pk3 = pp3(BBC(1)%idx)
      
          pk3 = BBC(1)%K*(pk1 - BBC(1)%h)/BBC(1)%D
          DP2 = sum(pk3)/size(h)
      
          DEALLOCATE(pk1,pk2,pk3,pk4,pk5)
          else
          DP2 = 0D0
          endif
      
          meanH = sum((HH-E)*S*Mapp)/size(HH)
          meanE = meanH-(meanH0 + meanW*dt - DP2*dt)
          if (abs(meanE)>5e-7) THEN
          call display("Aquifer2D MBE:", meanE)
          endif
      ENDIF
      !TEMP = (hNew - h(N1:NM,M1:MM))*S(N1:NM,M1:MM)/dt
      ! Percolation to lower layer is computed outside
      END SUBROUTINE Aquifer2D


      SUBROUTINE Aquifer3D(ny,nx,nz,nG,h,E,D,K,KZ,Wm,Ss,MAPP,BBC,NBC,                       &
     &             dd,dt,hN,DR)
      ! Dimension of the equation is now [T^-1] (not the same as aquifer2D)
      ! Units of W and S are also [T^-1] 
      ! This code follows the same format for Aquifer2D
      ! Assuming nG([1,2,3])>=1
      ! Assuming DX=DY, but D (it is DZ) can be spatially varying
      ! DXHalf = DX(i+1/2). Using T(i+1/2)/DX(i+1/2) to get K(i+1/2) to take thickness into account
      ! Assuming that KX = KY but generally not equal to KZ
      
      ! Recharge is implemented as a source term
      ! 
      
      USE VData, Only : oneD_ptr, pcg, elem, BC_type, gdG
      USE gd_Mod !  DebuggingTool
      IMPLICIT NONE
      INTEGER ny,nx,nz,NBC,MTD,nG(3),nCond ! NBC is the size arrays for BC specifications
      INTEGER ATYPE 
      REAL*8,DIMENSION(ny,nx)::MAPP
      REAL*8,DIMENSION(ny,nx,nz)::WTD,K,KZ,WTDNEW,D,hN,Ss
      REAL*8,DIMENSION(ny,nx,nz),Target::E,HH,DR,h,Wm ! E for a 2D aquifer is its bottom
      !REAL*8,DIMENSION(nz,nt) :: THES,FC ! Database items. For computing Specific Yield
      TYPE(BC_type):: BBC(NBC)
      REAL*8 DD(2),DX, DY,dt       ! Assuming DX, DY is the same (uniform grid)
      
      REAL*8,DIMENSION(ny,nx,nz),target:: T,Thalfx,Thalfy                       &
     &       ,Thalfz,Tempx,Tempy,Tempz,Dperc2,DZZHalf
      INTEGER I,J,nGx,nGy,nGz,MD(7)
      REAL*8,DIMENSION(ny-2*nG(1),nx-2*nG(2),nz-2*nG(3)),target::                       &
     &           A1b,A1,A2,A3,A4,A5,A6,A7,Rm,TEMP
      REAL*8,DIMENSION((ny-2*nG(1))*(nx-2*nG(2))*(nz-2*nG(3)))::RHS
      REAL*8,DIMENSION((ny-2*nG(1))*(nx-2*nG(2))*(nz-2*nG(3)),7)::M
      REAL*8,DIMENSION(:,:,:),Pointer:: hNew=>NULL(),hE=>NULL(),hW=>NULL()
      REAL*8,DIMENSION(:,:,:),Pointer:: hU=>NULL(),hD=>NULL()
      REAL*8,DIMENSION(:,:,:),Pointer:: hNN=>NULL(),hS=>NULL(),DRR=>NULL()
      INTEGER N1,NM,M1,MM,NNY,NNX,NN,NNP,n,iter,itmax,np,IND(2),K1,KM
      REAL*8 err,meanW,meanH0,meanH,meanE,H0(ny,nx,nz),DP2,PK(2)
      REAL*8,DIMENSION(:),POINTER::pp1=>NULL(),pp2=>NULL(),pp3=>NULL()
      REAL*8,DIMENSION(:),POINTER::pp4=>NULL(),pp5=>NULL() ! print 
      REAL*8,DIMENSION(:),ALLOCATABLE::pk1,pk2,pk3,pk4,pk5
      REAL*8,DIMENSION(3,3):: TEMPI1,TEMPI2
      REAL*8,DIMENSION(:,:),pointer:: output
      !meanW = sum(W*Mapp)/size(h); meanH0=sum((H-E)*S*Mapp)/size(h) ! print *, debugging
      output => gdG .G. 'GW' .FF. 'output'
      
      H0 = H
      
      nGy = int(nG(1)); nGx = int(nG(2)); nGz = int(nG(3))
      DY  = DD(1); DX = DD(2)
      N1=nGy+1;NM=ny-nGy;M1=nGx+1;MM=nx-nGx; K1=nGz + 1;KM=nz - nGz !@@@@@@@@@@@
      NNY = ny-2*nG(1); NNX=nx-2*nG(2); NN = NNY*NNX; n=NN*(KM-K1+1) !@@@@@@@@@@@
      Dperc2 = 0D0 !@@@@@@@@@@@
      ! This solver itself is indiscriminate about the method used for the modeling. 
      ! The actually method used is determined from discretization and Boundary Condition
      DO I = 1,nz
        WHERE (Mapp<=0) K(:,:,i)=0.0D0
      ENDDO
      KZ(:,:,1) = 0D0
      KZ(:,:,nz)= 0D0
      
      !IF (ATYPE .eq. 0) THEN
      !ELSEIF (ATYPE .eq. 1) THEN
      !  D = h - E  ! Linearize Transmissivity Term
      !ENDIF
      HH = h ! Make a copy
      
      ! DXHalf = DZ(i+1/2); DYHalf = DZ(j+1/2); 
      !DZYHalf = (D + CSHIFT(D,dim=1,shift=1))/2D0
      !DZXHalf = (D + CSHIFT(D,dim=2,shift=1))/2D0
      DZZHalf = (D + CSHIFT(D,dim=3,shift=1))/2D0
      ! Apply Inactive Cell Condition, if any more than indicated by MAPP
      CALL applyCond(0,BBC,NBC,n,7,A1,A2,A3,A4,A5,M,RHS,K,HH)
      
      T = D * K
      
      ! Assuming DX is uniform
      ! K(i+1/2)=T(i+1/2)/DX(i+1/2) ! to account for varying DZ
      ! Tempx in the end is K(i+1/2) / (DX*DX(i+1/2)) == T(i+1/2) / (DX^2 * DZ(i+1/2)) = (2*T(i)*T(i+1))/[(T(i)+T(i+1)) * DX^2 * DZ(i+1/2)] 
      ! Tempz is just KZ(i+1/2) / (DZ * DZ(i+1/2)) because DZ(i,j,k)==DZ(i,j,k+1) for all (i,j)
      Tempx  = (T + CSHIFT(T,dim=2,shift=1)) * D
      Thalfx = (T*CSHIFT(T,dim=2,shift=1)/Tempx)
      
      Tempy  = (T + CSHIFT(T,dim=1,shift=1)) * D
      Thalfy = (T*CSHIFT(T,dim=1,shift=1)/Tempy)
      
      Tempz  = (KZ + CSHIFT(KZ,dim=3,shift=1)) * DZZHalf * D
      Thalfz = (KZ*CSHIFT(KZ,dim=3,shift=1)/Tempz) 
      
      WHERE (Tempx<=0) Thalfx = 0.0D0
      WHERE (Tempy<=0) Thalfy = 0.0D0
      WHERE (Tempz<=0) Thalfz = 0.0D0
      
      Tempx = Thalfx*(2D0/DX ** 2D0)
      Tempy = Thalfy*(2D0/DY ** 2D0)
      Tempz = Thalfz*(2D0)
      
      ! Calculate Coefficients
      ! Here A1,A2,A3,A4,A5 meaning have been updated and it's different from Matlab code GW_sol_fd
      ! A2 is relation to its SN, A3 to its NN, A4 to its WN, A5 to its EN
      ! This is consistent with Ghost Cell numbering
      ! Note that Tempx, etc are interface (i+1/2) values
      A2 = -Tempy(N1-1:NM-1,M1:MM,K1:KM)
      A3 = -Tempy(N1:NM,M1:MM,K1:KM)
      A4 = -Tempx(N1:NM,M1-1:MM-1,K1:KM)
      A5 = -Tempx(N1:NM,M1:MM,K1:KM)
      A6 = -Tempz(N1:NM,M1:MM,K1-1:KM-1)
      A7 = -Tempz(N1:NM,M1:MM,K1:KM)
      A1b= -(A2+A3+A4+A5+A6+A7)
      Temp = Ss(N1:NM,M1:MM,K1:KM)/dt
      A1 = Temp + A1b
      Rm  = Temp*HH(N1:NM,M1:MM,K1:KM) + Wm(N1:NM,M1:MM,K1:KM)
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
      CALL elem(N,A6,M(:,1))
      CALL elem(N,A4,M(:,2))
      CALL elem(N,A2,M(:,3))
      CALL elem(N,A1,M(:,4))
      CALL elem(N,A3,M(:,5))
      CALL elem(N,A5,M(:,6))
      CALL elem(N,A7,M(:,7))
      
      ! Flux BC should be converted to a source term and a no flow BC
      call applyCond(1,BBC,NBC,n,7,A1,A2,A3,A4,A5,M,Rm,K,HH)
      ! Apply Lower/Upper Boundary conditions
      call applyCond(2,BBC,NBC,n,7,A1,A2,A3,A4,A5,M,Rm,K,HH)
      ! Assembly Final Matrices, Save to rank 1 format
      
      CALL elem(N,Rm,RHS)
      
      ! Solve equation
      MD = (/-NN, -NNY,-1,0,1,NNY, NN/)
      output(:,1:7)=M(:,1:7); output(:,8)=RHS ! DebuggingTool
      itmax = 100
      hNew => HH(N1:NM,M1:MM,K1:KM) ! Use hNew as a pointer to reduce copying data
      !DO I = N1,NM
      !DO J = M1,MM
      ! if (ISNAN(W(I,J))) THEN
      ! ENDIF
      !ENDDO
      !ENDDO
      CALL pcg(n,size(MD),M,MD,RHS,hNew,4,1d-8,itmax,iter,err) !hNew is hydraulic head
      !hNN => HH(8:10,15:17) ! print debug
      !TEMPI1 =  hNN - hE
      !TEMPI2 = hW*dt/S(8:10,15:17) ! print debug
      !hN(N1:NM,M1:MM) = reshape(hNew,(/NNY,NNX/))-E(N1:NM,M1:MM)
      hN = HH
      ! Compute Lateral Fluxes
      hE => HH(N1:NM,M1+1:MM+1,K1:KM)
      hW => HH(N1:NM,M1-1:MM-1,K1:KM)
      hNN=> HH(N1+1:NM+1,M1:MM,K1:KM)
      hS => HH(N1-1:NM-1,M1:MM,K1:KM)
      hU => HH(N1:NM,M1:MM,K1-1:KM-1)
      hD => HH(N1:NM,M1:MM,K1+1:KM+1) ! Vertically lower layer has higher index
      hNew=> HH(N1:NM,M1:MM,K1:KM)
      DRR=> DR(N1:NM,M1:MM,K1:KM)
      DRR= -(A2*(hS-hNew)+A3*(hNN-hNew)+A4*(hW-hNew)+A5*(hE-hNew)+A6*(hU-hNew)+A7*(hD-hNew))
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
      END SUBROUTINE Aquifer3D

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
      REAL*8,DIMENSION(:,:),Pointer:: hNew=>NULL(),hE=>NULL(),hW=>NULL()
      REAL*8,DIMENSION(:,:),Pointer:: hNN=>NULL(),hS=>NULL(),DRR=>NULL()
      INTEGER N1,NM,M1,MM,NNY,NNX,NN,n,iter,itmax,np,IND(2),IIT
      INTEGER I1,J1
      REAL*8 err,meanW,meanH0,meanH,meanE,H0(ny,nx),DP2,PK(2),errP,DL
      REAL*8,DIMENSION(:),POINTER::pp1=>NULL(),pp2=>NULL(),pp3=>NULL()
      REAL*8,DIMENSION(:),POINTER::pp4=>NULL(),pp5=>NULL() ! print 
      REAL*8,DIMENSION(:),ALLOCATABLE::pk1,pk2,pk3,pk4,pk5
      REAL*8,DIMENSION(3,3):: TEMPI1,TEMPI2
      REAL*8,DIMENSION(:,:),POINTER::ES=>NULL(),RP=>NULL()
      REAL*8,DIMENSION(:),POINTER::NL=>NULL()
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
      
      SUBROUTINE fixedHeadGWupdate(ST,m,h,h0,flux)
      INTEGER m(2)
      REAL*8,DIMENSION(m(1),m(2)),intent(in)::ST,h0
      REAL*8,DIMENSION(m(1),m(2)),intent(inout)::flux,h
      LOGICAL*1,DIMENSION(m(1),m(2)) :: fixBCmap
      
      flux = 0
      WHERE(ST>100)
        flux = (h - h0)*ST
        h = h0
      END WHERE
      
      END SUBROUTINE

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

      SUBROUTINE applyCond(TYP,BBC,NBC,n,nC,A1,A2,A3,A4,A5,M,RHS,K,H)
      ! Apply boundary condition for the groundwater model
      ! DirichletBC does not support upper and lower boundaries so far therefore A6,A7 not passed in
      ! CP Shen
      USE VDATA
      integer TYP,n,NBC,nC
      TYPE(BC_type) :: BBC(NBC)
      REAL*8,DIMENSION(n)::RHS,A1,A2,A3,A4,A5,K,H,TEMP
      REAL*8 M(n,nC) ! nC: for 3D code, nC = 7, for 2D, nC = 5
      integer I,NI
      integer,DIMENSION(:),POINTER :: idxy=>NULL(),Cidx=>NULL()
      integer,DIMENSION(:),POINTER :: Nidx=>NULL(),Sidx=>NULL()
      integer,DIMENSION(:),POINTER :: Widx=>NULL(),Eidx=>NULL()
      REAL*8, POINTER:: P=>NULL()
      INTEGER :: iU,iW,iS,iC,iN,iE,iL
      ! M=[A4 A2 A1 A3 A5]
      ! M=[A6 A4 A2 A1 A3 A5 A7] in 3D
      IF (nC .eq. 7) THEN
        iU = 1; iW = 2; iS = 3; iC = 4; iN = 5; iE = 6; iL = 7
      ELSE 
        iU = 0; iW = 1; iS = 2; iC = 3; iN = 4; iE = 5; iL = 0
      ENDIF

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
             RHS(Sidx)=RHS(Sidx)-A3(Sidx)*BBC(I)%h;M(Sidx,iN)=0.D0 ! SN
             RHS(Nidx)=RHS(Nidx)-A2(Nidx)*BBC(I)%h;M(Nidx,iS)=0.D0 ! NN
             RHS(Widx)=RHS(Widx)-A5(Widx)*BBC(I)%h;M(Widx,iE)=0.D0 ! WN
             RHS(Eidx)=RHS(Eidx)-A4(Eidx)*BBC(I)%h;M(Eidx,iW)=0.D0 ! EN
             M(Cidx,:)=0.0D0; M(Cidx,3)=1.0D0
             RHS(Cidx)=BBC(I)%h
            ! Boundary Cells:
            !! CURRENTLY SOME PROBLEMS
            !! IF GHOST CELL IS MARKED 0 in MAPP, their K will be set 0
            !! THEN the A values are all wrong!
            case (1) ! This is a South Boundary, Correct its North Neighbor (To NN's S)
             RHS(Nidx)=RHS(Nidx)-A2(Nidx)*BBC(I)%h
             !! BUG CATCH! The following was not added in original code
             ! in 2D, that is ok, because E/W boundaries automatically goes out of matrix bound
             ! in 3D, this will cause trouble as linkage to the East neighbor can be mistaken as linkage to Lower neighbor 
             !    as the start of low layer follows from the end of the upper layer
             ! Actually, even in 2D, the missing part is a bug for S,N BCs
             M(Nidx,iS)=0D0
             !A2(Nidx) = 0D0 ! Modifying M instead of A's in order to calculate DR in the end
            case (2) ! This is a Nouth Boundary, Correct its South Neighbor (To SN's N)
             RHS(Sidx)=RHS(Sidx)-A3(Sidx)*BBC(I)%h
             M(Sidx,iN)=0D0
             !A3(Sidx) = 0D0
            case (3) ! This is a West Boundary, Correct its East Neighbor (To EN's W)
             RHS(Eidx)=RHS(Eidx)-A4(Eidx)*BBC(I)%h
             M(Eidx,iW)=0D0
             !A4(Eidx) = 0D0
            case (4) ! This is a East Boundary, Correct its West Neighbor (To WN's E)
             RHS(Widx)=RHS(Widx)-A5(Widx)*BBC(I)%h
             M(Widx,iE)=0D0
             !A5(Widx) = 0D0
           end select
          CASE (2) ! Upper/Lower Layers Implicit Condition
          ! CONSIDER ARRAY SIZES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !TEMP = 0D0
          NI = size(BBC(I)%K)
             idxy=>BBC(I)%Cidx
             TEMP(1:NI) = BBC(I)%K/BBC(I)%D
             RHS(idxy)=RHS(idxy)+TEMP(1:NI)*BBC(I)%h
             !A1(idxy)=A1(idxy)+TEMP(1:NI)
             M(idxy,iC)=M(idxy,iC)+TEMP(1:NI)
         END SELECT
       ENDIF
      ENDDO
      END SUBROUTINE applyCond

      SUBROUTINE vdz1c(NZ,NA,E,DF,Dperc,h,THE,K,CI,SRC,ddt,dtP,                &
     &   MAXIT,DRAIN,WTL,THER,THES,Lambda,alpha,n,KS,DZ,DDZ,BBCI,                  &
     &   Rm,STATE,GWL,CC,PRatio,ho,SN,Rff,m_err,ICE,TN,TNstate,ice_alpha,      &
     &   latOptionIn,LatData)
      ! Vadose Zone Model using Celia Picard iteration scheme
      ! CP Shen. CEE, MSU
      ! PDE Unit: T^-1
      ! The unit doesn't have to be meter-day, as long as they are consistent with themselves that is fine-
      ! Boundary Condition is implemented with a ghost cell on the top or in the bottom
      ! For Equation BC the top layer is included in the computation.
      ! modified h, THE, K (unimportant), WTL, GWL, Rff, DF, Dperc
      ! Added 11/02/10 freezing algorithm

       USE VData, Only: BC_type,rtnobd,THOMAS,DEBUGG
       use displayMod, Only: display
       IMPLICIT NONE
       !
       ! !ARGUMENTS:
       ! ========= Input ========================> 
       integer NZ,NA                                        ! NZ: number of cells, NA: number of active cells
       REAL*8, intent(in) :: SN                            ! runoff coefficient, =(S0^0.5)/(ManningN*ldist). If turning off runoff is desired, set this to 0.
       REAL*8, intent(in) :: ho                            ! Ground surface interception (minimum depth for runoff) (m)
       REAL*8, intent(in) :: dtP(3)                        ! [minDT maxDT macroDT] (day)
       REAL*8, intent(in) :: PRatio                        ! Permeable fraction of the cell (1 minus impervious fraction)
       REAL*8, DIMENSION(NZ), intent(inout) :: ICE           ! ice content [kg/m2] ! ICE may be temporarily modified inside but will recover its value upon return. Overall will not change
       REAL*8, DIMENSION(NZ), intent(in) :: E             ! Cell elevation (m)
       REAL*8, DIMENSION(NZ), intent(in) :: THER          ! Residual water content (-)
       REAL*8, DIMENSION(NZ), intent(inout) :: THES       ! Saturated water content (-) ! THES may be temporarily modified inside but will recover its value upon return. Overall will not change
       REAL*8, DIMENSION(NZ), intent(in) :: DZ            ! Cell size (m)
       REAL*8, DIMENSION(NZ), intent(in) :: DDZ           ! Distance between i and i+1 cell centers (m)
       REAL*8, DIMENSION(NZ), intent(in) :: Lambda        ! van Genuchten tortuosity parameter (-)
       REAL*8, DIMENSION(NZ), intent(in) :: alpha         ! van Genuchten parameter (m^-1)
       REAL*8, DIMENSION(NZ), intent(in) :: n             ! van Genuchten parameter (-)
       REAL*8, intent(in)                 ::  ice_alpha     ! ice impermeable fraction parameter
       TYPE(BC_type), intent(in) :: BBCI(2)                 ! a derived type containing boundary conditions; 
       REAL*8, DIMENSION(NZ,5), intent(in),target,optional:: LatData  ! Lateral flow data. See below for explanations
       Integer, intent(in),target,optional:: LatOptionIn  ! Lateral flow data. See below for explanations
       
       !BBCI(2)%BV: boundary condition values. explained: 
       ! [1 ST,2 zl (=DZ(uc)/2+D(aquitard)),3 DR, 4, Kl, 5, hl(pressure head), 6 E(center) uc aquifer, 7,GW.DZ, 8 EB(1)]
       ! <======== InOut ========================>
       REAL*8, DIMENSION(NZ), intent(inout) :: KS         ! Saturated hydraulic conductivity (m/day) -- maybe modified during calculation due to freezing-- but changed back
       REAL*8, intent(inout) :: DF                         ! Accumulated Infiltration (m)
       REAL*8, intent(inout) :: Dperc                      ! Accumulated percolation (to unconfined aquifer) (m)
       REAL*8, intent(inout) :: ddt                        ! an initial guess of time step (recorded from last time step) (day)
       INTEGER, intent(inout):: WTL                        ! the first layer whose top is below water table. This is an initial guess, may be modified in this subroutine
       REAL*8, DIMENSION(NZ), intent(inout) :: SRC        ! Source term (positive for inflow, for the first cell: m/day; for other cells, day^-1)
       REAL*8, DIMENSION(NZ), intent(inout) :: h          ! Pressure head (m)
       REAL*8, intent(inout) :: m_err                      ! mass balance error (m)
       REAL*8, intent(inout) :: Rff                        ! Accumulated Runoff (m)
       REAL*8, DIMENSION(NZ,4), intent(inout),target :: TN    ! Transient numbers -- some temporary calculations saved to cut down computation
       INTEGER TNState                                      ! the state of these calculations --- whether to do them or not
       ! <======== Out      ======================
       REAL*8, intent(out) :: Rm(*)                        ! use R to record steps solved by the botdrain and saturation status
       REAL*8, intent(out) :: GWL                          ! Groundwater level (elevation, m) analyzed by soil moisture profile
       INTEGER, intent(out) :: State                       ! Initialization flag, no real meaning outside
       INTEGER, intent(out) :: MAXIT                       ! Maximum iteration (to reduce time step or to use botdrain), no real use outside
       REAL*8, DIMENSION(NZ), intent(out) :: K            ! Unsaturated Conductivity (m/day)
       REAL*8, DIMENSION(NZ), intent(out) :: THE          ! Water content (-)
       REAL*8, DIMENSION(NZ), intent(out) :: CI           ! Differential water capacity (m^-1)
       ! ========= NOT USED ======================
       REAL*8 :: DRAIN,CC
       ! ========= Local =========================
       LOGICAL:: SATH,BC_NC=.TRUE.,DEBUG                   ! SATH: whether column is close to saturation, BC_NC: whether BC is not changed in last fractional time step
       LOGICAL:: CVG,ICED                                  ! CVG: whether convergence is reached last time step; ICED: whether there is ice in the column
       REAL*8,DIMENSION(nz):: THESSAVE,KSSAVE             ! Local saves for recovery after freezing adjustment
       REAL*8,DIMENSION(nz):: HSAVE2                      ! HTEMP for the sake of deciding if flip-flopping is happening
       REAL*8,DIMENSION(nz):: fracice                     ! fractional impermeable area due to freezing
       REAL*8,DIMENSION(nz):: theta_ice                   ! ice volume content scaled by water density
       REAL*8,DIMENSION(nz):: ice_extra                   ! is it possible for combined the to be greater than thes?
       INTEGER,DIMENSION(nz):: filter_frz                 ! filter with freezing layer indices
       INTEGER num_frz                                     ! number of frozen layers
       REAL*8 Dperc2,hh4,denice
       PARAMETER(DEBUG = .TRUE.)                           ! CHANGE THIS TO ACTIVATE MASS CONSERVATION DEBUGGING
       !PARAMETER(ice_alpha = 3D0)                           ! scale-dependent freezing fraction coefficient
       PARAMETER(denice = 917D0)                           ! ice density (kg/m3)
       REAL*8 :: DR                                         ! Regional groundwater drainage term (positive inflow day^-1)
       REAL*8, DIMENSION(NZ) :: THESAVE                   ! A local copy of original THE profile
       REAL*8  err, h_err, t, GWL0, RR, Rf                  ! error indices; Rff is passed in, Rf is local
       REAL*8, DIMENSION(NZ) :: Khalf                     ! Unsat hydraulic conductivity at interface i+1/2
       REAL*8, DIMENSION(NZ) :: temp                      ! =Khalf/DDZ
       REAL*8, DIMENSION(NZ) :: a,b,c,d                   ! Tri-diagonal matrix coefficients. a(i)X(i-1)+b(i)X(i)+c(i)X(i+1)=d(i)
       REAL*8, DIMENSION(NZ) :: THEP,CP,HP                ! Values saved from last iteration
       REAL*8,POINTER:: BV(:)=>NULL()                              ! Points to Lower boundary values. 
       REAL*8 :: x,xl                                       ! x:distance from upper boundary of bottom cell to water table, xl: saturated thickness of the soil column (m)
       REAL*8, DIMENSION(:),POINTER:: tn1=>NULL(),tn2=>NULL(),tn4=>NULL()        ! temporary calculations, tn1 = n/(n-1) = 1/m; tn2 = m = (n-1)/n= 1-1/n; tn4 = (2n-1)/n = 1+ m
       integer :: ITER, NBC, N1,NM, JType1, JType2,NMAX
       REAL*8,DIMENSION(nz) :: HTEMP,THETEMP, temp1, temp2           &
     &                         ,tn3,tn5,DTHES,f,THEOLD      ! temporaries
       REAL*8 :: PREC , temp3(NZ-2), temp4, tempDF, DTHE, THE0, THE2,m
       TYPE(BC_type) :: BBC(2)  ! Boundary Condition!
       !ddt(1):current dt, ddt(2):dtmin, ddt(3):dtmax, ddt(4):MacroTimeStep(before going out of program)
       REAL*8 :: dt, RDT, dtmin, dtmax, SDT,tf,S(NZ),WTD,WTDL,SS(NZ)
       REAL*8 :: dtPREC,HH,hh1,hh2,hh3,t1,t2,t3,q1,FULL,FULL1,km,EUC
       REAL*8 :: ERRP,RS1,RS2,SS1,m0,qbot,qtop,qm,DFMX,Imx,P(6)
       integer I,J,IMAX,ITOTAL,MAXIT2,NM1 !,IWTL,NM2,WTL0,W1,W2,EN1
       integer ROUTE(15),sigBC,BC_OLD,NAA,WTL0,totalK,ki
       REAL*8 ::m_err2,THE1,DFF,h0(nz),SF,Dpercc,tt,Qin,GST,sumICE
       REAL*8 :: evapStemp, layerEvapTemp
       REAL*8 :: BC1,BC2,BC3,BC4,THEB,hNZ,DZB,gm(10),sumS,ZT,KT ! g is temporary storage variable
       LOGICAL :: PRE_FLUX,HnegS_F
       REAL*8 :: pflux,m_errf(200),m_err2f(200),hnegSum,hnegS
       REAL*8,PARAMETER :: THER_MIN = 1D-4  ! minimum distance to THER, sometimes due to freezing, THE may drop below THER, which may cause program explode
       REAL*8 :: THE_ICE_fix(NZ)
       ! Below added for new lateral flow coupling
       REAL*8, DIMENSION(NZ) :: latDN1      !  d^{n+1}
       REAL*8,POINTER:: latQn(:),latBeta(:),latHn(:),latTao1(:),latTao2(:) ! qn
       INTEGER :: latOption = 1 ! == 1, 2D Dupuit Assumption; ==2, 3D
       
    !debugg = .true.
       ! Lateral flow data. [LatOption; 1. beta; 2. tao1; 3. tao2; 4. q^n; 5. h^{star,n};]
       IF (PRESENT(LatOptionIn)) THEN
           LatOption = LatOptionIn
           if (LatOption .eq. 2) then
             latBeta   => LatData(:,1)
             latTao1   => LatData(:,2)
             latTao2   => LatData(:,3)
             latQn   => LatData(:,4)
             !latHn   => LatData(:,5) 
           endif
       ELSE
           LatOption = 1
       ENDIF
       
       ! f is the fraction of the Dupuit flow

       !if (h(2)<-30 .OR. h(2)>0 .OR. SRC(1)>1) THEN
       !print *,3
       !endif

       !IF (SRC(2)<0) THEN
       !PRINT*,3
       !endif

       ERRP = 1.D-4
       BBC = BBCI ! Do not modify input BBC
       t = 0.0D0
       dt = ddt;dtmin=dtP(1); dtmax=dtP(2);tf = dtP(3);
       dtPREC=0.01*dtmin
       RDT = tf-t;
       MAXIT = 12
       MAXIT2= 13
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
       num_frz = 0 
       !THE_SAVE=THE
       !IWTL = 0;
       !hh1 = 0.D0; EN1 = 1;
       ICED = .FALSE.
       SATH = .FALSE.
       PRE_FLUX = .FALSE.
       CVG = .TRUE.
       Rf =0D0
       BC_OLD = 0
       sumICE = 0D0
       THEB = THE(NA)
       hNZ  = h(NZ)
       totalK = 0
       hnegSum = 0D0
       hnegS = 0D0
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
       IF (DEBUGG) THEN
       THE0 = sum(THE(2:NA-1)*DZ(2:NA-1))  !
       WTL0 = WTL
       GWL0 = GWL
       THESAVE = THE
       ENDIF

       tn1 => TN(:,1)
       tn2 => TN(:,2)
       tn4 => TN(:,4)
       DO I = 2,NZ-1
         IF (ICE(I)>0D0) THEN
           ICED = .TRUE.
           ! rescale THES and K
           ! Niu, Effects of Frozen soil on snowmelt runoff ..., 2006, JOURNAL OF HYDROMETEOROLOGY
           ! did not follow CLM4.0 here, used original parameterization
           ! What is passed in is THE_Water
           theta_ice(I) = ICE(I)/(1000D0*DZ(I)) !*1.09051254089422D0 ! 1.09 is water:ice density ratio
           ! formulation 1: add ice to theta
           !THE(I) = THE(I) + theta_ice(I) ! the= the_w + the_ice * denice/denw, so, just use weight to be divided by denw
           
            
           IF (THE(I)>THES(I)) THEN
             ICE_EXTRA(I) = THE(I) - THES(I)
             THE(I) = THES(I)
           ELSE
             ICE_EXTRA(I) = 0D0
           ENDIF
           
!           THE_ICE_fix(I) = max(0D0, THER_MIN+THER(I)-THE(I))
!           IF (THE_ICE_fix(I)>0D0) THEN
!             THE(I) = THE(I) + THE_ICE_FIX(I)
!           ENDIF
           
           ! In formulation 2, below need not to be done
!           SS(I)=(THE(I)-THER(I))/DTHES(I)
!           h(I)   = ((SS(I)**(-1.D0/tn2(I))-1D0)**(1.D0/n(I)))/(-abs(alpha(I)))
           num_frz = num_frz + 1
           filter_frz(num_frz) = I
           KSSAVE(I) = KS(I)
           fracice(I) = (exp(-ice_alpha*(1D0-theta_ice(I)/THES(I)))-exp(-ice_alpha))/(1D0-exp(-ice_alpha)) ! rescale it to one.
           IF (THE(I) <= THER(I)+THER_MIN) THEN
             KS(I) = 0D0
           ELSE
             KS(I)  = max(0D0,(1D0-fracice(I))*KS(I))
           ENDIF
           sumICE = sumICE + (theta_ice(I)-ice_extra(I))*DZ(I)
           
           ! formulation 2: do not add theta_ice(I)
           THESSAVE(I) = THES(I)
           THES(I)     = max(THE(I),THES(I) - theta_ice(I),THER(I)+PREC) 
         ENDIF
       ENDDO
       
       
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
!       
!       IF (TNState .eq. 0) then
!         tn1 = n/(n-1D0)      ! 
!         tn2 = 1D0/tn1        ! 1-1/n, sometimes called m
!         tn4 = (2D0*n-1D0)/n
!         TNState = 1
!       ENDIF
       
       DTHES = THES-THER
       S = SRC
       NAA = NA
       
       ! ########11111111111111111111111111111111111111111111111111111111111111111111111111
       ! time loop, input h, SRC, mass balance must be conserved
       DO WHILE (t < tf-1e-12)
       totalK = totalK + 1 ! number of fractional steps
       RDT = tf-t;
       dt = min(dt,RDT,dtmax)
       dt = max(dt,min(tf-t,dtmin))
       !THEOLD = THE ! save for computing recharge
       !S(1:NA) = SRC(1:NA);  !temporary source term storage for the fractional time steps
       S(1) = SRC(1)  ! only the first one may get changed
       
       PRE_FLUX = .FALSE.
       HnegS_F  = .FALSE.
       pflux = 0D0    ! pflux is the flux that has already been used as infiltration in last step (when choosing route(5)), but it could be unsuccessful.
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
         !tn1 = n/(n-1) = 1/m; tn2 = m = (n-1)/n= 1-1/n; tn4 = (2n-1)/n = 1+ m
         SS = (THE - THER)/(THES - THER)
         !m = 1 - 1/n
      CALL vGh(h(2:NM1),ss(2:NM1),n(2:NM1),alpha(2:NM1),tn1(2:NM1),NA-2)
         DO I = 2,NM
           IF (THE(I)<THES(I)) THEN
           !SS=(THE(I)-THER(I))/(THES(I)-THER(I))
           !m = 1D0 - 1D0/n(I);
           !h(I)=((SS**(-1.D0/m)-1D0)**(1.D0/n(I)))/(-abs(alpha(I)))
           ELSEIF (h(I-1)>=0D0) THEN
            h(I)=h(I-1)+DDZ(I-1)
           ELSE
            h(I)=0D0
           ENDIF
         ENDDO
         ! BELOW IS WITH h initialization
         temp1 = abs(alpha*h)**n
         temp2 = 1+temp1
         tn3 = min(max((THE-THER)/(DTHES),PREC),1D0)
         !THE = THER + (DTHES)/(temp2**(tn2))
         CI = -DTHES/(temp2 **tn4)                                          &
     &    * temp1 *(n-1D0)/h
         !CI = -(THES-THER)/(temp2 **((2*n-1)/n))                                   &
         !&    * temp1 *(n-1)/h
       else
         STATE = 0 ! unless made 1 in the convergence exit
       ENDIF
       K = KS*(tn3)**Lambda*                                                &
     &    (1D0-(1D0-(tn3)**tn1)**tn2)**2.0D0
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
       
       ! print debug
       if (DEBUGG .and. any(K<0D0)) THEN
         call display( 'VDZ1C negative K')
         PAUSE
         call display( 'VDZ1C negative K')
         PAUSE
       endif

       THEP = THE
       CP = CI
       
       ! for debug
       !THE2 = sum(THE(2:NA-1)*DZ(2:NA-1))

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
       Khalf(1) = sqrt(K(1) * K(2))

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
         BV => BBC(2)%BV
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
       qtop = RS1;
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
          IF (BBC(2)%Jtype .eq. 2) then
            qbot = BBC(2)%BData
          else 
            qbot = Khalf(NM1)*((h(NM1)-h(NM1+1))/DZB+1D0)
          endif
          qbot = min(qbot,qbot/20D0) ! flux limiter for positive flux
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
              !IF (t .eq. 0D0 .and. qtop>-2D-1) THEN
              !IF (t .eq. 0D0 .and. qtop>-1) THEN  ! old code, why am I doing t .eq. 0D0??
              IF (qtop>-1 .or. qtop*dt>-0.005) THEN
                 BBC(1)%Bdata = qtop
              ELSE
              ! once it gets inside there would be mass imbalance
              ! some times it is caused by last fractional step incorrectly using BBC=3, thus resulting in negative flux
                 hnegS    = qtop*dt
                 hnegSum  = hnegSum + hnegS
                 HnegS_F  = .TRUE.
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
                if (d(2)>0 .and. rr>500) then
                 ! may lead to convergence issue, now we first infiltrate that water
                 THE(2)=THE(2)+qtop*dt/dz(2)
                 t2 = (THE(2)-THER(2))/(THES(2)-THER(2))
                 t3 = -N(2)/(N(2)-1D0)
                 h(2)=((t2**t3-1D0)**(1D0/n(2)))/(-abs(alpha(2)))
                 BBC(1)%Bdata = 0D0
                 pflux = qtop*dt
                 DF = DF + pflux
                 !DF = DF + BBC(1)%BData*dt ! This is an accounting error 11/05/10
                 ROUTE(5)=ROUTE(5)+1
                 PRE_FLUX = .TRUE.
                else
                 BBC(1)%Bdata = qtop
                endif
                BBC(1)%JType = 2
                ROUTE(8)=ROUTE(8)+1
              ENDIF
            ENDIF
          ENDIF
       ENDIF
       !IF (BBC(1)%JType .eq. 2 .and. h(1)>0D0) THEN ! This is the code before CLM
       ! With CLM, I am allowing negative fluxes en masse, therefore, h(1)>0D0 must be removed to keep mass balance
       ! If not converged, h(1) will be recovered from hh1
       IF (BBC(1)%JType .eq. 2) THEN
         S(1)=0D0
         BBC(1)%BData2 = h(1)
         HP(1) = 0D0; h(1) = 0D0
         
         IF (BBC(1)%BData < 0D0 .and. abs(THE(2)-THES(2))<1D-4) THEN
            DO I = 3,NM1
              IF (abs(THE(I)-THES(I))>1D-3) EXIT
            ENDDO
            KI = I
            
            IF (K(KI)<1D-12) THEN
            ! This is the really difficult scenario that may cause NaN in the solver
            ! It happens very rarely. can be fixed with specific storage at saturation later on
            
               evapStemp = -BBC(1)%BData * dt
               DO I = 2,KI - 1
                  IF (((E(1) - E(I))< 0.2D0 .OR. KI .eq. 3) .AND. evapStemp > 1D-15) THEN
                     layerEvapTemp = max(min(evapStemp, (THES(I)- THER(I))*DZ(I)*0.3D0 ),0D0)
                     THE(I) = THE(I) - layerEvapTemp/DZ(I)
                     SS(I) = (THE(I) - THER(I))/(THES(I) - THER(I))
                      !m = 1 - 1/n
                     CALL vGh(h(I),ss(I),n(I),alpha(I),tn1(I))
                     temp1(I) = abs(alpha(I)*h(I))**n(I)
                     temp2(I) = 1+temp1(I)
                     tn3(I) = min((THE(I)-THER(I))/(DTHES(I)),1D0)
                     CP(I) = -DTHES(I)/(temp2(I) **tn4(I)) * temp1(I) *(n(I)-1D0)/h(I)
                     evapStemp = evapStemp - layerEvapTemp
                  ENDIF
               ENDDO
               BBC(1)%BData = evapStemp/ dt
               
            ENDIF
         ENDIF
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
             IF (BBC(1)%BData2>=1 .and. BBC(1)%BData2<=NZ-1) THEN
             N1 = FLOOR(BBC(1)%BData2)+1 ! can be cells other than 2, to work with GGA
              IF (N1 > 2) THEN
                t1 = sum(S(2:N1-1)*DZ(2:N1-1))
                S(N1) = S(N1)+t1/DZ(N1)
                S(2:N1-1) = 0D0
              ENDIF
             ELSE
             N1 = 2
             call display( 'Bug with GGA?')
             ENDIF
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
         
         ! LatOption == 2, 3D lateral flow: Alternating Direction Implicit, expressed as weights on diagonals and RHS
         IF (LatOption .eq. 2) THEN
             b = b + latBeta * (1-latTao1)
             d = d + LatQn + LatBeta*LatTao2
         ENDIF
         

         ! Solve EQ
         CALL THOMAS(NZ,N1,NM,a,b,c,d,HTEMP)
         
         IF (ITER >= 5) THEN
!         ! suspected of flip-flopping
!           IF (BBC(1)%JType .eq. 2) THEN
!             ! Check for upper boundary flip-flopping
!             IF (HTEMP(2)>100D0 .OR. HP(2)>100D0) THEN
!!               DO I = 2,4
!!                 IF (HTEMP(I)>100D0) THEN
!!                   HTEMP = min(HTEMP/1000D0, 0.1D0)
!!                 ELSEIF (HTEMP(I)<-1000D0) THEN
!!                   HTEMP(I) = max(-0.04D0,HTEMP(I)/30D0)
!!                 ENDIF
!!               ENDDO
!             ENDIF
!           ENDIF
           
           ! lower boundary helping
           DO I = 2,NM
             IF (HP(I)>0D0 .AND. HTEMP(I)<=0D0.AND. HSAVE2(I)<=0D0) THEN   !                    &
      !&         .AND. ABS(HTEMP(I)-HSAVE2(I))<1D-3) THEN
             ! prime suspects of flip-flopping. Try to help it
             ! normally, this result from a difficulty with large cells near saturation,
             ! that a reflection from saturation always points back to too large negative head
               !H(1)=H(1)
               HTEMP(I) = max(-0.04D0,HTEMP(I)/30D0)
             ENDIF
           ENDDO
         ENDIF
         
         ! Error Convergence
         temp1 = abs(alpha*HTEMP)**n
         temp2 = 1+temp1
         THETEMP = THER + DTHES/(temp2**tn2)
         WHERE (HTEMP > -1e-5) THETEMP = THES
         !temp3=THETEMP(N1:NM) - THEP(N1:NM)
         err = sum(abs(THETEMP - THEP))
        !temp = HTEMP - HP
         h_err = sum(abs(HTEMP(N1:NM) - HP(N1:NM)))
         IF (ITER >=4) HSAVE2 = HP
         HP = HTEMP
         THEP = THETEMP
         CP = - DTHES/(temp2 **tn4) * temp1 *(n-1)/HP
         !CP = -(THES-THER)/(temp2 **((2.D0*n-1.D0)/n))                              &
         !&      * temp1 *(n-1)/HP
         if (err<1D-14) h_err = err  ! when it is very dry and err approaches this small value. round-off error can throw h_err off. No need to pursue that
         
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
           !Khalf(1)=0.5D0*(K(1)+K(2))
           Khalf(1) = sqrt(K(1) * K(2))
         ENDIF

       ENDDO
       ! ##############222222222222222222222222222222222222222222222222222222222222

       IF (h_err>errp .AND. ITER >= MAXIT .AND. dt>dtmin) THEN ! Adjust Time Step***********
         dt = max(dtmin,dt/3.0D0)
         SDT = dt
         !IF (BBC(1)%JType .eq. 2 .AND. MMASS) h(1)=BBC(1)%BData2
         ROUTE(1)=ROUTE(1)+1
         CVG = .FALSE.
         !IF (MMass) h(1:NA) = h0(1:NA) ! return to old states 'coz may have been changed by DFMX step
         IF (PRE_FLUX) THEN
         ! then this pflux need to undone
           THE(2) = THE(2) - pflux/DZ(2)
           DF     = DF - pflux
         ENDIF 
         IF (HnegS_F) THEN
           hnegSum= hnegSum - hnegS
         ENDIF
       ELSEIF (h_err < errp) THEN
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
         state = 1 ! no need to do h->THE or THE-> next fractional time step
         SELECT CASE(BBC(1)%JType)
         CASE(2)
           ! GGA BBC should be differently treated
           IF (JType1 .eq. 5) THEN
             !h(1) = 0D0 ! If marched a step using flux BC, this ponding depth should be 0
             DFF= BBC(1)%BData*dt
             DF = DF + DFF
             BBC(1)%BData= 0.0D0 ! Free surface type BC, will drain in one time step
           ENDIF
         CASE(1)
           HH = BBC(1)%BData - E(1)
           DFF = Khalf(1)*t1*dt
           DF = DF + DFF  ! Don't forget gravity ! Don't forget the BC
         CASE(3)
           DFF= Khalf(1)*((h(1)-h(2))/DDZ(1)+1.0D0)*dt
           DF = DF + DFF
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
         Dpercc = Khalf(NA-1)*                                                &
     &  ((h(NA-1)-h(NA))/DZB+1.0D0)*dt-(BV(3)-DR*BV(7))*dt ! NA-1 is to accomodate BBC=3
         Dperc=Dperc+ Dpercc
         !   Dperc=DF-sum((THE(2:NMAX)-THEOLD(2:NMAX)-                            &
         !&           S(2:NMAX)*dt)*DZ(2:NMAX))+Dperc
         CASE(2)
         !ELSEIF (BBC(2)%JType .eq. 2) THEN
           Dpercc= BBC(2)%BData*dt
           Dperc = Dperc + Dpercc
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
     &          LAMBDA,S,DZ,DDZ,BBC,DFF,DPercc,dt,E,DZB,BV,tn1,tn2,tn4)
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
       ! Below are mass balance checks for each fractional step:
       if (DEBUGG) then
       THE2 = sum(THE(2:NA-1)*DZ(2:NA-1))
       DTHE = THE2 - (THE0+sumICE)
       t1 = sum(SRC(2:NAA-1)*DZ(2:NAA-1))*t
       m_errf(totalK) = DTHE-(DF-Dperc+t1)
       m_err2f(totalK) = h0(1) - DF + SRC(1)*t - Rf - h(1)
       if (abs(m_errf(totalK))>1D-6) then
        !write(*,*) 'm_errf',m_errf,'totalK',totalK
        m_errf(totalK) = m_errf(totalK)+0D0
       endif
       endif
       
 100   CONTINUE
       ENDDO !  TIME LOOP

       Rff = Rff + Rf
       h(1) = h(1)+hh3+ hnegSum
       THE(NA) = THEB
       h(NZ)   = hNZ ! This final state doesn't matter. It shouldn't be modified (occasionally it can be negative)
       if (h(1)<-1D-11 .AND. JType1 .eq. 5) THEN
         CALL h1negative(nz,na-1,h,THE,THER,THES,DZ,DF,ALPHA,LAMBDA,n)
         ROUTE(6)=ROUTE(6)+1
         K(1) = 0D0 ! debugging info
         K(2) = 1D0
       ENDIF
       CALL waterTable(na,h,THE,THES,DZ,0.95D0,WTL,x,BV(2)) !!! incorrect when WTL = 21
       GWL = E(WTL)+DZ(WTL)/2D0+x*DZ(WTL-1)
       
       
       DO J = 1,num_frz
         ! get back water and ice content
           I      = filter_frz(J)
           KS(I)  = KSSAVE(I)
           ! formulation 1: subtract ice from THE after done with flow
           !THE(I) = THE(I) + ICE_EXTRA(I) - theta_ice(I)
           
           ! formulation 2: do not add theta_ice
           THE(I) = THE(I) + ICE_EXTRA(I) !- THE_ICE_FIX(I)
           THES(I)= THESSAVE(I)
           
           IF (THE(I)< 0D0) THEN
           ! need a fix-up
             ICE(I) = ICE(I) + THE(I)*DZ(I)*1000D0
             THE(I) = 0D0
             print *, 'negative water content with freezing soil'
           ENDIF
       ENDDO
       
       !K
       !error indices
       !mask this portion out after debugging is complete
       ! DEBUGG is global 
       IF (DEBUGG) THEN
       if (JType1 .eq. 5) then
         dff = df
       else
         dff = 0D0 ! only for GGA sealed top
         hh2 = BBC(1)%JIndex ! h(1)_n in fact
       endif
       THE2 = sum(THE(2:NA-1)*DZ(2:NA-1))
       DTHE = THE2 - THE0
       t1 = sum(SRC(2:NAA-1)*DZ(2:NAA-1))*dtP(3)
       m_err = DTHE-(DFF-Dperc+t1)
       if (BBCI(1)%JType .eq. 2) THEN
         t2 = h(1)
         m_err2 = 0D0
       else
         t2 = h(1)
          m_err2 = hh2 - DF + ss1*dtP(3) - t2 - Rf
       endif
       ! IF BBC(1)=2 prescribed this need to add that as well.
       K(5) = m_err ! print debug
       K(6) = m_err2
       IF (abs(m_err)>1D-6 .OR. abs(m_err2)>1e-9                             &
       !IF (abs(m_err)>ERRP*0.1                                               &
     &       .OR. isnan(m_err)) THEN
       call display('VDZ m_err =',m_err,'m_err2 =',m_err2)
       endif
       
       if (any(THE<THER)) then
        call display("VDZ1c:: THE<THER generated")
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
      SUBROUTINE Runoff(RF,frac,ho,hc,dt,hI,h,VParm,nzl,hback,imp,typ,RfI,RfF)
      !  Computes runoff from impervious region
      ! calculates runoff from ponding to flow only if last argument > 0
       implicit none
       integer nzl
       real*8 RF,fk,ho,h(nzl),dt,hI,hc,hback,imp
       real*8 VParm(*)
       REAL*8 frac,f,k,sep,f2,fI,h0,typ
       real*8 RfI, RfF
       ! typ<0 -- only calculate run-on, not going to do runoff from h(1) to OVN.h (do it in vdz1c)
       !parameter(sep = 0.5D0)
       ! function to calculate runoff
       ! hc is h in flow domain, hI is h in impervious ponding, h is vadose zone h
       ! ho is initial abstraction of ponding layer
       ! hback is the barrier height for flow domain to flow back (an offset)
       ! RfI: runoff from Impervious fraction
       ! RfF: runoff from ponding domain to flow domain solved in this function
       ! which should be the difference between average elevation and bottom elevation of the cell
       !RETURN
       ! imp is the impervious cover (fraction). higher the imp, higher direct runoff (higher f2)
       !Rf = 0D0
       !k = SN*(VParm(10)**(2D0/3D0))
       !frac = exp(-k*dt)

       IF (frac <= 0D0) RETURN

       if (hI > 0D0) then
         fI = max((1D0 - frac)*max(hI-0.000,0D0),0D0) ! temporarily adding an arbitrary urban runoff abstraction
         !f2 = fI*imp ! fraction that goes to channel filaments
         ! I want to have less water going to nearby permeable zone (ponding zone), just make f2 larger:
         !f2 = FI * (imp**(1D0/3D0))
         f2 = FI !* (imp ** VParm(13))
         !hc= hc + f2
         h(1)=h(1)+(fI-f2)
         hI = hI - fI
         Rf = Rf + f2
         RfI= RfI + fI
       endif

       if (h(1)> ho .and. h(1)>hc-hback .and. typ>0D0) then
          h0 = max(h(1)-max(hc-hback,0D0),0D0)
          f = max((1D0 - frac)*h(1),0D0)
          f = min(f,(h(1)-hc+hback)/2D0,h(1)-ho)
          !hc= hc + f  ! g.OVN.h uses Rf as a source term, so, not adding it here
          h(1)=h(1)-f
          Rf = Rf + f
          RFF = Rff + f
       elseif (h(1)<hc-hback .and. hc>hback) then
          ! allow flowing back from flow domain to ponding domain??
          ! comment out if not allowed
          ! run-on
          f = 0.5D0*(hc - hback + h(1))
          h(1) = h(1) + f
          !hc = hc - f
          Rf = Rf - f
          RFF = Rff - f
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
       REAL*8,DIMENSION(:),POINTER :: gg=>NULL()
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
       !print *, 'lowland m_err' ! debug
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
       !DO WHILE (WTL<nz .AND. (.NOT. isnan(DZ(WTL+1))) ) ! bottom layer is groundwaters
       DO WHILE (WTL<nz)
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
     &            LAMBDA,S,DZ,DDZ,BBC,DDF,DPerc,dt,E,DZB,BV,tn1,tn2,tn4)
       USE Vdata, Only: BC_type,DEBUGG
       IMPLICIT NONE
       integer nz,na
       REAL*8,DIMENSION(nz):: h,THE,THES,KS,ALPHA,LAMBDA,KHALF,S,DDZ
       REAL*8,DIMENSION(nz):: THER,DZ,E,N,tn1,tn2,tn4
       REAL*8 DDF,Dperc,dt,m,RF,DZB,BV(*),C
       TYPE(BC_type) :: BBC(2)
       ! This is a linearized fail-safe bottom drainage model that RE can fall back on
       ! it is simplified as it only means to help the RE solver (occasionally) in case of convergence problems
       ! It computes flux from bottom up, the flux is the min of
       ! 1.explicit flux (limited) 2.saturation of receiving cell, 3.available moisture of providing cell, 4. uni-drection equilibrium state flux
       ! less accurate, but more stable
       REAL*8 f1,f2,f3,f4,KU,Q,FU,FL,THE0,THE1,test,SS,test2,SSS,hh2,h_eq,t1
       real*8 tt(5),LIMITER1,LIMITER2, LIMITER3
       PARAMETER(LIMITER1 = 10D0, LIMITER2 = 1D0, LIMITER3 = 2D0)
       integer I,J,NM,RC
       NM = NA-1
       hh2 = h(1)
       DDF = 0D0
       Dperc  = 0D0; RF= 0D0
       C = 0.2D0
       
       IF (DEBUGG) THEN
       !h(NA) = BBC(2)%Bdata - E(NA)
       THE0 = sum(THE(2:NA)*DZ(2:NA))
       ENDIF

       ! now we use an implicit formulation for botdrain for stability
       ! we assume h(NM-1) not changing, use implicit method
       ! [1 ST,2 dzl (=DZ(uc)/2+D(aquitard)),3 DR, 4, Kl, 5, hl(pressure head), 6 E(center) uc aquifer]
       
       IF (BBC(1)%JType .eq. 2 .and. BBC(1)%Bdata < 0D0) THEN
        ! This has not been handled before
        FU = min(abs(BBC(1)%BData*dt),sum(DZ(2:5)*(THE(2:5)-THER(2:5))*C))
        THE(2:5) = THE(2:5) - FU/sum(DZ(2:5))
        DDF = -FU
        IF (FU < abs(BBC(1)%BData*dt)) then
         ! there will be mass loss
          BBC(1)%BData = 0D0
        ELSE
          BBC(1)%BData = 0D0
        ENDIF
       ENDIF

       tt(1)=BV(1)/dt; tt(2)=Khalf(NA-1)/DZB;
       tt(3)=tt(1)+tt(2)
       tt(4)=tt(2)*h(NA-1)+BV(3)+tt(1)*h(NA)+Khalf(NA-1)
       h(NA)=tt(4)/tt(3)

       ! Bottom Flux
       I = NM
       f1 = Khalf(NM)*((h(NM)-h(NA))/DZB+1)
       !f1 = f1; ! empirical scaling parameter
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
         t1 = the(I)
         THE(I)=max(THE(I)-FL*dt/DZ(I)+S(I)*dt, THER(I)+1D-10)
         IF (THE(I)>THES(I)) THEN
           Q = (THE(I)-THES(I))*DZ(I)/dt
           THE(I) = THES(I)
           FU = -Q ! flux should go up
         ELSE
           !IF (FL > 0D0) THEN
           ! cell already gives water to below, be careful, don't provide up too much
             SS=(THE(I)-THER(I))/(THES(I)-THER(I))
             h(I)=((SS**(-tn1(I))-1D0)**(1.D0/n(I)))/(-abs(alpha(I)))
           !ENDIF
           f1 = Khalf(I-1)*((h(I-1)-h(I))/DDZ(I-1)+1)/LIMITER2 ! when alpha is large this can be erraneous
           !IF (abs(ALPHA(I))>2D0)
           
           if (f1<0D0) THEN
             ss = 0.5D0*(THE(I)+THE(I-1))
             f1 = max(f1,-(ss - THE(I-1))*DZ(I-1)/dt)
             h_eq = h(I-1)+DDZ(I-1) ! h_eq for this cell
             call vGTHE(h_eq,ss,n(I),alpha(I),tn2(I))
             ss = (ss * (THES(I)-THER(I)) + THER(I))
             f4 = (ss - THE(I))*DZ(I)/dt/LIMITER3
             f1 = max(f1,min(0D0,f4))
           else
             h_eq = h(I) - DDZ(I-1) ! h_eq for upper cell
             call vGTHE(h_eq,ss,n(I-1),alpha(I-1),tn2(I-1))
             ss = (ss * (THES(I)-THER(I)) + THER(I))
             f4 = (THE(I-1) - ss)*DZ(I-1)/dt/LIMITER3
             f1 = min(f1,max(f4,0D0))
           endif
           if (f1>0D0) THEN ! .or. FL>0D0) THEN
           ! FL>0 condition is to prevent cell from giving water to both directions
            !if (F1<0D0) f1 = 10000
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
         SS=(THE(I)-THER(I))/(THES(I)-THER(I))
         call vGh(h(I),ss,n(I),alpha(I),tn1(I))
       ENDDO

       ! Top Cell
       FL = FU; I=2;
       THE(2) = max(THE(2)-FL*dt/DZ(2)+S(I)*dt,  THER(2)+1D-10)

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
             h(1) = h(1)+(THE(2)-THES(2))*DZ(2) !+ S(1)*dt ! IF BBC(1)%Type==2, shouldn't use S(1)
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
           f1 = Khalf(1)*((h(1)-h(2))/DDZ(2)+1)/LIMITER2
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
       SS=(THE(I)-THER(I))/(THES(I)-THER(I))
       call vGh(h(I),ss,n(I),alpha(I),tn1(I))
       
       ! Update h using THE
       DO I = 2,NM
         IF (THE(I)<THES(I)) THEN
           !SS=(THE(I)-THER(I))/(THES(I)-THER(I))
           !m = 1D0 - 1D0/n(I);
           !h(I)=((SS**(-1.D0/m)-1D0)**(1.D0/n(I)))/(-abs(alpha(I)))
         ELSEIF (h(I-1)>=0D0) THEN
           h(I)=h(I-1)+DDZ(I-1)
         ELSE
           h(I)=0D0
         ENDIF
       ENDDO
       IF (BBC(1)%Jtype .eq. 2) RF = BBC(1)%Bdata*dt
       IF (DEBUGG) THEN
       THE1 = sum(THE(2:NA)*DZ(2:NA))
       SSS =  sum(S(2:NA)*DZ(2:NA))*dt
       ENDIF
       !test = THE1 + h(1)- (THE0 + hh2+S(1)*dt+SSS                        &
       !&       - Dperc + RF)
       !test2 = THE1 -(THE0+SSS+DDF-Dperc)
       !IF (max(abs(test),abs(test2))>1d-7) THEN
       !print*,3
       !ENDIF
       if (any(abs(h)>20D0)) THEN
        h(1)=h(1)
       endif
       
       END SUBROUTINE
       
       subroutine vGh_scalar(h,s,n,alpha,tn1)
       ! vG head from relative saturation s
       ! tn1 = n/(n-1) = 1/m;
       real*8 h,s,n,alpha,tn1
       h=((S**(-tn1)-1D0)**(1.D0/n))/(-abs(alpha))
       end subroutine
       
       subroutine vGh_array(h,s,n,alpha,tn1,nn)
       ! vG head from relative saturation s
       integer nn
       real*8 h(nn),s(nn),n(nn),alpha(nn),tn1(nn)
       WHERE(S<0D0) S = 1D-12
       h=((S**(-tn1)-1D0)**(1.D0/n))/(-abs(alpha))
       end subroutine
       
       subroutine vGTHE_scalar(h,s,n,alpha,tn2)
       ! vG head from relative saturation s
       ! tn2 = m = (n-1)/n= 1-1/n;
       real*8 h,s,n,alpha,tn2
       real*8 temp1,temp2
       if (h > 0D0) THEN
         s = 1D0
       else
         temp1 = abs(alpha*h)**n
         temp2 = 1D0+temp1
         s = (temp2**(-tn2))
       endif
       end subroutine
       
       subroutine vGTHE_array(h,s,n,alpha,tn2,nn)
       ! vG head from relative saturation s
       ! tn2 = m = (n-1)/n= 1-1/n;
       integer nn
       real*8 h(nn),s(nn),n(nn),alpha(nn),tn2(nn)
       real*8 temp1(nn),temp2(nn)
       WHERE (h>0D0)
         s = 1D0
       ELSEWHERE
         temp1 = abs(alpha*h)**n
         temp2 = 1D0+temp1
         s = (temp2**(-tn2))
       ENDWHERE
       end subroutine

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

       SUBROUTINE ovn_cv(m,h,U,V,S,d,E,Ex,Ey,Mann,Mask,ddt,lastDt,NNBC,BBC)
       USE Vdata, Only : debugg, BC_type,num2str
       Use Vdata, Only: w
       Use exportASCII, Only: writeDimMat
       USE displayMod
       IMPLICIT NONE
       integer nStencil, NNBC ! nStencil<=5: 5-point scheme, ==9, 9-point scheme
       TYPE(BC_type):: BBC(NNBC)
       integer, save:: status=0, nBC=0
       integer m(2),nG(2)
       logical :: setBoundZero = .true.
       REAL*8 d(2),ddt,lastDt
       REAL*8,DIMENSION(m(1),m(2)),TARGET::h,E,Ex,Ey,Mann,U,V,S
       REAL*8,DIMENSION(m(1),m(2)) :: slopeX,slopeY
       REAL*8,DIMENSION(m(1),m(2),0:2),TARGET::dh
       REAL*8,DIMENSION(m(1),m(2)),TARGET::hx,hy,hSE,hNE,Eta,hnew,h1,hN
       REAL*8,DIMENSION(m(1),m(2)) :: fhuX,fhvY,fhNE,fhSE,UU,VV,VNE,VSE
       LOGICAL*1,DIMENSION(m(1),m(2)),INTENT(in) :: mask
       LOGICAL*1,DIMENSION(m(1),m(2)) :: MM, LIMMAP, Mask2
       REAL*8,POINTER,DIMENSION(:,:):: dhh=>NULL()
       REAL*8 t,rdt,dt,ft,PREC,minDt
       REAL*8 temp,temp2,sgn,dx,dy,t23,t12,limitF,sqdxy
       REAL*8 hm, um, thrshd, diffH, HOvN
       real*8 h0, s_tot, mass_error
       LOGICAL:: DRY, LIMITED, BC,tmpDbg=.false.
       integer,pointer,save:: filterX(:),filterY(:)
       integer i,j,k,n,nL,kF
       REAL*8 :: grav= 9.81D0, cflCond=0.90D0
       real*8 :: umax,hmin,cfl,dtminMxU
       ! hN is input, hNew is output (after RK3)
       nStencil = 9
       t = 0D0
       rdt = ddt
       ft = ddt  ! final time
       if (abs(lastDt)<1D-14 .or. abs(lastDt)>1D5) then
           dt = ddt*86400D0
       else
           dt = lastDt
       endif    
       minDt = ddt*86400D0/8D0
       PREC = 1D-10
       dy = d(1); dx = d(2)
       DRY = .TRUE.
       t23 = 2D0/3D0
       t12 = 0.5D0
       
       !debugg = .true.
       if (debugg) then
         h0 = sum(h)/size(h)
         s_tot = sum(S)/size(h) * dt
       endif

       if (status .eq. 0) then
         ! mask out the faces that should be 0
          allocate(filterX(m(1)*m(2)))
          allocate(filterY(m(1)*m(2)))
          k = 0
          do i=2,m(1)-1
          do j=2,m(2)-1
            if (.Not. Mask(i,j)) then
              BC = .FALSE.
              if (j .eq. 1) then
                if (Mask(i,j+1) .or. Mask(i-1,j+1) .or. Mask(i+1,j+1)) BC = .TRUE.
              elseif (j .eq. m(2)) then
                if (Mask(i,j-1) .or. Mask(i-1,j-1) .or. Mask(i+1,j-1)) BC = .TRUE.
              elseif (i .eq. 1) then
                if (Mask(i+1,j) .or. Mask(i+1,j-1) .or. Mask(i+1,j+1)) BC = .TRUE.
              elseif (i .eq. m(1)) then
                if (Mask(i-1,j) .or. Mask(i-1,j-1) .or. Mask(i-1,j+1)) BC = .TRUE.
              else
                if (Mask(i+1,j) .or. Mask(i-1,j) .or. Mask(i,j-1) .or. Mask(i,j+1) .or.  &
             &   Mask(i-1,j-1) .or. Mask(i-1,j+1) .or. Mask(i+1,j-1) .or. Mask(i+1,j+1))  BC = .TRUE.
              endif
              
              IF (BC) then
                k = k + 1
                filterY(k)=I
                filterX(k)=J
              endif
            endif
          enddo
          enddo
          nBC = k
          status = 1
       endif
       
!       DO i = 2,m(1)-1
!       DO j = 2,m(2)-1 ! negative source term implemented here
!         IF (S(I,J)<0D0) THEN
!           h(I,J) = h(I,J) - S(I,J)*dt
!           S(I,J) = 0D0
!           IF (h(I,J)< 0D0) THEN
!             write(*,*) 'what??'
!             h(I,J)=0D0
!           ENDIF
!         ENDIF
!       ENDDO
!       ENDDO
       kF = 0
       do while (t<ft-PREC)
        kF = kF + 1
        dh = 0D0; umax = 0
        rdt = ft - t; lastDt = dt;
        dt = min(rdt*86400D0,dt); 
        dtminMxU = cflCond * dx/minDt
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
              if (nStencil .eq. 9)  then
                MM(I-1,J-1) = .TRUE.;MM(I+1,J-1) = .TRUE.;
                MM(I-1,J+1) = .TRUE.;MM(I+1,J+1) = .TRUE.;
              endif
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
          sqdxy = sqrt(dx*dy)
          if (nStencil .eq. 9)  then
           hSE = 0D0; hNE=0D0; fhNE=0D0; fhSE=0D0
          endif 
          Eta = hN + E
          call applyCondOVN(3,BBC,NNBC,m,h,hx,hy,Eta,E,U,V,dx,dt,0,n)
          
          !UU = 0D0; VV = 0D0
          CALL interfaceH(m(1),m(2),hN,Eta,Ex,Ey,hx,hy,MM)
          if (nStencil .eq. 9) CALL interfaceHC(m(1),m(2),hN,Eta,Ex,Ey,hNE,hSE,MM)
          call applyCondOVN(0,BBC,NNBC,m,h,hx,hy,Eta,E,U,V,dx,dt,0,n)

          ! set boundary h to 0: this means no outflow BC on the overland will be permitted.
          if (setBoundZero) then
          do k=1,nBC
             i = filterY(k); j=filterX(k)
             hx(i,j)=0; hx(i,j-1)=0; hy(i,j)=0; hy(i-1,j)=0;
             hNE(i,j)=0; hNE(i-1,j-1)=0; hSE(i,j)=0; hSE(i+1,j-1)=0;
          ENDDO
          endif

          ! assume no boundary condition needs to be applied
          ! use diffusive wave
          DO i = 2,m(1)-1
          DO j = 2,m(2)-1
            IF (MM(I,J)) THEN
              diffH = (Eta(i,j)-Eta(i,j+1))
              temp = diffH/dx
              ! should limit slope if deep, flat areas.
              if (temp>0) then
                sgn = 1D0
              else
                sgn = -1D0
              endif
              if (hx(i,j)>0.3D0) then
                HOvN = (hx(i,j)**t23)/mann(i,j)
                temp2 = min(abs(diffH),hx(i,j))/2D0*dx/dt/hOvN
                temp = sgn*max(abs(temp),temp2)
              endif
              
              UU(i,j)=sgn*HOvN*(abs(temp)**t12)
              ! in Diffusive wave, inertia is neglected. Therefore, there can be no super-critical conditions. We restrict velocity to be less than critical velocity
              ! limit outflow only
              hm = hx(i,j); um=sqrt(grav*hm)
              UU(i,j)=sgn*min(um,abs(UU(i,j)),abs(dtminMxU))

              diffH = (Eta(i,j)-Eta(i+1,j))
              temp = (Eta(i,j)-Eta(i+1,j))/dy
              if (temp>0) then
                sgn = 1D0
              else
                sgn = -1D0
              endif
              if (hy(i,j)>0.3D0) then
                HOvN = (hy(i,j)**t23)/mann(i,j)
                temp2 = diffH/2D0*dx/dt/hOvN
                temp = sgn*max(abs(temp),temp2)
              endif
              VV(i,j)=sgn*HOvN*(abs(temp)**t12)
              !VV(i,j)=(sgn/mann(i,j))*(hy(i,j)**t23)*(abs(temp)**t12)
              hm = hy(i,j); um=sqrt(grav*hm)
              VV(i,j)=sgn*min(um,abs(VV(i,j)),abs(dtminMxU))
          
              fhux(I,J)=UU(i,j)*hx(i,j)
              fhvy(I,J)=VV(i,j)*hy(i,j)
              
            umax = max(umax, abs(UU(i,j)), abs(VV(i,j)))

            IF (nStencil .eq. 9) then
              temp = (Eta(i,j)-Eta(i+1,j+1))/(sqrt(2D0)*dy)
              if (temp>0) then
                sgn = 1D0
              else
                sgn = -1D0
              endif
          VNE(i,j)=(sgn/mann(i,j))*(hNE(i,j)**t23)*(abs(temp)**t12)
          !hm = min(sqrt(max(h(i,j)*h(i+1,j+1),0D0)),hNE(i,j)); 
          hm= hNE(i,j); um=sqrt(grav*hm)
          VNE(i,j)=sgn*min(um,abs(VNE(i,j)))

              temp = (Eta(i,j)-Eta(i-1,j+1))/(sqrt(2D0)*dy)
              if (temp>0) then
                sgn = 1D0
              else
                sgn = -1D0
              endif
          VSE(i,j)=(sgn/mann(i,j))*(hSE(i,j)**t23)*(abs(temp)**t12)
          !hm = min(sqrt(max(h(i,j)*h(i-1,j+1),0D0)),hSE(i,j)); 
          hm=hSE(i,j); um=sqrt(grav*hm)
          VSE(i,j)=sgn*min(um,abs(VSE(i,j)))

              fhNE(I,J)=VNE(i,j)*hNE(i,j)
              fhSE(I,J)=VSE(i,j)*hSE(i,j)
 
            ENDIF
              
              !if (isnan(UU(i,j)) .or. isnan(VV(i,j)) ) then
              !call display( 'UU VV nan')
              !    endif
            else 
              UU(i,j) = 0D0
              VV(i,j) = 0D0
              VNE(i,j) = 0D0
              VSE(i,j) = 0D0
            ENDIF
          ENDDO
          ENDDO
          fhux(:,m(2)-1) = 0D0;
          fhvy(m(1)-1,:) = 0D0;
          if (nStencil == 9) then 
            fhNE(:,m(2)-1) = 0D0
            fhNE(m(1)-1,:) = 0D0
            fhSE(:,m(2)-1) = 0D0
            fhSE(2,:) = 0D0
          endif
          LIMITED = .FALSE.
          !
          !if (tmpDbg) then
          !     call  writeDimMat(h,'h',w%t)
          !     call  writeDimMat(UU,'UU',w%t)
          !     call  writeDimMat(VV,'VV',w%t)
          !     call  writeDimMat(Eta,'Eta',w%t)
          !     call  writeDimMat(VNE,'VNE',w%t)
          !     call  writeDimMat(VSE,'VSE',w%t)
          !     call  writeDimMat(fhux,'fhux',w%t)
          !     call  writeDimMat(fhvy,'fhvy',w%t)
          !     call  writeDimMat(h,'fhne',w%t)
          !     call  writeDimMat(h,'fhse',w%t)
          !endif
          
          
          DO i = 2,m(1)-1
          DO j = 2,m(2)-1
            IF (MM(I,J)) THEN

            if (nStencil <= 5) then 
         dhh(I,J)=-(fhuX(I,J)-fhuX(I,J-1))/dx-(fhvY(I,J)-fhvY(I-1,J))/dy
            elseif (nStencil == 9) then 
         dhh(I,J)=-(fhuX(I,J)-fhuX(I,J-1))/dx-(fhvY(I,J)-fhvY(I-1,J))/dy
         dhh(I,J)=dhh(I,J)*t23-((fhNE(I,J)-fhNE(I-1,J-1))/sqdxy+(fhSE(I,J)-fhSE(I+1,J-1))/sqdxy)/6D0
         !dhh(I,J)=-((fhNE(I,J)-fhNE(I-1,J-1))/sqdxy+(fhSE(I,J)-fhSE(I+1,J-1))/sqdxy) ! for debug
            endif
              
              nL = 0
              DO WHILE (h(I,J)>0D0 .and. h(I,J)+dhh(I,J)*dt<-0.1D0 .and. dhh(I,J)<0D0 .and. nL<3)
                ! very occassional over-estimation of outflow flux
                ! newly-added (ugly) flux-limiter
                limitF = -dhh(i,j)*dt/h(i,j) ! this formulation is from LISFLOOD-FP, Bates (2001) doi:10.1016/S0022-1694(00)00278-X
                IF (UU(i,j)>0D0) fhux(i,j)=fhux(i,j)/limitF
                IF (UU(i,j-1)<0D0) fhux(i,j-1)=fhux(i,j-1)/limitF
                IF (VV(i,j)>0D0) fhvy(i,j)=fhvy(i,j)/limitF
                IF (VV(i-1,j)<0D0) fhvy(i-1,j)=fhvy(i-1,j)/limitF
                !fhux(I,J)=UU(i,j)*hx(i,j)
                !fhux(I,J-1)=UU(i,j-1)*hx(i,j-1)
                !fhvy(I,J)=VV(i,j)*hy(i,j)
                !fhvy(I-1,J)=VV(i-1,j)*hy(i,j)
                dhh(I,J)=-(fhuX(I,J)-fhuX(I,J-1))/dx-(fhvY(I,J)-fhvY(I-1,J))/dy

                IF (nStencil == 9) THEN
                    IF (fhNE(i,j)>0D0) fhNE(i,j)=fhNE(i,j)/limitF
                    IF (fhNE(i-1,j-1)<0D0) fhNE(i-1,j-1)=fhNE(i-1,j-1)/limitF
                    IF (fhSE(i,j)>0D0) fhSE(i,j)=fhSE(i,j)/limitF
                    IF (fhSE(i+1,j-1)<0D0) fhSE(i+1,j-1)=fhSE(i+1,j-1)/limitF
                    dhh(I,J)=dhh(I,J)*t23-((fhNE(I,J)-fhNE(I-1,J-1))/sqdxy+(fhSE(I,J)-fhSE(I+1,J-1))/sqdxy)/6D0
                ENDIF
                nL = nL + 1
                LIMITED = .TRUE.
                LIMMAP(I-1:I+1,J-1:J+1)=.TRUE.
              ENDDO
              dhh(I,J)=dhh(I,J)+S(I,J)
              
            ENDIF
          ENDDO
          ENDDO
          
          IF (LIMITED) THEN
          DO i = 2,m(1)-1
          DO j = 2,m(2)-1
            IF (LIMMAP(I,J) .AND. MM(I,J)) THEN
              dhh(I,J)=-(fhuX(I,J)-fhuX(I,J-1))/dx-(fhvY(I,J)-fhvY(I-1,J))/dy
              dhh(I,J)=dhh(I,J)*t23-((fhNE(I,J)-fhNE(I-1,J-1))/sqdxy+(fhSE(I,J)-fhSE(I+1,J-1))/sqdxy)/6D0
              dhh(I,J)=dhh(I,J)+S(I,J) ! positive implemented here
              LIMMAP(I,J) = .FALSE.
            ENDIF
          ENDDO
          ENDDO
          ENDIF
          ! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        enddo
        ! #################################
        hNew = h +dt/6D0*(dh(:,:,0)+4D0*dh(:,:,2)+dh(:,:,1))
        U = UU
        V = VV
        thrshd = 0.3D0
        hmin = minval(hnew)
        cfl = umax*dt/dx
        if ((hmin<-thrshd .or. cfl>1.5D0) .and. dt>mindt) THEN
         dt = min(rdt*86400D0,max(dt/2d0,mindt))
          !DO i = 2,m(1)-1
          !DO j = 2,m(2)-1
          !  IF (hNew(i,j)<-thrshd) THEN
          !    call display('hNew Neg Found',num2str(i),num2str(j),num2str(t))
          !  ENDIF
          !ENDDO
          !ENDDO
        else
         t = t + dt/86400D0
         !h1 = h ! print debug
         h = hnew
         if (cfl<cflCond) dt = min(cflCond*dx/umax,ddt*86400D0)
        endif
       enddo
       call applyCondOVN(3,BBC,NNBC,m,h,hx,hy,Eta,E,U,V,dx,dt,1,n)
       
       t = ft ! avoid runoff error
       NULLIFY(dhh)
       
       if (debugg) then
         mass_error = (sum(h)/size(h)-h0) - s_tot
         if (abs(mass_error)>1D-12) then
           call display("ovn_cv: MBE")
         endif
       endif
       
       END SUBROUTINE ovn_cv

       subroutine applyCondOVN(TYP,BBC,NBC,m,h,hx,hy,Eta,E,U,V,dx,dt,dir,rkStage)
       ! dir 
      USE VDATA
      integer TYP,m(2),NBC,nC,dir,rkStage 
      TYPE(BC_type) :: BBC(NBC)
      REAL*8,DIMENSION(m(1),m(2)),target::h,hx,hy,Eta,E,U,V
      REAL*8,DIMENSION(m(1),m(2))::hN,bBnd,denom ! local variables
      integer :: I,j,NI,nk,nV=11,state=0
      integer,DIMENSION(:),POINTER :: idxy=>NULL(),Cidx=>NULL()
      integer,DIMENSION(:),POINTER :: Nidx=>NULL(),Sidx=>NULL()
      integer,DIMENSION(:),POINTER :: Widx=>NULL(),Eidx=>NULL()
      REAL*8, POINTER:: P=>NULL(), p2D(:,:) =>null()
      REAL*8, DIMENSION(:),POINTER:: E1d, h1d, el1d
      REAL*8, DIMENSION(:),POINTER:: EBM2N1,EBM1N2,EBM1N1,EBMN2,EBMN1,EBM1N,EBMN,EBMNP1S,EBM2NS,EBM1NS,EBMNS
      REAL*8 :: dx,dt
      DO I=1,NBC
       IF (BBC(I)%JType .eq. TYP) THEN
         SELECT CASE (BBC(I)%JType)
          CASE (0) ! Wall
           idxy=>BBC(I)%idx
           h1d => oneD_ptr(m(1)*m(2), hx)
           h1d(BBC(I)%idx)=0; h1d(BBC(I)%BIndexIn1)=0;
           h1d => oneD_ptr(m(1)*m(2), hy)
           h1d(BBC(I)%idx)=0; h1d(BBC(I)%BIndexIn2)=0;
          CASE (1) ! Dirichilet
           idxy=>BBC(I)%idx ! In the original matrix
          CASE (2) ! 
          CASE (3,33) ! Orlandski outflow condition: one key assumption: no timestep resizing is occuring (dt = min(rdt*86400D0,max(dt/2d0,mindt)) NOT executed)
              nk = size(BBC(I)%idx)
              idxy=>BBC(I)%idx
              e1d => oneD_ptr(m(1)*m(2), Eta)
              h1d => oneD_ptr(m(1)*m(2), h)
              el1d => oneD_ptr(m(1)*m(2), E)
              if (.not. associated(BBC(I)%BV)) then
                  allocate(BBC(I)%BV(nk*nV)) ! temporary values
                  p2d => twoD_Ptr(nV,nK,BBC(I)%BV)
                  ! each rows are, respectively:
                  ! EBM2N1,EBM1N2,EBM1N1,EBMN2,EBMN1
                  state = 0
              endif
              p2d => twoD_Ptr(nk,nv,BBC(I)%BV)
              J = 1
              EBM2N1 => p2d(j,:); j=j+1
              EBM1N2 => p2d(j,:); j=j+1
              EBM1N1 => p2d(j,:); j=j+1
              EBMN2 => p2d(j,:); j=j+1
              EBMN1 => p2d(j,:); j=j+1
              EBMN  => p2d(j,:); j=j+1
              EBM1N => p2d(j,:); j=j+1
              
              EBM2NS => p2d(j,:); j=j+1
              EBM1NS => p2d(j,:); j=j+1
              EBMNS  => p2d(j,:); j=j+1
              EBMNP1S  => p2d(j,:); j=j+1
              if (state .eq. 0) then
                  EBM2N1  = e1d(int(BBC(I)%BIndexIn2))
                  EBM1N1  = e1d(int(BBC(I)%BIndexIn1))
                  EBM1N2  = e1d(int(BBC(I)%BIndexIn1))
                  EBMN1   = e1d(idxy) ! output
                  EBMN2   = e1d(idxy) ! output
                  state = 1
              endif    
              if (dir .eq. 0) then 
                  ! compute boundary value
                  EBMN   = e1d(idxy) ! output
                  EBM1N  = e1d(int(BBC(I)%BIndexIn1))
                  if (rKStage .eq. 0) then
                    call orlanskiOutflowPredict(nk,dx,dt,EBM2N1,EBM1N2,EBM1N1,EBMN2,EBM1N,EBMN)
                    EBM2NS = e1d(int(BBC(I)%BIndexIn2))
                    EBM1NS = EBM1N
                    EBMNS = EBMN
                  elseif (rKStage .eq. 1) then ! Eta coming in is Eta^{N+1}
                    ! upgrade on time level and use Saved values
                    call orlanskiOutflowPredict(nk,dx,dt,EBM2NS,EBM1N1,EBM1NS,EBMN1,EBM1N,EBMN)
                    EBMNP1S = EBMN
                  elseif (rKStage .eq. 2)  then
                    EBMN = (EBMNS+EBMNP1S)/2D0
                  endif    
                  h1d(idxy) = EBMN - El1d(idxy)
                  ! final implementation into BC point values
                  do j=1,size(idxy)
                      h1d(idxy(j)) = min(h1d(idxy(j)),0D0)
                  enddo
              else
                  ! save bv for next time step!
                  EBM1N2  = EBM1N1
                  EBMN2   = EBMN1
                  EBM2N1  = EBM2NS
                  EBM1N1  = EBM1NS
                  EBMN1   = EBMNS ! output
              endif    
              
              ! [1.EBm1n EBM2n eBm1n1 ebn
              !Demon=EtaBmin+EtaBmin2+EtaBM2N1
         END SELECT
       ENDIF
      ENDDO
       
       end subroutine applyCondOVN
       
       subroutine orlanskiOutflowPredict(n,dx,dt,BM2N1,BM1N2,BM1N1,BMN2,BM1N,BMN)
       ! need complete savings of M1 M2 and N1 N2
       integer n
       real*8,dimension(n),intent(in) :: BM1N,BM2N1,BM1N1,BM1N2,BMN2
       real*8,dimension(n),intent(out) :: BMN
       real*8,dimension(n) :: C, denom, temp
       real*8 dx, dt, dtOdx
       
       denom = BM1N + BM1N2 - 2D0*BM2N1
       C = (BM1N-BM1N2)/denom * dx/dt
       where (C<0D0) C=0D0
       temp = (1D0+(dt/dx)*C)
       BMN = (temp*BMN2+2D0*dt/dx*(C*BM1N1))/temp
       
       end subroutine
       
       SUBROUTINE h1negative(nz,na,h,THE,THER,THES,DZ,DF,                 &
     &             ALPHA,LAMBDA,n)
       ! this subroutine tries to fill up upper cell hh3<0
       integer nz,na
       REAL*8,DIMENSION(nz)::h,THE,THES,DZ,THER,ALPHA,                      &
     &             LAMBDA,n
       integer I,J,K
       REAL*8 TH,PREC,mm,SS,DF,C,C2

       I = 2
       C = 0.15D0 ! how much you can extract
       C2 = 0.07D0
       DO WHILE (I <= na)
         TH = min((THE(I)-THER(I))*DZ(I)*C,-h(1),(THES(I)-THER(I))*C2)
         h(1) = h(1)+TH
         DF = DF - TH ! too much infiltration, reduce correspondingly
         THE(I) = THE(I)-TH/DZ(I)
         mm = 1D0-1/(n(I))
         SS = (THE(I)-THER(I))/(THES(I)-THER(I))
         h(I) = ((SS**(-1D0/mm)-1)**(1D0/n(I)))/(-abs(ALPHA(I)))
         IF (h(1)>=0D0) EXIT
         if (h(I)<-100) then
         !print *, 'h1negative?'
         endif
         I = I+1
       ENDDO


       END SUBROUTINE h1negative

       !!!!! GENERALIZED GREEN AND AMPT METHOD with source term
       SUBROUTINE GGA(NZ,WFP,DF,THE,THES,THER,KS,h,DZ,P,                         &
     &            S,E,KBar,FM,dt,alpha,N,LAMBDA,FL,Code,THEW)
       ! Solve 1D GGA
       USE Vdata, Only: rtsafe, rtnewt, DEBUGG 
       use displayMod, Only: display
       USE locfunc, Only: sumDZ
       IMPLICIT NONE
       integer :: NZ
       REAL*8, DIMENSION(NZ) :: THE,THES,THER,KS,DZ,E,KBar,FM,h
       REAL*8, DIMENSION(NZ) :: ALPHA,N,LAMBDA,S,THESAVE
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
       REAL*8 :: L, IFC, Ft, Fex, IFSlow, temp1,temp2,err,S0(nz)
       ! Wetting Front Length, Infiltration Capacity, Cumulative Infiltration
       integer :: IWFP, ICASE, K, IU,I,J,SCODE
       LOGICAL :: T_F
       REAL*8 :: PREC,PREC2,FMM,t1,t2,t3,ss,x1,THE0,M1,ME,M2,M3,ME2,hh0
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
       THE0 = sum(THE*DZ); hh0 = h(1); S0 = S
       M2 = DF
       if (debugg) THESAVE = THE

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
       PREC = 1.0D-10
       PREC2 = 1.0D-8
       K = 0
       !DF = 0.0D0
       DPerc = 0.0D0
       T_F = .TRUE.

       ! Note that first layer is the ponding layer! Soil starts from the second one!
       ! due to the addition of CLM, the surface cells may have their mass changed outside of UNSAT
       ! therefore, these cells may already be saturated, need to move to next cells in this case
       IWFP = min(max(FLOOR(WFP+PREC)+1,2),nz-1) ! WFP = 2.2 means it is in the third layer, 0.2 cell from the top
       IF (Code .eq. 0 .OR. THE(IWFP)>=THES(IWFP)-PREC) THEN ! last time step was not GGA, so WFP cannot be used
       !IF (Code .eq. 0) THEN ! old code before CLM
         DO I = 2,NZ
           IF (THE(I)>=THES(I)-PREC2) THEN
           ELSE
             EXIT
           ENDIF
         ENDDO
         IF (I<nz .OR. (I .eq. nz .and. THE(I)<=THES(I)-PREC2)) THEN
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
         IF (SCODE .eq. 1) THEN
           goto 1234
         ENDIF
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
             do while (abs(dTHE) < PREC .and. IWFP <=NZ)
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
             IF (K > 1) THEN
             CALL GGA_Source(IWFP,IWFP,DZ,S,h,WFP,prcp,dt,tr,DF,SCODE)
             ENDIF
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
             ! This should not happen, for some reason, a cell is filled up but IWFP was not moved to the next cell
             call display('GGA, for some reason, a cell is filled up but IWFP was not moved to the next cell')
             tr = 0D0
             tf = 0D0
             tp = 0D0
             GOTO 2000
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
           if (tp <= 0D0) tp = 1D20
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
              err = 1D-6

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

 2000    t = DMIN1(tr,tp,tf)  ! Whichever time comes sooner
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
         IF (t1<-PREC) then
         ! maybe the 'ponding time' tp iteration did not converge. 
         ! as a result, too much infiltration
           Ft = h(1) + t*prcp
           T_F = .FALSE.
         endif
         t1 = max(t1,0D0)

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
!       IF (DEBUGG) THEN
!         print*, 'here3',k,dTHE,THES(IWFP),THE(IWFP),IWFP,nz
!         print*, 'hh',h
!       ENDIF
       ENDDO

       ! Update state variables (h, THE, WPF)
       IF (IWFP>2) THEN
         h(IWFP-1)=LEN + DZ(IWFP-1)/2D0
         DO I=IWFP-2,2,-1
            h(I) = h(I+1) + E(I) - E(I+1)
         ENDDO
         THE(2:IWFP-1)=THES(2:IWFP-1)
       ENDIF


 1234  IF (IWFP>1) THEN
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
       IF (DEBUGG) THEN
       M3 = sum(THE*DZ)
       M1 = sum(DZ(2:IWFP)*(S0(2:IWFP)-S(2:IWFP)) )*dt
       ME = M3 - (THE0 + DF - M2 + M1)
       ME2= M3+h(1)-(THE0+hh0+M1+P*dt)
       if (abs(ME)>PREC2 .or. abs(ME2)>PREC2) THEN
       Print *,'GGA m_err'
       endif
       ENDIF
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

       SUBROUTINE surfaceDer(h,u,dudh,h0,h1,h2,u1,u2,dx,w)
       ! for two-sided outflow case
       ! compute inner ^1/2 term and derivative
       REAL*8 u,dudh,h,h0,h1,h2,u1,u2,dx,w
       u = max(h0,0D0)
       
       END SUBROUTINE
       
       SUBROUTINE GGA_Source(N1,NM,DZ,S,h,WFP,prcp,dt,tr,DF,SCODE)
       IMPLICIT NONE
       integer N1,NM,SCODE
       REAL*8,DIMENSION(NM)::S,h,DZ
       REAL*8 :: WFP,prcp,dt,DF,tr
       REAL*8 t1,t2,t3,ss,s2
       ! dt is the full time step
       ! tr is remaining
       SCODE = 0
       !ss = SUM(DZ(2:IWFP)*S(2:IWFP)) * tr
       ss = SUM(DZ(N1:NM)*S(N1:NM)) * dt
       t1 = h(1) + ss
!       if (h(1)<0D0) then
!         s2   = min(prcp*tr,-h(1))
!         prcp = prcp - s2/tr
!         h(1) = h(1)+s2
!       endif
       
       if (t1 < 0) THEN
         t3 = ss
         ss = ss + max(0D0,h(1))
         DF = DF + max(0D0,h(1))
         h(1) = min(h(1),0D0) ! pond gone
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
      
      
!=====================================================================================


      SUBROUTINE SURFACE(m,n,nz,nt,ns,Veg,V,H,dt,OVH,OVF)
      ! This function takes care of ET, CANOPY, DEPRESSION STORAGE
      ! Rule of Thumb: Should not use pointer for output unless it is very clear given the context
      ! C H A O P E N G    S H E N
      ! Env. Engr.
      ! Code Unit: m

      USE Vdata
      IMPLICIT NONE
      integer m,n,nt,ns,nz,H ! H is hour of day
      integer,DIMENSION(:,:),POINTER:: nRPT=>NULL()
      integer,dimension(:,:),pointer:: SMAP=>NULL()
      integer,DIMENSION(:,:,:),POINTER::RPT=>NULL()
      REAL*8 dt
      REAL*8,DIMENSION(:,:),POINTER::ET=>NULL(),Cpc=>NULL(),hh=>NULL()
      REAL*8,DIMENSION(:,:,:),POINTER::FRAC=>NULL(),LAI=>NULL(),HI=>NULL(),ROOT=>NULL(),KC=>NULL()
      REAL*8,DIMENSION(:,:),POINTER::MAPP=>NULL(),Imp=>NULL()
      REAL*8,DIMENSION(:,:),POINTER::prcp=>NULL(),rad=>NULL(),wnd=>NULL(),hmd=>NULL(),ret=>NULL()
      REAL*8,DIMENSION(:,:),POINTER::rets=>NULL(),tempe=>NULL()
      REAL*8,DIMENSION(:,:,:),POINTER::EB=>NULL()
      REAL*8,DIMENSION(:),POINTER::EE=>NULL(),THE=>NULL(),THES=>NULL(),THER=>NULL()
      REAL*8,DIMENSION(:),POINTER::VS=>NULL(),FC=>NULL(),WP=>NULL(),hhh=>NULL()
      REAL*8,DIMENSION(m,n)::OVH,OVF
      TYPE(Vege_type)::Veg
      !TYPE(STA)::DAY(2)
      TYPE(VDZ_type)::V
      TYPE(THERMAL_type),POINTER:: th=>NULL()
      integer I,J,K,L,P,RDN,EDN,ETN,II,RMX
      REAL*8 RPET,PREC,NewS,RPrcp,FPET,EvapF
      REAL*8 PT(nt),RDF(nz),RFrac(nt)
      REAL*8,DIMENSION(nz):: PET,TP,EVAP,EDF,F
      REAL*8 THRSD, SS, RATIO, TEMP, temp2,temp3,temp4,temp5,temp6
      REAL*8,DIMENSION(:),POINTER:: TEMPP=>NULL(), tempp2=>NULL(),tempp3=>NULL(),tempp4=>NULL()
      integer NR,NCODE(3)
      REAL*8::EU=0D0,EL=0D0,RSnow=0D0,RS=0D0,RadS,Sa,TT,albedo,SNOMP
      REAL*8 ETS,ST,ETA,test,CPET,VA,GAM,OM,test2,OCOV,Newsnow,tt2

      REAL*8,DIMENSION(7) :: inpt
      REAL*8,DIMENSION(:),POINTER:: sitev=>NULL(),statev=>NULL(),outv=>NULL()
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
#ifndef MEX        
         CALL SnowLSub(DT*24,1,input,th%Ssitev(I,J,:),th%Sstatev(I,J,:)         &
     &  ,th%SParm,iflag,nstepday, cump,cume,cummr                                &
     &  ,outv,th%tss(I,J,:),th%tsa(I,J,:),th%nd,th%dfc)
#endif
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
         temp =  Veg%BSoil(I,J,1)*RPET  ! still not sure whether we should do the re-scaling ! changed to BSoil(I,J,1) on May 1st 2013
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
     &   ,nGho,Code,CType,rtr,UPST,UPSV,DST,DSVIn,STATE,umax,hg,mr,KB,fG,               &
     &   tt,vrv)
      ! RK3 river flow
      ! Chaopeng Shen
      ! PhD Environmental Engineering
      ! Michigan State Univ.
      IMPLICIT NONE
      integer nx, nGho,EQQ !EQQ: 1--KW, 2--DW, 3--DyW
       ! ========= Input ========================> 
      REAL*8,DIMENSION(nx), intent(In)   :: E    ! River bed elevation at mid-point of reach (m)
      REAL*8,DIMENSION(nx), intent(In)   :: Mann ! Manning's coefficient
      REAL*8,DIMENSION(nx), intent(In)   :: S    ! Source term (m2/s), mainly triburary input because ovn input applied directly, qgc solved inside here, others are 0
      REAL*8,DIMENSION(nx), intent(In)   :: W    ! River width (assuming rectangular channel) (m)
      REAL*8,DIMENSION(nx), intent(In)   :: SM   ! Not really used
      REAL*8,DIMENSION(nx), intent(In)   :: SE   ! Slope of eta-- river stage 
      REAL*8,DIMENSION(nx), intent(In)   :: hg   ! Groundwater head --> not changing. In new implicit scheme, this is no longer used
      REAL*8,DIMENSION(nx), intent(In)   :: KB   ! river bed conductivity (m/day)
      REAL*8,DIMENSION(*), intent(In)   :: DSVIn ! downstream boundary condition values [h,E] of downstream river junction cell. Now vrv also contains a part of the info
       ! <======== InOut ========================>
      REAL*8,DIMENSION(nx), intent(InOut):: h    ! river depth, (r(rid)%Riv%h) (m)
      REAL*8,DIMENSION(nx), intent(InOut):: A    ! cross-sectional area. The true state variable in river flow module  (m2)
      REAL*8,DIMENSION(nx), intent(InOut):: u    ! velocity at downstream end of each cell (m/s)
      REAL*8,DIMENSION(nx), intent(InOut):: fG   ! baseflow flux (positive to channel) (m3) -- Accumulated
      REAL*8,DIMENSION(nx), intent(InOut):: wF   ! Width at 
      REAL*8,DIMENSION(nx), intent(InOut):: mr   ! connection thickess between river bed and groundwater (maximum of bed sediment thickness or E - Hg) (m)
      REAL*8,DIMENSION(30,nx), intent(InOut),OPTIONAL,target::vrv ! virtual river variables [1. hg_accu (m); 2. Qgc_n (m/day); 3. gamma (-); 4. beta (m/s) ; 5 .ST; 6. zeta=beta*dt/mr/(1+K*dt/mr)] -- see notes
       ! <======== Out      ======================
      REAL*8,DIMENSION(nx), intent(Out)  :: Qx   ! Discharge at downstream end of each cell (m3/s)
      REAL*8,DIMENSION(nx), intent(Out)  :: R    ! Hydraulic Radius (m)
       ! ========= NOT USED ======================
      REAL*8,DIMENSION(nx), intent(In)   :: S0   ! Average slope (m/m) of reach. Not really used currently (because not doing Kinematic wave)
      REAL*8 dx,ddt, rtr, code
      integer CType(nx) ! Boundary Condition Types
      REAL*8 UPSV,UPST,DST
      INTEGER STATE, ISTATE,ICODE
       ! ========= Local =========================
      !REAL*8,DIMENSION(nx), OPTIONAL:: Fact, phi
      INTEGER hgImplicit
      
      integer T(nx),I,K,KK,NN, MAXNN !SHIFT
      REAL*8 G, tt, tt1, tt2, sff, Neg, RDT, dt,EL,ER, PREC
      REAL*8,DIMENSION(nx):: dAdt0,dAdt1,dAdt2,fG0,fG1,fG2,hN
      REAL*8,DIMENSION(nx):: QQ,Q0,Q1,Q2,U0,U1,U2,hx,ht,rt,hgg,fgg,hNtemp,QQgctemp
      REAL*8 UMAX,COU,sgn(nx),temp(nx),temp2(nx),hgNtemp(nx),alpha(nx)
      REAL*8 TQ(2),SUMQ,SUMQ2,T1,T2,TA(nx),QQgc(nx)
      REAL*8 AA_old,AAA,QQgc_old,ww,test,test2,Qout,DSV(10)
      REAL*8,POINTER:: hgU(:), qgc_n(:), gam(:), beta(:), st(:), zeta(:)
      
      G = 9.81D0      
      PREC = 1.0D-11
      sgn = 1.0D0
      Neg = 1e-6
      QQ = 0.0D0
      Qx = 0.0D0
      RDT = DDT
      DT = DDT
      DSV(1:2)=DSVIn(1:2) ![h, Qx, predictedH]
      DSV(3)=0D0; DSV(4)=w(nx); DSV(5)=dt/w(nx)/wF(nx)
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
        DSV(3)=DSV(1)
        CALL one_step(A,dAdt0,Q0,u,hx,h,R,fG0,dx,EQQ) ! TESTING WHAT IS WRONG
        STATE = 0
        
        DSV(3)=DSV(1)+DSV(5)*(Q0(nx-1)-DSV(2))
        CALL one_step(A+dt*dAdt0,dAdt1,Q1,u1,hx,h,R,fG1,dx,EQQ)
        
        DSV(3)=DSV(1)+DSV(5)/2D0*((Q0(nx-1)+Q1(nx-1))/2D0-DSV(2)) ! wF(nx) is downstream river dx
        CALL one_step(A+(dt/4D0)*(dAdt0+dAdt1),dAdt2,Q2,u2,hx,h,R,                     &
     &               fG2,dx,EQQ)
        
        umax = max(maxval(abs(u(2:nx))),maxval(abs(u1)))*(3D0/2D0)
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

      hgImplicit = 0
      IF (PRESENT(vrv)) hgImplicit = 1

      IF (hgImplicit .eq. 0) THEN
        hN = (TEMP*hg+h+E)/(1D0+TEMP) - E
      ELSE
!        hNtemp = (TEMP*hg+h+E)/(1D0+TEMP) - E
!        hN = (h+TEMP*(fact*h+Hg+phi-E) )/(1D0+TEMP+TEMP*fact)
!        
!        ! for comparison
!        WHERE (hNtemp<0.0D0) hNtemp=0.0D0
!        QQgctemp = (hNtemp-h)*w*dx
        hgU   => vrv(1,:)
        qgc_n => vrv(2,:)
        gam   => vrv(3,:)
        beta  => vrv(4,:)
        st    => vrv(5,:)
        zeta  => vrv(6,:)
        alpha = ddt/mr  ! s/m
        zeta  = beta*alpha/(1D0+KB/86400D0*alpha)
        hgNTemp = hgU
        hgU = (ST*hgU+qgc_n*(ddt/86400D0)+zeta*(h+E))/(ST+zeta)
        hN = (hgNTemp - hgU)*ST/gam + h

        !temp = hgU - hg
        !temp2= (1D0+KB/86400D0*alpha)
        !ta   = hN - h
      ENDIF
      !write(444,*) 'Riv internal QQgctemp', QQgctemp(37:38)/w/dx
      !TEMP2 = (hN - h)/(ddt/86400)

      WHERE (hN<0.0D0) hN=0.0D0
      hN(1)=hN(2); hN(nx)=hN(nx-1)

      QQgc = (hN-h)*w*dx
      fG = fG + QQgc ! Unit :m3
      A = hN * w      
      ENDIF
      
      DO I=2,nx
          if (hg(I) - E(I)>0) vrv(11,I) = vrv(11,I)+ddt
      ENDDO
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
      
      
      DO I=2,nx
          if (hg(I) - E(I)>0) vrv(11,I) = vrv(11,I)+ddt
      ENDDO
      
      CONTAINS
! BEGIN CONTAINS ##########################################################
      SUBROUTINE one_step(AA,dAdt,Q,UU,hhx,hh,RR,fGG,dx,EQ)
      ! Input: AA,(UU,wF,hhx). Output: dAdt,Q,UU,hhx,hh,RR
      ! Does not modify AA
      IMPLICIT NONE
      REAL*8,DIMENSION(nx):: AA,dAdt,Q,hN,hgg,TEMPC
      INTEGER OPT,EQ, nnx, I
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
       
      ! upwinding UU
      nnx = nx
      DO I=1,nnx-1
        !IF (UU(I)<-1D-15) THEN
        !  wF(I) = min(w(I)*2D0,w(I+1))
        !ELSE
        !  wF(I) = min(w(I),2D0*w(I+1))
        !ENDIF
        wF(I) = min(w(I), w(I+1))
      ENDDO
      Q  = wF * hhx * UU ! This Q is Q at this moment, this is different from the output Qx, which is a mass-balance Q (Q in a period of time)
      ! Now compute interchange with groundwater, use an implicit update (operator splitting)
      ! for rectangular channel:
      
      IF (int(UPST) .eq. 2) Q(1)=UPSV  ! UPST type 2 is not working inside BC. it will be overwritten outside
      
      dAdt =  (CSHIFT(Q,-1)-Q)/dx + S
      
      dAdt(1) =0.0D0
      dAdt(nx)=0.0D0
      
      DO I=2,nnx-1
         IF (dAdt(I)*dt + A(I) < 0) THEN ! too much outflow. 
            IF (Q(I-1)<0D0 .AND. Q(I)>0D0) THEN
                ! flow out on both sides and apparently hx will not be right
                !print*, 'double-sized outflow'
            ELSE
                !print*, 'not a double-sided outflow problem'
            ENDIF
         ENDIF
      ENDDO
      
      
      ! As long as channel geometry is longitidinally homogeneous, this relation holds
      END SUBROUTINE one_step
! END CONTAINS ##############################
      end subroutine riv_cv
      
      
      
      

      !CALL BC(nx,hh,AA,u,w,UPST,UPSV,DST,DSV)
      SUBROUTINE BC(nx,h,A,u,w,UPST,UPSV,DST,DSV,dx)
      integer nx
      REAL*8,DIMENSION(nx):: h,A,u,w
      REAL*8 UPSV, DSV(*), UPST, DST
      INTEGER UPS, DS
      REAL*8 sgn, EL,ER,AAA,dx
      ! BC
      
      UPS = int(UPST)
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
       h(nx)=DSV(3)
       w(nx)=min(w(nx-1),DSV(4))
       
       A(nx)=h(nx)*w(nx)
       
       u(nx)=u(nx-1)
      END SELECT
      
      END SUBROUTINE BC

      !=================================================================================
       Subroutine interfaceHC(ny,nx,h,E,Ex,Ey,hx,hy,M)
       ! Corner version of this code.
       ! hx is h_{i,j}^{NE}, pointing to (i+1,j+1)
       ! hy is h_{i,j}^{SE}, pointing to (i-1,j+1)
       implicit none
       integer ny, nx
       real*8, DIMENSION(ny,nx):: h, E, Ex,Ey,hx,hy,t1,t2
       logical*1:: M(ny,nx)
      !!! THIS ONE DOES NOT LOOK AT Ex or Ey, DISREGARD EDGE EFFECT

       integer i,j
       real*8 ER(ny,nx)
       REAL*8 PREC,tt1,tt2,phi

       hx = 0.0D0
       hy = 0.0D0
       ER = 0.0D0
       PREC = 1.0D-10
      ! hx!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       IF (nx>1) THEN
       DO I = 1,ny
         DO J = 1,nx-1
          IF (M(I,J)) THEN
           if (E(I,J)>=E(I+1,J+1)) THEN ! Positive Wind
             if (h(I,J)<PREC) THEN
               hx(I,J)=0.0D0
             elseif (h(I+1,J+1)>PREC) THEN
               if (h(I,J)<h(I+1,J+1)) THEN ! changed to 10/10/2009
                 hx(I,J)=h(I,J)
               ELSE
                hx(I,J)=0.5D0*(h(I,J)+h(I+1,J+1)) ! modification 05/06/2009
               endif
             else ! half-wet
               hx(I,J)=h(I,J)
             endif

           else
             if (h(I+1,J+1)<PREC) THEN
               hx(I,J) = 0.0D0
             elseif (h(I,J)>PREC) THEN
              !hx(I,J) = 0.5D0*(h(I,J)+h(I+1,J+1))

               if (h(I+1,J+1)<h(I,J)) THEN
                 hx(I,J)=h(I+1,J+1)
               ELSE
                 hx(I,J)=0.5D0*(h(I,J)+h(I+1,J+1)) ! modification 05/06/2009
               endif
             else ! half-wet
               hx(I,J) = h(I+1,J+1)
             endif

           endif
         ENDIF ! mask if
         ENDDO
       ENDDO
       hx = max(hx,0.0D0)
       ENDIF

       ! hy!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       IF (ny>1) THEN
       DO I = 1,ny-1
         DO J = 1,nx
        !IF ((I .eq. 29) .AND. (J .eq. 45)) THEN
        !print *,3
        !ENDIF

           if (E(I,J)>=E(I-1,J+1)) THEN
             if (h(I,J)<PREC) THEN
               hy(I,J)=0.0D0
             elseif (h(I-1,J+1)>PREC) THEN
              !hy(I,J)=0.5D0*(h(I,J)+h(I-1,J+1))
              !tt1=h(I,J); tt2=h(I-1,J+1)
               if (h(I,J)<h(I-1,J+1)) THEN
                 hy(I,J)=h(I,J)
               ELSE
                 hy(I,J)=0.5D0*(h(I,J)+h(I-1,J+1)) ! modification 05/06/2009
               endif
             else ! half-wet
               hy(I,J)=h(I,J)
             endif
           else
             if (h(I-1,J+1)<PREC) THEN
               hy(I,J) = 0.0D0
             elseif (h(I,J)>PREC) THEN
              !hy(I,J) = 0.5D0*(h(I,J)+h(I-1,J+1))

               if (h(I-1,J+1)<h(I,J)) THEN
                 hy(I,J)=h(I-1,J+1)
               ELSE
                 hy(I,J)=0.5D0*(h(I,J)+h(I-1,J+1)) ! modification 05/06/2009
               endif
             else ! half-wet
               hy(I,J) = h(I-1,J+1)
             endif
           endif
         ENDDO
       ENDDO
       hy = max(hy,0.0D0)
       ENDIF

       END SUBROUTINE interfaceHC
      
!=================================================================================
       Subroutine interfaceH(ny,nx,h,E,Ex,Ey,hx,hy,M)
       implicit none
       integer ny, nx
       real*8, DIMENSION(ny,nx):: h, E, Ex,Ey,hx,hy,t1,t2
       logical*1, optional:: M(ny,nx)
       logical IM, IDM
      !!! THIS ONE DOES NOT LOOK AT Ex or Ey, DISREGARD EDGE EFFECT

       integer i,j
       real*8 ER(ny,nx)
       REAL*8 PREC,tt1,tt2,phi

       hx = 0.0D0
       hy = 0.0D0
       ER = 0.0D0
       PREC = 1.0D-10
      ! hx!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       if (PRESENT(M)) then
          IM = .true.
       else
          IM = .false.
       endif

       IF (nx>1) THEN
       DO I = 1,ny
         DO J = 1,nx-1
          IF (IM) then 
            IDM = M(I,J)
          ELSEIF (J>1 .AND. J<nx) THEN
          ! Called by river flow subroutine.
            IDM = .TRUE.
          ELSE
            IDM = .FALSE.
          ENDIF
          IF (IDM) THEN
           if (E(I,J)>=E(I,J+1)) THEN ! Positive Wind
             if (h(I,J)<PREC) THEN
               hx(I,J)=0.0D0
             elseif (h(I,J+1)>PREC) THEN
               if (h(I,J)<h(I,J+1)) THEN ! changed to 10/10/2009
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

               if (h(I,J+1)<h(I,J)) THEN
                 hx(I,J)=h(I,J+1)
               ELSE
                 hx(I,J)=0.5D0*(h(I,J)+h(I,J+1)) ! modification 05/06/2009
               endif
             else ! half-wet
               hx(I,J) = h(I,J+1)
             endif

           endif
         ENDIF ! mask if
         ENDDO
       ENDDO
       hx = max(hx,0.0D0)
       ENDIF

       ! hy!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       IF (ny>1) THEN
       DO I = 1,ny-1
         DO J = 1,nx
        !IF ((I .eq. 29) .AND. (J .eq. 45)) THEN
        !print *,3
        !ENDIF
        IF (IM) then 
            IDM = M(I,J)
          ELSEIF (I>1 .AND. I<ny) THEN
          ! Called by river flow subroutine.
            IDM = .TRUE.
          ELSE
            IDM = .FALSE.
          ENDIF
          IF (IDM) THEN
           if (E(I,J)>=E(I+1,J)) THEN
             if (h(I,J)<PREC) THEN
               hy(I,J)=0.0D0
             elseif (h(I+1,J)>PREC) THEN
              !hy(I,J)=0.5D0*(h(I,J)+h(I+1,J))
              !tt1=h(I,J); tt2=h(I+1,J)
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
          ENDIF
         ENDDO
       ENDDO
       hy = max(hy,0.0D0)
       ENDIF

       END SUBROUTINE interfaceH
      
      END MODULE FLOW

