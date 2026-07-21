#include <preproc.h> 

      subroutine prism()
      ! This function is the gateway to the Fortran programs for surface and subsurface processes.
      ! Depending on the value of the first input, mode, the gateway function calls different subroutines-
      ! the subroutines named link***, bridge*** are all interface functions between Fortran and Matlab.-
      ! (has been modified, moved to matRW.F90)
      ! In Fortran, the structures (derived types) are defined in vdata.F90, which contain similarly-named fields in pointer format.-
      ! Essentially, a memory is allocated while loading the .mat file (or its address passed in from matlab).-
      ! The interface function then obtains the addresses of the fields in Bridge*** subroutines and then-
      ! associate these addresses to Fortran pointers in link*** subroutines
      ! Each set of Bridge*** and Link*** subroutines work on a sub-object. As a result the program is modularized and object-based.
      ! To understand these linking codes, here is the explanation:
      ! Remember, anything that can be viewed in Matlab is an mxArray.
      ! Each matlab data (even double array) is data structure (mxArray) that contains slightly more info than its data like size, type and a-
      ! pointer which points to the starting address of the actual data.-
      ! mxArray cannot be directly used in C or Fortran due to this wrapping.-
      ! Rather, we need to specifically get the addresses to the data.-
      ! pr=mxGetPr gets you the address of the real data. This step is done in Bridge*** routines.-
      ! However, Fortran still cannot understand pr as it still needs to know what is the size and rank of the array.-
      ! This is done in Link*** subroutines with declare the sizes and ranks.-
      ! The calling function pass in %val(pr) that is the starting address of the data. (It will not work if you pass in pr without %val()!)
      ! Alternatively, you can use mxCopyPtrToReal8 to get the data to an array declared in Fortran. This will get you the value in this function itself without %val()
      ! However this is a hard-copy (data will be duplicated).-
      ! Therefore we only use it when there is
      ! In addition:
      ! mxGetField get you the pointer to a certain field (any field of a mxArray is also a mxArray) to a mxArray
      ! mxGetDimensions gets the address of the size properties of the mxArray.-
      ! It needs to work with mxCopyPtrToInteger4 to get the numbers into your fortran workspace
      ! In general, mwPointer defines appropriate size for the integers so no worry about 32bit or 64bit machines


      USE Vdata
      USE Flow
      IMPLICIT NONE
      integer :: mode, H
      real*8 :: tf, prec
      TYPE(VDZ_type),POINTER::V

      PREC = 1D-10
      !I = 1; pm = mxGetPr(Prhs(I))
      !CALL mxCopyPtrToReal8(pm, VV, 1); mode=int(VV)
      SELECT CASE(mode)
      CASE(1,101) ! Run surface process together with vadose zone process
        !I = 2; pm = mxGetPr(Prhs(I))
        !CALL mxCopyPtrToReal8(pm, tf, 1); ! tf is final time
        !I = 3; pm = mxGetPr(Prhs(I))
        !CALL mxCopyPtrToReal8(pm, H, 1);
        !I=4; CALL bridgeVeg(prhs(I),1)
        !I=5; CALL bridgeVDZ(prhs(I),1)
        !thM = mxGetField(prhs(I),1,'th')
        !CALL bridgeTH(thM,1)
        !gaM = mxGetField(prhs(I),1,'GA')
        !CALL bridgeGA(gaM,1)
        !I=6; CALL bridgeDay(prhs(I),1)
        !I=7; CALL bridgeGW(prhs(I),1)
        !I=7; CALL bridgeOVN(prhs(I),1)
        !I=8; HB=mxGetPr(prhs(I))
        !CondM = mxGetField(prhs(I),1,'Cond')
        !CALL bridgeCond(BBC,CondM,NBC)
        if (mode .eq. 1) then
          CALL land(tf,int(H),g%GW(1)%h)
        else
          V => g%VDZ
          CALL UNSAT(V%MAPP,V,g%OVN%H,g%GW(1)%h) !!!!!!!! TEMP!
        endif
      CASE(2,102) ! Call the 1 column VDZ1c solver
       ! 101 is testUnsat2
       ! Boundary Conditions are passed in here
       ! not considering runoff here
       !I = 2; CALL bridgeVDZ(prhs(I),1)
       !I = 3; pm = mxGetPr(Prhs(I))
       !CALL mxCopyPtrToReal8(pm, TA(1), 1); ! Upper BC Type
       !I = 4; pm = mxGetPr(Prhs(I))
       !CALL mxCopyPtrToReal8(pm, TA(2), 1); ! Upper BC Bdata
       !I = 5; pm = mxGetPr(Prhs(I))
       !CALL mxCopyPtrToReal8(pm, TA(3), 1); ! Lower BC Type
       !I = 6; pm = mxGetPr(Prhs(I))
       !CALL mxCopyPtrToReal8(pm, TA(4), 1); ! Lower BC BData
       !I = 7; pm = mxGetPr(Prhs(I))
       !CALL mxCopyPtrToReal8(pm, BV, 7); ! Upper BC BV, if don't know what it is, put /1 1 0 0 0 0/
       !V => VD(1)

       !TE = (/0D0,0D0,0D0,0D0,1D0,0D0,0D0,0D0,0D0,0D0/)
       !BBC(1)%JType = int(TA(1))
       !BBC(1)%BData = TA(2)
       !BBC(2)%JType = int(TA(3))
       !BBC(2)%BData = TA(4)
       !BBC(2)%BV => BV
       !MAXIT = 20
       !TI(1) = 0
       !CCode = .true.

!       if (mode .eq. 2) then
!       CALL vdz1c(V%m(3),V%m(3),V%E,V%DF(1,1),V%Dperc(1,1),V%h,V%THE,              &
!     &   V%K,V%CI,V%SRC,V%dt(1,1),V%dtP,MAXIT,0D0,V%m(3),V%THER,V%THES,             &
!     &   V%Lambda,V%alpha,V%n,V%KS,V%DZ,V%DDZ,BBC,                                  &
!     &   TE(1),TI(1),TE(3),TE(4),TE(5),TE(6),TE(7),TE(8),TE(9))
!       else
!         I = 2
!         gaM = mxGetField(prhs(I),1,'GA')
!         CALL bridgeGA(gaM,1)
!
!         I = 1; J=1; nzl = V%NZC(I,J)
!         P => V%GA%GCODE(2:3); II = 0;
!         RS1 = V%SRC(I,J,1)+V%h(I,J,1)/V%dtp(3)
!         if (RS1>P(1)*V%KS(I,J,2) .and. V%GA%GCODE(1)>0.5) THEN ! potentially use GGA
!           DO K = 2,nzl
!             IF (V%THE(I,J,K)<=V%THES(I,J,K)-P(2)) THEN ! simpler decision of suction conditions
!               II = II + 1
!               IF (II>=3) EXIT
!             ENDIF
!           ENDDO
!           IF (II>=3) THEN
!             ! Going into GGA
!             IF (V%SRC(I,J,1)<0) THEN
!               Prcp = 0D0
!               V%h(I,J,1) = V%h(I,J,1) - V%SRC(I,J,1)*V%dtp(3)
!             ELSE
!               Prcp = V%SRC(I,J,1)
!             ENDIF
!             CALL GGA(nzl,GA%WFP(I,J),V%DF(I,J),V%THE(I,J,:),V%THES(I,J,:),    &
!         &        V%THER(I,J,:),V%KS(I,J,:),V%h(I,J,:),V%DZ(I,J,:),Prcp,         &
!         &    V%SRC(I,J,:), V%E(I,J,:),GA%KBar(I,J,:),GA%FM(I,J,:),V%dtP(3),               &
!         &        V%alpha(I,J,:),V%N(I,J,:),V%LAMBDA(I,J,:),FL,GA%Code(I,J),      &
!         &        GA%THEW(I,J))
!             ! Wrapping up boundary condition for Lower columns
!             GA%Code(I,J) = 1D0
!             BBC(1)%JType = 2
!             !BBCI(1)%BData = FL ! now the output from GGA is THEavg, in terms of operator splitting, the flux should be 0
!             BBC(1)%BData = 0D0
!             BBC(1)%BData2= GA%WFP(I,J) ! for N1
!           ELSE
!             GA%Code(I,J) = 0D0 ! If not modified, then it would be 0
!           ENDIF ! second condition if
!         endif ! first condition if
!
!         if (CCode) then
!         CALL vdz1c(V%m(3),V%m(3),V%E,V%DF(1,1),V%Dperc(1,1),V%h,V%THE,              &
!     &   V%K,V%CI,V%SRC,V%dt(1,1),V%dtP,MAXIT,0D0,V%m(3),V%THER,V%THES,             &
!     &   V%Lambda,V%alpha,V%n,V%KS,V%DZ,V%DDZ,BBC,                                  &
!     &   TE(1),TI(1),TE(3),TE(4),TE(5),TE(6),TE(7),TE(8),TE(9))
!         endif
!       endif ! mode = 2 if
!
!       !SUBROUTINE vdz1c(NZ,NA,E,DF,Dperc,h,THE,K,CI,SRC,ddt,dtP,                &
!       !&   MAXIT,DRAIN,WTL,THER,THES,Lambda,alpha,n,KS,DZ,DDZ,BBCI,                  &
!       !&   R,STATE,GWL,CC,PRatio,ho,SN,Rf,m_err)
!
!      CASE(3,5) ! Call Groundwater Solver
!        I = 2
!        Mapp = mxGetPr(prhs(I))
!
!        I = 3;
!        pm=mxGetPr(mxGetField(prhs(I),1,'STATE'))
!        GWM = prhs(I)
!        CALL mxCopyPtrToReal8(pm, STATE1, 1)
!        IF (abs(STATE1)<PREC) CALL bridgeGW(prhs(I),1)
!        CondM = mxGetField(prhs(I),1,'Cond')
!        CALL bridgeCond(BBC,CondM,NBC)
!        G => GW(1)
!        CALL Aquifer(G,%val(MAPP),BBC,NBC,mode)
!
!     CASE(4,11,12,34) ! General unsaturated/unconfined problem, Vauclin test problem and Hilberts example, Comp3D
!        I = 2; Mapp = mxGetPr(Prhs(I))
!        I = 3; CALL bridgeVDZ(prhs(I),1)
!        I = 4; CALL bridgeGW(prhs(I),1)
!        CondM = mxGetField(prhs(I),1,'Cond')
!        CALL bridgeCond(BBC,CondM,NBC)
!        I = 5; UBCT = mxGetPr(Prhs(I))
!        I = 6; UBCV = mxGetPr(Prhs(I))
!        G => GW(1)
!        IF (mode .eq. 34) THEN
!          I = 3
!          gaM = mxGetField(prhs(I),1,'GA')
!          CALL bridgeGA(gaM,1)
!        ENDIF
!
!        CALL unconfined(VD(1)%m(1),VD(1)%m(2),VD(1)%m(3),%val(Mapp),               &
!     &        VD,G,mode,NBC,BBC,%val(UBCT),%val(UBCV))
!        !UNCONFINED(ny,nx,nz,MAPP,V,G,MODE)
!      CASE(6) ! overland flow
!        I=2; CALL bridgeOVN(prhs(I),1)
!        I=3; CALL bridgeTopo(prhs(I),1)
!        I=4; pm = mxGetPr(prhs(I)) ! d
!        CALL mxCopyPtrToReal8(pm, t2, 2);
!        O => OVN(1)
        !CALL ovn_cv(O%ms,O%h,O%U,O%V,O%S,%val(d),Topo%E,Topo%Ex,Topo%Ey,                 &
        !&            %val(Mapp),O%Mann,O%dt)
        ! can't do this directly because functions inside module flow cannot use %val(Mapp) as input
!        CALL ovn_cv(O%ms,O%h,O%U,O%V,O%S,t2,Topo%E,Topo%Ex,Topo%Ey                 &
!     &            ,O%Mann,O%dt)
!     CASE(99) ! test bench for linking subroutines, modify freely
!        I=2; CALL bridgeGA(prhs(I),1)
      END SELECT
      END SUBROUTINE prism
!=====================================================================================================
      SUBROUTINE timeStepping(modeBC, T)
      ! time stepping with only surface and subsurface flow processes
      ! used for test cases 
      use VData
      use memoryChunk, Only: tempStorage,mm,mn,nn
      use gd_Mod
      use displayMod, Only: display
      use flow
      use surface_proc_mod
      use budgeter_mod, Only: budgeter, nitemw
      use displayMod, Only: printCode
      use exportASCII
      implicit none
      integer modeBC
      integer, save :: status = 0
      integer nGho(3), n
      TYPE(GW2D_type),pointer::GW
      REAL*8 T, tt
      REAL*8, DIMENSION(:,:), POINTER:: CC, MAPD, WTE
      REAL*8, DIMENSION(:,:,:), POINTER, SAVE::h,E,D,K,WW,ST,hNew,DR,Kz,beta,S
      REAL*8, DIMENSION(:,:,:,:), POINTER, SAVE::tau
      Integer digit
      Integer :: test3D = 0 ! a choice between 0 == legacy 2D code and 1 == 3D code.
      character*32 :: rMB, rSig, wSig
      integer :: scalar(2), kRiv, c, rid, ntStep = 0
      integer, dimension(size(w%CR)) :: rL
      real*8, dimension(:,:), pointer :: mask=>NULL()
      integer, dimension(size(g%OVN%h,1)*size(g%OVN%h,2)):: filter
      logical*1 :: rivp, ovnp
      real*8, dimension(size(w%CR)) :: DTS
      real*8 :: minST = 5D-4, rdt =0D0,dtnew,dtOld
      integer i,j,sig
      real*8,save:: lastDt = 0D0
      real*8, parameter :: PREC = 1D-6
      
      mask => w%DM%mask2
      
      ! specific to groundwater flow problem:
      mm = (/g%VDZ%m(1),g%VDZ%m(2),g%VDZ%m(3),g%Veg%nt/)
      mn = mm(1)*mm(2); nn = mm(1)*mm(2)*mm(3)  ! in memoryChunk
      CALL tempStorage('S',S,size(g%VDZ%h,3)) ! not used
      
      call distWea(0,1) ! continuous forcing. H=0, JD whatever. load forcing to g.OVN.prcp
      
      IF (status .eq. 0) then
        open(UNIT=272,FILE='problemLogs.txt',STATUS='UNKNOWN',ACCESS='APPEND')
        call allocVData()
        if (modeBC .eq. 0) then
            gd_ptr => gd_base .G. 'entrance' .G. 'g' .G. 'GW'
        else
            gd_ptr => gdG .G. 'GW'
        endif    
        h => gd_ptr .FFF. 'h3D'
        E => gd_ptr .FFF. 'E3D'
        D => gd_ptr .FFF. 'D3D'
        K => gd_ptr .FFF. 'K3D'
        WW => gd_ptr .FFF. 'W3D'
        ST => gd_ptr .FFF. 'ST3D'
        hNew => gd_ptr .FFF. 'hNew3D'
        DR => gd_ptr .FFF. 'DR3D'
        n = size(K,3);  CALL tempStorage('Kz',Kz,n) ! not used
        status = 1
        Kz = K
        CALL tempStorage('MAPD',MAPD,1) ! not used
        MAPD = w%DM%MAP
        printCode = 1
        nitemw = 51
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
      ENDIF
        
      
      GW => g%GW(1)
      nGho(1:2)=g%DM%nGho(1:2); nGho(3)=1
      
      if (modeBC .eq. 2) then
        !MAPD(1,:)=1;MAPD(g%DM%msize(1),1)=1;
        MAPD(g%DM%Nxy(1):g%DM%Nxy(2),1)=1;
        MAPD(g%DM%Nxy(1):g%DM%Nxy(2),g%DM%msize(2))=1
      endif
      
      ! initialization of beta, tau, h3D and DR3D
      IF (w%t .eq. 0D0) then
      h => g%VDZ%h
      WTE => g%GW(1)%h
      do i = 2, g%VDZ%m(3)
        h(:,:,i) = WTE - g%VDZ%E(:,:,i)
      end do
      h(:,:,1) = 0
      call vanGenuchten(g%VDZ%h,g%VDZ%THES,g%VDZ%THER,g%VDZ%ALPHA,g%VDZ%LAMBDA,g%VDZ%N,g%VDZ%KS,int(0),g%VDZ%THE,S,g%VDZ%K,g%VDZ%h)
      endif
      ! end initialization
            ! riv_cv is called in Matlab: because not all BCs were implemented in the Fortran version
      ! to connect in the right order, riv_cv should be called first
      rL = w%CR
      
      do i = 1, size(rL)
        DTS(i) = r(rL(i))%Riv%dt
      enddo
          kRiv = 0
          RDT = g%OVN%dt
          dtNew = g%OVN%dt
          do while (RDT > PREC)   
            kRiv = kRiv+1
            call rMinDtSet(rL, RDT, dtNew, dtOld)
            
            call rAvgStages()
            if(.not. associated(g%Riv%Qoc)) allocate(g%Riv%Qoc(size(rL)))
            do i = 1, size(rL)
                call F_oc_dw(i)
            enddo
            do i = 1, size(rL)
                c = rL(i)
                call rLateral(c)
                call r_MacCormack(c)
                call recAcc(rMB,c)
                call rClearTFlux(c)
            enddo
            call recAcc(rSig,0)
            rivp = .false.
            do j = 1, size(w%CR)
                !if (~isempty(find(isnan(r(w.CR(j)).Riv.h), 1)) || ~isempty(find(r(w.CR(j)).Riv.h<0, 1)))
                if(any(isnan(r(w%CR(j))%Riv%h)) .or. any(r(w%CR(j))%Riv%h<0)) then
                    rivp = .true.
                    sig = 1
                    call find1DL(isnan(r(w%CR(j))%Riv%h),size(r(w%CR(j))%Riv%h),filter,k)
                    !call display( 'River NaN!!!, j:', j, ' k: ',k, ' location(1):',filter(1))
                    return
                endif
            enddo
            if (any(isnan(g%OVN%h))) then
                ovnp = .true.
                call find1DL(isnan(g%OVN%h),w%m(1)*w%m(2),filter,k)
                !call display( 'OVN Exchange NaN!!! k:', k,' location(1): ', filter(1))
            endif
            
            !if rivp  || ~isempty(find(isnan(g(gid).OVN.h), 1)) || subp || ovnp
            if(rivp .or. ovnp) then
                !warning('4');
                call ch_assert('wtrshd_day :: NaN detected')
                return
            endif
            RDT = RDT - dtNew
          enddo
          do i = 1, size(rL)
            r(rL(i))%Riv%dt = DTS(i)
          enddo   !% TEMPORARY
      
      DO WHILE (w%t < T)
      
        ! source term ! might need to be modified as we go
        ntStep = ntStep + 1
        
        CALL gClearTFlux()
        
        WHERE(w%DM%MASK)
        ! g%VDZ%SRC(:,:,1) = g%Wea%Day%Prcp(1,1) ! comment out previous adhoc code
        g%OVN%S = g%OVN%prcp / 86400.0D0
        g%Veg%Gprcp= g%OVN%prcp * w%dt  ! for budgeter
        ENDWHERE
        
        if (test3D .eq. 1) then 
          !CALL UNSAT(g%VDZ%MAPP,g%VDZ,CC, g%GW(1)%H, 2)
          CALL surface_proc(1,2) ! first H doesn't make a difference
          
          CALL Aquifer3D(GW%ny,GW%nx,size(h,3),nGho, h,E,D,K,Kz,WW,                           &
     &          ST, MAPD, GW%Cond,SIZE(GW%Cond),GW%dd,GW%dt,hNew,DR)
          h =  hNew
        else
          !CALL UNSAT(g%VDZ%MAPP,g%VDZ,CC, g%GW(1)%H, 1)
          CALL surface_proc(1,1)
          
          CALL GW_source(1)
          CALL GW_sol(1)
          
        endif
        
        scalar = (/1, 1/)
        
        wSig = 'w'
        call budgeter(wSig,scalar(1),mask,1,1)
        
        g%OVN%h = g%VDZ%h(:,:,1)
        CALL ovn_cv(g%OVN%ms,g%OVN%h,g%OVN%U,g%OVN%V,g%OVN%S,g%DM%d,g%Topo%E,g%Topo%Ex,g%Topo%Ey    &
     &    , g%OVN%Mann, w%DM%Mask,g%OVN%dt,lastDt,size(g%OVN%Cond),g%OVN%Cond) ! Does this handle BC?
        g%VDZ%h(:,:,1) = g%OVN%h
     ! SUBROUTINE Aquifer3D(ny,nx,nz,nG,h,E,D,K,KZ,Wm,Ss,MAPP,BBC,NBC,                       &
     !&             dd,dt,hN,DR)
        ! now just trying GW3D!
        
        w%t = w%t + w%dt
        if (mod(ntStep,5) .eq. 1) then
          !call writeDimMat(g%OVN%Qx,'g.OVN.Qx',w%t)
          !call writeDimMat(g%OVN%Qy,'g.OVN.Qy',w%t)
          call writeDimMat(g%OVN%h,'g.OVN.h',w%t)
          call writeDimMat(g%GW(1)%h,'g.GW(1).h',w%t)
          call writeDimMat(g%VDZ%THE,'g.VDZ.THE',w%t)
        endif
      ENDDO
      
      close(UNIT=272)

      END SUBROUTINE
!=====================================================================================================
      SUBROUTINE AquiferSystem(GW,Mapp,BBC,NBC,OPT)
      ! calls the 2D groundwater solver, either static or nonlinear storage coefficient
      ! OPT = 3 (static  -- this should be used in most cases), else (nonlinear storage)
      USE Vdata
      USE Flow

      !GW%h = GW%hNew
      !call find1DL((abs(GW%h-GW%hNew)>15),w%m(1)*w%m(2),filter,k1)
      END SUBROUTINE AquiferSystem
!=====================================================================================================
!=====================================================================================================
      SUBROUTINE Aquifer(GW,Mapp,BBC,NBC,OPT)
      ! calls the 2D groundwater solver, either static or nonlinear storage coefficient
      ! OPT = 3 (static  -- this should be used in most cases), else (nonlinear storage)
      USE Vdata
      USE Flow
      INTEGER NBC,nG(2),OPT,k1
      TYPE(GW2D_type)::GW
      TYPE(BC_type)::BBC(NBC)
      REAL*8,DIMENSION(GW%ny,GW%nx)::MAPP,WTENEW
      integer, dimension(size(g%OVN%h,1)*size(g%OVN%h,2)):: filter
      nG = (/1, 1/)
      WTENEW = 0.D0
      IF (OPT .eq. 3) THEN  ! linear 
         CALL Aquifer2D(GW%ny,GW%nx,nG, GW%h,GW%EB,GW%D,GW%K,GW%W,                           &
     &          GW%ST, Mapp, BBC,NBC,GW%dd,GW%dt,GW%hNew,GW%DR)
      ELSE                  ! nonlinear
         CALL Aquifer2DNL(GW%ny,GW%nx,nG, GW%h,GW%EB,GW%D,GW%K,GW%W,                         &
     &          GW%ST, Mapp, BBC,NBC,GW%dd,GW%dt,GW%hNew,GW%DR,                              &
     &          1D0,GW%NLC,GW%STO,GW%nc)
      ENDIF
      GW%h = GW%hNew
      !call find1DL((abs(GW%h-GW%hNew)>15),w%m(1)*w%m(2),filter,k1)
      END SUBROUTINE Aquifer
!=====================================================================================================

      SUBROUTINE LAND(tf,H,HB)
       ! This file calculates land surface processes
       ! it calls surface (hourly), richards every 15 min (or its dt) and groundwater layer every day
       USE VData
       USE FLOW
#if (defined CLM)
       USE PAWS_CLM_Mod, Only: PAWS_CLM,SURFACE_URB
#endif
       IMPLICIT NONE
       INTEGER, SAVE:: state=0  ! start from scratch, obsolete option
       !INTEGER, SAVE:: state=2 ! CLM.mat read in but values need initialization
       !INTEGER, SAVE:: state=3 ! restart run
       integer dir
       REAL*8 DTS(4),tf, PREC
       REAL*8 DTW,DTV,DTL
       integer nx,ny,nz,nt,H,nG(2),k,k1
       TYPE(BC_type)::BBC(size(g%GW))
       TYPE(VDZ_type),POINTER::V
       TYPE(VEGE_type),POINTER::Veg
       TYPE(OV_type),POINTER::OVV
       REAL*8,DIMENSION(g%VEG%ny,g%VEG%nx)::CC,GWL,DR
       REAL*8,DIMENSION(g%VEG%ny,g%VEG%nx)::HB
       integer, dimension(size(g%OVN%h,1)*size(g%OVN%h,2)):: filterDBG
       INTEGER I,J
       REAL*8,DIMENSION(:),POINTER:: tt2=>NULL()
       logical::dbgsave=.false.

       PREC = 1d-10

       V => g%VDZ
       Veg => g%VEG
       OVV => g%OVN

       IF (w%t-floor(w%t)<w%dt) THEN
         CALL UpdateVeg(Veg,V%MAPP,V,g%Wea%Day) ! Update Delta, Cpc
       ENDIF

       CALL prepSnow(g%VDZ%th%SParm)
#if (defined CLM)
       DO WHILE (veg%t < tf - PREC)
         if (w%t .eq. w%ModelStartTime) then
           state = 2
         elseif (state .eq. 0) then
           state = 3
         else
           state = 1
         endif
         CALL PAWS_CLM(1,state,H)
         veg%t = veg%t + veg%dt
         
#ifndef CLM_Only
         CALL SURFACE_URB(V%m(1),V%m(2),V%m(3),veg%nt,g%Wea%DAY%ns,                          &
     &                Veg,V,H,veg%dt,OVV%h,OVV%OVF)
#endif

       ENDDO
#endif
#if (defined PAWS_ONLY)
       DO WHILE (veg%t < tf - PREC)
         CALL SURFACE(V%m(1),V%m(2),V%m(3),veg%nt,g%Wea%DAY%ns,                          &
     &                Veg,V,H,veg%dt,OVV%h,OVV%OVF)
         veg%t = veg%t + veg%dt
       ENDDO
#endif
#ifndef CLM_Only
       tt2=>V%h(18,18,:)
       !SURFACE(m,n,nz,nt,ns,Veg,V,DAY,MAPP,H,dt)
       ! Convert ET into source term for the VDZ

       nG=(/1,1/)
       DO WHILE (V%t < tf - PREC)
       !  CALL UNCONFINED(G%ny,G%nx,V%m(3),nG,V%MAPP,                           &
       !&        V,G,BBC,NBC,1)
       GWL = HB
       CALL UNSAT(V%MAPP,V,CC,HB,0)
         V%t = V%t + V%dtP(3)
         !GW(1).tt = GW(1).tt + GW(1).dt
       ENDDO
#endif 
#if (defined CLM)
       CALL PAWS_CLM(2,state,H)
#endif
       
      END SUBROUTINE LAND


      SUBROUTINE updateVeg(Veg,MAPP,V,DD)
      ! Computes the layer of rooting depth, evaporation extinction depth, canopy interception, bare soil fraction and albedo
      Use Vdata
      IMPLICIT NONE
      TYPE(Vege_type):: Veg
      TYPE(VDZ_TYPE)::V
      TYPE(STA_type)::DD
      REAL*8::MAPP(veg%ny,veg%nx)
      REAL*8:: Rfrac(veg%nt)
      REAL*8,DIMENSION(:,:,:),POINTER::LAI=>NULL(),Delta=>NULL()
      REAL*8,DIMENSION(:,:),POINTER::Cpc=>NULL()
      REAL*8,POINTER::pp=>NULL()
      REAL*8,DIMENSION(:),POINTER::EE=>NULL()
      REAL*8 TMAX,TMIN,TAVG,SFTMP,PREC
      REAL*8 tt1,tt2,tt3
      INTEGER II,IT1,RMX
      REAL*8,POINTER::tt=>NULL(),pp1(:)=>NULL(),pp2(:)=>NULL(),pp3(:)=>NULL()
      INTEGER,POINTER::KK(:)=>NULL()

      INTEGER I,J,K,KKK
      PREC = 1D-6
      LAI => veg%LAI
      Delta => veg%Delta
      Cpc => veg%Cpc
      K = 0; KKK = 0;

      DO I = 1,veg%ny
      DO J = 1,veg%nx
        IF (MAPP(I,J)>PREC) THEN
          RMX = V%NZC(I,J)-1 ! correction: NZC-th cell is GW table, should not have roots
          EE => V%EB(I,J,:)
          DO K = 1,veg%nt
            Delta(I,J,K)=(1D0-dexp(-0.5D0*LAI(I,J,K)))
            IF (veg%RPT(I,J,K) .eq. 2) THEN
              Delta(I,J,K) = 1
            ENDIF
            !pp=>LAI(I,J,K)
            !pp=>veg%Frac(I,J,K)
            !pp=>veg%Root(I,J,K)
            IF (veg%Root(I,J,K)<=0) THEN ! class has not root extraction
               veg%RDFN(I,J,K)=1
            ELSE
              DO II = max(2,veg%RDFN(I,J,K)),RMX
              IF (veg%ROOT(I,J,K) <= (EE(1)-EE(II))) EXIT
              ENDDO
              veg%RDFN(I,J,K) = min(II,RMX)
              DO II = veg%RDFN(I,J,K),2,-1
              IF (veg%ROOT(I,J,K) >= (EE(1)-EE(II-1))) EXIT
              ENDDO
              veg%RDFN(I,J,K) = II
            ENDIF
           ! IF (veg%RDFN(I,J,K)>22) THEN
           !tt1=veg%ROOT(I,J,K); IT1=V%NZC(I,J)
           !IT1 = veg%RDFN(I,J,K) ! see what it is
           !PRINT*,3
           !endif
          ENDDO
          Rfrac=veg%Frac(I,J,:)*Delta(I,J,:)
          veg%BSoil(I,J,1)=1.D0-sum(RFrac) ! impervious and water already considered
          !pp1=>veg%Frac(I,J,:)
          !pp2=>veg%Delta(I,J,:)
      veg%albedo(I,J)=veg%BSoil(I,J,1)*0.17D0+0.23D0*(1D0-veg%BSoil(I,J,1))
          !Cpc(I,J)=(2.0D-4)*sum(Rfrac*LAI(I,J,:))
          Cpc(I,J)=(2.0D-3)*sum(Rfrac*LAI(I,J,:))
          !pp3 => LAI(I,J,:)
          !tt => Cpc(I,J)
          !KK => veg%RPT(I,J,:)
        ENDIF

      ENDDO
      ENDDO

      ! decide if snow or not
      SFTMP = 0.75D0 ! snow threshold. Reference Brubaker (1996, Hydrological processes)
      DO I = 1,DD%ns
        TMAX = MAXVAL(DD%temp(:,I))
        TMIN = MINVAL(DD%temp(:,I))
        DD%TAVG(I) = 0.5D0*(TMAX + TMIN)
        IF (DD%TAVG(I) < SFTMP) THEN
          DD%SNOW(I) = .true.
        ELSE
          DD%SNOW(I) = .false.
        ENDIF
      ENDDO


      END subroutine updateVeg

      SUBROUTINE prepSnow(Param)
      ! initialize some snow parameters in the module
      USE COMM
      REAL*8,DIMENSION(*):: Param
      tr = Param(1)     !  Temperature above which all is rain (3 C)
      ts = Param(2)     !  Temperature below which all is snow (-1 C)
      es = Param(3)     !  emmissivity of snow (nominally 0.99)
      cg  = Param(4)    !  Ground heat capacity (nominally 2.09 KJ/kg/C)
      z = Param(5)      !  Nominal meas. height for air temp. and humidity (2m)
      zo = Param(6)     !  Surface aerodynamic roughness (m)
      rho = Param(7)    !  Snow Density (Nominally 450 kg/m^3)
      rhog = Param(8)   !  Soil Density (nominally 1700 kg/m^3)
      lc = Param(9)     !  Liquid holding capacity of snow (0.05)
      ks = Param(10)    !  Snow Saturated hydraulic conductivity (20 !160 m/hr)
      de = Param(11)    !  Thermally active depth of soil (0.1 m)
      abg = Param(12)   !  Bare ground albedo  (0.25)
      avo = Param(13)   !  Visual new snow albedo (0.95)
      anir0 = Param(14) !  NIR new snow albedo (0.65)
	  lans =  Param(15) !  the thermal conductivity of fresh (dry) snow (0.0576 kJ/m/k/hr)
	  lang =  Param(16) !  the thermal conductivity of soil (9.68 kJ/m/k/hr)
	  wlf =  Param(17)  !  Low frequency fluctuation in deep snow/soil layer (1/4 w1 = 0.0654 radian/hr)
	  rd1 =  Param(18)  !  Amplitude correction coefficient of heat conduction (1)
      fstab = Param(19) !  Stability correction control parameter 0 = no corrections, 1 = full corrections
	  Tref = Param(20)  !  Reference temperature of soil layer in ground heat calculation input
	  dNewS = Param(21) !  The threshold depth of for new snow (0.001 m)
	  gsurf = Param(22) !  The fraction of surface melt that runs off (e.g. from a glacier)

      ComState = 1

      ! think about tempAvg!!!

      END subroutine prepSnow




      SUBROUTINE UNSAT(MAPP,V,C,HB, OPT)
      ! Unsaturated zone subsurface processes
      ! Loop through the valid cells to do:
      ! (1) update lowland storage
      ! (2) if applicable (large enough prcp, dry portion exists), use GGA to compute infiltration
      ! (3) use Picard iteration scheme with variable top boundary condition for soil water flow
      USE Vdata
      USE Flow
      USE gd_Mod
      use displayMod, Only: display
      USE memoryChunk, Only: TEMPZ,ITZ,TEMPL,ITL,FINDEMPTY,tempStorage
      IMPLICIT NONE
      TYPE(VDZ_type),target ::V
      TYPE(GW2D_type),POINTER :: GW
      TYPE(BC_type):: BBCI(2)
      TYPE(OV_type),POINTER :: O
      INTEGER, OPTIONAL :: OPT
      integer nzl,nz,TP,I,J,ny,nx,MAXIT,n0,N,N1
      integer, POINTER:: nzp=>NULL()
      integer:: MODE,K,II
      REAL*8:: MAPP(V%m(1),V%m(2)),Rm,GWL(V%m(1),V%m(2)),m_err
      REAL*8,DIMENSION(V%m(1),V%m(2)):: C,HB,SN,Rf
      REAL*8,DIMENSION(:,:),pointer::FFrac_Urb=>NULL(),RfI=>NULL(),RfF=>NULL(),EvapSoil=>null()
      REAL*8 z,h0,ddr,drh,FL,x,th,t1,t2,t3
      real*8, save :: tk1,tk2,tk3,tk4,tk5
      integer:: WTL(V%m(1),V%m(2)),Co,IP
      REAL*8,DIMENSION(:),POINTER::RL=>NULL(),RL2=>NULL(),hh=>NULL()
      REAL*8 RS1,Prcp,ho,VGW,Buf
      REAL*8,POINTER::P(:)=>NULL()
      REAL*8 TEMPS,TEMPM,tt
      integer,POINTER:: nzi=>NULL()
      LOGICAL GCODE,MOVE,SealedBottom
      logical::dbgsave =.false.
      REAL*8,POINTER  :: TN(:,:,:,:)=>NULL()
      INTEGER,POINTER :: TNState(:,:)=>NULL()
      INTEGER,SAVE:: TNI = 0,TNSTATEI = 0, TNURBI = 0
      ! 06/16 added
      ! testCode: ==0, watershed simulation with full processes
      ! == 1, skip all other processes, run original 2D version GW flow;
      ! == 2, skip all other processes, run new 3D version GW flow; 
      INTEGER:: testCode = 0; ! May be changed by optional input. 
      REAL*8, DIMENSION(:,:,:,:),POINTER :: LatData
      REAL*8, DIMENSION(:,:,:),POINTER :: ST3D,h3D,DR3D
      REAL*8, DIMENSION(:),pointer :: beta,tau
      REAL*8, DIMENSION(V%m(3)) ::cc
      TYPE(gd), pointer:: gdGW,gdVDZ
      MODE = 1
      MAXIT = 100; n0=0
      ny = V%m(1)
      nx = V%m(2)
      nz = V%m(3)
 
      !BBCI(2)%JType=5;
      K = 0
      GW => g%GW(1)
      O => g%OVN
      P => g%VDZ%GA%GCODE(2:4);

      IF (PRESENT(OPT)) testCode = OPT
      
      
      
      IF (TNSTATEI .eq. 0) THEN
        TNSTATEI = FINDEMPTY(ITL,1,size(ITL));
        TNI      = FINDEMPTY(ITZ,4,size(ITZ));
        TNSTATE => TEMPL(:,:,TNSTATEI);    ITL(TNSTATEI) = 1;
        TNSTATE  = 0
        TN      => TEMPZ(:,:,:,TNI:TNI+3); ITZ(TNI:TNI+3) = 1;
        TN       = -999999D0
        DO I = 1,ny
        DO J = 1,nx
        IF (MAPP(I,J)>0) THEN
        DO K = 1,nz
        TN(I,J,K,1) = g%VDZ%n(I,J,K)/(g%VDZ%n(I,J,K)-1D0)      ! 
        TN(I,J,K,2) = 1D0/TN(I,J,K,1)        ! 1-1/n, sometimes called n
        TN(I,J,K,4) = (2D0*g%VDZ%n(I,J,K)-1D0)/g%VDZ%n(I,J,K)
        ENDDO
        ENDIF
        ENDDO
        ENDDO
        
      ELSE
        TNSTATE => TEMPL(:,:,TNSTATEI);    
        TN      => TEMPZ(:,:,:,TNI:TNI+3);
      ENDIF
      
      
      call tempStorage('ffrac_urb',FFrac_Urb,1)
      call tempStorage('EvapSoil',EvapSoil,1)
      call tempStorage('RfF',RfF,1)
      call tempStorage('RfI',RfI,1)
      ! TEMPORARY code:
      IF (testCode .eq. 2) then
          gdGW => gdG .G. 'GW'  ! will only get the first one
          call getPtr(gdGW, 'ST3D', ST3D)
          call getPtr(gdGW, 'h3D', h3D)
          call getPtr(gdGW, 'DR3D', DR3D)
          gdVDZ => gdG .G. 'VDZ'
          call getPtr(gdVDZ, 'LatData', LatData)
          !FFrac_Urb(20,30)=1001;
      else
          LatData => TN ! not really used
      endif
      
      DO I = 1,ny
      DO J = 1,nx
      IF (MAPP(I,J)>0) THEN
      n0 = 0
      BBCI(1)%JType=5; !automatically decide
      ! print debug
      K=k+1
      !IF (J<=5) THEN
      ! BBCI(1)%JType=2
      ! BBCI(1)%BData=0.14791D0
      !ELSE
      ! BBCI(1)%JType=2
      ! BBCI(1)%BData=0
      !ENDIF
      nzl = V%NZC(I,J)
      hh => V%h(I,J,:)
      ! from the head of the bottom layer, compute THE(nzl). If less than a threshold
      ! then move the interface downward

      !nzl = min(nz,V%WTL(I,J)+4)
      !IF ((I .eq. 13) .AND. (J .eq. 54) .AND. V%t>=731218.41666D0) THEN
      IF ((I .eq. 77) .AND. (J .eq. 102)) THEN
      !RL => V%h(I,J,:)
      m_err = 1
      !write(*,*) V%t  !, RL
      !DEBUGG = .TRUE.
      endif
        
      IF (testCode > 0) then !! test code 00000000000000000000000000000000000000000000000
          ! testing 3D code, skipping all land surface procs. The IF block goes all the way until before vdz1c and 3D assembly
          ! This IF block should be commented out if development is finished
          nzi => V%NZC(I,J)
          ho = 0D0
      else
      
      RfI(I,J)=0D0
      RfF(I,J)=0D0
      
      ! groundwater exchange depression storage
      CALL Lowland(O%RF(I,J),O%hd(I,J),O%h(I,J),g%VEG%hI(I,J),                &
     &    HB(I,J),GW%dt,V%VParm,O%dP(I,J,:),g%VEG%imp(I,J),O%Ex(I,J))

      ! Test if going into GGA
      ! Note how source terms are handled, what if only available for infiltrating over a short time span
      ! Also note inheritance issue from last time step
      ! How h0, WFP and DF are updated, how mass balance is kept
      ! GCode is a signal whether to consider GGA at all.
      ! GA%Code is a signal whether last time step used GGA
      II = 0;
      RS1 = V%SRC(I,J,1)+V%h(I,J,1)/V%dtp(3)+min(0D0,sum(V%SRC(I,J,2:3)*V%DZ(I,J,2:3)) )
      if (RS1>P(1)*V%KS(I,J,2) .and. g%VDZ%GA%GCODE(1)>0.5) THEN ! potentially use GGA
        DO K = 2,nzl
          IF (V%THE(I,J,K)<=V%THES(I,J,K)-P(2)) THEN ! simpler decision of suction conditions
            II = II + 1
            IF (II>=3) EXIT
          ENDIF
        ENDDO
        IF (II>=2 .AND. (.NOT. ANY(V%ICE(I,J,2:nzl)>0D0))) THEN
          ! Going into GGA
          IF (V%SRC(I,J,1)<0) THEN
            Prcp = 0D0
            V%h(I,J,1) = V%h(I,J,1) - V%SRC(I,J,1)*V%dtp(3)
          ELSE
            Prcp = V%SRC(I,J,1)
          ENDIF
          BBCI(1)%JINDEX = V%h(I,J,1) ! save for debugging
          !if (debugg) tk1 = sum(V%THE(I,J,:)*V%DZ(I,J,:))
          CALL GGA(nzl,g%VDZ%GA%WFP(I,J),V%DF(I,J),V%THE(I,J,:),V%THES(I,J,:),    &
     &       V%THER(I,J,:),P(3)*V%KS(I,J,:),V%h(I,J,:),V%DZ(I,J,:),Prcp,     &
     &       V%SRC(I,J,:),V%E(I,J,:),g%VDZ%GA%KBar(I,J,:),g%VDZ%GA%FM(I,J,:),            &
     &       V%dtP(3),V%alpha(I,J,:),V%N(I,J,:),V%LAMBDA(I,J,:),FL,          &
     &       g%VDZ%GA%Code(I,J),g%VDZ%GA%THEW(I,J))
          !if (debugg) tk2 = sum(V%THE(I,J,:)*V%DZ(I,J,:))
          ! Wrapping up boundary condition for Lower columns
          g%VDZ%GA%Code(I,J) = 1D0
          BBCI(1)%JType = 2
          !BBCI(1)%BData = FL ! now the output from GGA is THEavg, in terms of operator splitting, the flux should be 0
          BBCI(1)%BData = 0D0
          BBCI(1)%BData2= g%VDZ%GA%WFP(I,J) ! for N1

          ! Runoff for unfiltrated water
         CALL Runoff(O%RF(I,J),O%ffrac(I,J),O%ho(I,J),O%h(I,J),V%dtP(3),       &
     &            g%VEG%hI(I,J),V%h(I,J,:),V%VParm,nzl,O%hback(I,J),          &
     &            g%VEG%Imp(I,J),1D0,RfI(I,J),RfF(I,J))
          !Runoff(RF,SN,ho,hc,dtP,hI,h,VParm,nzl,hback)
          V%R(I,J,2) =0D0
        ELSE
          g%VDZ%GA%Code(I,J) = 0D0
          BBCI(1)%JINDEX = 0D0 ! debug
        ENDIF
      else
        g%VDZ%GA%Code(I,J) = 0D0
        !BBCI(1)%JINDEX = 0D0 ! debug
      endif
      ENDIF ! testCode 00000000000000000000000000000000000000000000000

      if (g%VDZ%GA%Code(I,J) .eq. 1D0 .and. g%VDZ%GA%WFP(I,J)>=nzl) then
      else

        IF (g%VDZ%GA%Code(I,J) .eq. 0 .and. (g%VEG%hI(I,J)>1D-8 .or.                 &
     &        V%h(I,J,1)<O%h(I,J)-O%hback(I,J))) THEN ! hI or flow-in
        ! calculate run-on and urban runoff
        ! does not do ordinary runoff here because it will be calculated in vdz1c
        ! taking advantage of this fact (large argument < 0)
        ! sending in FFRAC_URB(I,J)
        ! Hard-wiring urban N to 0.011 (from TR-55)
        tt = (86400.0D0*g%Topo%S0(I,J)**0.5D0)/(0.011*g%OVN%dist(I,J))
        tt = tt*(g%VDZ%VParm(10)**(2.0D0/3.0D0))
        FFRAC_URB(I,J)  = dexp(-tt*g%VDZ%dtP(3))
        IF (testCode .eq. 0) CALL Runoff(O%RF(I,J),FFRAC_URB(I,J),      &
     &       O%ho(I,J),O%h(I,J),V%dtP(3),g%VEG%hI(I,J),V%h(I,J,:),          &
     &       V%VParm,nzl,O%hback(I,J),g%VEG%Imp(I,J),-1D0,RfI(I,J),RfF(I,J))
        ENDIF
        ho = max(O%ho(I,J),O%h(I,J)-O%hback(I,J))

      ! [1 ST,2 zl (=DZ(uc)/2+D(aquitard)),3 DR, 4, Kl, 5, hl(pressure head), 6 E(center) uc aquifer, 7,GW.DZ, 8 EB(1)]
      ! floating last cell
      ! general requirement:
      IF (testCode .eq. 0) then !!!!!!22222222222222222222222222
      VGW = min(V%GWL(I,J),HB(I,J))
      nzi => V%NZC(I,J)
      ! WTL = min(insertLoc(nz,V%EB(I,J,:),VGW,V%WTL(I,J),1)+1,nz)
      ! where water table is right now (between EB(WTL-1) & EB(WTL)), meaning different from inside vdz1c
      TH = (V%EB(I,J,1)-V%EB(I,J,nz))
      IF (TH > 10) THEN
        Buf = 0.5D0
      ELSE
        Buf = max(TH/20,0.1D0)
      ENDIF
      !TOBE = min(V%EB(I,J,nzi-1)+0.5D0,V%EB(I,J,nzi-3))
      IF (VGW > V%EB(I,J,nzi-2)+Buf) THEN
        nzi = max(nzi -1,5)
        MOVE = .TRUE.
      !ELSEIF (VGW < V%EB(I,J,nzi)+Buf .AND. nzi<nz-1) THEN
      ELSEIF (VGW < V%EB(I,J,nzi-1)+Buf .AND. nzi<nz) THEN
        nzi = min(nz,nzi +1) ! move down. nzi>=nz-1 EB(nz-1) is GW(1).EB
        V%h(I,J,nzi-1) = 0D0 ! mass balance, previous negative value is fictional
        IF (V%THE(I,J,nzi-1)<V%THES(I,J,nzi-1)) THEN
         tempm = 1.0D0 - 1.0D0/V%N(I,J,nzi-1)
         tempS = (V%THE(I,J,nzi-1)-V%THER(I,J,nzi-1))/                    &
     &          (V%THES(I,J,nzi-1)-V%THER(I,J,nzi-1))
         V%h(I,J,nzi-1) = ((TEMPS**(-1.0D0/TEMPM)-1.0D0)**                &
     &     (1.0D0/V%N(I,J,nzi-1)))/(-abs(V%alpha(I,J,nzi-1)))
        ENDIF
        MOVE = .TRUE.
      ELSE
        MOVE = .FALSE.
      ENDIF
      ! If changed, some coefficients need to be modified
      IF (MOVE) THEN
        !V%Cond(I,J,2) = V%Cond(I,J,2)+(V%DZ(I,J,nzi)/2)*(nzl-nzi) ! for some reason this generates negative value for some watershed. need to think why
        V%Cond(I,J,2) = max((g%VDZ%EB(I,J,nzl)+g%GW(1)%EB(I,J))/2-g%GW(2)%E(I,J),0.5D0)
        V%Cond(I,J,6) = 0.5D0*(V%EB(I,J,nzi-1)+V%EB(I,J,nz-1))
        V%Cond(I,J,7) = V%EB(I,J,nzi-1)-V%EB(I,J,nz-1)
        ! there is one ghost cell with VDZ.EB
      ENDIF
      IF (V%E(I,J,nzi-1)<V%Cond(I,J,6) .or. V%Cond(I,J,2)<0D0) then
      call display( 'unsat: shouldn''t be here!i:'//num2str(i)//'j:'//num2str(j))
      PAUSE
      STOP
      ENDIF
      
      else
         nzi = g%VDZ%m(3)-1
      endif !!!!!!!!! testCode 222222222222222222222

      IF (V%h(I,J,nzl)<=-9000 .OR. testCode .eq. 2) THEN  ! bottom is not saturated
        !IF (HB(I,J)<-9000) THEN
          BBCI(2)%JType = 2
          BBCI(2)%BData = 0.D0
          BBCI(2)%BV => V%Cond(I,J,:)
          IF (testCode .eq. 2) THEN
            ! ADI3D assemble
            ![1. beta; 2. tao1; 3. tao2; 4. q^n; 5. h^{*,n}]
            CC   = ST3D(I,J,:)/V%dtP(3)
            LatData(I,J,2:nzi,1) = 0.01D0 ! trial -- beta
            LatData(I,J,2:nzi,2) = LatData(I,J,2:nzi,1)/(CC+LatData(I,J,2:nzi,1))  ! tau1
            LatData(I,J,2:nzi,3) = (CC*h3D(I,J,2:nzi)+DR3D(I,J,2:nzi)) /(CC + LatData(I,J,2:nzi,1)) ! tau2
            LatData(I,J,2:nzi,4) = DR3D(I,J,2:nzi) ! qn
          ENDIF
        !ELSE
          !BBCI(2)%JType = 1
          !BBCI(2)%BData = HB(I,J)
        !ENDIF
      ELSE
        !BBCI(2)%JType=1;
        !BBCI(2)%BData= HB(I,J)
        SealedBottom = .false.
        if (nzi .EQ. nz .AND. testCode .eq. 0 .AND. max(V%GWL(I,J),HB(I,J))<max(V%E(I,J,nzi-2),V%EB(I,J,nzi)+0.5D0)) then
        ! last condition added 10/19/2013: Allow saturated cells (but for some reason HB is low) interact to break the insolvability issue with all Neumann conditions
        ! test if the column is very wet (and therefore may have solvability issues)
        ! further test
          tt = sum( (V%THES(I,J,2:nzi) - V%THE(I,J,2:nzi))*V%DZ(I,J,2:nzi) )/(V%E(I,J,1)-V%E(I,J,nzi))
          if (tt > 2D-4) THEN
            SealedBottom =  .TRUE.
          else
            SealedBottom =  .FALSE.
          endif
        endif
        if (SealedBottom) then
        ! Already at last cell (top of bedrock) and still drying. no more suction
          BBCI(2)%JType = 2
          BBCI(2)%BData = 0.D0
          BBCI(2)%BV => V%Cond(I,J,:)
        else
        
        if (testCode .eq. 1) then
          V%Cond(I,J,2) = V%DZ(I,J,nzi)
          V%Cond(I,J,6) = V%E(I,J,nzi)
          V%Cond(I,J,7) = V%DZ(I,J,nzi)
        endif

        BBCI(2)%JType=3
        BBCI(2)%BV => V%Cond(I,J,:)
        z = BBCI(2)%BV(6)-BBCI(2)%BV(8) ! 8 is groundwater EB
        !h0 = max(HB(I,J)-BBCI(2)%BV(8),1D-6)
        h0 = HB(I,J)-BBCI(2)%BV(8)
        if (h0<=1D-6) THEN
          BBCI(2)%BData = BBCI(2)%BV(8)+h0
        else
        ddr = V%Cond(I,J,3)/h0
        drH = 0.5D0*ddr/V%KS(I,J,nzi)*(h0**2D0-z**2D0)
        drH = max(drH,-h0+0.1D0) ! head loss cannot make the pressure head there <0
        BBCI(2)%BData = BBCI(2)%BV(8)+h0+drH
        BBCI(2)%BData = min(BBCI(2)%BData,V%E(I,J,1)+2)
        endif
        !BBCI(2)%Bdata = HB(I,J)
        
        endif ! end bottom drying if
      ENDIF
      
      IF (MOVE .AND. abs(V%SRC(I,J,nzi))>1D-13) THEN
      ! In this case, we won't be able to apply the source term in the last cell
      ! And it won't be greater than nz
        V%SRC(I,J,nzi-1) = V%SRC(I,J,nzi-1)+V%SRC(I,J,nzi)*V%DZ(I,J,nzi)/V%DZ(I,J,nzi-1)
        V%SRC(I,J,nzi) = 0D0
      ENDIF


          !if (debugg) tk3 = sum(V%THE(I,J,:)*V%DZ(I,J,:))
        CALL vdz1c(nzi,nzi,V%E(I,J,:),V%DF(I,J),V%Dperc(I,J),                 &
     &      hh,V%THE(I,J,:),V%K(I,J,:),V%CI(I,J,:),V%SRC(I,J,:),              &
     &      V%dt(I,J),V%dtp,MAXIT,  0.D0,V%WTL(I,J),V%THER(I,J,:),           &
     &      V%THES(I,J,:),V%Lambda(I,J,:),V%alpha(I,J,:),                    &
     &      V%N(I,J,:),V%KS(I,J,:),V%DZ(I,J,:),V%DDZ(I,J,:),                 &
     &         BBCI,V%R(I,J,:),n0,V%GWL(I,J),C(I,J),V%Pratio(I,J),               &
     &      ho,O%SN(I,J),O%Rf(I,J),m_err,V%ICE(I,J,:),TN(I,J,1:nzi,:),        &
     &     TNSTATE(I,J),g%VDZ%VParm(12),mode,LatData(I,J,:,:)) !R,STATE,GWL,CC) ! mode ==1, conventional; mode == 2, 3D
     
          !if (debugg) tk4 = sum(V%THE(I,J,:)*V%DZ(I,J,:))
        IF (V%DF(I,J)<-1D-13) THEN
          EvapSoil(I,J) = -V%DF(I,J)
          V%DF(I,J)=0D0
        ELSE
          EvapSoil(I,J) = 0D0
        ENDIF
        
        !IF (ISNAN(m_err) .or. (m_err>1d-7)) THEN
        IF (ANY(ISNAN(V%h(I,J,:))) .or. any(abs(hh)>1D4)) THEN
           WRITE(1001,*) I,J,V%t,nzl
        ENDIF
        
      endif ! (g%VDZ%GA%Code(I,J) .eq. 1D0 .and. g%VDZ%GA%WFP(I,J)>=nzl) if

      ENDIF  ! 'Mask' IF
      ENDDO
      ENDDO
      !CLOSE(3)
    END SUBROUTINE UNSAT
