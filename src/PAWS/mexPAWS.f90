#include "fintrf.h"
      
      
      ! This subroutine maintains the interactive interface
      ! this function should NOT be added into PAWS VS project because it requires mex lib
      subroutine mexFunction(nlhs, plhs, nrhs, prhs)
      ! prism(option, T, w, g, r)
      ! Input: VDZ, Soil
      ! Now use the same interface, but different procedures depending call mode (first input)
      USE Vdata
      use gd_Mod
      use matRW
      use MatTable, Only: connectField
      IMPLICIT NONE
      mwpointer :: plhs(*), prhs(*)
      mwpointer :: thM
      mwPointer :: pm, Cond, pm2,DR,C,ho,SN,UBCT,UBCV
      mwPointer :: mxGetField, mxGetDimensions, mxGetM, mxGetN
      Integer :: nlhs, nrhs, I,J, nz, m(3), ntype,mc(3),NBC,nG(2)
      INTEGER:: mode,MAXIT,TI(10), mode_in, dig ! mode_in, the first digit of the input mode, is the program control, what follows are the case number for BC
      mwpointer :: GWM,VDZM,CondM,SoilsM,Mapp,HB
      REAL*8:: STATE1,STATE2,STATE3,PREC,VV,tf,H,TE(10),TA(10)
      REAL*8,TARGET::BV(7)
      REAL*8,ALLOCATABLE::WTENEW(:,:)
      TYPE(BC_type)::BBC(NSTRUCT)
      TYPE(VDZ_type),POINTER::V
      integer, save :: status = 0
      integer digit, modeBC
      REAL*8 T
      type(gd), pointer:: gdt
      real*8, pointer:: temp1d(:),temp2d(:,:),temp3d(:,:,:),            &
     &   temp4d(:,:,:,:)
      integer, pointer:: temp1L(:),temp2L(:,:),temp3L(:,:,:),           &
     &   temp4L(:,:,:,:)

      !mwPointer :: pGlobalG, pGlobalW, pGlobalR
      
      PREC = 1D-10
      I = 1; pm = mxGetPr(Prhs(I))
      CALL mxCopyPtrToReal8(pm, VV, 1); mode=int(VV)
            
      mode_in = int(mode/(10** digit(mode)) ) ! first digit of the number
      if (digit(mode) .eq. 0) then
        modeBC  = mode
      else
        modeBC  = mode - mode_in * 10 ** digit(mode)
      endif
      
      SELECT CASE(mode_in)
      CASE(1,1002) ! Run surface and subsurface flow process
      
      if (status .eq. 0) then
        I = 2; pm = mxGetPr(Prhs(I))
        CALL mxCopyPtrToReal8(pm, tf, 1); ! tf is final time
        I = 3; pGlobalW = Prhs(I)
        I = 4; pGlobalG = Prhs(I)
        I = 5; pGlobalR = Prhs(I)
        call WGRentrance()
!      else ! will waste some time but we are doing ovn_cv outside
!        I = 4; pGlobalG = Prhs(I)
!        call connectField(gd_base%g(1),pGlobalG,'g','entrance')
!        gd_ptr => gdG .G. 'OVN'        
!        call getptr(gd_ptr,'h',g%OVN%h)
!        call getptr(gd_ptr,'U',g%OVN%U)
!        call getptr(gd_ptr,'V',g%OVN%V)
!        call getptr(gd_ptr,'S',g%OVN%S)
!        call getptr(gd_ptr,'ho',g%OVN%ho)
!        call getptr(gd_ptr,'hd',g%OVN%hd)
!        call getptr(gd_ptr,'dP',g%OVN%dP)
!        call getptr(gd_ptr,'SN',g%OVN%SN)
!        call getptr(gd_ptr,'SNSAVE',g%OVN%SNSAVE)
!        call getptr(gd_ptr,'prcp',g%OVN%prcp)
!        call getptr(gd_ptr,'frac',g%OVN%frac)
!        call getptr(gd_ptr,'t',temp1d); g%OVN%t => temp1d(1)
!        call getptr(gd_ptr,'dt',temp1d); g%OVN%dt => temp1d(1)
!        !call getptr(gd_ptr,'Cond',g%OVN%Cond)
!        call getptr(gd_ptr,'Mann',g%OVN%Mann)
!        call getptr(gd_ptr,'code',temp1d); g%OVN%Code => temp1d(1)
!        call getptr(gd_ptr,'Rf',g%OVN%Rf)
!        call getptr(gd_ptr,'Ex',g%OVN%Ex)
!        call getptr(gd_ptr,'Qoc',g%OVN%Qoc)
!        !call getptr(gd_ptr,'dfrac',g%OVN%dfrac)
!        call getptr(gd_ptr,'dist',g%OVN%dist)
!        call getptr(gd_ptr,'hc',g%OVN%hc)
!        call getptr(gd_ptr,'hf',g%OVN%hf)
!        call getptr(gd_ptr,'SATH',g%OVN%SATH)
!        call getptr(gd_ptr,'OVF',g%OVN%OVF)
!        call getptr(gd_ptr,'ffrac',g%OVN%ffrac)
!        call getptr(gd_ptr,'hback',g%OVN%hback)
!        call getptr(gd_ptr,'hod',g%OVN%hod)
!        call getptr(gd_ptr,'FParm',g%OVN%FParm)
      endif

      END SELECT
      
      if (status .eq. 0) then
        call bridgeGlobalGd(gd_base%g(1))
        call allocVData()
        status = 1
      endif
      
      I = 2; pm = mxGetPr(Prhs(I))
      CALL mxCopyPtrToReal8(pm, VV, 1); T = VV
      

      SELECT CASE(mode_in)
      CASE(1) ! Run surface and subsurface flow process
        call timeStepping(modeBC,T)
      END SELECT
      
      END SUBROUTINE
      
      
      subroutine WGRentrance()
      use gd_Mod
      use matTable
      use matRW
      !mwPointer :: pGlobalG, pGlobalW, pGlobalR
      
      if (.not. associated(gd_base%p)) call allocate_gd(gd_base,10)
      ntemp = 100
      np = gd_base%np + 1
      call allocate_gd(gd_base%g(np),ntemp)
      
      gd_base%f(np) = 'entrance'
      gd_base%np    = np
      gd_base%name(1) = 'base'
      gd_base%g(np)%name(1) = 'entrance'
      
      if (.not. associated(gd_base%g)) then
         allocate(gd_base%g(nG))
      endif
      
      call connectField(gd_base%g(np),pGlobalW,'w','entrance')
      call connectField(gd_base%g(np),pGlobalG,'g','entrance')
      call connectField(gd_base%g(np),pGlobalR,'r','entrance')
      
      
      end subroutine
