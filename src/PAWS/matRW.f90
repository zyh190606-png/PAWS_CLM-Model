#include "fintrf.h"
 
    module matRW

    implicit none
    ! Functions
    public
    mwPointer :: pGlobalG, pGlobalW, pGlobalR
    
    mwPointer :: cR3, cR4, cL3, cL4, cPtr(20)
    
    mwPointer, external :: matOpen                   ! mwPointer matOpen(filename, mode)
                                                     ! character*(*) filename, mode

    mwPointer, external ::  mxCreateDoubleMatrix      ! mwPointer mxCreateDoubleMatrix(m, n, ComplexFlag)
                                                      ! mwSize m, n
                                                      ! integer*4 ComplexFlag
    mwPointer, external ::  mxCreateString            ! mwPointer mxCreateString(str)
                                                      ! character*(*) str

    mwPointer, external ::  matGetDir                 ! mwPointer matGetDir(mfp, num)
                                                      ! mwPointer mfp
                                                      ! integer*4 num

    mwPointer, external ::  matGetVariableInfo        ! mwPointer matGetVariableInfo(mfp, name);
                                                      ! mwPointer mfp
                                                      ! character*(*) name

    mwPointer, external ::  matGetNextVariableInfo    ! mwPointer matGetNextVariableInfo(mfp, name)
                                                      ! mwPointer mfp
                                                      ! character*(*) name

    mwPointer, external ::  matGetVariable            ! mwPointer matGetVariable(mfp, name)
                                                      ! mwPointer mfp
                                                      ! character*(*) name

    mwPointer, external ::  matGetNextVariable        ! mwPointer matGetNextVariable(mfp, name)
                                                      ! mwPointer mfp
                                                      ! character*(*) name

    mwPointer, external ::  mxGetNumberOfElements     ! mwPointer mxGetNumberOfElements(pm)
                                                      ! mwPointer pm

    mwPointer, external ::  mxGetDimensions           ! mwPointer mxGetDimensions(pm)
                                                      ! mwPointer pm

    mwSize, external ::  mxGetNumberOfDimensions      ! mwSize mxGetNumberOfDimensions(pm)
                                                      ! mwPointer pm

    mwPointer, external ::  mxGetPr                   ! mwPointer mxGetPr(pm)
                                                      ! mwPointer pm

    mwPointer, external ::  mxGetM                    ! mwPointer mxGetM(pm)
                                                      ! mwPointer pm

    mwPointer, external ::  mxGetN                    ! mwPointer mxGetN(pm)
                                                      ! mwPointer pm

    integer*4, external ::  matClose                  ! integer*4 matClose(mfp)
                                                      ! mwPointer mfp

    integer*4, external ::  matPutVariable            ! integer*4 matPutVariable(mfp, name, pm)
                                                      ! mwPointer mfp, pm
                                                      ! character*(*) name

    integer*4, external ::  matPutVariableAsGlobal    ! integer*4 matPutVariableAsGlobal(mfp, name, pm)
                                                      ! mwPointer mfp, pm
                                                      ! character*(*) name

    integer*4, external ::  matDeleteVariable         ! integer*4 matDeleteVariable(mfp, name)
                                                      ! mwPointer mfp
                                                      ! character*(*) name

    integer*4, external ::  mxIsEmpty                 ! integer*4 mxIsEmpty(pm)
                                                      ! mwPointer pm

    integer*4, external :: mxIsStruct                 ! integer*4 mxIsStruct(pm)
                                                      ! mwPointer pm

    integer*4, external :: mxIsCell                   ! integer*4 mxIsCell(pm)
                                                      ! mwPointer pm

    integer*4, external :: mxIsSparse                 ! integer*4 mxIsSparse(pm)
                                                      ! mwPointer pm

    mwSize, external :: mxGetNzmax                    ! mwSize mxGetNzmax(pm)
                                                      ! mwPointer pm

    mwPointer, external :: mxGetIr                    ! mwPointer mxGetIr(pm)
                                                      ! mwPointer pm

    mwPointer, external :: mxGetJc                    ! mwPointer mxGetJc(pm)
                                                      ! mwPointer pm

    mwPointer, external :: mxGetCell                  ! mwPointer mxGetCell(pm, index)
                                                      ! mwPointer pm
                                                      ! mwIndex index

    integer*4, external ::  mxGetNumberOfFields       ! integer*4 mxGetNumberOfFields(pm)
                                                      ! mwPointer pm

    mwPointer, external ::  mxGetField                ! mwPointer mxGetField(pm, index, fieldname)
                                                      ! mwPointer pm
                                                      ! mwIndex index
                                                      ! character*(*) fieldname

    character*32, external ::  mxGetFieldNameByNumber ! character*(*) mxGetFieldNameByNumber(pm, fieldnumber)
                                                      ! mwPointer pm
                                                      ! integer*4 fieldnumber

    integer*4, external :: mxGetFieldNumber           ! integer*4 mxGetFieldNumber(pm, fieldname)
                                                      ! mwPointer pm
                                                      ! character*(*) fieldname

    mwPointer, external :: mxGetFieldByNumber         ! mwPointer mxGetFieldByNumber(pm, index, fieldnumber)
                                                      ! mwPointer pm
                                                      ! mwIndex index
                                                      ! integer*4 fieldnumber
                                                      
    integer*4, external :: mxClassIDFromClassName     ! integer*4 mxClassIDFromClassName(classname)
                                                      ! character*(*) classname

    mwPointer, external :: mxCreateNumericArray       ! mwPointer, external :: mxCreateNumericArray(ndim, dims, classid,
                                                      ! ComplexFlag) 
                                                      ! mwSize ndim, dims
                                                      ! integer*4 classid, ComplexFlag                                                                                                            

    end module matRW


        
!==============================================================================================



    module pawsControl
    ! shared control variables (DEBUGG is placed in VData, that is an exception)
    ! these variables are for the convenience of run-time simulation control
    ! This is to allow us send in control at run time, so less compilation time
    ! this module will not replace preproc.h totally
    logical :: CLM_Only
    logical :: PAWS_Only
    logical :: NO_GWF
    logical :: NO_SURF
    real*8 CaliStartTime
    
    END MODULE pawsControl
    
    SUBROUTINE pawsControlDecide(control,cst)
    use matRW
    use pawsControl
    real*8 control(*), cst
    ! w%control: [1. NO_surf==1; 2. initializeTHE>0; 3. tsMapOption==1--saveVar]

    NO_SURF = (control(1) .eq. 1D0)
    CaliStartTime = cst
    !call display( 'CaliStartTime = '//num2str(CaliStartTime))
     
    END SUBROUTINE
    
    subroutine bridgeGlobal2(file_name)
    ! interface1 filename, buildMatTable, bridgeGlobalGd
    use matRW
    use Vdata
    use gd_Mod
    use MatTable
    implicit none
    character*(*) ::  file_name
    type(gd),pointer:: gd_tempLocal
    mwpointer, dimension(:), pointer :: arg_pa
    character*(100), dimension(:), pointer :: arg_names
    integer i
    
    call buildMatTable(file_name, arg_pa, arg_names) 
    
    ! legacy issue: whoever calls buildMatTable should links Global pointers.
    do i=1,size(arg_pa)
       if (arg_names(i) .eq. 'g') then
         pGlobalG = arg_pa(i)
       elseif (arg_names(i) .eq. 'w') then
         pGlobalW = arg_pa(i)
       elseif (arg_names(i) .eq. 'r') then
         pGlobalR = arg_pa(i)
       endif
    enddo
    
    gd_tempLocal => gd_base .G. file_name
    call bridgeGlobalGd(gd_tempLocal)
    
    end


    subroutine bridgeGlobalGd(gd_rt)
    ! link from gd to vdata. some remnant mex operations
    use matRW
    use VData
    use gd_Mod
    Use pawsControl
    mwPointer :: pG, pW, pR, pm, Rec
    mwIndex :: index = 1
    type(gd):: gd_rt
    real*8 :: r1
    
    call linkVdataFromGd(gd_rt)
    ! ad_hoc assemble 
    w%m(1:4) = (/g%VDZ%m, g%VEG%NT/)
    
    ! remnant -- too complex to handle currently
    pm = mxGetField(pGlobalG, index, 'Riv')
    call bridge_gRiv(pm)
    
    ! below is due to character
    Rec = mxGetField(pGlobalW, index, 'Rec')
    pm = mxGetField(Rec, index, 'robj')
    do i = 1, size(w%Rec%robj)
       index = i
       call bridge_robj(pm,index)
    end do
    
    ! this will be easy to do
    if (associated(w%CaliStartTime)) then
      call pawsControlDecide(w%control,w%CaliStartTime)
    else
      r1 = w%ModelStartTime
      call pawsControlDecide(w%control,r1)
    endif
    
    end subroutine bridgeGlobalGd
!==============================================================================================

        subroutine bridgeMChunk(file_name)
        use Vdata
        use matRW
        implicit none
        character(len=length_of_character) :: file_name
        character(len=300):: file_nameCLM
        integer*4 status
        mwPointer :: file_in, pm
        mwPointer,parameter :: nField = 11
        mwPointer :: m(5, nField), ptr(nField)
        mwSize :: number_of_dim
        character*32 :: Fields(nField)
        integer i
        i = index(file_name, '.')
        file_nameCLM =' '
        file_nameCLM(1:i-1) = file_name(1:i-1)
        file_nameCLM(i:i+7)='_CLM.mat'

        file_in = matOpen(trim(file_nameCLM), 'r')
        if (file_in .eq. 0) then
          call display( 'Can''t open '// trim(file_nameCLM)// ' for reading.')
          call display( '(Does the file exist?)')
          !pause
          stop
        end if

        Fields = (/ 'L1','L2','L3','L4','L5','R1','R2','R3','R4','R5','cLevelMap'/)
        do i = 1, nField
          pm = matGetVariable(file_in, trim(Fields(i)))
          cPtr(i) = pm
!          if (i .eq. 3) then
!            cL3 = pm
!          elseif (i .eq. 4) then
!            cL4 = pm
!          elseif (i .eq. 8) then
!            cR3 = pm
!          elseif (i .eq. 9) then
!            cR4 = pm
!          endif
          if (pm .ne. 0) then
          number_of_dim = mxGetNumberOfDimensions(pm)
          call mxCopyPtrToPtrArray(mxGetDimensions(pm),m(:,i),number_of_dim)
          ptr(i) = mxGetPr(pm)
          else 
          call CH_ASSERT('variable not found in bridgeMChunk, most likely cLevelMap which is added 09/12/12!!!!')
          endif
        enddo
        
        call link_MChunk(%val(ptr(1)),%val(ptr(2)),%val(ptr(3)),%val(ptr(4)),%val(ptr(5)),%val(ptr(6)), & 
       &            %val(ptr(7)),%val(ptr(8)),%val(ptr(9)),%val(ptr(10)),%val(ptr(11)),m,nField)
        
        end subroutine bridgeMChunk
        
!==============================================================================================        
        subroutine link_MChunk(L1m,L2m,L3m,L4m,L5m,R1m,R2m,R3m,R4m,R5m,cLevelMapm,m,nField)
        use memoryChunk
        use displayMod, Only: display
        implicit none        
        mwPointer :: nField
        mwPointer :: m(5, nField)
        integer,dimension(m(1,1)),target::L1m
        integer,dimension(m(1,2),m(2,2)),target::L2m
        integer,dimension(m(1,3),m(2,3),m(3,3)),target::L3m
        integer,dimension(m(1,4),m(2,4),m(3,4),m(4,4)),target::L4m
        integer,dimension(m(1,5),m(2,5),m(3,5),m(4,5),m(5,5)),target::L5m
        real*8,dimension(m(1,6)),target::R1m
        real*8,dimension(m(1,7),m(2,7)),target::R2m
        real*8,dimension(m(1,8),m(2,8),m(3,8)),target::R3m
        real*8,dimension(m(1,9),m(2,9),m(3,9),m(4,9)),target::R4m
        real*8,dimension(m(1,10),m(2,10),m(3,10),m(4,10),m(5,10)),target::R5m
        integer,dimension(m(1,11),m(2,11)),target:: cLevelMapm
        
        L1 => L1m
        L2 => L2m
        L3 => L3m
        L4 => L4m
        L5 => L5m
        R1 => R1m
        R2 => R2m
        R3 => R3m
        R4 => R4m
        R5 => R5m
        cLevelMap => cLevelMapm

        if ((Size(R3,1) == Size(R4,1)) .and. (Size(R3,2) == Size(R4,2)) .and. &  
       &        (Size(L3,1) == Size(L4,1)) .and. (Size(L3,2) == Size(L4,2)) ) then
            mn = size(R3,1)*size(R3,2)
            nn  = size(R4,1)*size(R4,2)*size(R4,3)            
        else
          call display( 'size not match!')
          stop
        endif
           
        end subroutine link_MChunk
        
        
!==============================================================================================
        subroutine bridgeGlobal(file_name)
        use matRW
        use Vdata
        implicit none
        character(len=length_of_character) ::  file_name
        integer*4 status
        mwPointer :: file_in,pGlobal,pm
        mwIndex :: index = 1
        mwIndex :: i
        mwSize :: number_of_dim
        mwPointer dims(3)
        
        !call test_mxCreateDoubleMatrix
        ! read data
        file_in = matOpen(trim(file_name), 'r')
        if (file_in .eq. 0) then
          call display( 'Can''t open '// trim(file_name)// ' for reading.')
          call display( '(Does the file exist?)')
          !pause
          stop
        end if
        pGlobal = matGetVariable(file_in, 'g')
        pGlobalG = pGlobal
        call bridge_globalG(pGlobal)
        
        pm = mxGetField(pGlobal, index, 'ET')
        call bridgeET(pm)

        pm = mxGetField(pGlobal, index, 'DM')
        call bridgeDM(pm, index, 'g')

        pm = mxGetField(pGlobal, index, 'VDZ')
        call bridgeVDZ(pm,1)

        pm = mxGetField(pGlobal, index, 'Soil')
        call bridgeSoils(pm)

        pm = mxGetField(pGlobal, index, 'GW')
        number_of_dim = mxGetNumberOfDimensions(pm)
        CALL mxCopyPtrToPtrArray(mxGetDimensions(pm),dims,number_of_dim)
        allocate(g%GW(dims(2)))
        do i = 1, dims(2)
          call bridgeGW(pm, i)
        end do

        pm = mxGetField(pGlobal, index, 'Veg')
        call bridgeVeg(pm, 1)

        pm = mxGetField(pGlobal, index, 'Topo')
        call bridgeTopo(pm, 1)

        pm = mxGetField(pGlobal, index, 'OVN')
        call bridgeOVN(pm, 1)

        pm = mxGetField(pGlobal, index, 'Wea')
        call bridgeWea(pm)

        pm = mxGetField(pGlobal, index, 'Riv')
        call bridge_gRiv(pm)

        pGlobal = matGetVariable(file_in, 'w')
        pGlobalW = pGlobal
        call bridge_w(pGlobal)
        
        pm = mxGetField(pGlobal, index, 'Stats')
        call bridge_wStats(pm)
        
        pm = mxGetField(pGlobal, index, 'Tolerances')
        call bridge_Tol(pm)
        
        pm = mxGetField(pGlobal, index, 'Stations')
        number_of_dim = mxGetNumberOfDimensions(pm)
        CALL mxCopyPtrToPtrArray(mxGetDimensions(pm),dims,number_of_dim)
        allocate(w%Stations(dims(2)))
        do i = 1, dims(2)
          call bridgeStations(pm,i)
        end do

        pm = mxGetField(pGlobal, index, 'tData')
        call bridge_tData(pm)

        pm = mxGetField(pGlobal, index, 'Rec')
        call bridge_Rec(pm)
        
        pm = mxGetField(pm, index, 'robj')
        number_of_dim = mxGetNumberOfDimensions(pm)
        call mxCopyPtrToPtrArray(mxGetDimensions(pm),dims,number_of_dim)
        allocate(w%Rec%robj(dims(2)))
        do i = 1, dims(2)
          call bridge_robj(pm,i)
        end do
        
        pm = mxGetField(pGlobal, index, 'Rec')
        pm = mxGetField(pm, index, 'acc') 
        call bridge_acc(pm) 
                
        pm = mxGetField(pGlobal, index, 'DM')
        call bridgeDM(pm, index, 'w')

        pGlobal = matGetVariable(file_in, 'r')
        pGlobalR = pGlobal
        number_of_dim = mxGetNumberOfDimensions(pGlobal)
        call mxCopyPtrToPtrArray(mxGetDimensions(pGlobal),dims,number_of_dim)
        allocate(r(dims(2)))

        do i = 1, dims(2)
          pm = mxGetField(pGlobal, i, 'status')
          if(pm/=0) then
            if(mxIsEmpty(pm)==0) then
              call link_rStatus(%val(mxGetPr(pm)),i)
            endif
            if(r(i)%status /= 0) then
              pm = mxGetField(pGlobal, i, 'Riv')
              if(pm/=0) then
                if(mxIsEmpty(pm)==0) then
                  call bridge_rRiv(pm, i)
                endif
              endif

              pm = mxGetField(pGlobal, i, 'Net')
              if(pm/=0) then
                if(mxIsEmpty(pm)==0) then
                  call bridge_Net(pm, i)
                endif
              endif

              pm = mxGetField(pGlobal, i ,'DM')
              if(pm/=0) then
                if(mxIsEmpty(pm)==0) then
                  call bridgeDM(pm, i, 'r')
                endif
              endif    
            else
              call nullify_r_pointers(i)           
            endif             
          else
            call nullify_r_pointers(i)
          endif
                   
        end do


        ! close file
        status = matClose(file_in)
        if (status .ne. 0) then
          call display( 'Error closing MAT-file')
          !pause
          stop
        end if

        end subroutine bridgeGlobal


!====================================================================================================
        subroutine bridge_wStats(Stats)
        use matRW
        implicit none
        mwPointer :: Stats
        mwPointer :: THES, THE_FULL, THED, TA, TAV
        mwIndex :: index = 1
        
        THES = mxGetPr(mxGetField(Stats, index, 'THES'))
        THE_FULL = mxGetPr(mxGetField(Stats, index, 'THE_FULL'))
        THED = mxGetPr(mxGetField(Stats, index, 'THED'))
        TA = mxGetPr(mxGetField(Stats, index, 'TA'))
        TAV = mxGetPr(mxGetField(Stats, index, 'TAV'))
        
        call link_wStats(%val(THES),%val(THE_FULL),%val(THED),%val(TA),%val(TAV))
        
        end subroutine bridge_wStats
        
!====================================================================================================
        subroutine link_wStats(THES,THE_FULL,THED,TA,TAV)
        use Vdata
        implicit none
        real*8,target::THES,THE_FULL,THED,TA,TAV
        
        w%Stats%THES => THES
        w%Stats%THE_FULL => THE_FULL
        w%Stats%THED => THED
        w%Stats%TA => TA
        w%Stats%TAV => TAV

        end subroutine link_wStats        
!====================================================================================================
        subroutine bridgeET(ET)
        use matRW
        implicit none
        mwPointer :: ET, t, dt
        mwIndex :: index = 1
        
        t = mxGetPr(mxGetField(ET, index, 't'))
        dt = mxGetPr(mxGetField(ET, index, 'dt'))
        
        call linkET(%val(t), %val(dt))
        
        end subroutine bridgeET        
!====================================================================================================
        subroutine linkET(t, dt)
        use Vdata
        implicit none
        real*8, target :: t, dt
        g%ET%t => t
        g%ET%dt => dt
        end subroutine linkET
!====================================================================================================
        subroutine bridge_globalG(g)
        use matRW
        implicit none
        mwPointer :: g, t, dt
        mwIndex :: index = 1
        
        t = mxGetPr(mxGetField(g, index, 't'))
        dt = mxGetPr(mxGetField(g, index, 'dt'))
        
        call link_globalG(%val(t), %val(dt))
        
        end subroutine bridge_globalG        
!====================================================================================================
        subroutine link_globalG(t, dt)
        use Vdata
        implicit none
        real*8, target :: t, dt
        g%t => t
        g%dt => dt
        end subroutine link_globalG        
!====================================================================================================
        subroutine nullify_r_pointers(i)
        use Vdata
        implicit none
        !mwPointer i
        integer i
        type(rRiv_type) Riv
        type(Net_type) Net
        type(DM_TYPE) :: DM
        
        if (associated(r(i)%Riv%h)) then
        nullify(r(i)%Riv%h,r(i)%Riv%hx,r(i)%Riv%Ac,r(i)%Riv%U,r(i)%Riv%w,r(i)%Riv%wF,r(i)%Riv%R)
        nullify(r(i)%Riv%E,r(i)%Riv%Eta,r(i)%Riv%Ex,r(i)%Riv%ZBank,r(i)%Riv%Qx,r(i)%Riv%EstE)
        nullify(r(i)%Riv%Mann,r(i)%Riv%S,r(i)%Riv%SE,r(i)%Riv%S0,r(i)%Riv%SM,r(i)%Riv%DSV,r(i)%Riv%K)
        nullify(r(i)%Riv%hg,r(i)%Riv%mr,r(i)%Riv%fG,r(i)%Riv%Qoc,r(i)%Riv%Qgc,r(i)%Riv%Ql)
        nullify(r(i)%Riv%prcp,r(i)%Riv%evap,r(i)%Riv%trib,r(i)%Riv%Cond,r(i)%Riv%Ab)
        nullify(r(i)%Riv%rtr,r(i)%Riv%X,r(i)%Riv%Y,r(i)%Riv%DDist,r(i)%Riv%CType)
        nullify(r(i)%Riv%UPST,r(i)%Riv%UPSV,r(i)%Riv%DST,r(i)%Riv%Dx,r(i)%Riv%code)
        nullify(r(i)%Riv%l,r(i)%Riv%t,r(i)%Riv%dt,r(i)%Riv%umax,r(i)%Riv%State)

        nullify(r(i)%Net%level,r(i)%Net%UPSC,r(i)%Net%INSERT,r(i)%Net%Parent,r(i)%Net%Trib,r(i)%Net%TribC)
        nullify(r(i)%Net%DSQ, r(i)%Net%UPST, r(i)%Net%UPS, r(i)%Net%DS, r(i)%Net%DSC, r(i)%Net%UPSQ)
      
        nullify(r(i)%DM%Nxy, r(i)%DM%dim, r(i)%DM%Bindex, r(i)%DM%origin_idx)
        nullify(r(i)%DM%origin, r(i)%DM%x, r(i)%DM%y, r(i)%DM%Pbound, r(i)%DM%Cbound, r(i)%DM%Ebound)
        nullify(r(i)%DM%Map, r(i)%DM%MapL)
        nullify(r(i)%DM%mask, r(i)%DM%mask2)
        nullify(r(i)%DM%d, r(i)%DM%msize, r(i)%DM%nGho)
        endif

        end subroutine nullify_r_pointers

!====================================================================================================
        subroutine bridge_Rec(Rec)
        use matRW
        implicit none
        mwPointer :: Rec, pm
        mwPointer,parameter :: nField = 6
        mwPointer :: m(3, nField), ptr(nField)
        mwSize :: number_of_dim
        mwIndex :: index = 1
        character*32 :: Fields(nField)
        integer i
        
        Fields = (/ 'hG    ', 'Mem   ', 'MemR  ', 'MemQgc', 'hSave ', 'tss   '/)
        do i = 1, nField
          pm = mxGetField(Rec, index, trim(Fields(i)))
          number_of_dim = mxGetNumberOfDimensions(pm)
          call mxCopyPtrToPtrArray(mxGetDimensions(pm),m(:,i),number_of_dim)
          ptr(i) = mxGetPr(pm)
        enddo
        
        call link_Rec(%val(ptr(1)),%val(ptr(2)),%val(ptr(3)),%val(ptr(4)),%val(ptr(5)),%val(ptr(6)),m,nField)
        
        end subroutine bridge_Rec
!====================================================================================================
        subroutine link_Rec(hG, Mem, MemR, MemQgc, hSave, tss, m, nField)
        use Vdata
        implicit none
        mwPointer :: nField
        mwPointer :: m(3, nField)
        real*8, target :: hG(m(1,1),m(2,1)), Mem(m(2,2)), MemR(m(2,3)), MemQgc(m(2,4)), hSave(m(1,5),m(2,5)), tss(m(2,6))
        
        w%Rec%hG => hG
        w%Rec%Mem => Mem
        w%Rec%MemR => MemR
        w%Rec%MemQgc => MemQgc
        w%Rec%hSave => hSave
        w%Rec%tss => tss
        
        end subroutine link_Rec        
!====================================================================================================
        subroutine bridge_acc(acc)        
        use matRW
        implicit none
        mwPointer :: acc
        mwIndex :: index = 1
        mwPointer :: Qoc, Qout, trib
        mwPointer :: m(3)
        mwSize :: number_of_dim
        
        Qoc = mxGetField(acc, index, 'Qoc')
        number_of_dim = mxGetNumberOfDimensions(Qoc) 
        call mxCopyPtrToPtrArray(mxGetDimensions(Qoc),m,number_of_dim)
        Qoc = mxGetPr(Qoc)
        
        Qout = mxGetPr(mxGetField(acc, index, 'Qout'))
        trib = mxGetPr(mxGetField(acc, index, 'trib'))
        
        call link_acc(m, %val(Qoc),%val(Qout), %val(trib))      
         
        endsubroutine bridge_acc  
        
!====================================================================================================
        subroutine link_acc(m, Qoc, Qout, trib)
        use Vdata
        implicit none
        mwPointer :: m(3)
        real*8, dimension(m(2)), target :: Qoc, Qout, trib
        
        w%Rec%acc%Qoc => Qoc
        w%Rec%acc%Qout => Qout
        w%Rec%acc%trib => trib
        
        end subroutine link_acc        
!====================================================================================================
        subroutine link_rStatus(status, IS)
        use Vdata
        implicit none
        mwPointer :: IS
        real*8,target :: status
        
        r(IS)%status => status
        
        end subroutine link_rStatus        

!====================================================================================================
        subroutine bridge_Tol(Tolerances)
        use matRW
        implicit none
        mwPointer :: Tolerances, COU, Qmax
        mwIndex :: index = 1
        
        COU = mxGetPr(mxGetField(Tolerances, index, 'COU'))
        Qmax = mxGetPr(mxGetField(Tolerances, index, 'Qmax'))
        
        call link_Tol(%val(COU),%val(Qmax))
        
        end subroutine bridge_Tol

!====================================================================================================
        subroutine link_Tol(COU, Qmax)
        use Vdata
        implicit none
        real*8, target :: COU, Qmax
        
        w%Tolerances%COU => COU
        w%Tolerances%Qmax => Qmax
        
        end subroutine link_Tol
!====================================================================================================
        subroutine bridgeWea(Wea)
        use matRW
        implicit none
        mwPointer :: Wea, pm, pdt, pt
        mwIndex :: index = 1

        pm = mxGetField(Wea, index, 'day')
        call bridgeDay(pm, 1)
        
        pm = mxGetField(Wea, index, 'month')
        call bridgeMonth(pm)

        pdt = mxGetPr(mxGetField(Wea, index, 'dt'))
        pt = mxGetPr(mxGetField(Wea, index, 't'))
        call linkWea(%val(pdt),%val(pt))

        end subroutine bridgeWea

!====================================================================================================
        subroutine bridgeMonth(month)
        use matRW
        implicit none
        mwPointer :: month
        mwPointer :: prcp, temp, tmx, tmn, rad, wnd, hmd, dew
        mwIndex :: index = 1
        
        prcp = mxGetPr(mxGetField(month, index, 'prcp'))
        temp = mxGetPr(mxGetField(month, index, 'temp'))
        tmx = mxGetPr(mxGetField(month, index, 'tmx'))
        tmn = mxGetPr(mxGetField(month, index, 'tmn'))
        rad = mxGetPr(mxGetField(month, index, 'rad'))
        wnd = mxGetPr(mxGetField(month, index, 'wnd'))
        hmd = mxGetPr(mxGetField(month, index, 'hmd'))
        dew = mxGetPr(mxGetField(month, index, 'dew'))
        
        call linkMonth(%val(prcp),%val(temp),%val(tmx),%val(tmn),%val(rad),%val(wnd),%val(hmd),%val(dew))
        
        end subroutine bridgeMonth
        
!====================================================================================================
        subroutine linkMonth(prcp, temp, tmx, tmn, rad, wnd, hmd, dew)
        use Vdata
        implicit none
        real*8, target :: prcp, temp, tmx, tmn, rad, wnd, hmd, dew
        
        g%Wea%month%prcp => prcp
        g%Wea%month%temp => temp
        g%Wea%month%tmx => tmx
        g%Wea%month%tmn => tmn
        g%Wea%month%rad => rad
        g%Wea%month%wnd => wnd
        g%Wea%month%hmd => hmd
        g%Wea%month%dew => dew
        
        end subroutine linkMonth        
!====================================================================================================
        subroutine linkWea(dt,t)
        use Vdata
        implicit none
        real*8,target :: dt,t

        g%Wea%dt => dt
        g%Wea%t => t

        end subroutine linkWea
!====================================================================================================
        subroutine bridge_w(w)
        use matRW
        implicit none
        mwIndex :: index = 1
        mwPointer w
        mwPointer, parameter :: nFields = 9
        mwPointer i
        character*32 :: Fields(nFields)
        mwPointer :: ptr(nFields), pm, p2
        mwPointer :: m(3, nFields)
        mwSize number_of_dim
        integer*4 fn
        
        m = 0
        Fields = (/ 't','dt','ModelStartTime','ModelEndTime','plotT','g','CR','ER','control'/)

        do i = 1, nFields
          pm = mxGetField(w,index,TRIM(Fields(i)))
          if(pm /= 0) then
            if(mxISEmpty(pm)==0) then
              number_of_dim = mxGetNumberOfDimensions(pm)
              call mxCopyPtrToPtrArray(mxGetDimensions(pm),m(:,i),number_of_dim)
            end if
            ptr(i) = mxGetPr(pm)
          else
            ptr(i) = 0
          end if
        end do

        call link_w(%val(ptr(1)),%val(ptr(2)),%val(ptr(3)),%val(ptr(4)),%val(ptr(5)),%val(ptr(6)), &
       & %val(ptr(7)),%val(ptr(8)),m,nFields)

        fn = mxGetFieldNumber(w, 'CaliStartTime')
        if (fn>0) then
            pm = mxGetField(w,index,'CaliStartTime')
        else
            pm = mxGetField(w,index,'ModelStartTime')
        endif
        p2 = mxGetPr(pm)

        call pawsControlDecide(%val(ptr(9)),%val(p2))
        end subroutine bridge_w

!====================================================================================================
        subroutine link_w(t,dt,ModelStartTime,ModelEndTime,plotT,gm,CR,ER,dims,nFields)
        use Vdata
        use pawsControl
        implicit none
        mwPointer nFields
        mwPointer dims(3,nFields)
        real*8,target:: t, dt, ModelStartTime, ModelEndTime, plotT, gm
        integer, dimension(dims(2,7)),target :: CR
        integer, dimension(dims(2,8)),target :: ER
        integer i

        w%msg = ''
        w%t => t
        w%dt => dt
        w%ModelStartTime => ModelStartTime
        w%ModelEndTime => ModelEndTime
        w%plotT => plotT
        w%g => gm
        if(dims(2,7)/=0) w%CR => CR
        if(dims(2,8)/=0) w%ER => ER
        w%m(1:4) = (/g%VDZ%m, g%VEG%NT/)
        
        ! control variables
        
        end subroutine link_w
!====================================================================================================
        subroutine bridge_Net(Net, IS)
        use matRW
        implicit none
        mwPointer Net
        mwIndex :: IS, index = 1
        mwPointer, parameter :: nFields = 11
        mwPointer i
        character*32 :: Fields(nFields)
        mwPointer :: ptr(nFields), pm
        mwPointer :: m(3, nFields)
        mwSize number_of_dim
        
        m = 0
        Fields = (/ 'level','UPS','UPSC','UPSQ','DS','DSC','INSERT','DSQ','Parent','Trib','TribC'/)

        do i = 1, nFields
          pm = mxGetField(Net,index,TRIM(Fields(i)))
          if(pm /= 0) then
            if(mxISEmpty(pm)==0) then
              number_of_dim = mxGetNumberOfDimensions(pm)
              call mxCopyPtrToPtrArray(mxGetDimensions(pm),m(:,i),number_of_dim)
            end if
            ptr(i) = mxGetPr(pm)
          else
            ptr(i) = 0
          end if
        end do

        call link_Net(%val(ptr(1)),%val(ptr(2)),%val(ptr(3)),%val(ptr(4)),%val(ptr(5)),%val(ptr(6)), &
       & %val(ptr(7)),%val(ptr(8)),%val(ptr(9)),%val(ptr(10)),%val(ptr(11)),m,nFields,IS)

        end subroutine bridge_Net

!====================================================================================================
        subroutine link_Net(level,UPS,UPSC,UPSQ,DS,DSC,INSERT,DSQ,Parent,Trib,TribC,dims,nFields,IS)
        use VData
        implicit none
        mwPointer nFields
        mwPointer dims(3,nFields)
        mwIndex IS
        integer, dimension(dims(2,1)), target :: level
        real*8, target :: UPS
        integer, dimension(dims(2,3)), target :: UPSC
        real*8, target :: UPSQ
        real*8, target :: DS
        real*8, target :: DSC
        integer, dimension(dims(2,7)), target :: INSERT
        real*8, target :: DSQ
        integer, dimension(dims(2,9)), target :: Parent
        integer, dimension(dims(2,10)), target :: Trib
        integer, dimension(dims(2,11)), target :: TribC

        if(dims(2,1)/=0) r(IS)%Net%level => level
        r(IS)%Net%UPS => UPS
        if(dims(2,3)/=0) r(IS)%Net%UPSC => UPSC
        r(IS)%Net%UPSQ => UPSQ
        r(IS)%Net%DS => DS
        r(IS)%Net%DSC => DSC
        if(dims(2,7)/=0) r(IS)%Net%INSERT => INSERT
        r(IS)%Net%DSQ => DSQ
        if(dims(2,9)/=0) r(IS)%Net%Parent => Parent
        if(dims(2,10)/=0) r(IS)%Net%Trib => Trib
        if(dims(2,11)/=0) r(IS)%Net%TribC => TribC

        end subroutine link_Net
!====================================================================================================
        subroutine bridge_rRiv(Riv, IS)
        use matRW
        implicit none
        mwPointer Riv
        mwIndex :: IS, index = 1
        mwPointer, parameter :: nFields = 46
        mwPointer i
        character*32 :: Fields(nFields)
        mwPointer :: ptr(nFields), pm
        mwPointer :: m(3, nFields)
        mwSize number_of_dim
        
        m = 0
        Fields = (/ 'h','hx','Ac','U','w','wF','R','E','Eta','Ex','ZBank','Qx','EstE','DX', &
       &    'Mann','S','SE','S0','SM','code','l','CType','UPST','UPSV','DST','DSV','K', &
       &    'hg','mr','fG','Qoc','Qgc','Ql','prcp','evap','trib','Cond','Ab','rtr','t','dt', &
       &    'State','umax','X','Y','DDist'/)

        do i = 1, nFields
          pm = mxGetField(Riv,index,TRIM(Fields(i)))
          if(pm /= 0) then
            if(mxISEmpty(pm)==0) then
              number_of_dim = mxGetNumberOfDimensions(pm)
              call mxCopyPtrToPtrArray(mxGetDimensions(pm),m(:,i),number_of_dim)
            end if
            ptr(i) = mxGetPr(pm)
          else
            ptr(i) = 0
          end if
        end do

        call link_rRiv(%val(ptr(1)),%val(ptr(2)),%val(ptr(3)),%val(ptr(4)),%val(ptr(5)),%val(ptr(6)), &
       & %val(ptr(7)),%val(ptr(8)),%val(ptr(9)),%val(ptr(10)),%val(ptr(11)),%val(ptr(12)), &
       & %val(ptr(13)),%val(ptr(14)),%val(ptr(15)),%val(ptr(16)),%val(ptr(17)),%val(ptr(18)), &
       & %val(ptr(19)),%val(ptr(20)),%val(ptr(21)),%val(ptr(22)),%val(ptr(23)),%val(ptr(24)), &
       & %val(ptr(25)),%val(ptr(26)),%val(ptr(27)),%val(ptr(28)),%val(ptr(29)),%val(ptr(30)), &
       & %val(ptr(31)),%val(ptr(32)),%val(ptr(33)),%val(ptr(34)),%val(ptr(35)),%val(ptr(36)), &
       & %val(ptr(37)),%val(ptr(38)),%val(ptr(39)),%val(ptr(40)),%val(ptr(41)),%val(ptr(42)), &
       & %val(ptr(43)),%val(ptr(44)),%val(ptr(45)),%val(ptr(46)),m,nFields,IS)

        end subroutine bridge_rRiv

!====================================================================================================
        subroutine link_rRiv(h,hx,Ac,U,wm,wF,Rm,E,Eta,Ex,ZBank,Qx,EstE,Dx, &
       &    Mann,S,SE,S0,SM,code,l,CType,UPST,UPSV,DST,DSV,K, &
       &    hg,mr,fG,Qoc,Qgc,Ql,prcp,evap,trib,Cond,Ab,rtr,t,dt, &
       &    State,umax,X,Y,DDist,dims,nFields,IS)
        use VData
        implicit none
        mwPointer nFields
        mwPointer dims(3,nFields)
        mwIndex IS
        real*8, dimension(dims(2,1)), target :: h
        real*8, dimension(dims(2,2)), target :: hx
        real*8, dimension(dims(2,3)), target :: Ac
        real*8, dimension(dims(2,4)), target :: U
        real*8, dimension(dims(2,5)), target :: wm
        real*8, dimension(dims(2,6)), target :: wF
        real*8, dimension(dims(2,7)), target :: Rm
        real*8, dimension(dims(2,8)), target :: E
        real*8, dimension(dims(2,9)), target :: Eta
        real*8, dimension(dims(2,10)), target :: Ex
        real*8, dimension(dims(2,11)), target :: ZBank
        real*8, dimension(dims(2,12)), target :: Qx
        real*8, dimension(dims(2,13)), target :: EstE
        real*8, target :: DX
        real*8, dimension(dims(2,15)), target :: Mann
        real*8, dimension(dims(2,16)), target :: S
        real*8, dimension(dims(2,17)), target :: SE
        real*8, dimension(dims(2,18)), target :: S0
        real*8, dimension(dims(2,19)), target :: SM
        real*8, target :: code
        real*8, target :: l
        integer, dimension(dims(2,22)), target :: CType
        real*8, target :: UPST
        real*8, target :: UPSV
        real*8, target :: DST
        real*8, dimension(dims(2,26)), target :: DSV
        real*8, dimension(dims(2,27)), target :: K
        real*8, dimension(dims(2,28)), target :: hg
        real*8, dimension(dims(2,29)), target :: mr
        real*8, dimension(dims(2,30)), target :: fG
        real*8, dimension(dims(2,31)), target :: Qoc
        real*8, dimension(dims(2,32)), target :: Qgc
        real*8, dimension(dims(2,33)), target :: Ql
        real*8, dimension(dims(2,34)), target :: prcp
        real*8, dimension(dims(2,35)), target :: evap
        real*8, dimension(dims(2,36)), target :: trib
        real*8, dimension(dims(2,37)), target :: Cond
        real*8, dimension(dims(2,38)), target :: Ab
        real*8, dimension(dims(2,39)), target :: rtr
        real*8, target :: t
        real*8, target :: dt
        integer, target :: State
        real*8, target :: umax
        real*8, dimension(dims(2,44)), target :: X
        real*8, dimension(dims(2,45)), target :: Y
        real*8, dimension(dims(2,46)), target :: DDist

        if(dims(2,1)/=0) r(IS)%Riv%h => h
        if(dims(2,2)/=0) r(IS)%Riv%hx => hx
        if(dims(2,3)/=0) r(IS)%Riv%Ac => Ac
        if(dims(2,4)/=0) r(IS)%Riv%U => U
        if(dims(2,5)/=0) r(IS)%Riv%w => wm
        if(dims(2,6)/=0) r(IS)%Riv%wF => wF
        if(dims(2,7)/=0) r(IS)%Riv%R => Rm
        if(dims(2,8)/=0) r(IS)%Riv%E => E
        if(dims(2,9)/=0) r(IS)%Riv%Eta => Eta
        if(dims(2,10)/=0) r(IS)%Riv%Ex => Ex
        if(dims(2,11)/=0) r(IS)%Riv%ZBank => ZBank
        if(dims(2,12)/=0) r(IS)%Riv%Qx => Qx
        if(dims(2,13)/=0) r(IS)%Riv%EstE => EstE
        r(IS)%Riv%DX => DX
        if(dims(2,15)/=0) r(IS)%Riv%Mann => Mann
        if(dims(2,16)/=0) r(IS)%Riv%S => S
        if(dims(2,17)/=0) r(IS)%Riv%SE => SE
        if(dims(2,18)/=0) r(IS)%Riv%S0 => S0
        if(dims(2,19)/=0) r(IS)%Riv%SM => SM
        r(IS)%Riv%code => code
        r(IS)%Riv%l => l
        if(dims(2,22)/=0) r(IS)%Riv%CType => CType
        r(IS)%Riv%UPST => UPST
        r(IS)%Riv%UPSV => UPSV
        r(IS)%Riv%DST => DST
        if(dims(2,26)/=0) r(IS)%Riv%DSV => DSV
        r(IS)%Riv%DSV(1) = DSV(1)
        if(dims(2,27)/=0) r(IS)%Riv%K => K
        if(dims(2,28)/=0) r(IS)%Riv%hg => hg
        if(dims(2,29)/=0) r(IS)%Riv%mr => mr
        if(dims(2,30)/=0) r(IS)%Riv%fG => fG
        if(dims(2,31)/=0) r(IS)%Riv%Qoc => Qoc
        if(dims(2,32)/=0) r(IS)%Riv%Qgc => Qgc
        if(dims(2,33)/=0) r(IS)%Riv%Ql => Ql
        if(dims(2,34)/=0) r(IS)%Riv%prcp => prcp
        if(dims(2,35)/=0) r(IS)%Riv%evap => evap
        if(dims(2,36)/=0) r(IS)%Riv%trib => trib
        if(dims(2,37)/=0) r(IS)%Riv%Cond => Cond
        if(dims(2,38)/=0) r(IS)%Riv%Ab => Ab
        if(dims(2,39)/=0) r(IS)%Riv%rtr => rtr
        r(IS)%Riv%t => t
        r(IS)%Riv%dt => dt
        r(IS)%Riv%State => State
        r(IS)%Riv%umax => umax
        if(dims(2,44)/=0) r(IS)%Riv%X => X
        if(dims(2,45)/=0) r(IS)%Riv%Y => Y
        if(dims(2,46)/=0) r(IS)%Riv%DDist => DDist

        end subroutine link_rRiv


!====================================================================================================
        subroutine bridge_gRiv(Riv)
        use matRW
        implicit none
        mwPointer Riv, pm, ptr, pIr, pJc
        mwIndex :: index = 1, i, n
        mwSize :: number_of_dim, Nzmax, columns
        mwPointer :: m(3), dims(3)
        mwPointer :: ID, sides, K
        mwPointer, parameter :: nFields = 11
        character*32 :: Fields(nFields), Fields2(2)
        integer*4 iscell, issparse

        m = 0
        Fields = (/ 'ID','sides','Lidx', 'Cidx', 'ZBank', 'Len', 'Stages', 'W', 'RepLidx', 'K', 'Qoc'/)

        do n = 1, 2
            ptr = mxGetField(Riv, index, Fields(n))
            number_of_dim = mxGetNumberOfDimensions(ptr)
            call mxCopyPtrToPtrArray(mxGetDimensions(ptr),dims,number_of_dim)
            if (dims(1)*dims(2) .ne. 0) then
            if(Fields(n)=='ID') then
                call link_gRivInt(%val(mxGetPr(ptr)),m,dims,index,Fields(n))
            else
                call link_gRiv(%val(mxGetPr(ptr)),m,dims,index,Fields(n))
            endif
            endif
        end do

        do n = 3, nFields
            pm = mxGetField(Riv, index, Fields(n))
            number_of_dim = mxGetNumberOfDimensions(pm)
            call mxCopyPtrToPtrArray(mxGetDimensions(pm),m,number_of_dim)
            do i = 1, m(2)
              ptr = mxGetCell(pm, i)
              number_of_dim = mxGetNumberOfDimensions(ptr)
              call mxCopyPtrToPtrArray(mxGetDimensions(ptr),dims,number_of_dim)
              if(Fields(n)=='Lidx' .or. Fields(n)=='Cidx' .or. Fields(n)=='RepLidx') then
                call link_gRivInt(%val(mxGetPr(ptr)),m,dims,i,Fields(n))
              else
                call link_gRiv(%val(mxGetPr(ptr)),m,dims,i,Fields(n))
              endif
            end do
        end do

        Fields2 = (/'TM','TN'/)
        do n = 1, 2
            pm = mxGetField(Riv, index, Fields2(n))
            number_of_dim = mxGetNumberOfDimensions(pm)
            call mxCopyPtrToPtrArray(mxGetDimensions(pm),m,number_of_dim)
            do i = 1, m(2)
              ptr = mxGetCell(pm, i)
              Nzmax = mxGetNzmax(ptr)
              columns = mxGetN(ptr)
              pIr = mxGetIr(ptr)
              pJc = mxGetJc(ptr)
              number_of_dim = mxGetNumberOfDimensions(ptr)
              call mxCopyPtrToPtrArray(mxGetDimensions(ptr),dims,number_of_dim)
              call link_gRiv_TMTN(%val(mxGetPr(ptr)),%val(pIr),Nzmax,%val(pJc),columns+1,m,dims,i,Fields2(n))
            end do
        end do

        end subroutine bridge_gRiv

!====================================================================================================
        subroutine SparseToFull(values, Ir, Nzmax, Jc, columns, dims, ptr)
        implicit none
        mwSize Nzmax, columns
        real*8 values(Nzmax)
        mwPointer :: Ir(Nzmax), Jc(columns)
        mwPointer :: dims(3), n, i, j, row_start, row_end
        real*8, dimension(dims(1),dims(2)) :: ptr

        n = 0
        do j = 1, columns - 1
          row_start = Ir(Jc(j)+1)+1
          row_end = Ir(Jc(j+1))+1
          do i = row_start, row_end
            n = n + 1
            ptr(i,j) = values(n)
          end do
        end do

        end subroutine


!====================================================================================================
        subroutine link_gRiv_TMTN(values, Ir, Nzmax, Jc, columns, m, dims, index, Field)
        use Vdata
        implicit none
        mwSize Nzmax, columns
        mwIndex index
        real*8 values(Nzmax)
        mwPointer :: Ir(Nzmax), Jc(columns)
        mwPointer :: m(3), dims(3)
        character(len=*) Field

        select case (Field)
          case ('TM')
            if(.not. associated(g%Riv%TM)) then
              allocate(g%Riv%TM(m(2)))
            end if
            allocate(g%Riv%TM(index)%ptr(dims(1),dims(2)))
            g%Riv%TM(index)%ptr(:,:) = 0
            call SparseToFull(values, Ir, Nzmax, Jc, columns, dims, g%Riv%TM(index)%ptr)
          case ('TN')
            if(.not. associated(g%Riv%TN)) then
              allocate(g%Riv%TN(m(2)))
            end if
            allocate(g%Riv%TN(index)%ptr(dims(1),dims(2)))
            g%Riv%TN(index)%ptr(:,:) = 0
            call SparseToFull(values, Ir, Nzmax, Jc, columns, dims, g%Riv%TN(index)%ptr)
        end select

        end subroutine

!====================================================================================================
        subroutine link_gRiv(field_values, m, dims, index, field)
        use Vdata
        implicit none
        mwIndex :: index
        mwPointer :: m(3), dims(3)
        real*8, dimension(dims(1),dims(2)), target :: field_values
        character(len=*) field

        select case (field)
          case ('sides')
            g%Riv%sides => oneD_Ptr(size(field_values),field_values)
          case ('ZBank')
            if(.not. associated(g%Riv%ZBank)) then
              allocate(g%Riv%ZBank(m(2)))
            end if
            g%Riv%ZBank(index)%ptr => oneD_Ptr(size(field_values),field_values)
          case ('Len')
            if(.not. associated(g%Riv%Len)) then
              allocate(g%Riv%Len(m(2)))
            end if
            g%Riv%Len(index)%ptr => oneD_Ptr(size(field_values),field_values)
          case ('Stages')
            if(.not. associated(g%Riv%Stages)) then
              allocate(g%Riv%Stages(m(2)))
            end if
            g%Riv%Stages(index)%ptr => oneD_Ptr(size(field_values),field_values)
          case ('W')
            if(.not. associated(g%Riv%W)) then
              allocate(g%Riv%W(m(2)))
            end if
            g%Riv%W(index)%ptr => oneD_Ptr(size(field_values),field_values)
          case ('K')
            if(.not. associated(g%Riv%K)) then
              allocate(g%Riv%K(m(2)))
            end if
            g%Riv%K(index)%ptr => oneD_Ptr(size(field_values),field_values)
          case ('Qoc')
            if(.not. associated(g%Riv%Qoc)) then
              allocate(g%Riv%Qoc(m(2)))
            end if
            g%Riv%Qoc(index)%ptr => oneD_Ptr(size(field_values),field_values)            
        end select

        end subroutine link_gRiv


!====================================================================================================
        subroutine link_gRivInt(field_values, m, dims, index, field)
        use Vdata
        implicit none
        mwIndex :: index
        mwPointer :: m(3), dims(3)
        integer, dimension(dims(1),dims(2)), target :: field_values
        character(len=*) field

        select case (field)
          case ('ID')
            g%Riv%ID => oneD_ptrI(size(field_values),field_values)
          case ('Lidx')
            if(.not. associated(g%Riv%Lidx)) then
              allocate(g%Riv%Lidx(m(2)))
            end if
            g%Riv%Lidx(index)%ptr => oneD_ptrI(size(field_values),field_values)
          case ('Cidx')
            if(.not. associated(g%Riv%Cidx)) then
              allocate(g%Riv%Cidx(m(2)))
            end if
            g%Riv%Cidx(index)%ptr => oneD_ptrI(size(field_values),field_values)
          case ('RepLidx')
            if(.not. associated(g%Riv%RepLidx)) then
              allocate(g%Riv%RepLidx(m(2)))
            end if
            g%Riv%RepLidx(index)%ptr => oneD_ptrI(size(field_values),field_values)     
        end select

        end subroutine link_gRivInt
!====================================================================================================
        subroutine bridge_tData(tData)
        use matRW
        use Vdata
        implicit none
        mwPointer tData
        mwPointer pm,dc,h1save,h2save,objfun,m3R,m3E
        mwIndex :: index = 1, i
        mwPointer, parameter :: nField = 6
        mwPointer :: m(3,nField)
        mwSize :: number_of_dim
        mwPointer dims(3)
        
        pm = mxGetField(tData, index, 'usgs')
        if(pm/=0) then
            number_of_dim = mxGetNumberOfDimensions(pm)
            CALL mxCopyPtrToPtrArray(mxGetDimensions(pm),dims,number_of_dim)
            allocate(w%tData%usgs(dims(2)))
            do i = 1, dims(2)
              call bridge_usgs(pm,i)
            enddo
        endif

        pm = mxGetField(tData, index, 'usgs_gw')
        if(pm/=0) then
            number_of_dim = mxGetNumberOfDimensions(pm)
            CALL mxCopyPtrToPtrArray(mxGetDimensions(pm),dims,number_of_dim)
            allocate(w%tData%usgs_gw(dims(2)))
            do i = 1, dims(2)
              call bridge_usgs_gw(pm,i)
            enddo
        endif
        
        pm = mxGetField(tData, index, 'maps')
        call bridge_maps(pm)
        
        pm = mxGetField(tData, index, 'par')
        number_of_dim = mxGetNumberOfDimensions(pm)
        CALL mxCopyPtrToPtrArray(mxGetDimensions(pm),dims,number_of_dim)
        allocate(w%tData%par(dims(2)))
        do i = 1, dims(2)
          call bridge_par(pm,i)
        enddo
        
        m = 0
        pm = mxGetField(tData, index, 'dc')
        number_of_dim = mxGetNumberOfDimensions(pm)
        call mxCopyPtrToPtrArray(mxGetDimensions(pm),m(:,1),number_of_dim)
        dc = mxGetPr(pm)
        
        pm = mxGetField(tData, index, 'h1save')
        number_of_dim = mxGetNumberOfDimensions(pm)
        call mxCopyPtrToPtrArray(mxGetDimensions(pm),m(:,2),number_of_dim)
        h1save = mxGetPr(pm)
        
        pm = mxGetField(tData, index, 'h2save')
        number_of_dim = mxGetNumberOfDimensions(pm)
        call mxCopyPtrToPtrArray(mxGetDimensions(pm),m(:,3),number_of_dim)
        h2save = mxGetPr(pm)
        
        objfun = mxGetPr(mxGetField(tData, index, 'objfun'))
           
        pm = mxGetField(tData, index, 'm3')
        pm = mxGetField(pm, index, 'R')
        number_of_dim = mxGetNumberOfDimensions(pm)
        call mxCopyPtrToPtrArray(mxGetDimensions(pm),m(:,5),number_of_dim)
        m3R = mxGetPr(pm)
        
        pm = mxGetField(tData, index, 'm3')
        pm = mxGetField(pm, index, 'E')
        number_of_dim = mxGetNumberOfDimensions(pm)
        call mxCopyPtrToPtrArray(mxGetDimensions(pm),m(:,6),number_of_dim)
        m3E = mxGetPr(pm)

        call link_tData(%val(dc),%val(h1save),%val(h2save),%val(objfun),%val(m3R),%val(m3E),m,nField)

        end subroutine bridge_tData

!====================================================================================================
        subroutine link_tData(dc,h1save,h2save,objfun,m3R,m3E,dims,nField)
        use Vdata
        implicit none
        mwPointer :: nField
        mwPointer :: dims(3,nField)
        real*8, dimension(dims(1,1),dims(2,1)), target :: dc
        real*8, dimension(dims(1,2),dims(2,2)), target :: h1save
        real*8, dimension(dims(1,3),dims(2,3)), target :: h2save
        real*8, target :: objfun
        real*8, dimension(dims(1,5),dims(2,5)), target :: m3R
        real*8, dimension(dims(1,6),dims(2,6)), target :: m3E

        if(dims(1,1)/=0 .and. dims(2,1)/=0) w%tData%dc => dc
        if(dims(1,2)/=0 .and. dims(2,2)/=0) w%tData%h1save => h1save
        if(dims(1,3)/=0 .and. dims(2,3)/=0) w%tData%h2save => h2save
        w%tData%objfun => objfun
        if(dims(1,5)/=0 .and. dims(2,5)/=0) w%tData%m3%R => m3R
        if(dims(1,6)/=0 .and. dims(2,6)/=0) w%tData%m3%E => m3E

        end subroutine link_tData

!====================================================================================================
        subroutine bridge_usgs(usgs,index)
        use matRW
        implicit none
        mwPointer usgs
        mwPointer, parameter :: nFields = 8
        mwPointer i
        character*32 :: Fields(nFields)
        mwPointer :: ptr(nFields), pm
        mwPointer :: m(3, nFields)
        mwIndex :: index
        mwSize number_of_dim
        
        m = 0
        Fields = (/ 'id','name','dates','v','datenums','t','sd','ed'/)

        do i = 1, nFields
          pm = mxGetField(usgs,index,TRIM(Fields(i)))
          if(pm /= 0) then
            if(mxISEmpty(pm)==0) then
              number_of_dim = mxGetNumberOfDimensions(pm)
              call mxCopyPtrToPtrArray(mxGetDimensions(pm),m(:,i),number_of_dim)
            end if
            ptr(i) = mxGetPr(pm)
          else
            ptr(i) = 0
          end if
        end do

        call link_usgs(%val(ptr(1)),%val(ptr(3)),%val(ptr(4)),%val(ptr(5)), &
       & %val(ptr(6)),%val(ptr(7)),%val(ptr(8)),m,nFields,index)

        end subroutine bridge_usgs

!====================================================================================================
        subroutine link_usgs(id,dates,v,datenums,t,sd,ed,dims,nFields,index)
        use Vdata
        implicit none
        mwPointer nFields
        mwIndex :: index
        mwPointer dims(3,nFields)
        real*8, target :: id
        !character*32, dimension(dims(2,2)), target :: name
        real*8, dimension(dims(1,3)), target :: dates
        real*8, dimension(dims(1,4)), target :: v
        real*8, dimension(dims(1,5)), target :: datenums
        real*8, dimension(dims(1,6)), target :: t
        real*8, target :: sd
        real*8, target :: ed

        if(dims(1,1)/=0) w%tData%usgs(index)%id => id
        !w%tData%usgs%name => name
        if(dims(1,3)/=0) w%tData%usgs(index)%dates => dates
        if(dims(1,4)/=0) w%tData%usgs(index)%v => v
        if(dims(1,5)/=0) w%tData%usgs(index)%datenums => datenums
        if(dims(1,6)/=0) w%tData%usgs(index)%t => t
        if(dims(1,7)/=0) w%tData%usgs(index)%sd => sd
        if(dims(1,8)/=0) w%tData%usgs(index)%ed => ed

        end subroutine link_usgs
        
!====================================================================================================
        subroutine bridge_usgs_gw(usgs_gw,index)
        use matRW
        implicit none
        mwPointer usgs_gw
        mwPointer, parameter :: nFields = 10
        mwPointer i
        character*32 :: Fields(nFields)
        mwPointer :: ptr(nFields), pm
        mwPointer :: m(3, nFields)
        mwIndex :: index
        mwSize number_of_dim
        
        m = 0
        Fields = (/ 'dates','dtwt','datenums','v','t','XY','indices','elev','sd','ed'/)

        do i = 1, nFields
          pm = mxGetField(usgs_gw,index,TRIM(Fields(i)))
          if(pm /= 0) then
            if(mxISEmpty(pm)==0) then
              number_of_dim = mxGetNumberOfDimensions(pm)
              call mxCopyPtrToPtrArray(mxGetDimensions(pm),m(:,i),number_of_dim)
            end if
            ptr(i) = mxGetPr(pm)
          else
            ptr(i) = 0
          end if
        end do

        call link_usgs_gw(%val(ptr(1)),%val(ptr(2)),%val(ptr(3)),%val(ptr(4)),%val(ptr(5)), &
       & %val(ptr(6)),%val(ptr(7)),%val(ptr(8)),%val(ptr(9)),%val(ptr(10)),m,nFields,index)

        end subroutine bridge_usgs_gw

!====================================================================================================
        subroutine link_usgs_gw(dates,dtwt,datenums,v,t,XY,indices,elev,sd,ed,dims,nFields,index)
        use Vdata
        implicit none
        mwPointer nFields
        mwIndex :: index
        mwPointer dims(3,nFields)
        real*8, dimension(dims(1,1)), target :: dates
        real*8, dimension(dims(1,2)), target :: dtwt
        real*8, dimension(dims(1,3)), target :: datenums
        real*8, dimension(dims(1,4)), target :: v
        real*8, dimension(dims(1,5)), target :: t
        real*8, dimension(dims(2,6)), target :: XY
        real*8, dimension(dims(2,7)), target :: indices
        real*8, target :: elev, sd, ed

        if(dims(1,1)/=0) w%tData%usgs_gw(index)%dates => dates
        if(dims(1,2)/=0) w%tData%usgs_gw(index)%dtwt => dtwt
        if(dims(1,3)/=0) w%tData%usgs_gw(index)%datenums => datenums
        if(dims(1,4)/=0) w%tData%usgs_gw(index)%v => v
        if(dims(1,5)/=0) w%tData%usgs_gw(index)%t => t
        if(dims(1,6)/=0) w%tData%usgs_gw(index)%XY => XY
        if(dims(1,7)/=0) w%tData%usgs_gw(index)%indices => indices
        w%tData%usgs_gw(index)%elev => elev
        w%tData%usgs_gw(index)%sd => sd
        w%tData%usgs_gw(index)%ed => ed

        end subroutine link_usgs_gw      
        

!====================================================================================================
        subroutine bridge_maps(maps)
        use matRW
        use vData
        implicit none
        mwPointer maps
        mwPointer, parameter :: nFields = 18
        mwPointer i
        character*32 :: Fields(nFields)
        mwPointer :: ptr(nFields), pm
        mwPointer :: m(3, nFields)
        mwIndex :: index
        mwSize number_of_dim, nMaps(2)
        mwPointer :: nMdim = 2
        integer j
        
        m = 0
        Fields = (/ 'Prcp','Inf','PET','ET','EvapG','Dperc','Dperc2','SatE','InfE','botT',  &
                 &  'Qoc','Qgc','SUBM','SMELT','GW','ovnH','n','clmmap'/)

        number_of_dim = 2
        call mxCopyPtrToPtrArray(mxGetDimensions(maps),nMaps,number_of_dim)
        j = nMaps(1)* nMaps(2)
        allocate(w%tData%maps(j))
        
        do index = 1, nMaps(1)* nMaps(2) ! do not know whether if this is an array, if it is, which dimension it is on
        
        do i = 1, nFields
          pm = mxGetField(maps,index,TRIM(Fields(i)))
          if(pm /= 0) then
            if(mxISEmpty(pm)==0) then
              number_of_dim = mxGetNumberOfDimensions(pm)
              call mxCopyPtrToPtrArray(mxGetDimensions(pm),m(:,i),number_of_dim)
            end if
            ptr(i) = mxGetPr(pm)
          else
            ptr(i) = 0
          end if
        end do
        
        j = index
        call link_maps(w%tData%maps(j),%val(ptr(1)),%val(ptr(2)),%val(ptr(3)),%val(ptr(4)),%val(ptr(5)), &
       & %val(ptr(6)),%val(ptr(7)),%val(ptr(8)),%val(ptr(9)),%val(ptr(10)),%val(ptr(11)),&
       & %val(ptr(12)),%val(ptr(13)),%val(ptr(14)),%val(ptr(15)),%val(ptr(16)),%val(ptr(17))&
       & ,%val(ptr(18)),m,nFields)
        
        enddo

        end subroutine bridge_maps

!====================================================================================================
        subroutine link_maps(maps,Prcp,Inf,PET,ET,EvapG,Dperc,Dperc2,SatE,InfE,botT,Qoc,Qgc,SUBM,SMELT,GW,ovnH,n,clmmap,dims,nFields)
        use Vdata
        implicit none
        mwPointer nFields
        mwIndex :: index
        mwPointer dims(3,nFields)
        real*8, dimension(dims(1,1),dims(2,1)), target :: Prcp
        real*8, dimension(dims(1,2),dims(2,2)), target :: Inf
        real*8, dimension(dims(1,3),dims(2,3)), target :: PET
        real*8, dimension(dims(1,4),dims(2,4)), target :: ET
        real*8, dimension(dims(1,5),dims(2,5)), target :: EvapG
        real*8, dimension(dims(1,6),dims(2,6)), target :: Dperc
        real*8, dimension(dims(1,7),dims(2,7)), target :: Dperc2
        real*8, dimension(dims(1,8),dims(2,8)), target :: SatE
        real*8, dimension(dims(1,9),dims(2,9)), target :: InfE
        real*8, dimension(dims(1,10),dims(2,10)), target :: botT
        real*8, dimension(dims(1,11),dims(2,11)), target :: Qoc
        real*8, dimension(dims(1,12),dims(2,12)), target :: Qgc
        real*8, dimension(dims(1,13),dims(2,13)), target :: SUBM
        real*8, dimension(dims(1,14),dims(2,14)), target :: SMELT
        real*8, dimension(dims(1,15),dims(2,15)), target :: GW
        real*8, dimension(dims(1,16),dims(2,16)), target :: ovnH
        real*8, target :: n
        real*8, dimension(dims(1,18),dims(2,18),dims(3,18)), target :: clmmap
        type(maps_type) maps

        if(dims(1,1)/=0) maps%Prcp => Prcp
        if(dims(1,2)/=0) maps%Inf => Inf
        if(dims(1,3)/=0) maps%PET => PET
        if(dims(1,4)/=0) maps%ET => ET
        if(dims(1,5)/=0) maps%EvapG => EvapG 
        if(dims(1,6)/=0) maps%Dperc => Dperc
        if(dims(1,7)/=0) maps%Dperc2 => Dperc2
        if(dims(1,8)/=0) maps%SatE => SatE
        if(dims(1,9)/=0) maps%InfE => InfE
        if(dims(1,10)/=0) maps%botT => botT
        if(dims(1,11)/=0) maps%Qoc => Qoc
        if(dims(1,12)/=0) maps%Qgc => Qgc
        if(dims(1,13)/=0) maps%SUBM => SUBM
        if(dims(1,14)/=0) maps%SMELT => SMELT
        if(dims(1,15)/=0) maps%GW => GW
        if(dims(1,16)/=0) maps%ovnH => ovnH
        maps%n => n
        if(dims(1,18)/=0) maps%clmmap => clmmap

        end subroutine link_maps        
!====================================================================================================
        subroutine bridge_par(parm,index)
        use matRW
        implicit none
        mwPointer parm
        mwPointer, parameter :: nFields = 7
        mwPointer i
        character*32 :: Fields(nFields)
        mwPointer :: ptr(nFields), pm
        mwPointer :: m(3, nFields)
        mwIndex :: index
        mwSize number_of_dim
        
        m = 0
        Fields = (/ 'id','imet','code','mapp','v','ul','ll'/)

        do i = 1, nFields
          pm = mxGetField(parm,index,TRIM(Fields(i)))
          if(pm /= 0) then
            if(mxISEmpty(pm)==0) then
              number_of_dim = mxGetNumberOfDimensions(pm)
              call mxCopyPtrToPtrArray(mxGetDimensions(pm),m(:,i),number_of_dim)
            end if
            ptr(i) = mxGetPr(pm)
          else
            ptr(i) = 0
          end if
        end do

        call link_par(%val(ptr(1)),%val(ptr(2)),%val(ptr(3)),%val(ptr(4)),%val(ptr(5)), &
       & %val(ptr(6)),%val(ptr(7)),m,nFields,index)

        end subroutine bridge_par

!====================================================================================================
        subroutine link_par(id,imet,code,mapp,v,ul,ll,dims,nFields,index)
        use Vdata
        implicit none
        mwPointer nFields
        mwIndex :: index
        mwPointer dims(3,nFields)
        integer,target::id,imet,code
        integer,dimension(dims(1,4)),target::mapp
        real*8,target::v,ul,ll

        w%tData%par(index)%id => id
        w%tData%par(index)%imet => imet
        w%tData%par(index)%code => code
        if(dims(1,4)/=0) w%tData%par(index)%mapp => mapp
        w%tData%par(index)%v => v
        w%tData%par(index)%ul => ul
        w%tData%par(index)%ll => ll

        end subroutine link_par          
!====================================================================================================
        subroutine bridge_robj(robj,index)
        use matRW
        use Vdata
        implicit none
        mwPointer robj,index
        mwPointer, parameter :: nFields = 13
        mwPointer i
        character*32 :: Fields(nFields)
        mwPointer :: ptr(nFields), pm, fp
        mwPointer :: m(3, nFields), mf(3), dimf(3)
        mwSize number_of_dim
        
        m = 0
        Fields = (/ 'type','id','dist','ddist','xy','ix','sidx','objnum','acc','accT','v','cc','weight'/)

        do i = 1, nFields
          pm = mxGetField(robj,index,TRIM(Fields(i)))
          if(pm /= 0) then
            if(mxISEmpty(pm)==0) then
              number_of_dim = mxGetNumberOfDimensions(pm)
              call mxCopyPtrToPtrArray(mxGetDimensions(pm),m(:,i),number_of_dim)
            end if
            ptr(i) = mxGetPr(pm)
          else
            ptr(i) = 0
          end if
        end do

        call link_robj(%val(ptr(1)),%val(ptr(2)),%val(ptr(3)),%val(ptr(4)),%val(ptr(5)), &
       & %val(ptr(6)),%val(ptr(7)),%val(ptr(8)),%val(ptr(9)),%val(ptr(10)),%val(ptr(11)), &
       & %val(ptr(12)),%val(ptr(13)),m,nFields,index)
       
        pm = mxGetField(robj, index, 'Field')
        if(pm /= 0) then
          if(mxIsEmpty(pm) == 0) then
            number_of_dim = mxGetNumberOfDimensions(pm)
            call mxCopyPtrToPtrArray(mxGetDimensions(pm),mf,number_of_dim)
            allocate(w%Rec%robj(index)%Field(mf(2)))
            do i = 1, mf(2)
              fp = mxGetCell(pm, i)
              number_of_dim = mxGetNumberOfDimensions(fp)
              call mxCopyPtrToPtrArray(mxGetDimensions(fp),dimf,number_of_dim)
              call link_robjField(%val(mxGetPr(fp)),dimf(2),i,index)
            enddo
          endif
        endif

        end subroutine bridge_robj

!====================================================================================================
        subroutine link_robjField(Field, dim, index, IS)
        use Vdata
        implicit none
        mwPointer :: dim, index, IS, i
        integer(kind=2),dimension(dim),target :: Field
        
        w%Rec%robj(IS)%Field(index) = ''
        do i = 1, dim
          w%Rec%robj(IS)%Field(index) = trim(w%Rec%robj(IS)%Field(index))//char(Field(i))
        enddo
        end subroutine link_robjField
        
        
!====================================================================================================
        subroutine link_robj(type,id,dist,ddist,xy,ix,sidx,objnum,acc,accT,v,cc,weight,dims,nFields,IS)
        use Vdata
        implicit none
        mwPointer nFields,IS
        mwPointer :: dims(3,nFields)
        character, target :: type
        integer, target :: id
        real*8, target :: dist
        real*8, target :: ddist
        real*8, dimension(dims(2,5)), target :: xy
        integer, dimension(dims(2,6)), target :: ix
        integer, target :: sidx
        integer, target :: objnum
        real*8, dimension(dims(2,9)), target :: acc
        real*8, target :: accT
        real*8, dimension(dims(2,11)), target :: v
        real*8, target :: cc
        real*8, target :: weight

        w%Rec%robj(IS)%type => type
        if(dims(2,2)==0) then
          !allocate(w%Rec%robj(IS)%id)
          !w%Rec%robj(IS)%id = 0.0D0
        else
          w%Rec%robj(IS)%id => id
        endif
        if(dims(2,3)==0) then
          !allocate(w%Rec%robj(IS)%dist)
          !w%Rec%robj(IS)%dist = 0.0D0
        else
          w%Rec%robj(IS)%dist => dist
        endif
        if(dims(2,4)==0) then
          !allocate(w%Rec%robj(IS)%ddist)
          !w%Rec%robj(IS)%ddist = 0.0D0
        else
          w%Rec%robj(IS)%ddist => ddist
        endif
        if(dims(2,5)/=0) w%Rec%robj(IS)%xy => xy
        if(dims(2,6)/=0) w%Rec%robj(IS)%ix => ix
        if(dims(2,7)==0) then
          !allocate(w%Rec%robj(IS)%sidx)
          !w%Rec%robj(IS)%sidx = 0.0D0
        else
          w%Rec%robj(IS)%sidx => sidx
        endif
        if(dims(2,8)==0) then
          !allocate(w%Rec%robj(IS)%objnum)
          !w%Rec%robj(IS)%objnum = 0.0D0
        else
          w%Rec%robj(IS)%objnum => objnum
        endif
        if(dims(2,9)==0) then
          !allocate(w%Rec%robj(IS)%acc(1))
          !w%Rec%robj(IS)%acc = 0.0D0
        else
          w%Rec%robj(IS)%acc => acc
        endif
        if(dims(2,10)==0) then
          !allocate(w%Rec%robj(IS)%accT)
          !w%Rec%robj(IS)%accT = 0.0D0
        else
          w%Rec%robj(IS)%accT => accT
        endif        
        if(dims(2,11)/=0) w%Rec%robj(IS)%v => v
        if(dims(2,12)==0) then
          !allocate(w%Rec%robj(IS)%cc)
          !w%Rec%robj(IS)%cc = 0.0D0
        else
          w%Rec%robj(IS)%cc => cc
        endif
        if(dims(2,13)==0) then
          !allocate(w%Rec%robj(IS)%weight)
          !w%Rec%robj(IS)%weight = 0.0D0
        else
          w%Rec%robj(IS)%weight => weight
        endif

        end subroutine link_robj
!====================================================================================================
        subroutine bridgeStations(sta,index)
        use matRW
        implicit none
        mwPointer sta
        mwIndex index
        mwPointer, parameter :: nFields = 36
        mwPointer i
        character*32 :: Fields(nFields)
        mwPointer :: ptr(nFields), pm
        mwPointer :: m(3, nFields)
        mwSize number_of_dim
        
        m = 0
        Fields = (/ 'id','XYElev','LatLong','gIndices','wPct','dates','datenums', &
       &    'vratio','prcp','tmin','tmax','sol','wnd','hmd','dptp','awnd','evap','Rad','RD','RI', &
       &    'tau','rrad','s1','s0','snow','mnrh','mntp','mxrh','neighbor','ref','ptr','DData','pct', &
       &    'dat_time','cosz','decl'/)
        do i = 1, nFields
          pm = mxGetField(sta,index,TRIM(Fields(i)))
          if(pm /= 0) then
              if(mxISEmpty(pm)==0) then
                number_of_dim = mxGetNumberOfDimensions(pm)
                call mxCopyPtrToPtrArray(mxGetDimensions(pm),m(:,i),number_of_dim)
              else
                m(:,i) = 0
              end if
              ptr(i) = mxGetPr(pm)
          else
              ptr(i) = 0
              m(:,i) = 0
          end if
        end do

        call linkStations(%val(ptr(1)),%val(ptr(2)),%val(ptr(3)),%val(ptr(4)),%val(ptr(5)), &
       & %val(ptr(6)),%val(ptr(7)),%val(ptr(8)),%val(ptr(9)),%val(ptr(10)),%val(ptr(11)), &
       & %val(ptr(12)),%val(ptr(13)),%val(ptr(14)),%val(ptr(15)),%val(ptr(16)),%val(ptr(17)), &
       & %val(ptr(18)),%val(ptr(19)),%val(ptr(20)),%val(ptr(21)),%val(ptr(22)),%val(ptr(23)), &
       & %val(ptr(24)),%val(ptr(25)),%val(ptr(26)),%val(ptr(27)),%val(ptr(28)),%val(ptr(29)), &
       & %val(ptr(30)),%val(ptr(31)),%val(ptr(32)),%val(ptr(33)),%val(ptr(34)),%val(ptr(35)), &
       & %val(ptr(36)),m,nFields,index)

        end subroutine bridgeStations

!=====================================================================================================

        subroutine linkStations(id,XYElev,LatLong,gIndices,wPct,dates,datenums,vratio,prcp, &
       &    tmin,tmax,sol,wnd,hmd,dptp,awnd,evap,Rad,RD,RI,tau,rrad,s1,s0,snow,mnrh,mntp,mxrh,neighbor, &
       &    ref,ptr,DData,pct,dat_time,cosz,decl,dims,nFields,IS)
        use Vdata
        implicit none
        mwPointer nFields, IS
        mwPointer :: dims(3,nFields)
        integer, dimension(dims(2,1)),target :: id
        real*8, dimension(dims(2,2)),target :: XYElev
        real*8, dimension(dims(2,3)),target :: LatLong
        !character*32, dimension(dims(2,4)),target :: Name
        integer, dimension(dims(1,4)),target :: gIndices
        real*8, dimension(dims(1,5)),target :: wPct
        real*8, dimension(dims(1,6)),target :: dates
        real*8, dimension(dims(1,7)),target :: datenums
        real*8, dimension(dims(1,8)),target :: vratio
        real*8, dimension(dims(1,9)),target :: prcp
        real*8, dimension(dims(1,10)),target :: tmin
        real*8, dimension(dims(1,11)),target :: tmax
        real*8, dimension(dims(1,12)),target :: sol
        real*8, dimension(dims(1,13)),target :: wnd
        real*8, dimension(dims(1,14)),target :: hmd
        real*8, dimension(dims(1,15)),target :: dptp
        real*8, dimension(dims(1,16)),target :: awnd
        real*8, dimension(dims(1,17)),target :: evap
        real*8, dimension(dims(1,18)),target :: Rad
        real*8, dimension(dims(1,19)),target :: RD
        real*8, dimension(dims(1,20)),target :: RI        
        real*8, dimension(dims(1,21)),target :: tau
        real*8, dimension(dims(1,22)),target :: rrad
        real*8, dimension(dims(1,23)),target :: s1
        real*8, dimension(dims(1,24)),target :: s0
        real*8, dimension(dims(1,25)),target :: snow
        real*8, dimension(dims(1,26)),target :: mnrh
        real*8, dimension(dims(1,27)),target :: mntp
        real*8, dimension(dims(1,28)),target :: mxrh
        integer, dimension(dims(1,29),dims(2,29)),target :: neighbor
        real*8, dimension(dims(1,30)),target :: ref
        real*8, dimension(dims(1,31)),target :: ptr
        real*8, dimension(dims(1,32)),target :: DData
        real*8, dimension(dims(1,33)),target :: pct
        real*8, dimension(dims(1,34)),target :: dat_time
        real*8, dimension(dims(1,35)),target :: cosz
        real*8, dimension(dims(1,35)),target :: decl

        if(dims(2,1)/=0) w%Stations(IS)%id => id
        if(dims(2,2)/=0) w%Stations(IS)%XYElev => XYElev
        if(dims(2,3)/=0) w%Stations(IS)%LatLong => LatLong
        !w%Stations(IS)%Name => Name
        if(dims(1,4)/=0) w%Stations(IS)%gIndices => gIndices
        if(dims(1,5)/=0) w%Stations(IS)%wPct => wPct
        if(dims(1,6)/=0) w%Stations(IS)%dates => dates
        if(dims(1,7)/=0) w%Stations(IS)%datenums => datenums
        if(dims(1,8)/=0) w%Stations(IS)%vratio => vratio
        if(dims(1,9)/=0) w%Stations(IS)%prcp => prcp
        if(dims(1,10)/=0) w%Stations(IS)%tmin => tmin
        if(dims(1,11)/=0) w%Stations(IS)%tmax => tmax
        if(dims(1,12)/=0) w%Stations(IS)%sol => sol
        if(dims(1,13)/=0) w%Stations(IS)%wnd => wnd
        if(dims(1,14)/=0) w%Stations(IS)%hmd => hmd
        if(dims(1,15)/=0) w%Stations(IS)%dptp => dptp
        if(dims(1,16)/=0) w%Stations(IS)%awnd => awnd
        if(dims(1,17)/=0) w%Stations(IS)%evap => evap
        if(dims(1,18)/=0) w%Stations(IS)%Rad => Rad
        if(dims(1,19)/=0) w%Stations(IS)%RD => RD
        if(dims(1,20)/=0) w%Stations(IS)%RI => RI
        if(dims(1,21)/=0) w%Stations(IS)%tau => tau
        if(dims(1,22)/=0) w%Stations(IS)%rrad => rrad
        if(dims(1,23)/=0) w%Stations(IS)%s1 => s1
        if(dims(1,24)/=0) w%Stations(IS)%s0 => s0
        if(dims(1,25)/=0) w%Stations(IS)%snow => snow
        if(dims(1,26)/=0) w%Stations(IS)%mnrh => mnrh
        if(dims(1,27)/=0) w%Stations(IS)%mntp => mntp
        if(dims(1,28)/=0) w%Stations(IS)%mxrh => mxrh
        if(dims(1,29)/=0) w%Stations(IS)%neighbor => neighbor
        if(dims(1,30)/=0) w%Stations(IS)%ref => ref
        if(dims(1,31)/=0) w%Stations(IS)%ptr => ptr
        if(dims(1,32)/=0) w%Stations(IS)%DData => DData
        if(dims(1,33)/=0) w%Stations(IS)%pct => pct
        if(dims(1,34)/=0) w%Stations(IS)%dat_time => dat_time
        if(dims(1,35)/=0) w%Stations(IS)%cosz => cosz
        if(dims(1,35)/=0) w%Stations(IS)%decl => decl

        end subroutine linkStations

!=====================================================================================================
        subroutine bridgeDM(DM, IS, OBJ)
        use matRW
        implicit none
        mwPointer DM
        mwIndex IS
        character(len=*) :: OBJ
        mwPointer, parameter :: nFields = 17
        mwPointer :: i
        character*32 :: Fields(nFields)
        mwpointer :: ptr(nFields), pm
        mwPointer :: m(3, nFields)
        mwIndex :: index = 1
        mwSize number_of_dim
        
        m = 0
        Fields = (/'d','msize','nGho','Map','MapL','Nxy','dim','origin','origin_idx','x', &
       &     'y','Pbound','Cbound','Ebound','BIndex','mask','mask2'/)

        do i = 1, nFields
          pm = mxGetField(DM,index,TRIM(Fields(i)))
          if(pm /= 0) then
            if(mxISEmpty(pm)==0) then
              number_of_dim = mxGetNumberOfDimensions(pm)
              call mxCopyPtrToPtrArray(mxGetDimensions(pm),m(:,i),number_of_dim)
            end if
            ptr(i) = mxGetPr(pm)
          else
            ptr(i) = 0
          end if
        end do

        call linkDM(%val(ptr(1)),%val(ptr(2)),%val(ptr(3)),%val(ptr(4)),%val(ptr(5)), &
       & %val(ptr(6)),%val(ptr(7)),%val(ptr(8)),%val(ptr(9)),%val(ptr(10)),%val(ptr(11)), &
       & %val(ptr(12)),%val(ptr(13)),%val(ptr(14)),%val(ptr(15)),%val(ptr(16)),%val(ptr(17)),m,nFields,IS,OBJ)

        end subroutine bridgeDM

!======================================================================================================================
        subroutine linkDM(d,msize,nGho,Map,MapL,Nxy,dim,origin,origin_idx,x,y,pBound,Cbound, &
       &        Ebound,Bindex,mask,mask2,dims,nFields,IS,OBJ)
          use Vdata
          implicit none
          mwPointer nFields
          mwIndex IS
          character(len=*) :: OBJ
          mwPointer :: dims(3,nFields)
          real*8, dimension(dims(2,1)), target :: d
          real*8, dimension(dims(2,2)), target :: msize
          real*8, dimension(dims(2,3)), target :: nGho
          real*8, dimension(dims(1,4),dims(2,4)),target :: Map
          integer, dimension(dims(1,5),dims(2,5)),target :: MapL
          real*8, dimension(dims(2,6)),target :: Nxy
          real*8, dimension(dims(2,7)),target :: dim
          real*8, dimension(dims(2,8)),target :: origin
          real*8, dimension(dims(2,9)),target :: origin_idx
          real*8, dimension(dims(2,10)),target :: x
          real*8, dimension(dims(2,11)),target :: y
          real*8, dimension(dims(2,12)),target :: Pbound
          real*8, dimension(dims(2,13)),target :: Cbound
          real*8, dimension(dims(2,14)),target :: Ebound
          real*8, dimension(dims(2,15)),target :: Bindex
          logical*1, dimension(dims(1,16),dims(2,16)),target :: mask
          real*8, dimension(dims(1,17),dims(2,17)),target :: mask2

          select case (OBJ)
            case('g')
              if(dims(2,1)/=0) g%DM%d => d
              if(dims(2,2)/=0) g%DM%msize => msize
              if(dims(2,3)/=0) g%DM%nGho => nGho
              if(dims(1,4)/=0 .and. dims(2,4)/=0) g%DM%Map => Map
              if(dims(1,5)/=0 .and. dims(2,5)/=0) g%DM%MapL => MapL
              if(dims(2,6)/=0) g%DM%Nxy => Nxy
              if(dims(2,7)/=0) g%DM%dim => dim
              if(dims(2,8)/=0) g%DM%origin => origin
              if(dims(2,9)/=0) g%DM%origin_idx => origin_idx
              if(dims(2,10)/=0) g%DM%x => x
              if(dims(2,11)/=0) g%DM%y => y
              if(dims(2,12)/=0) g%DM%Pbound => Pbound
              if(dims(2,13)/=0) g%DM%Cbound => Cbound
              if(dims(2,14)/=0) g%DM%Ebound => Ebound
              if(dims(2,15)/=0) g%DM%Bindex => Bindex
            case('r')
              if(dims(2,1)/=0) r(IS)%DM%d => d
              if(dims(2,2)/=0) r(IS)%DM%msize => msize
              if(dims(2,3)/=0) r(IS)%DM%nGho => nGho
              if(dims(1,4)/=0 .and. dims(2,4)/=0) r(IS)%DM%Map => Map
              if(dims(1,5)/=0 .and. dims(2,5)/=0) r(IS)%DM%MapL => MapL
              if(dims(2,6)/=0) r(IS)%DM%Nxy => Nxy
              if(dims(2,7)/=0) r(IS)%DM%dim => dim
              if(dims(2,8)/=0) r(IS)%DM%origin => origin
              if(dims(2,9)/=0) r(IS)%DM%origin_idx => origin_idx
              if(dims(2,10)/=0) r(IS)%DM%x => x
              if(dims(2,11)/=0) r(IS)%DM%y => y
              if(dims(2,12)/=0) r(IS)%DM%Pbound => Pbound
              if(dims(2,13)/=0) r(IS)%DM%Cbound => Cbound
              if(dims(2,14)/=0) r(IS)%DM%Ebound => Ebound
              if(dims(2,15)/=0) r(IS)%DM%Bindex => Bindex
            case('w')
              if(dims(1,16)/=0 .and. dims(2,16)/=0) w%DM%mask => mask
              if(dims(1,17)/=0 .and. dims(2,17)/=0) w%DM%mask2 => mask2
          end select

        end subroutine linkDM


!========================================================================================================
        SUBROUTINE bridgeVDZ(VDZM,IS)
        use matRW
        IMPLICIT NONE
        !subroutine mexFunction(nlhs, plhs, nrhs, prhs)
        mwPointer IS
        mwPointer m(3),np(3),nD(3),nr(3),nrep(3)
        mwpointer :: pm, pm2,VDZM,VParm,GWL
        mwpointer :: NZC,DZ,DDZ,h,SRC,JType,K,C,WTL,ET,Cond
        mwpointer :: THE,DF,Dperc,t,dt,Pratio,EB,E,dtP,R
        mwPointer :: THER,THES,KS,Ko,nn,Lambda,alpha,wp,fc,MAPP,REP,STATE,SMR,KSSAVE
        mwPointer :: Sand,Clay,OM,ICE,soic,gti
        mwIndex :: index = 1
        mwSize :: number_of_dim
        ! linkage subroutine from matlab to fortran struct
        pm = mxGetField(VDZM,index,'h')
        number_of_dim = mxGetNumberOfDimensions(pm)
        !CALL mxCopyPtrToInteger4(mxGetDimensions(pm),m,number_of_dim)
        CALL mxCopyPtrToPtrArray(mxGetDimensions(pm),m,number_of_dim)
        ! VDZ object!
        h = mxGetPr(pm)
        EB = mxGetPr(mxGetField(VDZM,index,'EB'))
        E = mxGetPr(mxGetField(VDZM,index,'E'))
        NZC = mxGetPr(mxGetField(VDZM,index,'NZC'))
        DZ = mxGetPr(mxGetField(VDZM,index,'DZ'))
        DDZ = mxGetPr(mxGetField(VDZM,index,'DDZ'))
        !h = mxGetPr(mxGetField(VDZM,index,'h'))
        SRC = mxGetPr(mxGetField(VDZM,index,'SRC'))
        JType = mxGetPr(mxGetField(VDZM,index,'JType'))
        WTL = mxGetPr(mxGetField(VDZM,index,'WTL'))
        K = mxGetPr(mxGetField(VDZM,index,'K'))
        C = mxGetPr(mxGetField(VDZM,index,'C'))
        THE = mxGetPr(mxGetField(VDZM,index,'THE'))
        DF = mxGetPr(mxGetField(VDZM,index,'DF'))
        Dperc = mxGetPr(mxGetField(VDZM,index,'Dperc'))
        t = mxGetPr(mxGetField(VDZM,index,'t'))
        dt = mxGetPr(mxGetField(VDZM,index,'dt'))
        dtP = mxGetPr(mxGetField(VDZM,index,'dtP'))
        PRatio = mxGetPr(mxGetField(VDZM,index,'PRatio'))
        !R = mxGetPr(mxGetField(VDZM,index,'R'))
        pm = mxGetField(VDZM,index,'R')
        R = mxGetPr(pm)
        number_of_dim = mxGetNumberOfDimensions(pm)
        !CALL mxCopyPtrToInteger4(mxGetDimensions(pm),nr,number_of_dim)
        call mxCopyPtrToPtrArray(mxGetDimensions(pm),nr,number_of_dim)

        THER = mxGetPr(mxGetField(VDZM,index,'THER'))
        THES = mxGetPr(mxGetField(VDZM,index,'THES'))
        KS = mxGetPr(mxGetField(VDZM,index,'KS'))
        Ko = mxGetPr(mxGetField(VDZM,index,'Ko'))
        nn = mxGetPr(mxGetField(VDZM,index,'N'))
        Lambda = mxGetPr(mxGetField(VDZM,index,'LAMBDA'))
        alpha = mxGetPr(mxGetField(VDZM,index,'ALPHA'))
        fc = mxGetPr(mxGetField(VDZM,index,'FC'))
        wp = mxGetPr(mxGetField(VDZM,index,'WP'))
        MAPP = mxGetPr(mxGetField(VDZM,index,'Mapp'))
        ET = mxGetPr(mxGetField(VDZM,index,'ET'))
        GWL = mxGetPr(mxGetField(VDZM,index,'GWL'))
        STATE = mxGetPr(mxGetField(VDZM,index,'STATE'))
        SMR = mxGetPr(mxGetField(VDZM,index,'SMR'))
        KSSAVE = mxGetPr(mxGetField(VDZM,index,'KSSAVE'))
        SAND = mxGetPr(mxGetField(VDZM,index,'Sand'))
        Clay = mxGetPr(mxGetField(VDZM,index,'Clay'))
        OM = mxGetPr(mxGetField(VDZM,index,'OM'))
        ICE = mxGetPr(mxGetField(VDZM,index,'ICE'))
        soic = mxGetPr(mxGetField(VDZM,index,'soic'))
        gti = mxGetPr(mxGetField(VDZM,index,'gti'))

        pm = mxGetField(VDZM,index,'Cond')
        Cond = mxGetPr(pm)
        number_of_dim = mxGetNumberOfDimensions(pm)
        !CALL mxCopyPtrToInteger4(mxGetDimensions(pm),nD,number_of_dim)
        call mxCopyPtrToPtrArray(mxGetDimensions(pm),nD,number_of_dim)

        pm = mxGetField(VDZM,index,'VParm')
        VParm = mxGetPr(pm)
        number_of_dim = mxGetNumberOfDimensions(pm)
        !CALL mxCopyPtrToInteger4(mxGetDimensions(pm),np,number_of_dim)
        CALL mxCopyPtrToPtrArray(mxGetDimensions(pm),np,number_of_dim)

        pm = mxGetField(VDZM,index,'REP')
        REP = mxGetPr(pm)
        number_of_dim = mxGetNumberOfDimensions(pm)
        CALL mxCopyPtrToPtrArray(mxGetDimensions(pm),nrep,number_of_dim)

        pm = mxGetField(VDZM,index,'GA')
        call bridgeGA(pm,1)

        pm = mxGetField(VDZM,index,'th')
        call bridgeTH(pm,1)

        CALL linkVDZ(IS,m(1),m(2),m(3),np(1)*np(2),%val(NZC),%val(DZ)                  &
     &   ,%val(DDZ),%val(h),%val(SRC),%val(JType),%val(WTL),%val(K)          &
     &   ,%val(C),%val(THE),%val(DF),%val(Dperc),%val(t),%val(dt)          &
     &   ,%val(dtP),%val(EB),%val(E),%val(R),%val(THER),%val(THES)         &
     &   ,%val(KS),%val(Ko),%val(nn),%val(Lambda),%val(alpha),%val(ET)                      &
     &   ,%val(fc),%val(wp),%val(MAPP),%val(Cond),%val(PRatio)             &
     &   ,%val(VParm),nD(3),%val(GWL),%val(STATE),%val(SMR),%val(KSSAVE),nr(3),%val(REP),nrep(3) &
     &   ,%val(Sand),%val(Clay),%val(OM),%val(ICE),%val(soic),%val(gti))

        END SUBROUTINE bridgeVDZ

!================================================================================================================
        SUBROUTINE linkVDZ(IS,m,n,nz,np,NZC,DZ,DDZ,h,SRC,                     &
     &            JType,WTL,K,CI,THE,DF,Dperc,t,dt,dtP,EB,E,Rm,              &
     &   THER,THES,KS,Ko,nn,Lambda,alpha,ET,fc,wp,mapp,cond,PRatio,            &
     &   VParm,nD,GWL,STATE,SMR,KSSAVE,nr,REP,nrep,Sand,Clay,OM,ICE,soic,gti)
        USE Vdata
        mwPointer IS,n,nz,ntype
        mwPointer m,np,nD,nr,nrep
        REAL*8,DIMENSION(m,n,nz),target::h,THE,K,CI,SRC,DZ,DDZ,fc,wp
        REAL*8,DIMENSION(m,n,nz),target::Sand,Clay,OM,ICE
        integer,DIMENSION(m,n),target::JType,NZC,WTL,soic
        REAL*8,DIMENSION(m,n),target::DF,DPerc,dt,Mapp,PRatio,GWL,SMR,gti
        REAL*8,DIMENSION(m,n,nz),target::EB,E,ET ! ET is evapotranspiration
        REAL*8,DIMENSION(m,n,nz),target::THER,THES,KS,Ko,nn,Lambda,alpha,KSSAVE
        REAL*8,DIMENSION(m,n,nD),target::Cond
        REAL*8,DIMENSION(np),Target::VParm
        REAL*8,DIMENSION(m,n,nr),Target::Rm
        REAL*8,DIMENSION(m,n,nrep),Target::REP
        REAL*8,Target:: t,dtP(3),STATE
        g%VDZ%NZC => NZC
        g%VDZ%m(1)=m; g%VDZ%m(2)=n; g%VDZ%m(3)=nz;
        g%VDZ%h => h
        g%VDZ%THE => THE
        g%VDZ%K => K
        g%VDZ%CI => CI
        g%VDZ%SRC => SRC
        g%VDZ%DZ => DZ
        g%VDZ%DDZ => DDZ
        g%VDZ%EB => EB
        g%VDZ%E => E
        g%VDZ%ET => ET
        g%VDZ%DF => DF
        g%VDZ%DPerc => DPerc
        g%VDZ%JType => JType
        g%VDZ%WTL => WTL
        g%VDZ%t => t
        g%VDZ%dt => dt
        g%VDZ%R => Rm
        g%VDZ%dtP => dtP
        g%VDZ%STATE => STATE
        g%VDZ%THER=> THER
        g%VDZ%THES=> THES
        g%VDZ%KS=> KS
        g%VDZ%Ko=> Ko
        g%VDZ%N=> nn
        g%VDZ%Lambda=> Lambda
        g%VDZ%alpha=> alpha
        g%VDZ%fc=> fc
        g%VDZ%wp=> wp
        g%VDZ%MAPP=> MAPP
        g%VDZ%Cond => Cond
        g%VDZ%PRatio => PRatio
        g%VDZ%VParm => VParm
        if (size(VParm)<12) then
          call CH_ASSERT('Not enough VParm parameters!!!')
          ! otherwise the program actually won't be able to detect it, and use some weird number!!
        endif
        g%VDZ%GWL => GWL
        g%VDZ%REP => REP
        g%VDZ%SMR => SMR
        g%VDZ%KSSAVE => KSSAVE
        g%VDZ%Sand => Sand
        g%VDZ%Clay => Clay
        g%VDZ%OM => OM
        g%VDZ%ICE => ICE
        g%VDZ%gti => gti
        g%VDZ%soic => soic
        END SUBROUTINE linkVDZ


!==========================================================================================================
      SUBROUTINE bridgeSoils(SOILSM)
      use matRW
      mwpointer :: KS,THER,THES,alpha,lambda,nn,pm,SOILSM,FC
      mwpointer :: DP,DPN,DPUL,DPUU,STATE
      mwPointer mc(3),n(3)
      mwIndex :: index = 1
      mwSize :: number_of_dim

      pm = mxGetField(SOILSM,index,'KS')
      number_of_dim = mxGetNumberOfDimensions(pm)
      !CALL mxCopyPtrToInteger4(mxGetDimensions(pm),mc,3)
      call mxCopyPtrToPtrArray(mxGetDimensions(pm),mc,number_of_dim)
      KS = mxGetPr(mxGetField(SOILSM,index,'KS'))
      THER = mxGetPr(mxGetField(SOILSM,index,'THER'))
      FC = mxGetPr(mxGetField(SOILSM,index,'FC'))
      THES = mxGetPr(mxGetField(SOILSM,index,'THES'))
      alpha = mxGetPr(mxGetField(SOILSM,index,'alpha'))
      lambda = mxGetPr(mxGetField(SOILSM,index,'lambda'))
      nn = mxGetPr(mxGetField(SOILSM,index,'n'))
      pm = mxGetField(SOILSM,index,'DP')
      DP = mxGetPr(pm)
      number_of_dim = mxGetNumberOfDimensions(pm)
      !CALL mxCopyPtrToInteger4(mxGetDimensions(pm),n,3)
      CALL mxCopyPtrToPtrArray(mxGetDimensions(pm),n,number_of_dim)
      DPN = mxGetPr(mxGetField(SOILSM,index,'DPN'))
      DPUL = mxGetPr(mxGetField(SOILSM,index,'DPUL'))
      DPUU = mxGetPr(mxGetField(SOILSM,index,'DPUU'))
      STATE = mxGetPr(mxGetField(SOILSM,index,'STATE'))

      CALL linkSoils(mc(1),mc(2),n(1),%val(KS),%val(THER),%val(FC),                          &
     &     %val(THES),%val(alpha),%val(lambda),%val(nn),%val(DP)                        &
     &     ,%val(DPN),%val(DPUL),%val(DPUU),%val(STATE))

      END SUBROUTINE bridgeSoils

!=====================================================================================================
      SUBROUTINE linkSoils(m,n,nD,KS,THER,FC,THES,alpha,lambda,nn,                  &
     &          DP,DPN,DPUL,DPUU,STATE)
      USE Vdata
      mwPointer J
      mwPointer :: m,n,nD
      REAL*8,DIMENSION(m,n),target::KS,THER,THES,alpha,lambda,nn,FC
      REAL*8,DIMENSION(nD,n),target::DP
      REAL*8,DIMENSION(n),target::DPUL,DPUU
      integer,DIMENSION(n),target::DPN
      real*8,target::STATE
      g%SOILS%nl = m
      g%SOILS%nt = n
      g%SOILS%KS => KS
      g%SOILS%FC => FC
      g%SOILS%THER => THER
      g%SOILS%THES => THES
      g%SOILS%alpha => alpha
      g%SOILS%lambda => lambda
      g%SOILS%nn => nn
      g%SOILS%STATE => STATE
      g%SOILS%DP => DP
      g%SOILS%DPUL => DPUL
      g%SOILS%DPUU => DPUU
      g%SOILS%DPN => DPN
      END SUBROUTINE linkSoils

!========================================================================================
      SUBROUTINE bridgeGW(GWM,index)
      use Vdata
      use matRW
      IMPLICIT none
      mwIndex index
      mwPointer :: nG(2), m(2),nn(3)
      mwpointer :: GWM,MAPP
      mwpointer :: h,EB,ES,E,D,K,Wm,W2,ST,T,DR,DLDR
      mwpointer :: hNEW,DZ,Seep,Qperc,WTL
      mwpointer :: A1B,A2,A3,A4,A5,pm,dd,dt,tt
      mwpointer :: topH,topK,topD,botH,botK,botD,topBound
      mwpointer :: STO,NLC,STATE
      mwSize number_of_dim, i

      pm = mxGetField(GWM,index,'Cond')
      number_of_dim = mxGetNumberOfDimensions(pm)
      call mxCopyPtrToPtrArray(mxGetDimensions(pm),m,number_of_dim)
      allocate(g%GW(index)%Cond(m(1)*m(2)))
      do i = 1, m(1)*m(2)
        call bridgeCond(pm,index,i)
      enddo

      pm = mxGetField(GWM,index,'h')
      number_of_dim = mxGetNumberOfDimensions(pm)
      !CALL mxCopyPtrToInteger4(mxGetDimensions(pm),m,2)
      call mxCopyPtrToPtrArray(mxGetDimensions(pm),m,number_of_dim)
      h = mxGetPr(pm)
      nG = (/1,1/)
      EB = mxGetPr(mxGetField(GWM,index,'EB'))
      MAPP = mxGetPr(mxGetField(GWM,index,'MAPP'))
      ES = mxGetPr(mxGetField(GWM,index,'ES'))
      E = mxGetPr(mxGetField(GWM,index,'E'))
      D = mxGetPr(mxGetField(GWM,index,'D'))
      DR = mxGetPr(mxGetField(GWM,index,'DR'))
      DLDR = mxGetPr(mxGetField(GWM,index,'DLDR'))
      K = mxGetPr(mxGetField(GWM,index,'K'))
      Wm = mxGetPr(mxGetField(GWM,index,'W'))
      W2 = mxGetPr(mxGetField(GWM,index,'W2'))
      ST = mxGetPr(mxGetField(GWM,index,'ST'))
      T = mxGetPr(mxGetField(GWM,index,'T'))
      hNew = mxGetPr(mxGetField(GWM,index,'hNew'))
      DZ = mxGetPr(mxGetField(GWM,index,'DZ')) ! Thickness of aquifer
      Seep = mxGetPr(mxGetField(GWM,index,'Seep'))
      Qperc = mxGetPr(mxGetField(GWM,index,'Qperc'))
      !A1B = mxGetPr(mxGetField(GWM,1,'A1B'))
      !A2 = mxGetPr(mxGetField(GWM,1,'A2'))
      !A3 = mxGetPr(mxGetField(GWM,1,'A3'))
      !A4 = mxGetPr(mxGetField(GWM,1,'A4'))
      !A5 = mxGetPr(mxGetField(GWM,1,'A5'))
      !topH = mxGetPr(mxGetField(GWM,1,'topH'))
      !topK = mxGetPr(mxGetField(GWM,1,'topK'))
      !topD = mxGetPr(mxGetField(GWM,1,'topD'))
      !botH = mxGetPr(mxGetField(GWM,1,'botH'))
      botK = mxGetPr(mxGetField(GWM,index,'botK'))
      botD = mxGetPr(mxGetField(GWM,index,'botD'))
      dd = mxGetPr(mxGetField(GWM,index,'dd'))
      dt = mxGetPr(mxGetField(GWM,index,'dt'))
      tt = mxGetPr(mxGetField(GWM,index,'tt'))
      STO = mxGetPr(mxGetField(GWM,index,'STO'))
      STATE = mxGetPr(mxGetField(GWM,index,'STATE'))
      topBound = mxGetPr(mxGetField(GWM,index,'topBound')) 
      pm = mxGetField(GWM,index,'NLC')
      number_of_dim = mxGetNumberOfDimensions(pm)
      !CALL mxCopyPtrToInteger4(mxGetDimensions(pm),m,2)
      call mxCopyPtrToPtrArray(mxGetDimensions(pm),nn,number_of_dim)
      NLC = mxGetPr(pm)

      CALL linkGW(index,m(1),m(2),nn,%val(dd),%val(h),%val(EB),%val(ES),%val(E)        &
     &     ,%val(D),%val(K),%val(Wm),%val(W2),%val(ST),%val(T)              &
     &     ,%val(hNew),%val(DZ),%val(Seep) ,%val(DLDR),%val(DR)                     &
     &     ,%val(dt),%val(Qperc),%val(A1B),%val(A2),%val(A3)               &
     &     ,%val(A4),%val(A5),%val(tt),%val(MAPP),%val(STO),                  &
     &     %val(NLC),%val(botK),%val(botD),%val(STATE),%val(topBound),nn) !,%val(topH),%val(topK),%val(topD)               &
      !&     ,%val(botH),%val(botK),%val(botD))


      END SUBROUTINE bridgeGW

!=======================================================================================================
      SUBROUTINE linkGW(J,ny,nx,nn,dd,h,EB,ES,E,D,K,W1,W2,ST,T,hNew,DZ,Seep,               &
     &        DLDR,DR,dt,Qperc,A1B,A2,A3,A4,A5,tt,MAPP,STO,NLC,botK,botD,STATE,topBound) !,topH,topK,topD,botH,botK,botD)
      USE Vdata
      mwIndex J
      mwPointer nc
      mwPointer :: ny,nx,nn(3)
      REAL*8,DIMENSION(ny,nx),Target::EB,ES,E,D,K,W1,W2,ST,T,hNew,DZ             &
     &      ,Seep,DLDR,DR,Qperc,A1B,A2,A3,A4,A5,h,STO,botK,botD,topBound !,topH,topK,topD,botH,botK,botD
      real*8,dimension(ny,nx),target:: MAPP
      !REAL*8,Target::NLC(ny,nx,nc)
      REAL*8,Target:: dd(2),dt,tt,STATE
      real*8,dimension(nn(1),nn(2),nn(3)),target:: NLC
      g%GW(J)%ny = ny
      g%GW(J)%nx = nx
      g%GW(J)%EB => EB
      g%GW(J)%ES => ES
      g%GW(J)%E => E
      g%GW(J)%D => D
      g%GW(J)%K => K
      g%GW(J)%W => W1
      g%GW(J)%W2 => W2
      g%GW(J)%ST => ST
      g%GW(J)%T => T
      g%GW(J)%hNew => hNew
      g%GW(J)%h => h
      g%GW(J)%DZ => DZ
      g%GW(J)%Seep => Seep
      g%GW(J)%DR => DR
      g%GW(J)%DLDR => DLDR
      g%GW(J)%Qperc => Qperc
      g%GW(J)%A1B => A1B
      g%GW(J)%A2 => A2
      g%GW(J)%A3 => A3
      g%GW(J)%A4 => A4
      g%GW(J)%A5 => A5
      g%GW(J)%dd => dd
      g%GW(J)%dt => dt
      g%GW(J)%tt => tt
      g%GW(J)%MAPP => MAPP
      g%GW(J)%STO => STO
      g%GW(J)%NLC => NLC
      g%GW(J)%nc = nc
      !g%GW(J)%topH => topH
      !g%GW(J)%topD => topD
      !g%GW(J)%topK => topK
      !g%GW(J)%botH => botH
      g%GW(J)%botK => botK
      g%GW(J)%botD => botD
      g%GW(J)%STATE => STATE
     g%GW(j)%topBound=>topBound 
      END SUBROUTINE linkGW

!==================================================================================================
      SUBROUTINE bridgeVeg(Veg,IS)
      use matRW
      IMPLICIT NONE
      !mwPointer :: mxGetField, mxGetDimensions, mxGetM, mxGetN,mxGetPr
      mwPointer IS
      mwPointer :: nG(2), m(3)
      mwpointer :: Veg,Frac,LAI,RPT,Evap,BSoil,RDF,RDFN,EDFN
      mwpointer :: Delta,hI,dt,Kc,hc,Cpc,CS,EvapCS,DEPC,h,Tp,ET,RPET
      mwpointer :: smap,t,Root,pm,PET,Gprcp,EvapG,SNOW,NewS,LUidx
      mwpointer :: ALBEDO,ASP,SWE,Salb,SNOCOV,SNOM,SNOCoef,Imp,VDB,day
      mwPointer :: index = 1
      mwSize number_of_dim

      pm = mxGetField(Veg,index,'RPT')
      number_of_dim = mxGetNumberOfDimensions(pm)
      !CALL mxCopyPtrToInteger4(mxGetDimensions(pm),m,3)
      CALL mxCopyPtrToPtrArray(mxGetDimensions(pm),m,number_of_dim)

      Frac = mxGetPr(mxGetField(Veg,index,'Frac'))
      SMAP = mxGetPr(mxGetField(Veg,index,'SMAP'))
      LAI = mxGetPr(mxGetField(Veg,index,'LAI'))
      RPT = mxGetPr(mxGetField(Veg,index,'RPT'))
      Evap = mxGetPr(mxGetField(Veg,index,'Evap'))
      EvapCS = mxGetPr(mxGetField(Veg,index,'EvapCS'))
      BSoil = mxGetPr(mxGetField(Veg,index,'BSoil'))
      Root = mxGetPr(mxGetField(Veg,index,'Root'))
      RDFN = mxGetPr(mxGetField(Veg,index,'RDFN'))
      EDFN = mxGetPr(mxGetField(Veg,index,'EDFN'))
      Delta = mxGetPr(mxGetField(Veg,index,'Delta'))
      hI = mxGetPr(mxGetField(Veg,index,'hI'))
      h = mxGetPr(mxGetField(Veg,index,'h'))
      Imp = mxGetPr(mxGetField(Veg,index,'Imp'))
      t = mxGetPr(mxGetField(Veg,index,'t'))
      dt = mxGetPr(mxGetField(Veg,index,'dt'))
      Kc = mxGetPr(mxGetField(Veg,index,'Kc'))
      hc = mxGetPr(mxGetField(Veg,index,'hc'))
      Cpc = mxGetPr(mxGetField(Veg,index,'Cpc'))
      CS = mxGetPr(mxGetField(Veg,index,'CS'))
      EvapCS = mxGetPr(mxGetField(Veg,index,'EvapCS'))
      h = mxGetPr(mxGetField(Veg,index,'h'))
      Evap = mxGetPr(mxGetField(Veg,index,'Evap'))
      Tp = mxGetPr(mxGetField(Veg,index,'Tp'))
      ET = mxGetPr(mxGetField(Veg,index,'ET'))
      RPET = mxGetPr(mxGetField(Veg,index,'RPET'))
      PET = mxGetPr(mxGetField(Veg,index,'PET'))
      SMAP = mxGetPr(mxGetField(Veg,index,'SMAP'))
      Albedo = mxGetPr(mxGetField(Veg,index,'albedo'))
      ASP = mxGetPr(mxGetField(Veg,index,'ASP'))
      SWE = mxGetPr(mxGetField(Veg,index,'SWE'))
      Salb = mxGetPr(mxGetField(Veg,index,'Salb'))
      SNOCOV = mxGetPr(mxGetField(Veg,index,'SNOCOV'))
      SNOM = mxGetPr(mxGetField(Veg,index,'SNOM'))
      Gprcp = mxGetPr(mxGetField(Veg,index,'Gprcp'))
      EvapG = mxGetPr(mxGetField(Veg,index,'EvapG'))
      SNOCoef = mxGetPr(mxGetField(Veg,index,'SNOCoef'))
      SNOW = mxGetPr(mxGetField(Veg,index,'SNOW'))
      NewS = mxGetPr(mxGetField(Veg,index,'NewS'))
      LUidx = mxGetPr(mxGetField(Veg,index,'LUidx'))

      VDB = mxGetField(Veg,index,'VDB')
      call bridgeVDB(VDB)
      
      day = mxGetField(Veg,index,'day')
      call bridgeVegDay(day)

      !ALBEDO,ASP,SWE,Salb,SNOCOV,SNOM,SNOCoef
      CALL linkVeg(IS,m,%val(Frac),%val(LAI),%val(RPT),%val(Evap)                   &
     &     ,%val(EvapCS),%val(BSoil),%val(Root),%val(RDFN),%val(EDFN)                     &
     &     ,%val(Delta),%val(hI),%val(t),%val(dt),%val(Kc),%val(hc),%val(Cpc)             &
     &     ,%val(CS),%val(Tp),%val(ET),%val(SMap),%val(RPET)                         &
     &     ,%val(PET),%val(Albedo),%val(ASP),%val(SWE),%val(Gprcp)                       &
     &     ,%val(EvapG),%val(Salb),%val(SNOCOV),%val(SNOM)                       &
     &     ,%val(SNOCoef),%val(SNOW),%val(NEWS),%val(h),%val(Imp),%val(LUidx))

      END SUBROUTINE bridgeVeg


!=========================================================================================================
      subroutine bridgeVegDay(day)
      use matRW
      implicit none
      mwPointer :: day
      mwPointer :: LAI, Root, hc, Kc
      mwPointer :: index = 1, m(3)
      mwSize number_of_dim      

      LAI = mxGetField(day,index,'LAI')
      number_of_dim = mxGetNumberOfDimensions(LAI)
      CALL mxCopyPtrToPtrArray(mxGetDimensions(LAI),m,number_of_dim)
      LAI = mxGetPr(LAI)
      
      Root = mxGetPr(mxGetField(day,index,'Root'))
      hc = mxGetPr(mxGetField(day,index,'hc'))
      Kc = mxGetPr(mxGetField(day,index,'Kc'))
      
      call linkVegDay(m, %val(LAI),%val(Root),%val(hc),%val(Kc))
      
      end subroutine bridgeVegDay
      
!=========================================================================================================
      subroutine linkVegDay(m, LAI, Root, hc, Kc)
      use Vdata
      implicit none
      mwPointer :: m(3)
      real*8, dimension(m(1)),target :: LAI, Root, hc, Kc
      
      g%Veg%day%LAI => LAI
      g%Veg%day%Root => Root
      g%Veg%day%hc => hc
      g%Veg%day%Kc => Kc
      end subroutine linkVegDay      
!=========================================================================================================
      SUBROUTINE linkVeg(J,m,Frac,LAI,RPT,Evap,EvapCS,BSoil,Root,             &
     &      RDFN,EDFN,Delta,hI,t,dt,Kc,hc,Cpc,CS,Tp,ET,SMAP,RPET,PET,            &
     &      ALBEDO,ASP,SWE,GPrcp,EvapG,Salb,SNOCOV,SNOM,SNOCoef,SNOW,          &
     &      NEWS,h,Imp,LUidx)
      USE Vdata
      IMPLICIT NONE
      mwPointer J
      mwPointer m(3)
      REAL*8,DIMENSION(m(1),m(2)),Target::EDF,Evap,                    &
     &        Tp,Cpc,CS,ET,RPET,EvapCS,PET,h,Imp,hI
      integer,DIMENSION(m(1),m(2),m(3)),Target::RPT,RDFN
      integer,DIMENSION(m(1),m(2)),Target::EDFN,SMAP
      REAL*8,DIMENSION(m(1),m(2),m(3)),Target::Delta,LAI,Root,              &
     &        BSoil,Frac,Kc,hc
      REAL*8,DIMENSION(m(1),m(2)),Target::ALBEDO,ASP,SWE,Salb,SNOCOV
      REAL*8,DIMENSION(m(1),m(2)),Target::SNOM,GPrcp,EvapG,SNOW,NewS
      REAL*8,DIMENSION(2),Target::SNOCoef
      REAL*8,DIMENSION(5),Target:: LUidx


      REAL*8,target:: dt,t

      g%VEG%ny = m(1)
      g%VEG%nx = m(2)
      if (m(3)<1 .or. m(3)>20) m(3)=1
      g%VEG%nt = m(3)
      g%VEG%Frac => Frac
      g%VEG%LAI => LAI
      g%VEG%RPT => RPT
      g%VEG%Frac => Frac
      g%VEG%BSoil => BSoil
      g%VEG%Root => Root
      g%VEG%EDFN => EDFN
      g%VEG%RDFN => RDFN
      g%VEG%Evap => Evap
      g%VEG%EvapCS => EvapCS
      g%VEG%Tp => Tp
      g%VEG%Cpc => Cpc
      g%VEG%CS => CS
      g%VEG%ET => ET
      g%VEG%PET => PET
      g%VEG%RPET => RPET
      g%VEG%SMAP => SMAP
      g%VEG%Delta => Delta
      g%VEG%t => t
      g%VEG%dt => dt
      g%VEG%Kc => Kc
      g%VEG%hc => hc
      g%VEG%Albedo => Albedo
      g%VEG%ASP => ASP
      g%VEG%SWE => SWE
      g%VEG%Salb => Salb
      g%VEG%SNOCOV => SNOCOV
      g%VEG%SNOM => SNOM
      g%VEG%SNOCoef => SNOCoef
      g%VEG%GPrcp => GPrcp
      g%VEG%EvapG => EvapG
      g%VEG%SNOW => SNOW
      g%VEG%NewS => NewS
      g%VEG%h => h
      g%VEG%hI => hI
      g%VEG%Imp => Imp
      g%VEG%LUidx => LUidx
      !ALBEDO,ASP,SWE,Salb,SNOCOV,SNOM,SNOCoef
      END SUBROUTINE linkVeg

!====================================================================================================
      SUBROUTINE bridgeVDB(VDB)
        use matRW
        implicit none
        mwPointer VDB
        mwPointer, parameter :: nFields = 10
        mwPointer :: i
        character*32 :: Fields(nFields)
        mwpointer :: ptr(nFields), pm
        mwPointer :: m(3, nFields)
        mwIndex :: index = 1
        mwSize number_of_dim
        
        m = 0
        Fields = (/ 'LAICoeff','hcCoeff','KcCoeff','RootCoeff','GrowthT','ho','LAI','Kc','Root','hc'/)

        do i = 1, nFields
          pm = mxGetField(VDB,index,TRIM(Fields(i)))
          if(pm /= 0) then
            if(mxISEmpty(pm)==0) then
              number_of_dim = mxGetNumberOfDimensions(pm)
              call mxCopyPtrToPtrArray(mxGetDimensions(pm),m(:,i),number_of_dim)
            endif
            ptr(i) = mxGetPr(pm)
          else
            ptr(i) = 0
          endif
        enddo

        call linkVDB(%val(ptr(1)),%val(ptr(2)),%val(ptr(3)),%val(ptr(4)),%val(ptr(5)), &
       & %val(ptr(6)),%val(ptr(7)),%val(ptr(8)),%val(ptr(9)),%val(ptr(10)),m,nFields)

      END SUBROUTINE bridgeVDB

!====================================================================================================
        subroutine linkVDB(LAICoeff,hcCoeff,KcCoeff,RootCoeff,GrowthT,ho,LAI,Kc,Root,hc,dims,nFields)
        use VData
        implicit none
        mwPointer nFields
        mwPointer dims(3,nFields)
        real*8, dimension(dims(1,1),dims(2,1)), target :: LAICoeff
        real*8, dimension(dims(1,2),dims(2,2)), target :: hcCoeff
        real*8, dimension(dims(1,3),dims(2,3)), target :: KcCoeff
        real*8, dimension(dims(1,4),dims(2,4)), target :: RootCoeff
        real*8, dimension(dims(1,5),dims(2,5)), target :: GrowthT
        real*8, dimension(dims(1,6),dims(2,6)), target :: ho
        real*8, dimension(dims(1,7),dims(2,7)), target :: LAI
        real*8, dimension(dims(1,8),dims(2,8)), target :: Kc
        real*8, dimension(dims(1,9),dims(2,9)), target :: Root
        real*8, dimension(dims(1,10),dims(2,10)), target :: hc

        if(dims(1,1)/=0 .and. dims(2,1)/=0) g%Veg%VDB%LAICoeff => LAICoeff
        if(dims(1,2)/=0 .and. dims(2,2)/=0) g%Veg%VDB%hcCoeff => hcCoeff
        if(dims(1,3)/=0 .and. dims(2,3)/=0) g%Veg%VDB%KcCoeff => KcCoeff
        if(dims(1,4)/=0 .and. dims(2,4)/=0) g%Veg%VDB%RootCoeff => RootCoeff
        if(dims(1,5)/=0 .and. dims(2,5)/=0) g%Veg%VDB%GrowthT => GrowthT
        if(dims(1,6)/=0 .and. dims(2,6)/=0) g%Veg%VDB%ho => ho
        if(dims(1,7)/=0 .and. dims(2,7)/=0) g%Veg%VDB%LAI => LAI
        if(dims(1,8)/=0 .and. dims(2,8)/=0) g%Veg%VDB%Kc => Kc
        if(dims(1,9)/=0 .and. dims(2,9)/=0) g%Veg%VDB%Root => Root
        if(dims(1,10)/=0 .and. dims(2,10)/=0) g%Veg%VDB%hc => hc

        end subroutine linkVDB


!====================================================================================================
      SUBROUTINE bridgeGA(GA,IS)
      use matRW
      IMPLICIT NONE
      !mwPointer :: mxGetField, mxGetDimensions, mxGetM, mxGetN,mxGetPr
      mwPointer NM,J
      PARAMETER (NM=16)
      mwPointer IS, i, nG(2), m(3)
      CHARACTER*100 :: F(NM)
      REAL*8 y
      mwpointer :: P(NM),GA,pm
      mwIndex :: index = 1
      mwSize :: number_of_dim

      F = (/ 'm','nspan','KBar','WFP','FM','Code','THEW','SW','MM','h',                      &
     &     'n','l','dn','dl','hspan','GCode'/)
      DO i = 1,NM
        P(i) = mxGetPr(mxGetField(GA,index,TRIM(F(i))))
      ENDDO

      pm = mxGetField(GA,index,TRIM(F(3)))
      number_of_dim = mxGetNumberOfDimensions(pm)
      !CALL mxCopyPtrToInteger4(mxGetDimensions(pm),m2,3) ! This is the WFP size
      call mxCopyPtrToPtrArray(mxGetDimensions(pm),m,number_of_dim)

      CALL linkGA(m,%val(P(1)),%val(P(2)),%val(P(3)),%val(P(4))                 &
     &  ,%val(P(5)),%val(P(6)),%val(P(7)),%val(P(8)),%val(P(9))                  &
     &  ,%val(P(10)),%val(P(11)),%val(P(12)),%val(P(13)),%val(P(14))             &
     &  ,%val(P(15)),%val(P(16)))

      END SUBROUTINE bridgeGA

!=======================================================================================================
      SUBROUTINE linkGA(ms,m,nspan,KBar,WFP,FM,Code,THEW,SW,MM,                         &
     &     h,n,l,dn,dl,hspan,GCode)
      USE Vdata
      IMPLICIT NONE
      integer,target:: m(3),nspan
      mwPointer ms(3)
      REAL*8,DIMENSION(ms(1),ms(2),ms(3)),Target::KBar,FM
      REAL*8,DIMENSION(ms(1),ms(2)),Target:: WFP,Code,THEW
      REAL*8,DIMENSION(m(1),m(2),m(3)),Target::SW
      REAL*8,Target:: n(m(1)),l(m(2)),h(m(3)),hspan(nspan,4)
      REAL*8,Target:: dn,dl,MM(m(1),m(2)),GCode(4)
      ! m2 is the watershed matrix size
      ! m is the GA matrix size (with n=m(1) and l=m(2))
      ! nspan is the number of spans
      ! SW is a lookup table for suction
      g%VDZ%GA%ms = ms
      g%VDZ%GA%m => m
      g%VDZ%GA%nspan => nspan
      g%VDZ%GA%hspan => hspan
      g%VDZ%GA%KBar => KBar
      g%VDZ%GA%WFP => WFP
      g%VDZ%GA%FM => FM
      g%VDZ%GA%Code => Code
      g%VDZ%GA%THEW => THEW
      g%VDZ%GA%SW => SW
      g%VDZ%GA%MM => MM
      g%VDZ%GA%h => h
      g%VDZ%GA%n => n
      g%VDZ%GA%l => l
      g%VDZ%GA%dn => dn
      g%VDZ%GA%dl => dl
      g%VDZ%GA%GCode => GCode
      END SUBROUTINE linkGA

!================================================================================================
      SUBROUTINE bridgeTopo(Topo,IS)
      use matRW
      IMPLICIT NONE
      !mwPointer :: mxGetField, mxGetDimensions, mxGetM, mxGetN,mxGetPr
      mwPointer NM,J
      PARAMETER (NM=6)
      mwPointer IS, nG(2), m(3),i
      mwPointer m2(2)
      CHARACTER*100 :: F(NM)
      REAL*8 y
      mwpointer :: P(NM),Topo,pm
      mwIndex :: index = 1
      mwSize :: number_of_dim

      F = (/ 'E','Ex','Ey','S0','S0x','S0y'/)
      DO i = 1,NM
        P(i) = mxGetPr(mxGetField(Topo,index,TRIM(F(i))))
      ENDDO

      pm = mxGetField(Topo,index,TRIM(F(1)))
      number_of_dim = mxGetNumberOfDimensions(pm)
      !CALL mxCopyPtrToInteger4(mxGetDimensions(pm),m2,2) ! This is the WFP size
      call mxCopyPtrToPtrArray(mxGetDimensions(pm),m2,number_of_dim)

      CALL linkTopo(m2,%val(P(1)),%val(P(2)),%val(P(3)),%val(P(4)),%val(P(5)),%val(P(6)))

      END SUBROUTINE bridgeTopo

!========================================================================================================
      SUBROUTINE linkTopo(ms,E,Ex,Ey,S0,S0x,S0y)
      USE Vdata
      IMPLICIT NONE
      mwPointer m(3),nspan,J
      mwPointer ms(2)
      REAL*8,DIMENSION(ms(1),ms(2)),Target:: E,Ex,Ey,S0,S0x,S0y
      ! m2 is the watershed matrix size
      ! m is the GA matrix size (with n=m(1) and l=m(2))
      ! nspan is the number of spans
      ! SW is a lookup table for suction
      g%Topo%ms = ms
      g%Topo%E => E
      g%Topo%Ex => Ex
      g%Topo%Ey => Ey
      g%Topo%S0 => S0
      g%Topo%S0x => S0x
      g%Topo%S0y => S0y
      END SUBROUTINE linkTopo

!====================================================================================================
      SUBROUTINE bridgeDay(Day,IS)
      use matRW
      implicit none
      !mwPointer :: mxGetField, mxGetDimensions, mxGetM, mxGetN,mxGetPr
      mwPointer IS, nG(2)
      mwPointer :: m(2)
      mwpointer:: Day,hli,cosz,pm
      mwpointer::prcp,temp,rad,ret,rets,wnd,hmd,hb,rn,dptp
      mwIndex :: index = 1
      mwSize :: number_of_dim

      m = 1
      pm = mxGetField(Day,index,'prcp')
      number_of_dim = mxGetNumberOfDimensions(pm)
      !CALL mxCopyPtrToInteger4(mxGetDimensions(pm),m,2)
      ! CALL mxCopyPtrToInteger4(mxGetDimensions(pm),m,3) (or make m(3)) and it will die on some machine!
      call mxCopyPtrToPtrArray(mxGetDimensions(pm),m,number_of_dim)
      prcp = mxGetPr(mxGetField(Day,index,'prcp'))
      temp = mxGetPr(mxGetField(Day,index,'temp'))
      rad = mxGetPr(mxGetField(Day,index,'rad'))
      ret = mxGetPr(mxGetField(Day,index,'ret'))
      rets = mxGetPr(mxGetField(Day,index,'rets'))
      wnd = mxGetPr(mxGetField(Day,index,'wnd'))
      hmd = mxGetPr(mxGetField(Day,index,'hmd'))
      hb = mxGetPr(mxGetField(Day,index,'hb'))
      hli = mxGetPr(mxGetField(Day,index,'hli'))
      cosz = mxGetPr(mxGetField(Day,index,'cosz'))
      rn = mxGetPr(mxGetField(Day,index,'rn'))
      dptp = mxGetPr(mxGetField(Day,index,'dptp'))

      CALL linkDay(1_8,m,%val(prcp),%val(temp),%val(rad),                              &
     &     %val(ret),%val(rets),%val(wnd),%val(hmd),%val(rn),%val(dptp),             &
     &     %val(hb),%val(hli),%val(cosz))

      END SUBROUTINE bridgeday

!======================================================================================================
      SUBROUTINE linkDay(J,m,prcp,temp,rad,ret,rets,wnd,hmd,               &
     &            rn,dptp,hb,hli,cosz)
      USE Vdata
      mwPointer J
      mwPointer m(2)
      REAL*8,DIMENSION(m(1),m(2)),Target::prcp,temp,rad,ret,                     &
     &        rets,wnd,hmd,rn,dptp,hb,hli,cosz
      g%Wea%DAY%ns = m(2)
      g%Wea%DAY%ny = m(1)
      g%Wea%DAY%prcp => prcp
      g%Wea%DAY%temp => temp
      g%Wea%DAY%rad => rad
      g%Wea%DAY%ret => ret
      g%Wea%DAY%rets => rets
      g%Wea%DAY%wnd => wnd
      g%Wea%DAY%hmd => hmd
      g%Wea%DAY%hb => hb
      g%Wea%DAY%hli => hli
      g%Wea%DAY%cosz => cosz
      g%Wea%DAY%rn => rn
      g%Wea%DAY%dptp => dptp
      END SUBROUTINE linkDay

!================================================================================================
      SUBROUTINE bridgeTH(TH,IS)
      use matRW
      !mwPointer :: mxGetField, mxGetDimensions, mxGetM, mxGetN,mxGetPr
      mwPointer IS, nG(2)
      mwPointer m(3), ndfc(2),np(2)
      mwpointer:: TH,pm
      mwpointer:: a,c,bbase,Temp,KT,CT
      mwpointer:: Sstatev,Ssitev,Soutv,tss,tsa
      mwpointer:: Sparm, dfc, SMR,SNOP,SNOE,Salb,Bsub,state
      mwIndex :: index = 1
      mwSize :: number_of_dim

      a = mxGetPr(mxGetField(TH,index,'a'))
      pm = mxGetField(TH,index,'a')
      number_of_dim = mxGetNumberOfDimensions(pm)
      !CALL mxCopyPtrToInteger4(mxGetDimensions(pm),m,3)
      call mxCopyPtrToPtrArray(mxGetDimensions(pm),m,number_of_dim)
      c = mxGetPr(mxGetField(TH,index,'c'))
      bbase = mxGetPr(mxGetField(TH,index,'bbase'))
      Temp = mxGetPr(mxGetField(TH,index,'Temp'))
      KT = mxGetPr(mxGetField(TH,index,'KT'))
      CT = mxGetPr(mxGetField(TH,index,'CT'))
      SStatev = mxGetPr(mxGetField(TH,index,'Sstatev'))
      Ssitev = mxGetPr(mxGetField(TH,index,'Ssitev'))
      Soutv = mxGetPr(mxGetField(TH,index,'Soutv'))
      tss = mxGetPr(mxGetField(TH,index,'tss'))
      tsa = mxGetPr(mxGetField(TH,index,'tsa'))
      Sparm = mxGetPr(mxGetField(TH,index,'SParm'))
      dfc = mxGetPr(mxGetField(TH,index,'dfc'))
      SMR = mxGetPr(mxGetField(TH,index,'SMR'))
      SNOP = mxGetPr(mxGetField(TH,index,'SNOP'))
      SNOE = mxGetPr(mxGetField(TH,index,'SNOE'))
      Salb = mxGetPr(mxGetField(TH,index,'Salb'))
      Bsub = mxGetPr(mxGetField(TH,index,'BSub'))
      state = mxGetPr(mxGetField(TH,index,'state'))

      pm = mxGetField(TH,index,'dfc')
      number_of_dim = mxGetNumberOfDimensions(pm)
      !CALL mxCopyPtrToInteger4(mxGetDimensions(pm),ndfc,2)
      call mxCopyPtrToPtrArray(mxGetDimensions(pm),ndfc,number_of_dim)

      pm = mxGetField(TH,index,'SParm')
      number_of_dim = mxGetNumberOfDimensions(pm)
      !CALL mxCopyPtrToInteger4(mxGetDimensions(pm),np,2)
      call mxCopyPtrToPtrArray(mxGetDimensions(pm),np,number_of_dim)

      CALL linkTH(IS,m,%val(a),%val(c),%val(bbase),%val(Temp),%val(KT)                &
     &,%val(CT),%val(Sstatev),%val(Ssitev),%val(Soutv),%val(tss)                      &
     &,%val(tsa),%val(SParm),%val(dfc),%val(SMR),%val(SNOP),%val(SNOE)                &
     &,%val(Salb),%val(BSub),ndfc,np(1)*np(2),%val(state))

      END SUBROUTINE bridgeTH

!==================================================================================================
      SUBROUTINE linkTH(J,m,a,c,bbase,Temp,KT,CT,Sstatev,Ssitev,Soutv                 &
     &,tss,tsa,SParm,dfc,SMR,SNOP,SNOE,Salb,BSub,nd,np,state)

      USE Vdata
      mwPointer J
      mwPointer m(3), np, nd(2)
      REAL*8,DIMENSION(m(1),m(2),m(3)),Target::a,c,bbase,Temp,KT,CT
      REAL*8,DIMENSION(m(1),m(2)),Target:: SMR,SNOP,SNOE,Salb,BSub
      REAL*8,TARGET::Sstatev(m(1),m(2),11),Ssitev(m(1),m(2),6)                          &
     &   ,Soutv(m(1),m(2),23),tss(m(1),m(2),24),tsa(m(1),m(2),24),state
      REAL*8,TARGET::Sparm(np),dfc(nd(1),nd(2))

      allocate(g%VDZ%th)
      g%VDZ%th%ny = m(1)
      g%VDZ%th%nx = m(2)
      g%VDZ%th%a => a
      g%VDZ%th%c => c
      g%VDZ%th%bbase => bbase
      g%VDZ%th%Temp => Temp
      g%VDZ%th%KT => KT
      g%VDZ%th%CT => CT
      g%VDZ%th%Sstatev => Sstatev
      g%VDZ%th%Ssitev => Ssitev
      g%VDZ%th%Soutv => Soutv
      g%VDZ%th%tss => tss
      g%VDZ%th%tsa => tsa
      g%VDZ%th%SParm => SParm
      g%VDZ%th%dfc => dfc
      g%VDZ%th%SMR => SMR
      g%VDZ%th%SNOP => SNOP
      g%VDZ%th%SNOE => SNOE
      g%VDZ%th%Salb => Salb
      g%VDZ%th%BSub => BSub
      g%VDZ%th%np = np
      g%VDZ%th%nd = nd(1)
      g%VDZ%th%state => state
      END SUBROUTINE linkTH

!===================================================================================================
      SUBROUTINE bridgeOVN(OVN,IS)
      use matRW
      !mwPointer :: mxGetField, mxGetDimensions, mxGetM, mxGetN,mxGetPr
      mwPointer IS, nG(2)
      mwPointer m(2), m3(3), nFParm(2)
      mwpointer:: OVN,pm,SATH,OVF,Ex,prcp,hd,dP,Qoc,hod
      mwpointer:: h,SN,SNSAVE,Mann,ho,hc,hf,Rf,hback,ffrac,dt,t,S,U,V,dist,FParm
      mwIndex :: index = 1
      mwSize :: number_of_dim
       ! if any more integers, add here. and do the continuation if necessary
      pm = mxGetField(OVN,index,'h')
      number_of_dim = mxGetNumberOfDimensions(pm)
      !CALL mxCopyPtrToInteger4(mxGetDimensions(pm),m,2)
      call mxCopyPtrToPtrArray(mxGetDimensions(pm),m,number_of_dim)
      ! first dimension must carry correct dimension info!!!!
      h = mxGetPr(mxGetField(OVN,index,'h'))
      SN = mxGetPr(mxGetField(OVN,index,'SN'))
      SNSAVE = mxGetPr(mxGetField(OVN,index,'SNSAVE'))
      Mann = mxGetPr(mxGetField(OVN,index,'Mann'))
      ho = mxGetPr(mxGetField(OVN,index,'ho'))
      hd = mxGetPr(mxGetField(OVN,index,'hd'))
      dP = mxGetPr(mxGetField(OVN,index,'dP'))
      pm = mxGetField(OVN,index,'dP')
      number_of_dim = mxGetNumberOfDimensions(pm)
      !CALL mxCopyPtrToInteger4(mxGetDimensions(pm),m3,3)
      call mxCopyPtrToPtrArray(mxGetDimensions(pm),m3,number_of_dim)
      hc = mxGetPr(mxGetField(OVN,index,'hc'))
      hf = mxGetPr(mxGetField(OVN,index,'hf'))
      Rf = mxGetPr(mxGetField(OVN,index,'Rf'))
      SATH = mxGetPr(mxGetField(OVN,index,'SATH'))
      OVF = mxGetPr(mxGetField(OVN,index,'OVF'))
      hback = mxGetPr(mxGetField(OVN,index,'hback'))
      ffrac = mxGetPr(mxGetField(OVN,index,'ffrac'))
      dt = mxGetPr(mxGetField(OVN,index,'dt'))
      U = mxGetPr(mxGetField(OVN,index,'U'))
      V = mxGetPr(mxGetField(OVN,index,'V'))
      S = mxGetPr(mxGetField(OVN,index,'S'))
      Ex = mxGetPr(mxGetField(OVN,index,'Ex'))
      prcp = mxGetPr(mxGetField(OVN,index,'prcp'))
      t = mxGetPr(mxGetField(OVN,index,'t'))
      dt = mxGetPr(mxGetField(OVN,index,'dt'))
      dist = mxGetPr(mxGetField(OVN,index,'dist'))
      Qoc = mxGetPr(mxGetField(OVN,index,'Qoc'))
      hod = mxGetPr(mxGetField(OVN,index,'hod'))
      
      FParm = mxGetField(OVN,index,'FParm')
      number_of_dim = mxGetNumberOfDimensions(FParm)
      call mxCopyPtrToPtrArray(mxGetDimensions(FParm),nFParm,number_of_dim)
      FParm = mxGetPr(FParm)
      
      CALL linkOVN(IS,m,%val(h),%val(SN),%val(SNSAVE),%val(Mann),%val(ho),%val(hc),%val(hf)  &
     & ,%val(Rf),%val(SATH),%val(OVF),%val(hback),%val(ffrac),%val(U),%val(V)   &
     & ,%val(S),%val(hd),%val(Ex),%val(prcp),%val(dP),m3(3),%val(t),%val(dt),%val(dist),%val(Qoc),%val(hod),%val(FParm),nFParm)

      END SUBROUTINE bridgeOVN

!========================================================================================================
      SUBROUTINE linkOVN(J,m,h,SN,SNSAVE,Mann,ho,hc,hf,Rf,SATH,OVF,hback,ffrac,U,V,    &
     &                  S,hd,Ex,prcp,dP,m3,t,dt,dist,Qoc,hod,FParm,nFParm)
      USE Vdata
      mwPointer J
      mwPointer m(2),m3,nFParm(2)
      REAL*8,DIMENSION(m(1),m(2)),Target::h,SN,SNSAVE,Mann,ho,hc,hf,U,V,S,Ex,prcp
      REAL*8,DIMENSION(m(1),m(2)),Target::Rf,SATH,OVF,hback,ffrac,hd,Qoc,hod,dist
      REAL*8,TARGET:: dP(m(1),m(2),m3)
      REAL*8,DIMENSION(nFParm(2)),TARGET::FParm
      REAL*8,Target :: t,dt
      ! OVN will not be used much as of now, so only link a few fields
      g%OVN%ms = m
      g%OVN%h => h
      g%OVN%SN => SN
      g%OVN%SNSAVE => SNSAVE
      g%OVN%Mann => Mann
      g%OVN%ho => ho
      g%OVN%hc => hc
      g%OVN%hf => hf
      g%OVN%Rf => Rf
      g%OVN%SATH => SATH
      g%OVN%OVF => OVF
      g%OVN%hback => hback
      g%OVN%ffrac => ffrac
      g%OVN%t => t
      g%OVN%dt => dt
      g%OVN%U => U
      g%OVN%V => V
      g%OVN%S => S
      g%OVN%hd => hd
      g%OVN%Ex => Ex
      g%OVN%prcp => prcp
      g%OVN%dP => dP
      g%OVN%dist => dist
      g%OVN%Qoc => Qoc
      g%OVN%hod => hod
      g%OVN%FParm => FParm
      END SUBROUTINE linkOVN

!=============================================================================================================
      SUBROUTINE bridgeCond(Cond,GWindex,Condindex)
      USE matRW
      mwPointer :: Cond,GWindex,Condindex,pm
      mwPointer :: mc(2)
      mwpointer :: idx,Cidx,Sidx,Nidx,Eidx,Widx,Dat,K,D,h
      mwSize :: number_of_dim, One = 1
      integer :: TYP,bnd

      pm = mxGetField(Cond,One,'type')
      if(pm==0) return

      CALL mxCopyPtrToInteger4(mxGetPr(mxGetField(Cond,One,'type')),TYP,One)
      IF (TYP .eq. 1) THEN
        CALL mxCopyPtrToInteger4(mxGetPr(mxGetField(Cond,One,'bnd')),bnd,1)
        pm = mxGetField(Cond,One,'indices')
        number_of_dim = mxGetNumberOfDimensions(pm)
        idx = mxGetPr(pm)
        !CALL mxCopyPtrToInteger4(mxGetDimensions(pm2),mc,2)
        CALL mxCopyPtrToPtrArray(mxGetDimensions(pm),mc,number_of_dim)
        Cidx= mxGetPr(mxGetField(Cond,One,'Cidx'))
        Sidx= mxGetPr(mxGetField(Cond,One,'Sidx'))
        Nidx= mxGetPr(mxGetField(Cond,One,'Nidx'))
        Eidx= mxGetPr(mxGetField(Cond,One,'Eidx'))
        Widx= mxGetPr(mxGetField(Cond,One,'Widx'))
        Dat = mxGetPr(mxGetField(Cond,One,'h'))
        CALL linkBC1(GWindex,Condindex,TYP,bnd,%val(idx),%val(Cidx),%val(Sidx),             &
       &     %val(Nidx),%val(Eidx),%val(Widx),%val(Dat),mc(1)*mc(2))
      ELSEIF (TYP .eq. 2) THEN
        pm = mxGetField(Cond,One,'indices')
        number_of_dim = mxGetNumberOfDimensions(pm)
        idx = mxGetPr(pm)
        !CALL mxCopyPtrToInteger4(mxGetDimensions(pm2),mc,2)
        CALL mxCopyPtrToPtrArray(mxGetDimensions(pm),mc,number_of_dim)
        Cidx= mxGetPr(mxGetField(Cond,One,'Cidx'))
        K= mxGetPr(mxGetField(Cond,One,'K'))
        h= mxGetPr(mxGetField(Cond,One,'h'))
        D= mxGetPr(mxGetField(Cond,One,'D'))
        !Q= mxGetPr(mxGetField(pm,J,'Q'))
        CALL linkBC2(GWindex,Condindex,TYP,%val(idx),%val(Cidx),                   &
       &     %val(K),%val(D),%val(h),mc(1)*mc(2))
      ELSEIF (TYP .eq. 0) THEN
        pm = mxGetField(Cond,One,'indices')
        number_of_dim = mxGetNumberOfDimensions(pm)
        idx = mxGetPr(pm)
        !CALL mxCopyPtrToInteger4(mxGetDimensions(pm2),mc,2)
        CALL mxCopyPtrToPtrArray(mxGetDimensions(pm),mc,number_of_dim)
        CALL linkBC0(GWindex,Condindex,TYP,%val(idx),mc(1)*mc(2))
      ENDIF


      END SUBROUTINE bridgeCond

!=============================================================================================================
      SUBROUTINE linkBC1(J,Cidx,Typ,bnd,idx,C,S,Nidx,E,Widx,h,n)
      ! Linking codes for type 1 boundary conditions for the aquifer
      ! Link Fortran BC struct to Matlab struct (use pointers)
      USE Vdata
      IMPLICIT NONE
      integer,target :: m,n,TYP,J,bnd,Cidx
      integer,Target,DIMENSION(n)::C,S,Nidx,E,Widx,idx
      REAL*8,Target,DIMENSION(n)::h

      g%GW(J)%Cond(Cidx)%JType=TYP
      g%GW(J)%Cond(Cidx)%bnd =bnd
      g%GW(J)%Cond(Cidx)%idx =>idx
      g%GW(J)%Cond(Cidx)%Cidx =>C
      g%GW(J)%Cond(Cidx)%Sidx =>S
      g%GW(J)%Cond(Cidx)%Nidx =>Nidx
      g%GW(J)%Cond(Cidx)%Eidx =>E
      g%GW(J)%Cond(Cidx)%Widx =>Widx
      g%GW(J)%Cond(Cidx)%h =>h
      END SUBROUTINE linkBC1

!==============================================================================================================
      SUBROUTINE linkBC2(J,Cidx,Typ,idx,C,K,D,h,n)
      ! Linking codes for type 2 boundary conditions for the aquifer (quasi-3D lower aquifer seepage)
      ! Link Fortran BC struct to Matlab struct (use pointers)
      USE Vdata
      IMPLICIT NONE
      integer,target :: m,n,TYP,J,bnd,Cidx
      integer,Target,DIMENSION(n)::C,idx
      REAL*8,Target,DIMENSION(n)::h,K,D

      g%GW(J)%Cond(Cidx)%JType=TYP
      g%GW(J)%Cond(Cidx)%bnd = bnd !bnd bnd does not exist, this may cause problem?
      g%GW(J)%Cond(Cidx)%idx =>idx
      g%GW(J)%Cond(Cidx)%Cidx =>C
      g%GW(J)%Cond(Cidx)%K =>K
      g%GW(J)%Cond(Cidx)%D =>D
      g%GW(J)%Cond(Cidx)%h =>h
      END SUBROUTINE linkBC2

!===========================================================================================================
      SUBROUTINE linkBC0(J,Cidx,Typ,idx,n)
      ! Linking codes for type 0 boundary conditions for the aquifer
      ! BCs with C,S,N,E indices
      USE Vdata
      IMPLICIT NONE
      integer,target :: m,n,TYP,J,bnd, Cidx
      integer,Target,DIMENSION(n)::idx

      g%GW(J)%Cond(Cidx)%Jtype=TYP
      !BBC(J)%bnd =bnd
      g%GW(J)%Cond(Cidx)%idx =>idx
      END SUBROUTINE linkBC0
