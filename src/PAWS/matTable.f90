#include "fintrf.h"
      MODULE matTable
      USE gd_mod
      ! mat file interface
      ! *****************************************************
      ! Chaopeng Shen
      ! Lawrence Berkeley National Lab
      ! *****************************************************
      ! build an intermediate table object which stores data read from a mat file. 
      ! It can be used later to conveniently link the data to computational objects
      ! syntax:
      ! ####################################################
      ! SOME WHERE IN THE PROGRAM YOU MUST FIRST DO:
      ! ++ use matTable
      ! ++ call buildMatTable(mat_file)
      ! ####################################################
      ! after this, you can access all variable stored in the mat_file via gd_base, which is a type(gd) variable in module matTable
      ! the general syntax is
      ! ++ call getPtr(gd_base,fieldname,ptr,num)  ! num is optional, default to 1
      ! this will link the data to the pointer, ptr, which may be a real*8,integer, logical array supported up to 6 dimension
      ! or a type(gd) array pointer.
      ! depending on what pointer is passed in, the calling interface will pick the corresponding type/dimension subroutine
      ! If the queried field is not there or if the field type/dimension does not match that of the ptr, an error is thrown.
      ! num is the index, if omitted, default to 1. type(gd) allows to access sub-structures.
      ! Overloaded operators as a fast access:
      ! Besides the calling interface getPtr, we have overloaded operator .G., .FF. and .FFF. serving as shorthand methods
      ! for extracting sub-structure, [real*8,dimension(:,:)] and [real*8,dimension(:,:,:)] data, respectively. 
      ! With these operators, num is always 1.
      !
      ! Examples below are stored in subroutine example, make it PROGRAM MAIN to run the test case.
      ! Example(1): in matlab, you have 
      ! -> a = rand(1,4); b = rand(1,2,3,4,5,6); c = int32(rand(1,2,3,4)*100); d = (rand(4,4)>0.5);
      ! to get these data into fortran:
      ! ++ use matTable
      ! ++ real*8,pointer:: a_ptr(:),b_ptr(:,:,:,:,:,:)
      ! ++ integer,pointer:: c_ptr(:,:,:,:)
      ! ++ logical*1,pointer:: d_ptr(:,:)
      ! ++ call getPtr(gd_base,'a',a_ptr)
      ! ++ call getPtr(gd_base,'b',b_ptr)
      ! ++ call getPtr(gd_base,'c',c_ptr)
      ! ++ call getPtr(gd_base,'d',d_ptr)
      
      ! Example(2): in matlab, you have
      ! -> a(1).f=rand(4); a(2).f=rand(3,4,5)
      ! in order to access a(2).f, you do:
      ! ++ use matTable
      ! ++ real*8,dimension(:,:),pointer:: f1
      ! ++ real*8,dimension(:,:,:),pointer:: f2
      ! ++ call getPtr(gd_base,'a',gd_ptr)    ! gd_ptr is a globally available temporary pointer of type(gd). This gets the structure a(1) to gd_ptr
      ! ++ call getPtr(gd_ptr,'f',f1)      ! if you put f2 in place of f1 here, it will be an error because a.f(1) is 2D. 
      ! ++ call getPtr(gd_base,'a',gd_ptr,2)  ! This gets the structure a(2) to gd_ptr
      ! ++ call getPtr(gd_ptr,'f',f2)      
      ! 
      ! Example(3): in matlab you have
      ! -> g.VDZ.h = rand(5,6,7)
      ! -> g.GW(2).h = rand(10,11)
      ! in Fortran:
      ! ++ use matTable
      ! ++ real*8,pointer:: vdzh(:,:),gw2h(:,:)
      ! ++ type(gd),pointer:: g
      ! ++ g => gd_base .G. 'g'
      ! ++ vdzh => (g .G. 'h') .FFF. 'h'
      ! ++ call getPtr(g, 'GW', gd_ptr, 2)
      ! ++ gw2h => gd_ptr.FF.h
      ! NOTE: when using these operators, thou shall not use parenthesis to enclose the RHS. some wierdo here.
      
      integer nG
      parameter (nG =300)
      PUBLIC :: buildMatTable, connectField
      
      CONTAINS
!      
!      function rptsearch(v1,str)
!      ! return pointer
!      real*8,dimension(:),pointer :: rptsearch
!      type(gd),intent(in)::v1
!      character(*),intent(in)::str
!      real*8,target:: z(1)
!      integer I,J
!        DO i = 1,v1%np
!          IF (TRIM(str) == TRIM(v1%fp(i))) THEN
!            rptsearch => v1%p(i)%p
!            RETURN
!          ENDIF
!        ENDDO
!        z = -99D0
!        rptsearch => z
!      end function rptsearch
      
      recursive subroutine connectField(main,dat,field,gdname)
      ! add "field" to the gd table "main" by connecting to "dat"
      ! so that main%f(np+1)=field and main%p(np+1)%ptr => dat (unwrapped from mx)
      ! if dat is a structure, main%g(np+1) => dat (unwrapped from mx)
!      use matTable
      implicit none
      type(gd),target:: main
      type(gd),pointer::nextLevel
      integer m,n,j
      character(*),intent(in)::field
      character(*),intent(in),optional::gdname
      
      integer*4 mxGetNumberOfFields,mxGetFieldNumber
      integer*4 mxIsDouble,mxIsInt32,mxIsLogical,mxIsStruct
      mwSize mxGetNumberOfDimensions
      mwPointer mxGetDimensions,mxGetFieldByNumber,mxGetPr
      character*(100) mxGetFieldNameByNumber
      
      integer*4 fieldnum,nfields,nallo
      integer*4 i
      integer di(6),ns(6)
      integer numel
      character*(100) fsub,str
      character*(100),dimension(nG):: fsubR
      character*(20) numstr
      character*(100), pointer :: namestr
      mwsize ndim
      mwPointer dims(6)
      
      mwPointer::  mxGetField
      mwIndex :: index
      mwPointer :: dat,pr
      mwPointer pm
      
      if (dat .eq. 0) return ! empty!!
      
      ndim = mxGetNumberOfDimensions(dat)
      CALL mxCopyPtrToPtrArray(mxGetDimensions(dat),dims,ndim)
      numel = 1; ns = 0
      DO j=1,ndim
        numel = numel*dims(j)
        if (dims(j) > 1) ns(j)=1
      ENDDO
      pr = mxGetPr(dat)
      di = dims
      index = 1
      
      if (numel .eq. 0) return ! empty!!
      
      if ((mxIsDouble(dat) .eq. 1)) then
        n = main%np + 1
        
        if (ndim .eq. 1 .or. (ndim<=2 .and. sum(ns)<=1)) then
          !main%p(n)%p1 => gen_Ptr(di,%val(pr))
          call assignR1(numel,%val(pr),main,n)
        elseif (ndim .eq. 2) then
          call assignR2(di,%val(pr),main,n)
        elseif (ndim .eq. 3) then
          call assignR3(di,%val(pr),main,n)
        elseif (ndim .eq. 4) then
          call assignR4(di,%val(pr),main,n)
        elseif (ndim .eq. 5) then
          call assignR5(di,%val(pr),main,n)
        elseif (ndim .eq. 6) then
          call assignR6(di,%val(pr),main,n)
        endif
        main%f(n) = field
        main%np = n
         
      elseif (mxIsInt32(dat) .eq. 1) then
        n = main%np + 1
     
        if (ndim .eq. 1 .or. (ndim<=2 .and. sum(ns)<=1)) then
          call assignI1(numel,%val(pr),main,n)  ! notice numel since arrays can be 1 x n or n x 1
        elseif (ndim .eq. 2) then
          call assignI2(di,%val(pr),main,n)
        elseif (ndim .eq. 3) then
          call assignI3(di,%val(pr),main,n)
        elseif (ndim .eq. 4) then
          call assignI4(di,%val(pr),main,n)
        elseif (ndim .eq. 5) then
          call assignI5(di,%val(pr),main,n)
        elseif (ndim .eq. 6) then
          call assignI6(di,%val(pr),main,n)
        endif
        main%f(n) = field
        main%np = n
     
      elseif (mxIsLogical(dat) .eq. 1) then
        n = main%np + 1
        
        if (ndim .eq. 1 .or. ((ndim <= 2) .and. sum(ns)<=1)) then
          call assignL1(numel,%val(pr),main,n)
        elseif (ndim .eq. 2) then
          call assignL2(di,%val(pr),main,n)
        elseif (ndim .eq. 3) then
          call assignL3(di,%val(pr),main,n)
        elseif (ndim .eq. 4) then
          call assignL4(di,%val(pr),main,n)
        elseif (ndim .eq. 5) then
          call assignL5(di,%val(pr),main,n)
        elseif (ndim .eq. 6) then
          call assignL6(di,%val(pr),main,n)
        endif
        main%f(n) = field
        main%np = n
      
      elseif (mxIsStruct(dat) .eq. 1) then
      
       nfields = mxGetNumberOfFields(dat)
       n = main%np
       do index = 1,numel
         nallo = countMxSubNumel(dat,index)
         IF (nallo > 0) THEN
         n = n + 1
       ! structure may be a structure array. e.g. g.GW(1:2) (GW is a struct itself)
       ! add each element as a new entry in the table
       ! name will ONLY be added to the sub-structure
         i = size(main%g); j = 0
         if (.not. associated(main%g)) then
           allocate(main%g(nG))
           j = 1
         endif
         call allocate_gd(main%g(n),nallo)
         
         allocate(namestr)
         if (numel> 1) then
             write(numstr,'(I16)') index; 
             str = trim(field)//trim('(')//trim(adjustl(numstr))//          &
     &               trim(')')
         else
             str = field
         endif
         if (present(gdname)) then
           namestr = trim(gdname) // trim('.')// trim(str)
         else
           namestr = trim(main%name(1)) // trim('.')// trim(str)
         endif
         main%g(n)%name(1) = namestr
         str = main%name(1)
         do i = 1,nfields
           pm = mxGetFieldByNumber(dat, index, i)
           fsub = mxGetFieldNameByNumber(dat,i)
           fsubR(i)=fsub
           !nextLevel => main%g(n)
           if (pm > 0) call connectField(main%g(n), pm, fsub)
         enddo
         
         main%f(n) = field
         main%np = n
         main%g(n)%parent => main
         main%g(n)%parent_fnum = n
         ELSE
         ENDIF
         enddo
     else
           ! we don't actually link textdata. 
           ! If connectField start linking character then countMxSubNumel needs to be changed too
     endif
      
      end subroutine
      
      function countMxSubNumel(dat,index)
      integer*4 mxGetNumberOfFields,mxGetFieldNumber
      integer*4 mxIsDouble,mxIsInt32,mxIsLogical,mxIsStruct
      mwSize mxGetNumberOfDimensions
      mwPointer mxGetDimensions,mxGetFieldByNumber,mxGetPr
      integer countMxSubNumel
      character*(100) mxGetFieldNameByNumber
      
      integer*4 fieldnum,nfields
      integer*4 i,j,k
      integer di(6),ns(6)
      integer numel
      character*(100) fsub,str
      character*(20) numstr
      character*(100), pointer :: namestr
      mwsize ndim
      mwPointer dims(6)
      
      mwPointer::  mxGetField
      mwIndex :: index
      mwPointer :: dat,pr
      mwPointer pm
      
      countMxSubNumel = 0
      nfields = mxGetNumberOfFields(dat)
      do i = 1,nfields
          
        pm = mxGetFieldByNumber(dat, index, i)
        fsub = mxGetFieldNameByNumber(dat,i)
        if (pm > 0) then
            ! otherwise empty pointer!
          if ((mxIsDouble(pm) .eq. 1)) then
            k = 1
            numel = 1
          elseif (mxIsInt32(pm) .eq. 1) then
            k = 2
            numel = 1
          elseif (mxIsLogical(pm) .eq. 1) then
            k = 3
            numel = 1
          elseif (mxIsStruct(pm) .eq. 1) then
            k = 4
            ndim = mxGetNumberOfDimensions(pm)
            CALL mxCopyPtrToPtrArray(mxGetDimensions(pm),dims,ndim)
            numel = 1; ns = 0
            DO j=1,ndim
              numel = numel*dims(j)
              if (dims(j) > 1) ns(j)=1
            ENDDO
          else
            k = 5
            numel = 0  ! we don't actually link textdata --- see connectField above. If connectField start linking character then here needs to be changed too
          endif
          countMxSubNumel = countMxSubNumel + numel
        endif
      enddo
      
      end function
      
      subroutine buildMatTable(mat_file, arg_pa, arg_names)
      ! read a workspace saved in mat_file into (gd_base .G. mat_file)
      !      mwPointer file_in, pGlobal
      integer, parameter :: cLen = 100
      mwpointer mp, dir, adir(100), pa
      mwpointer matOpen, matGetDir, matGetVariable
      mwpointer matGetVariableInfo, matGetNextVariableInfo
      mwpointer, dimension(:), optional, pointer :: arg_pa
      character*(cLen), dimension(:), optional, pointer :: arg_names
      integer*4 num, ndir2,nbase
      character*(cLen) names(100), name
      integer   i, stat, np, prepareEntryToGd_base
      integer ntemp
      character*(*) mat_file
      character*(cLen), dimension(10) :: skipVar
      integer nSkipVar, skip
      mwsize ndir
      !CP: I noticed that if ndir is set to integer, then problem may occur with openning multiple
      !mat files on 64-it platforms. This might have been neglected in matdemo2.F
      !with these mat****, mx**** commands, it is always wise to absolutely follow the data types
      !specified by Matlab help files

      nSkipVar = 1
      skipVar(1) = 'Prob'
      
      mp = 0
      dir = 0
      ndir = 0
      mp = matOpen(TRIM(mat_file), 'r')
      dir = matGetDir(mp, ndir2)
      ndir = ndir2
      
      call mxCopyPtrToPtrArray(dir, adir, ndir)
      do i=1,ndir
         call mxCopyPtrToCharacter(adir(i), names(i), 32)
      enddo
      
      !np = prepareEntryToGd_base(mat_file)
      nbase = 10
      if (.not. associated(gd_base%p)) call allocate_gd(gd_base,nbase)
      gd_base%name(1) = 'base'
      np = gd_base%np + 1
      gd_base%f(np) = mat_file
      gd_base%np    = np
      
      if (.not. associated(gd_base%g)) then
         allocate(gd_base%g(nbase))
      endif
      ntemp = nG
      call allocate_gd(gd_base%g(np),ntemp)
      gd_base%g(np)%name = mat_file
      
      if (present(arg_pa)) then
         allocate(arg_pa(ndir2))
         allocate(arg_names(ndir2))
      endif
      do i=1,ndir
         pa = matGetVariable(mp,names(i))
         if (present(arg_pa)) then
           arg_pa(i) = pa
           arg_names(i) = names(i)
         endif
         skip = 0
         do j=1, nSkipVar
           if (names(i) .eq. skipVar(j)) skip=1
         enddo
         if (skip .ne. 1) call connectField(gd_base%g(np),pa,names(i),mat_file)
      enddo
      stat = matClose(mp)
      
      end subroutine
      
      !#######################################################################
      !#######################################################################
      END MODULE matTable
      
      !PROGRAM MAIN  ! Uncomment this line to look at the example case
      ! I updated buildMatTable to put each mat file in its own workspace (under gd_base) so now this test doesn't work
      subroutine buildMatTable_example
      USE matTable
      
       real*8,pointer:: a_ptr(:)=>NULL(),b_ptr(:,:,:,:,:,:)=>NULL()
       integer,pointer:: c_ptr(:,:,:,:)=>NULL()
       logical*1,pointer:: d_ptr(:,:)=>NULL()
       
       real*8,dimension(:,:),pointer:: f1=>NULL()
       real*8,dimension(:,:,:),pointer:: f2=>NULL()
       
       real*8,pointer:: vdzh(:,:,:)=>NULL(),gw2h(:,:)=>NULL()
       type(gd),pointer:: g=>NULL()
       type(gd),pointer:: mat=>NULL()
      
       mwPointer :: file_in,pGlobal,pm
      
      
       call buildMatTable('matTable_example.mat') 
      ! Example(1): in matlab, you have 
      ! -> a = rand(1,4); b = rand(1,2,3,4,5,6); c = int32(rand(1,2,3,4)*100); d = (rand(4,4)>0.5);
       mat => gd_base .G. 'matTable_example.mat'
      
       call getPtr(mat,'a',a_ptr)
       call getPtr(mat,'b',b_ptr)
       call getPtr(mat,'c',c_ptr)
       call getPtr(mat,'d',d_ptr)
      
      ! Example(2): in matlab, you have
      ! -> aa(1).f=rand(4); aa(2).f=rand(3,4,5)
      ! in order to access a(2).f, you do:
       call getPtr(gd_base,'aa',gd_ptr)    ! gd_ptr is a globally available temporary pointer of type(gd). This gets the structure a(1) to gd_ptr
       call getPtr(gd_ptr,'f',f1)      ! if you put f2 in place of f1 here, it will be an error because a.f(1) is 2D. 
       call getPtr(gd_base,'aa',gd_ptr,2)  ! This gets the structure a(2) to gd_ptr
       call getPtr(gd_ptr,'f',f2)      
      ! 
      ! Example(3): in matlab you have
      ! -> g.VDZ.h = rand(5,6,7)
      ! -> g.GW(2).h = rand(10,11)
      ! in Fortran:
       g => gd_base .G. 'g'
       vdzh => g .G. 'VDZ' .FFF. 'h'
       call getPtr(g, 'GW', gd_ptr, 2)
       gw2h => gd_ptr.FF.'h'
       gd_ptr => gd_base%g(1)%g(2)%g(8)%g(8)%g(8)%g(8)%g(8)%g(8)
      END
      
