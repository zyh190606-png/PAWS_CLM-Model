
      MODULE gd_Mod
      ! *****************************************************
      ! Chaopeng Shen
      ! Pennsyvalnia State University
      ! *****************************************************
      ! the general syntax is
      ! ++ call getPtr(gd_base,fieldname,ptr,num)  ! num is optional, default to 1
      ! this will link the data to the pointer, ptr, which may be a real*8,integer, logical array supported up to 6 dimension
      ! or a type(gd) array pointer.
      ! depending on what pointer is passed in, the calling interface will pick the corresponding type/dimension subroutine
      ! If the queried field is not there or if the field type/dimension does not match that of the ptr, an error is thrown.
      ! num is the index, if omitted, default to 1. type(gd) allows to access sub-structures.
      ! Overloaded operators as a fast access:
      ! Besides the calling interface getPtr, we have overloaded operator .G.,.F., .FF. and .FFF. serving as shorthand methods
      ! for extracting sub-structure (.G.), [real*8,dimension(:,:)] and [real*8,dimension(:,:,:)] data, respectively. 
      ! With these operators, num is always 1.
      
      ! Look at matTable to see how this is used to create an interface to Matlab file
      ! How to directly put stuff into the structure?
      ! Since in Fortran, we cannot pass pointers around (they will just get re-directed)
      
      PRIVATE
      
      integer, parameter :: maxDim = 6
      integer, parameter :: maxStrLen = 100
      
      type pp
        real*8, dimension(:), pointer :: p1  => null()
        real*8, dimension(:,:), pointer :: p2 => null()
        real*8, dimension(:,:,:), pointer :: p3 => null()
        real*8, dimension(:,:,:,:), pointer :: p4 => null()
        real*8, dimension(:,:,:,:,:), pointer :: p5 => null()
        real*8, dimension(:,:,:,:,:,:), pointer :: p6 => null()
        integer,dimension(:), pointer :: j1 => null()
        integer,dimension(:,:), pointer :: j2 => null()
        integer,dimension(:,:,:), pointer :: j3 => null()
        integer,dimension(:,:,:,:), pointer :: j4 => null()
        integer,dimension(:,:,:,:,:), pointer :: j5 => null()
        integer,dimension(:,:,:,:,:,:), pointer :: j6 => null()
        logical*1,dimension(:), pointer :: l1 => null()
        logical*1,dimension(:,:), pointer :: l2 => null()
        logical*1,dimension(:,:,:), pointer :: l3 => null()
        logical*1,dimension(:,:,:,:), pointer :: l4 => null()
        logical*1,dimension(:,:,:,:,:), pointer :: l5 => null()
        logical*1,dimension(:,:,:,:,:,:), pointer :: l6 => null()
      end type pp
      
      type gd ! general data type
        ! it comes with a strs array, an integer array, a pointers array, and an integer nf, the number of valid fields
        integer np                                    ! number of associated fields
        type(pp),dimension(:),pointer::p  => null() ! for all pointers
        type(gd),dimension(:),pointer::g  => null() ! for sub-structures
        character(len=maxStrLen),dimension(:),pointer :: name => null()                 ! a name that traces from root to current object
        character(len=maxStrLen),dimension(:),pointer::f  => null() ! fieldnames
        type(gd),pointer:: parent=>null()           ! uplink, this link must be implemented when the hierarchy is created
        integer parent_fnum                         ! field number of this object in parent
      end type gd
      
      
      interface operator (.G.) ! pointer to field
      module procedure getPtrStrFunc
      end interface
      interface operator (.F.) ! pointer to field
      module procedure getPtrR1Func
      end interface
      interface operator (.FF.) ! pointer to field
      module procedure getPtrR2Func
      end interface
      interface operator (.FFF.) ! pointer to field
      module procedure getPtrR3Func
      end interface
      
!      interface connectField
!      module procedure connectFieldReal8,connectFieldInteger,connectFieldLogical
!      module procedure connectFieldStruct
!      end interface
      type(gd),save,target :: gd_base,gd_local,gd_temp
      type(gd),pointer::gd_ptr
      integer totalCounter
      
      ! get access to a field
      interface getPtr         ! call getPtr(main,field,ptr,[num])
      module procedure getPtrR1,getPtrR2,getPtrR3,getPtrR4,getPtrR5
      module procedure getPtrI1,getPtrI2,getPtrI3,getPtrI4,getPtrI5
      module procedure getPtrL1,getPtrL2,getPtrL3,getPtrL4,getPtrL5
      module procedure getPtrR6,getPtrI6,getPtrL6
      module procedure getPtrStr
      end interface
      
      ! assign a piece of data to a structure (n-th)
      interface assignPtr      ! call assignPtr(nn,datND,main,n)
      module procedure mAssignR1,mAssignR2,mAssignR3,mAssignR4,mAssignR5
      module procedure mAssignI1,mAssignI2,mAssignI3,mAssignI4,mAssignI5
      module procedure mAssignL1,mAssignL2,mAssignL3,mAssignL4,mAssignL5
      module procedure mAssignR6,mAssignI6,mAssignL6
      module procedure mAssignPtrStr
      end interface
      
      ! add a piece of data into structure, with its name also added
      interface addPtr         ! call addPtr(nn,field,datND,main)
      module procedure addR1,addR2,addR3,addR4,addR5,addR6
      module procedure addI1,addI2,addI3,addI4,addI5,addI6
      module procedure addL1,addL2,addL3,addL4,addL5,addL6
      module procedure addPtrStr
      end interface
      
      PUBLIC:: pp, gd, totalCounter, countRepNum
      PUBLIC:: gd_base,gd_local,gd_temp, gd_ptr, isempty
      PUBLIC:: getPtr, addPtr, assignPtr
      PUBLIC:: returnType, returnTypeN, addFieldName, setFieldName
      PUBLIC:: searchStr, searchStrF, findRepNum, returnSizes
      PUBLIC:: deallocate_p,  deallocate_gd, allocate_gd
      PUBLIC:: operator(.G.),operator(.F.),operator(.FF.),                &
     &  operator(.FFF.)
      
      ! ###########################################################################
      ! we cannot do similar pretty things for assignment because you cannot 
      !get a pointer (where in fact a local pointer is assigned to the memory) 
      !and then assign it to sth else, it will only redirect your local pointer without
      !touching the actual pointer in the structure
      
      CONTAINS
      
      function searchStr(Fs,str,nf,num)
      integer nf, searchStr
      integer, optional :: num
      character(*),intent(in)::str
      character(*),intent(in),dimension(nf)::Fs
      integer I,J,K
      
      K = 0
      searchStr = -99
      DO i = 1,nf
        IF (TRIM(str) == TRIM(Fs(i))) THEN
          K = K + 1
        !write(*,*) i,TRIM(str),K,TRIM(Fs(i)),num
          IF (K .eq. num) then
            searchStr = i
            RETURN
          ENDIF
        ENDIF
      ENDDO
      
      end function
      
      FUNCTION str2num(str)
       ! Converting a string to a double variable
       REAL*8 str2num
       CHARACTER*(*) str
       
       read(str,*) str2num
       
      END FUNCTION
      
      FUNCTION getdim(main,field)
      type(gd) main
      character*10 typeN
      character(*),intent(in)::field
      integer getdim
       ! Converting a string to a double variable
       typeN = returnType(main,field)
       getDim = floor(str2num(typeN(2:2))+1D-8)
      END FUNCTION
             
      function searchStrF(main,str,num)  ! wrapper
      ! searchStrF(main,str)
      type(gd) main
      character(*),intent(in)::str
      integer, optional :: num
      integer nf, searchStrF
      integer I,J,K,n
      
      n = 1
      if (present(num)) n=num
      searchStrF = searchStr(main%f,str,main%np,n)
      end function
      
      subroutine allocate_gd(main,nAllo)
      type(gd) main
      integer, optional :: nAllo
      integer n
      n = 500
      if (present(nAllo)) n=nAllo
      
      allocate(main%p(n))
      allocate(main%f(n))
      if (present(nAllo)) allocate(main%g(n))
      allocate(main%name(1))
      do i=1,n
        main%f(i) = ' '
      enddo
      main%np = 0
      main%name(1) = ' '
      totalCounter = totalCounter + 1
      end subroutine
   
      subroutine deallocate_gd(main)
      type(gd) main
      integer m,n
      character*10, rt, cur
      do i=1,main%np 
        call deallocate_p(main,i)
      enddo
      main%np = 0
      deallocate(main%f)
      deallocate(main%p)
      end subroutine

      
      function addFieldName(main,str)
      ! Caution, this must be immediately followed by field populating
      type(gd) main
      character(*),intent(in)::str
      integer addFieldName
      n = main%np + 1
      main%f(n) = str
      !main%np = n  ! don't do this in here
      addFieldName = n
      end function
      
      function setFieldName(main,str,n)
      ! return the same-name index
      type(gd) main
      character(*),intent(in)::str
      integer setFieldName,n
      main%f(n) = str
      setFieldName = findRepNum(main, str, n)
      !write(*,*) 'setFieldName = ', setFieldName
      end function
      
      function countRepNum(main,str)
      integer countRepNum, num, i
      type(gd) main
      character(*),intent(in)::str
      countRepNum = 0
      DO I=1,main%np
         IF (TRIM(str) == TRIM(main%f(i))) THEN
            countRepNum = countRepNum + 1  
         ENDIF
      ENDDO
      end function
      
      function findRepNum(main,str,n)
      ! some field names are replicated
      ! this function finds the same-name index
      type(gd) main
      character(*),intent(in)::str
      integer findRepNum, n, k
      num = 0; k = 0;
      do while (k >= 0)
        num = num + 1
        k = searchStr(main%f,str,main%np,num)
            !write(*,*) 'wtf ', num, k, n, TRIM(main%f(num))
        if (k .eq. n) then
          findRepNum = num
          return
        endif
      enddo
      findRepNum = -99
      end function
      
      
      
      function returnType(main,field, num)
      ! Give you a type of the argument
      ! order of else if arranged to improve performance on common types
      character*10 returnType
      character(*),intent(in)::field
      integer, optional::num
      integer n, i, rt, m
      ! get the type of the return
      type(gd) main
      n = 0; m = 1; 
      if (present(num)) m =  num
      
      n = searchStr(main%f,field,main%np,m)
      returnType = returnTypeN(main,n)
      
      end function
      
      function returnTypeN(main,n)
      integer n
      type(gd) main
      character*10 rt, returntypeN
      
      if (n < 0 .or. n>main%np) then
        rt = 'NULL'
      else
        if (associated(main%g) .AND. main%g(n)%np>0) then
          rt = 'gd'
        elseif (associated(main%p(n)%p1)) then
          rt = 'R1'
        elseif (associated(main%p(n)%p2)) then
          rt = 'R2'
        elseif (associated(main%p(n)%p3)) then
          rt = 'R3'
        elseif (associated(main%p(n)%p4)) then
          rt = 'R4'
        elseif (associated(main%p(n)%j1)) then
          rt = 'J1'
        elseif (associated(main%p(n)%j2)) then
          rt = 'J2'
        elseif (associated(main%p(n)%j3)) then
          rt = 'J3'
        elseif (associated(main%p(n)%j4)) then
          rt = 'J4'
        elseif (associated(main%p(n)%l1)) then
          rt = 'L1'
        elseif (associated(main%p(n)%l2)) then
          rt = 'L2'
        elseif (associated(main%p(n)%l3)) then
          rt = 'L3'
        elseif (associated(main%p(n)%l4)) then
          rt = 'L4'
        elseif (associated(main%p(n)%j5)) then
          rt = 'J5'
        elseif (associated(main%p(n)%j6)) then
          rt = 'J6'
        elseif (associated(main%p(n)%p5)) then
          rt = 'R5'
        elseif (associated(main%p(n)%p6)) then
          rt = 'R6'
        elseif (associated(main%p(n)%l5)) then
          rt = 'L5'
        elseif (associated(main%p(n)%l6)) then
          rt = 'L6'
        else
          rt = 'NULL'
        endif  
      endif
      returntypeN = rt
      end function
      
      
      
      function returnSizes(main,field, num)
      ! Give you a type of the argument
      ! order of else if arranged to improve performance on common types
      integer :: returnSizes(maxDim)
      character(*),intent(in)::field
      integer, optional::num
      integer n, i, rt, m
      ! get the type of the return
      type(gd) main
      n = 0; m = 1; 
      if (present(num)) m =  num
      
      n = searchStr(main%f,field,main%np,m)
      returnSizes = returnSizesN(main,n)
      
      end function
      
      function returnSizesN(main,n)
      integer :: n
      integer :: returnSizesN(maxDim)
      type(gd) main
      integer rt(maxDim)
      character(len=maxStrLen) field
      
      if (n < 0 .or. n>main%np) then
        rt(1) = nan
      else
        if (associated(main%g) .AND. main%g(n)%np>0) then
          field = main%f(n)
          rt = countRepNum(main,field)
        elseif (associated(main%p(n)%p1)) then
          rt = size(main%p(n)%p1,1)
        elseif (associated(main%p(n)%p2)) then
          rt(1) = size(main%p(n)%p2,1)
          rt(2) = size(main%p(n)%p2,2)
        elseif (associated(main%p(n)%p3)) then
          rt(1) = size(main%p(n)%p3,1)
          rt(2) = size(main%p(n)%p3,2)
          rt(3) = size(main%p(n)%p3,3)
        elseif (associated(main%p(n)%p4)) then
          rt(1) = size(main%p(n)%p4,1)
          rt(2) = size(main%p(n)%p4,2)
          rt(3) = size(main%p(n)%p4,3)
          rt(4) = size(main%p(n)%p4,4)
        elseif (associated(main%p(n)%j1)) then
          rt = size(main%p(n)%j1,1)
        elseif (associated(main%p(n)%j2)) then
          rt = size(main%p(n)%j2,1)
          rt = size(main%p(n)%j2,2)
        elseif (associated(main%p(n)%j3)) then
          rt(1) = size(main%p(n)%j3,1)
          rt(2) = size(main%p(n)%j3,2)
          rt(3) = size(main%p(n)%j3,3)
        elseif (associated(main%p(n)%j4)) then
          rt(1) = size(main%p(n)%j4,1)
          rt(2) = size(main%p(n)%j4,2)
          rt(3) = size(main%p(n)%j4,3)
          rt(4) = size(main%p(n)%j4,4)
        elseif (associated(main%p(n)%l1)) then
          rt = size(main%p(n)%l1,1)
        elseif (associated(main%p(n)%l2)) then
          rt(1) = size(main%p(n)%l2,1)
          rt(2) = size(main%p(n)%l2,2)
        elseif (associated(main%p(n)%l3)) then
          rt(1) = size(main%p(n)%l3,1)
          rt(2) = size(main%p(n)%l3,2)
          rt(3) = size(main%p(n)%l3,3)
        elseif (associated(main%p(n)%l4)) then
          rt(1) = size(main%p(n)%l4,1)
          rt(2) = size(main%p(n)%l4,2)
          rt(3) = size(main%p(n)%l4,3)
          rt(4) = size(main%p(n)%l4,4)
        elseif (associated(main%p(n)%j5)) then
          rt(1) = size(main%p(n)%j5,1)
          rt(2) = size(main%p(n)%j5,2)
          rt(3) = size(main%p(n)%j5,3)
          rt(4) = size(main%p(n)%j5,4)
          rt(5) = size(main%p(n)%j5,5)
        elseif (associated(main%p(n)%j6)) then
          rt(1) = size(main%p(n)%j6,1)
          rt(2) = size(main%p(n)%j6,2)
          rt(3) = size(main%p(n)%j6,3)
          rt(4) = size(main%p(n)%j6,4)
          rt(5) = size(main%p(n)%j6,5)
          rt(6) = size(main%p(n)%j6,6)
        elseif (associated(main%p(n)%p5)) then
          rt(1) = size(main%p(n)%p5,1)
          rt(2) = size(main%p(n)%p5,2)
          rt(3) = size(main%p(n)%p5,3)
          rt(4) = size(main%p(n)%p5,4)
          rt(5) = size(main%p(n)%p5,5)
        elseif (associated(main%p(n)%p6)) then
          rt(1) = size(main%p(n)%p6,1)
          rt(2) = size(main%p(n)%p6,2)
          rt(3) = size(main%p(n)%p6,3)
          rt(4) = size(main%p(n)%p6,4)
          rt(5) = size(main%p(n)%p6,5)
          rt(6) = size(main%p(n)%p6,6)
        elseif (associated(main%p(n)%l5)) then
          rt(1) = size(main%p(n)%l5,1)
          rt(2) = size(main%p(n)%l5,2)
          rt(3) = size(main%p(n)%l5,3)
          rt(4) = size(main%p(n)%l5,4)
          rt(5) = size(main%p(n)%l5,5)
        elseif (associated(main%p(n)%l6)) then
          rt(1) = size(main%p(n)%l6,1)
          rt(2) = size(main%p(n)%l6,2)
          rt(3) = size(main%p(n)%l6,3)
          rt(4) = size(main%p(n)%l6,4)
          rt(5) = size(main%p(n)%l6,5)
          rt(6) = size(main%p(n)%l6,6)
        else
          rt = nan
        endif  
      endif
      returnSizesN = rt
      end function
      
      recursive subroutine deallocate_p(main,n)
      ! (attempt to) deallocate dynamic memory
      ! we must use alloc_err because some memory cannot be deallocated
      integer n
      type(gd) main
      integer alloc_err 
      character*200, ev
      
      if (n < 0 .or. n>main%np) then
        ! nothing to be done
      else
        if (associated(main%p(n)%p1)) then
          deallocate(main%p(n)%p1, STAT = ALLOC_ERR)
        elseif (associated(main%p(n)%p2)) then
          deallocate(main%p(n)%p2, STAT = ALLOC_ERR)
        elseif (associated(main%p(n)%p3)) then
          deallocate(main%p(n)%p3, STAT = ALLOC_ERR)
        elseif (associated(main%p(n)%p4)) then
          deallocate(main%p(n)%p4, STAT = ALLOC_ERR)
        elseif (associated(main%p(n)%j1)) then
          deallocate(main%p(n)%j1, STAT = ALLOC_ERR)
        elseif (associated(main%p(n)%j2)) then
          deallocate(main%p(n)%j2, STAT = ALLOC_ERR)
        elseif (associated(main%p(n)%j3)) then
          deallocate(main%p(n)%j3, STAT = ALLOC_ERR)
        elseif (associated(main%p(n)%j4)) then
          deallocate(main%p(n)%j4, STAT = ALLOC_ERR)
        elseif (associated(main%p(n)%l1)) then
          deallocate(main%p(n)%l1, STAT = ALLOC_ERR)
        elseif (associated(main%p(n)%l2)) then
          deallocate(main%p(n)%l2, STAT = ALLOC_ERR)
        elseif (associated(main%p(n)%l3)) then
          deallocate(main%p(n)%l3, STAT = ALLOC_ERR)
        elseif (associated(main%p(n)%l4)) then
          deallocate(main%p(n)%l4, STAT = ALLOC_ERR)
        elseif (associated(main%p(n)%j5)) then
          deallocate(main%p(n)%l5, STAT = ALLOC_ERR)
        elseif (associated(main%p(n)%j6)) then
          deallocate(main%p(n)%l6, STAT = ALLOC_ERR)
        elseif (associated(main%p(n)%p5)) then
          deallocate(main%p(n)%p5, STAT = ALLOC_ERR)
        elseif (associated(main%p(n)%p6)) then
          deallocate(main%p(n)%p6, STAT = ALLOC_ERR)
        elseif (associated(main%p(n)%l5)) then
          deallocate(main%p(n)%l5, STAT = ALLOC_ERR)
        elseif (associated(main%p(n)%l6)) then
          deallocate(main%p(n)%l6, STAT = ALLOC_ERR)
        elseif (associated(main%g)) then
          if (main%g(n)%np>0) then
          do i = 1,main%g(n)%np
            call deallocate_p(main%g(n),i)
          enddo
          endif
        endif  
      endif
      end subroutine
      
      function isempty(main,field)
      logical isempty
      type(gd) main
      character*(*) field
      i = searchStr(main%f,field,main%np,1)
      if (i<=0) then
        isempty = .true.
      else
        isempty = .false.
      endif
      end function
      
      
      !! ###########  ALL THAT IS BELOW IS VERY NASTY ########################
      !########### But I've given you the interfaces above ###################
      
      !#############  accessor routines: dimension and type ##################
      !#######################################################################
      
      subroutine getPtrR1(main,field,ptr,num)    ! names change
      ! Get the "field" from "main" to "ptr"
      ! return <0 for trouble and =0 for success
      real*8,dimension(:),pointer:: ptr          ! dimension and type may change
      type(gd) main
      character*(*) field
      integer, optional :: num
      integer i,n,s
      n = 1
      if (present(num)) n=num
      i = searchStr(main%f,field,main%np,n)
      if (i<=0) then
        call CH_ASSERT("getPtrR1:: no field called "//field)
      elseif (.not. associated(main%p(i)%p1)) then ! dimension and type may change
        call CH_ASSERT("getPtrR1:: type/dimension between pointer and table    &
     & class do not match "//field)
      else
        ptr => main%p(i)%p1                          ! dimension and type may change
      endif
      end subroutine
      
      subroutine getPtrR2(main,field,ptr,num)    ! names change
      ! Get the "field" from "main" to "ptr"
      ! return <0 for trouble and =0 for success
      
      real*8,dimension(:,:),pointer:: ptr          ! dimension and type may change
      type(gd) main
      character*(*) field
      integer, optional :: num
      integer i,n,s
      n = 1
      if (present(num)) n=num
      i = searchStr(main%f,field,main%np,n)
      if (i<=0) then
        call CH_ASSERT("getPtrR2:: no field called "//field)
      elseif (.not. associated(main%p(i)%p2)) then ! dimension and type may change
        call CH_ASSERT("getPtrR2:: type/dimension between pointer and table    &
     & class do not match "//field)
      else
        ptr => main%p(i)%p2                          ! dimension and type may change
      endif
      end subroutine
      
      subroutine getPtrR3(main,field,ptr,num)    ! names change
      ! Get the "field" from "main" to "ptr"
      ! return <0 for trouble and =0 for success
      
      real*8,dimension(:,:,:),pointer:: ptr          ! dimension and type may change
      type(gd) main
      character*(*) field
      integer, optional :: num
      integer i,n,s
      n = 1
      if (present(num)) n=num
      i = searchStr(main%f,field,main%np,n)
      if (i<=0) then
        call CH_ASSERT("getPtrR3:: no field called "//field)
      elseif (.not. associated(main%p(i)%p3)) then ! dimension and type may change
        call CH_ASSERT("getPtrR3:: type/dimension between pointer and table    &
     & class do not match "//field)
      else
        ptr => main%p(i)%p3                          ! dimension and type may change
      endif
      end subroutine
      
      
      subroutine getPtrR4(main,field,ptr,num)    ! names change
      ! Get the "field" from "main" to "ptr"
      ! return <0 for trouble and =0 for success
      
      real*8,dimension(:,:,:,:),pointer:: ptr          ! dimension and type may change
      type(gd) main
      character*(*) field
      integer, optional :: num
      integer i,n,s
      n = 1
      if (present(num)) n=num
      i = searchStr(main%f,field,main%np,n)
      if (i<=0) then
        call CH_ASSERT("getPtrR4:: no field called "//field)
      elseif (.not. associated(main%p(i)%p4)) then ! dimension and type may change
        call CH_ASSERT("getPtrR4:: type/dimension between pointer and table    &
     & class do not match "//field)
      else
        ptr => main%p(i)%p4                          ! dimension and type may change
      endif
      end subroutine
      
      subroutine getPtrR5(main,field,ptr,num)    ! names change
      ! Get the "field" from "main" to "ptr"
      ! return <0 for trouble and =0 for success
      
      real*8,dimension(:,:,:,:,:),pointer:: ptr          ! dimension and type may change
      type(gd) main
      character*(*) field
      integer, optional :: num
      integer i,n,s
      n = 1
      if (present(num)) n=num
      i = searchStr(main%f,field,main%np,n)
      if (i<=0) then
        call CH_ASSERT("getPtrR5:: no field called "//field)
      elseif (.not. associated(main%p(i)%p5)) then ! dimension and type may change
        call CH_ASSERT("getPtrR5:: type/dimension between pointer and table    &
     & class do not match "//field)
      else
        ptr => main%p(i)%p5                          ! dimension and type may change
      endif
      end subroutine
      
      subroutine getPtrR6(main,field,ptr,num)    ! names change
      ! Get the "field" from "main" to "ptr"
      ! return <0 for trouble and =0 for success
      
      real*8,dimension(:,:,:,:,:,:),pointer:: ptr          ! dimension and type may change
      type(gd) main
      character*(*) field
      integer, optional :: num
      integer i,n,s
      n = 1
      if (present(num)) n=num
      i = searchStr(main%f,field,main%np,n)
      if (i<=0) then
        call CH_ASSERT("getPtrR6:: no field called "//field)
      elseif (.not. associated(main%p(i)%p6)) then ! dimension and type may change
        call CH_ASSERT("getPtrR6:: type/dimension between pointer and table    &
     & class do not match "//field)
      else
        ptr => main%p(i)%p6                          ! dimension and type may change
      endif
      end subroutine
      
      subroutine getPtrI1(main,field,ptr,num)    ! names change
      ! Get the "field" from "main" to "ptr"
      ! return <0 for trouble and =0 for success
      
      integer,dimension(:),pointer:: ptr          ! dimension and type may change
      type(gd) main
      character*(*) field
      integer, optional :: num
      integer i,n,s
      n = 1
      if (present(num)) n=num
      i = searchStr(main%f,field,main%np,n)
      if (i<=0) then
        call CH_ASSERT("getPtrI1:: no field called "//field)
      elseif (.not. associated(main%p(i)%j1)) then ! dimension and type may change
        call CH_ASSERT("getPtrI1:: type/dimension between pointer and table    &
     & class do not match "//field)
      else
        ptr => main%p(i)%j1                          ! dimension and type may change
      endif
      end subroutine
      
      
      subroutine getPtrI2(main,field,ptr,num)    ! names change
      ! Get the "field" from "main" to "ptr"
      ! return <0 for trouble and =0 for success
      
      integer,dimension(:,:),pointer:: ptr          ! dimension and type may change
      type(gd) main
      character*(*) field
      integer, optional :: num
      integer i,n,s
      n = 1
      if (present(num)) n=num
      i = searchStr(main%f,field,main%np,n)
      if (i<=0) then
        call CH_ASSERT("getPtrI2:: no field called "//field)
      elseif (.not. associated(main%p(i)%j2)) then ! dimension and type may change
        call CH_ASSERT("getPtrI2:: type/dimension between pointer and table    &
     & class do not match "//field)
      else
        ptr => main%p(i)%j2                          ! dimension and type may change
      endif
      end subroutine
      
      subroutine getPtrI3(main,field,ptr,num)    ! names change
      ! Get the "field" from "main" to "ptr"
      ! return <0 for trouble and =0 for success
      
      integer,dimension(:,:,:),pointer:: ptr          ! dimension and type may change
      type(gd) main
      character*(*) field
      integer, optional :: num
      integer i,n,s
      n = 1
      if (present(num)) n=num
      i = searchStr(main%f,field,main%np,n)
      if (i<=0) then
        call CH_ASSERT("getPtrI3:: no field called "//field)
      elseif (.not. associated(main%p(i)%j3)) then ! dimension and type may change
        call CH_ASSERT("getPtrI3:: type/dimension between pointer and table    &
     & class do not match "//field)
      else
        ptr => main%p(i)%j3                          ! dimension and type may change
      endif
      end subroutine
      
      
      subroutine getPtrI4(main,field,ptr,num)    ! names change
      ! Get the "field" from "main" to "ptr"
      ! return <0 for trouble and =0 for success
      
      integer,dimension(:,:,:,:),pointer:: ptr          ! dimension and type may change
      type(gd) main
      character*(*) field
      integer, optional :: num
      integer i,n,s
      n = 1
      if (present(num)) n=num
      i = searchStr(main%f,field,main%np,n)
      if (i<=0) then
        call CH_ASSERT("getPtrI4:: no field called "//field)
      elseif (.not. associated(main%p(i)%j4)) then ! dimension and type may change
        call CH_ASSERT("getPtrI4:: type/dimension between pointer and table    &
     & class do not match "//field)
      else
        ptr => main%p(i)%j4                          ! dimension and type may change
      endif
      end subroutine
      
      subroutine getPtrI5(main,field,ptr,num)    ! names change
      ! Get the "field" from "main" to "ptr"
      ! return <0 for trouble and =0 for success
      
      integer,dimension(:,:,:,:,:),pointer:: ptr          ! dimension and type may change
      type(gd) main
      character*(*) field
      integer, optional :: num
      integer i,n,s
      n = 1
      if (present(num)) n=num
      i = searchStr(main%f,field,main%np,n)
      if (i<=0) then
        call CH_ASSERT("getPtrI5:: no field called "//field)
      elseif (.not. associated(main%p(i)%j5)) then ! dimension and type may change
        call CH_ASSERT("getPtrI5:: type/dimension between pointer and table    &
     & class do not match "//field)
      else
        ptr => main%p(i)%j5                          ! dimension and type may change
      endif
      end subroutine
      
      subroutine getPtrI6(main,field,ptr,num)    ! names change
      ! Get the "field" from "main" to "ptr"
      ! return <0 for trouble and =0 for success
      
      integer,dimension(:,:,:,:,:,:),pointer:: ptr          ! dimension and type may change
      type(gd) main
      character*(*) field
      integer, optional :: num
      integer i,n,s
      n = 1
      if (present(num)) n=num
      i = searchStr(main%f,field,main%np,n)
      if (i<=0) then
        call CH_ASSERT("getPtrI6:: no field called "//field)
      elseif (.not. associated(main%p(i)%j6)) then ! dimension and type may change
        call CH_ASSERT("getPtrI6:: type/dimension between pointer and table    &
     & class do not match "//field)
      else
        ptr => main%p(i)%j6                          ! dimension and type may change
      endif
      end subroutine
            
      
      subroutine getPtrL1(main,field,ptr,num)    ! names change
      ! Get the "field" from "main" to "ptr"
      ! return <0 for trouble and =0 for success
      
      logical*1,dimension(:),pointer:: ptr          ! dimension and type may change
      type(gd) main
      character*(*) field
      integer, optional :: num
      integer i,n,s
      n = 1
      if (present(num)) n=num
      i = searchStr(main%f,field,main%np,n)
      if (i<=0) then
        call CH_ASSERT("getPtrL1:: no field called "//field)
      elseif (.not. associated(main%p(i)%l1)) then ! dimension and type may change
        call CH_ASSERT("getPtrL1:: type/dimension between pointer and table    &
     & class do not match "//field)
      else
        ptr => main%p(i)%l1                          ! dimension and type may change
      endif
      end subroutine
      
      
      subroutine getPtrL2(main,field,ptr,num)    ! names change
      ! Get the "field" from "main" to "ptr"
      ! return <0 for trouble and =0 for success
      
      logical*1,dimension(:,:),pointer:: ptr          ! dimension and type may change
      type(gd) main
      character*(*) field
      integer, optional :: num
      integer i,n,s
      n = 1
      if (present(num)) n=num
      i = searchStr(main%f,field,main%np,n)
      if (i<=0) then
        call CH_ASSERT("getPtrL2:: no field called "//field)
      elseif (.not. associated(main%p(i)%l2)) then ! dimension and type may change
        call CH_ASSERT("getPtrL2:: type/dimension between pointer and table    &
     & class do not match "//field)
      else
        ptr => main%p(i)%l2                          ! dimension and type may change
      endif
      end subroutine
      
      subroutine getPtrL3(main,field,ptr,num)    ! names change
      ! Get the "field" from "main" to "ptr"
      ! return <0 for trouble and =0 for success
      
      logical*1,dimension(:,:,:),pointer:: ptr          ! dimension and type may change
      type(gd) main
      character*(*) field
      integer, optional :: num
      integer i,n,s
      n = 1
      if (present(num)) n=num
      i = searchStr(main%f,field,main%np,n)
      if (i<=0) then
        call CH_ASSERT("getPtrL3:: no field called "//field)
      elseif (.not. associated(main%p(i)%l3)) then ! dimension and type may change
        call CH_ASSERT("getPtrL3:: type/dimension between pointer and table    &
     & class do not match "//field)
      else
        ptr => main%p(i)%l3                          ! dimension and type may change
      endif
      end subroutine
      
      
      subroutine getPtrL4(main,field,ptr,num)    ! names change
      ! Get the "field" from "main" to "ptr"
      ! return <0 for trouble and =0 for success
      
      logical*1,dimension(:,:,:,:),pointer:: ptr          ! dimension and type may change
      type(gd) main
      character*(*) field
      integer, optional :: num
      integer i,n,s
      n = 1
      if (present(num)) n=num
      i = searchStr(main%f,field,main%np,n)
      if (i<=0) then
        call CH_ASSERT("getPtrL4:: no field called "//field)
      elseif (.not. associated(main%p(i)%l4)) then ! dimension and type may change
        call CH_ASSERT("getPtrL4:: type/dimension between pointer and table    &
     & class do not match "//field)
      else
        ptr => main%p(i)%l4                          ! dimension and type may change
      endif
      end subroutine
      
      subroutine getPtrL5(main,field,ptr,num)    ! names change
      ! Get the "field" from "main" to "ptr"
      ! return <0 for trouble and =0 for success
      
      logical*1,dimension(:,:,:,:,:),pointer:: ptr          ! dimension and type may change
      type(gd) main
      character*(*) field
      integer, optional :: num
      integer i,n,s
      n = 1
      if (present(num)) n=num
      i = searchStr(main%f,field,main%np,n)
      if (i<=0) then
        call CH_ASSERT("getPtrL5:: no field called "//field)
      elseif (.not. associated(main%p(i)%l5)) then ! dimension and type may change
        call CH_ASSERT("getPtrL5:: type/dimension between pointer and table    &
     & class do not match "//field)
      else
        ptr => main%p(i)%l5                          ! dimension and type may change
      endif
      end subroutine
      
      subroutine getPtrL6(main,field,ptr,num)    ! names change
      ! Get the "field" from "main" to "ptr"
      ! return <0 for trouble and =0 for success
      
      logical*1,dimension(:,:,:,:,:,:),pointer:: ptr          ! dimension and type may change
      type(gd) main
      character*(*) field
      integer, optional :: num
      integer i,n,s
      n = 1
      if (present(num)) n=num
      i = searchStr(main%f,field,main%np,n)
      if (i<=0) then
        call CH_ASSERT("getPtrL6:: no field called "//field)
      elseif (.not. associated(main%p(i)%l6)) then ! dimension and type may change
        call CH_ASSERT("getPtrL6:: type/dimension between pointer and table    &
     & class do not match "//field)
      else
        ptr => main%p(i)%l6                          ! dimension and type may change
      endif
      end subroutine
      
      subroutine getPtrStr(main,field,ptr,num)    ! names change
      ! Get the "field" from "main" to "ptr"
      ! return <0 for trouble and =0 for success
      
      type(gd),pointer:: ptr          ! dimension and type may change
      type(gd) main
      character*(*) field
      integer, optional :: num
      integer i,n,s
      n = 1
      if (present(num)) n=num
      i = searchStr(main%f,field,main%np,n)
      if (i<=0) then
        call CH_ASSERT("getPtrStr:: no field called "//field)
      !elseif (.not. associated(main%g(i))) then ! dimension and type may change
      !  call CH_ASSERT("getPtrStr:: type/dimension between pointer and table class do not match "//field)
      else
        ptr => main%g(i)                            ! dimension and type may change
      endif
      end subroutine
      
      function getPtrStrFunc(main,field)
      ! same as getPtrStr, but in the form of a function, (in order to overload operator .G.
      ! also can only search for the first item
      type(gd),pointer:: getPtrStrFunc          ! dimension and type may change
      type(gd),intent(in) :: main
      character*(*),intent(in) :: field
      integer i,n,s
      n = 1
      i = searchStr(main%f,field,main%np,n)
      if (i<=0) then
        call CH_ASSERT("getPtrStr:: no field called "//field)
      else
        getPtrStrFunc => main%g(i)                ! dimension and type may change
      endif
      end function
      
      function getPtrR1Func(main,field)
      ! same as getPtrStr, but in the form of a function, (in order to overload operator .G.
      ! also can only search for the first item
      real*8,dimension(:),pointer:: getPtrR1Func          ! dimension and type may change
      type(gd),intent(in) :: main
      character*(*),intent(in) :: field
      integer i,n,s
      n = 1
      i = searchStr(main%f,field,main%np,n)
      if (i<=0) then
        call CH_ASSERT("getPtrStr:: no field called "//field)
      else
        getPtrR1Func => main%p(i)%p1                ! dimension and type may change
      endif
      end function
      
      function getPtrR2Func(main,field)
      ! same as getPtrStr, but in the form of a function, (in order to overload operator .G.
      ! also can only search for the first item
      real*8,dimension(:,:),pointer:: getPtrR2Func          ! dimension and type may change
      type(gd),intent(in) :: main
      character*(*),intent(in) :: field
      integer i,n,s
      n = 1
      i = searchStr(main%f,field,main%np,n)
      if (i<=0) then
        call CH_ASSERT("getPtrStr:: no field called "//field)
      else
        getPtrR2Func => main%p(i)%p2                ! dimension and type may change
      endif
      end function
      
      function getPtrR3Func(main,field)
      ! same as getPtrStr, but in the form of a function, (in order to overload operator .G.
      ! also can only search for the first item
      real*8,dimension(:,:,:),pointer:: getPtrR3Func          ! dimension and type may change
      type(gd),intent(in) :: main
      character*(*),intent(in) :: field
      integer i,n,s
      n = 1
      i = searchStr(main%f,field,main%np,n)
      if (i<=0) then
        call CH_ASSERT("getPtrStr:: no field called "//field)
      else
        getPtrR3Func => main%p(i)%p3                ! dimension and type may change
      endif
      
      end function
      
      !#############  assign pointer routines: dimension and type (in module) #############
      !####################################################################################
      ! the reason of existence for this module is the above dilemma of assign pointer:
      ! (1) you cannot put these assignPointers in a module due to %VAL (that will be used in matTable)
      ! (2) it already used gd_mod (therefore gd_mod cannot call them back)
      ! gd itself must be allocated (with allocate_gd) before it is used 
      ! Otherwise you get an error
      
      subroutine mAssignR1(nn,datND,main,n)
      INTEGER nn(1),n
      REAL*8,Target:: datND(nn(1))
      type(gd) main
      call deallocate_p(main,n)
      main%p(n)%p1 => datND
      END subroutine 
      
      subroutine mAssignR2(nn,datND,main,n)
      INTEGER nn(2),n
      REAL*8,Target:: datND(nn(1),nn(2))
      type(gd) main
      call deallocate_p(main,n)
      main%p(n)%p2 => datND
      END subroutine 
      
      subroutine mAssignR3(nn,datND,main,n)
      INTEGER nn(3),n
      REAL*8,Target:: datND(nn(1),nn(2),nn(3))
      type(gd) main
      call deallocate_p(main,n)
      main%p(n)%p3 => datND
      END subroutine 
      
      subroutine mAssignR4(nn,datND,main,n)
      INTEGER nn(4),n
      REAL*8,Target:: datND(nn(1),nn(2),nn(3),nn(4))
      type(gd) main
      call deallocate_p(main,n)
      main%p(n)%p4 => datND
      END subroutine 
      
      subroutine mAssignR5(nn,datND,main,n)
      INTEGER nn(5),n
      REAL*8,Target:: datND(nn(1),nn(2),nn(3),nn(4),nn(5))
      type(gd) main
      call deallocate_p(main,n)
      main%p(n)%p5 => datND
      END subroutine 
      
      subroutine mAssignR6(nn,datND,main,n)
      INTEGER nn(6),n
      REAL*8,Target:: datND(nn(1),nn(2),nn(3),nn(4),nn(5),nn(6))
      type(gd) main
      call deallocate_p(main,n)
      main%p(n)%p6 => datND
      END subroutine 
      
      subroutine mAssignI1(nn,datND,main,n)
      INTEGER nn(1),n
      integer,Target:: datND(nn(1))
      type(gd) main
      call deallocate_p(main,n)
      main%p(n)%j1 => datND
      END subroutine 
      
      subroutine mAssignI2(nn,datND,main,n)
      INTEGER nn(2),n
      integer,Target:: datND(nn(1),nn(2))
      type(gd) main
      call deallocate_p(main,n)
      main%p(n)%j2 => datND
      END subroutine 
      
      subroutine mAssignI3(nn,datND,main,n)
      INTEGER nn(3),n
      integer,Target:: datND(nn(1),nn(2),nn(3))
      type(gd) main
      call deallocate_p(main,n)
      main%p(n)%j3 => datND
      END subroutine 
      
      subroutine mAssignI4(nn,datND,main,n)
      INTEGER nn(4),n
      integer,Target:: datND(nn(1),nn(2),nn(3),nn(4))
      type(gd) main
      call deallocate_p(main,n)
      main%p(n)%j4 => datND
      END subroutine 
      
      subroutine mAssignI5(nn,datND,main,n)
      INTEGER nn(5),n
      integer,Target:: datND(nn(1),nn(2),nn(3),nn(4),nn(5))
      type(gd) main
      call deallocate_p(main,n)
      main%p(n)%j5 => datND
      END subroutine 
      
      subroutine mAssignI6(nn,datND,main,n)
      INTEGER nn(6),n
      integer,Target:: datND(nn(1),nn(2),nn(3),nn(4),nn(5),nn(6))
      type(gd) main
      call deallocate_p(main,n)
      main%p(n)%j6 => datND
      END subroutine 
      
      subroutine mAssignL1(nn,datND,main,n)
      INTEGER nn(1),n
      logical*1,Target:: datND(nn(1))
      type(gd) main
      call deallocate_p(main,n)
      main%p(n)%l1 => datND
      END subroutine 
      
      subroutine mAssignL2(nn,datND,main,n)
      INTEGER nn(2),n
      logical*1,Target:: datND(nn(1),nn(2))
      type(gd) main
      call deallocate_p(main,n)
      main%p(n)%l2 => datND
      END subroutine 
      
      subroutine mAssignL3(nn,datND,main,n)
      INTEGER nn(3),n
      logical*1,Target:: datND(nn(1),nn(2),nn(3))
      type(gd) main
      call deallocate_p(main,n)
      main%p(n)%l3 => datND
      END subroutine 
      
      subroutine mAssignL4(nn,datND,main,n)
      INTEGER nn(4),n
      logical*1,Target:: datND(nn(1),nn(2),nn(3),nn(4))
      type(gd) main
      call deallocate_p(main,n)
      main%p(n)%l4 => datND
      END subroutine 
      
      subroutine mAssignL5(nn,datND,main,n)
      
      INTEGER nn(5),n
      logical*1,Target:: datND(nn(1),nn(2),nn(3),nn(4),nn(5))
      type(gd) main
      call deallocate_p(main,n)
      main%p(n)%l5 => datND
      END subroutine 
      
      subroutine mAssignL6(nn,datND,main,n)
      INTEGER nn(6),n
      logical*1,Target:: datND(nn(1),nn(2),nn(3),nn(4),nn(5),nn(6))
      type(gd) main
      call deallocate_p(main,n)
      main%p(n)%l6 => datND
      END subroutine 
      
      subroutine mAssignPtrStr(ptr,main,n)    ! names change
      type(gd),pointer:: ptr                  ! dimension and type may change
      type(gd) :: main
      integer i,n,s
      !ptr => main%g(n)
      !main%g(n) => ptr                           ! dimension and type may change
      call CH_ASSERT('not coded yet')
      end subroutine
      
      !####################################################################################
      !####################################################################################
      ! Unfortunately there is nothing such as Code re-use in Fortran
      subroutine addR1(nn,field,datND,main)
      INTEGER nn(1),n
      REAL*8,Target:: datND(nn(1))
      type(gd) main
      character*(*),intent(in) :: field
      n = addFieldName(main,field)
      call assignPtr(nn,datND,main,n); main%np = n
      END subroutine 
      
      subroutine addR2(nn,field,datND,main)
      INTEGER nn(2),n
      REAL*8,Target:: datND(nn(1),nn(2))
      type(gd) main
      character*(*),intent(in) :: field
      n = addFieldName(main,field)
      call assignPtr(nn,datND,main,n); main%np = n
      END subroutine 
      
      subroutine addR3(nn,field,datND,main)
      INTEGER nn(3),n
      REAL*8,Target:: datND(nn(1),nn(2),nn(3))
      type(gd) main
      character*(*),intent(in) :: field
      n = addFieldName(main,field)
      call assignPtr(nn,datND,main,n); main%np = n
      END subroutine 
      
      subroutine addR4(nn,field,datND,main)
      INTEGER nn(4),n
      REAL*8,Target:: datND(nn(1),nn(2),nn(3),nn(4))
      type(gd) main
      character*(*),intent(in) :: field
      n = addFieldName(main,field)
      call assignPtr(nn,datND,main,n); main%np = n
      END subroutine 
      
      subroutine addR5(nn,field,datND,main)
      INTEGER nn(5),n
      REAL*8,Target:: datND(nn(1),nn(2),nn(3),nn(4),nn(5))
      type(gd) main
      character*(*),intent(in) :: field
      n = addFieldName(main,field)
      call assignPtr(nn,datND,main,n); main%np = n
      END subroutine 
      
      subroutine addR6(nn,field,datND,main)
      INTEGER nn(6),n
      REAL*8,Target:: datND(nn(1),nn(2),nn(3),nn(4),nn(5),nn(6))
      type(gd) main
      character*(*),intent(in) :: field
      n = addFieldName(main,field)
      call assignPtr(nn,datND,main,n); main%np = n
      END subroutine 
      
      subroutine addI1(nn,field,datND,main)
      INTEGER nn(1),n
      integer,Target:: datND(nn(1))
      type(gd) main
      character*(*),intent(in) :: field
      n = addFieldName(main,field)
      call assignPtr(nn,datND,main,n); main%np = n
      END subroutine 
      
      subroutine addI2(nn,field,datND,main)
      INTEGER nn(2),n
      integer,Target:: datND(nn(1),nn(2))
      type(gd) main
      character*(*),intent(in) :: field
      n = addFieldName(main,field)
      call assignPtr(nn,datND,main,n); main%np = n
      END subroutine 
      
      subroutine addI3(nn,field,datND,main)
      INTEGER nn(3),n
      integer,Target:: datND(nn(1),nn(2),nn(3))
      type(gd) main
      character*(*),intent(in) :: field
      n = addFieldName(main,field)
      call assignPtr(nn,datND,main,n); main%np = n
      END subroutine 
      
      subroutine addI4(nn,field,datND,main)
      INTEGER nn(4),n
      integer,Target:: datND(nn(1),nn(2),nn(3),nn(4))
      type(gd) main
      character*(*),intent(in) :: field
      n = addFieldName(main,field)
      call assignPtr(nn,datND,main,n); main%np = n
      END subroutine 
      
      subroutine addI5(nn,field,datND,main)
      INTEGER nn(5),n
      integer,Target:: datND(nn(1),nn(2),nn(3),nn(4),nn(5))
      type(gd) main
      character*(*),intent(in) :: field
      n = addFieldName(main,field)
      call assignPtr(nn,datND,main,n); main%np = n
      END subroutine 
      
      subroutine addI6(nn,field,datND,main)
      INTEGER nn(6),n
      integer,Target:: datND(nn(1),nn(2),nn(3),nn(4),nn(5),nn(6))
      type(gd) main
      character*(*),intent(in) :: field
      n = addFieldName(main,field)
      call assignPtr(nn,datND,main,n); main%np = n
      END subroutine 
      
      subroutine addL1(nn,field,datND,main)
      INTEGER nn(1),n
      logical*1,Target:: datND(nn(1))
      type(gd) main
      character*(*),intent(in) :: field
      n = addFieldName(main,field)
      call assignPtr(nn,datND,main,n); main%np = n
      END subroutine 
      
      subroutine addL2(nn,field,datND,main)
      INTEGER nn(2),n
      logical*1,Target:: datND(nn(1),nn(2))
      type(gd) main
      character*(*),intent(in) :: field
      n = addFieldName(main,field)
      call assignPtr(nn,datND,main,n); main%np = n
      END subroutine 
      
      subroutine addL3(nn,field,datND,main)
      INTEGER nn(3),n
      logical*1,Target:: datND(nn(1),nn(2),nn(3))
      type(gd) main
      character*(*),intent(in) :: field
      n = addFieldName(main,field)
      call assignPtr(nn,datND,main,n); main%np = n
      END subroutine 
      
      subroutine addL4(nn,field,datND,main)
      INTEGER nn(4),n
      logical*1,Target:: datND(nn(1),nn(2),nn(3),nn(4))
      type(gd) main
      character*(*),intent(in) :: field
      n = addFieldName(main,field)
      call assignPtr(nn,datND,main,n); main%np = n
      END subroutine 
      
      subroutine addL5(nn,field,datND,main)
      
      INTEGER nn(5),n
      logical*1,Target:: datND(nn(1),nn(2),nn(3),nn(4),nn(5))
      type(gd) main
      character*(*),intent(in) :: field
      n = addFieldName(main,field)
      call assignPtr(nn,datND,main,n); main%np = n
      END subroutine 
      
      subroutine addL6(nn,field,datND,main)
      INTEGER nn(6),n
      logical*1,Target:: datND(nn(1),nn(2),nn(3),nn(4),nn(5),nn(6))
      type(gd) main
      character*(*),intent(in) :: field
      n = addFieldName(main,field)
      call assignPtr(nn,datND,main,n); main%np = n
      END subroutine 
      
      subroutine addPtrStr(ptr,main,n)    ! names change
      type(gd),pointer:: ptr                  ! dimension and type may change
      type(gd) :: main
      integer i,n,s
      !ptr => main%g(n)
      !main%g(n) => ptr                           ! dimension and type may change
      call CH_ASSERT('not coded yet')
      end subroutine
      
      
      
      END MODULE gd_mod
      
      
      function prepareEntryToGd_base(name)
      use gd_Mod
      ! adds a field (sub-structure into base). This is likely mat_file, or 'entrance'
      integer, parameter :: maxStrLen = 100  !!Duplicated if declaration in gd_Mod is changed, this also needs to be changed
      character(len=maxStrLen) :: names(100), name
      integer prepareEntryToGd_base
      integer nbasePtr ! I found that if I declare integer gdbasePtr = 10, there will be retrieval problems, astonishing
      parameter(nbasePtr = 10) ! shouldn't be so many base level structures (mat_file, 'entrance')
      if (.not. associated(gd_base%p)) call allocate_gd(gd_base)
      gd_base%name(1) = 'base'
      
      np = gd_base%np + 1
      gd_base%f(np) = name
      gd_base%np    = np
      
      if (.not. associated(gd_base%g)) then
         allocate(gd_base%g(nbasePtr))
      endif
      call allocate_gd(gd_base%g(np))
      prepareEntryToGd_base = np
      gd_base%g(np)%name(1)= trim('base.') // trim(name)
      
      end function
      
      !#############  assign pointer routines: dimension and type ############
      !#######################################################################
      ! put this outside matTable because I may want to use it in mex environment
      ! if it is put inside a module, %VAL fails
!      subroutine assign1dR(nn,datND,ptr)
!      INTEGER nn(1)
!      REAL*8,Target:: datND(nn(1))
!      REAL*8,DIMENSION(:),POINTER::ptr
!      ptr => datND
!      END subroutine
      ! This subroutine above is a lost cause. pointers cannot be passed around when they are not associated
      ! Making this subroutine a function still doesn't work, because function cannot take %val construct
      ! therefore, below is our only choice: pass the damn structure itself in
      
      subroutine assignR1(nn,datND,main,n)
      use gd_mod
      INTEGER nn(1),n
      REAL*8,Target:: datND(nn(1))
      type(gd) main
      call deallocate_p(main,n)
      main%p(n)%p1 => datND
      END subroutine 
      
      subroutine assignR2(nn,datND,main,n)
      use gd_mod
      INTEGER nn(2),n
      REAL*8,Target:: datND(nn(1),nn(2))
      type(gd) main
      call deallocate_p(main,n)
      main%p(n)%p2 => datND
      END subroutine 
      
      subroutine assignR3(nn,datND,main,n)
      use gd_mod
      INTEGER nn(3),n
      REAL*8,Target:: datND(nn(1),nn(2),nn(3))
      type(gd) main
      call deallocate_p(main,n)
      main%p(n)%p3 => datND
      END subroutine 
      
      subroutine assignR4(nn,datND,main,n)
      use gd_mod
      INTEGER nn(4),n
      REAL*8,Target:: datND(nn(1),nn(2),nn(3),nn(4))
      type(gd) main
      call deallocate_p(main,n)
      main%p(n)%p4 => datND
      END subroutine 
      
      subroutine assignR5(nn,datND,main,n)
      use gd_mod
      INTEGER nn(5),n
      REAL*8,Target:: datND(nn(1),nn(2),nn(3),nn(4),nn(5))
      type(gd) main
      call deallocate_p(main,n)
      main%p(n)%p5 => datND
      END subroutine 
      
      subroutine assignR6(nn,datND,main,n)
      use gd_mod
      INTEGER nn(6),n
      REAL*8,Target:: datND(nn(1),nn(2),nn(3),nn(4),nn(5),nn(6))
      type(gd) main
      call deallocate_p(main,n)
      main%p(n)%p6 => datND
      END subroutine 
      
      subroutine assignI1(nn,datND,main,n)
      use gd_mod
      INTEGER nn(1),n
      integer,Target:: datND(nn(1))
      type(gd) main
      call deallocate_p(main,n)
      main%p(n)%j1 => datND
      END subroutine 
      
      subroutine assignI2(nn,datND,main,n)
      use gd_mod
      INTEGER nn(2),n
      integer,Target:: datND(nn(1),nn(2))
      type(gd) main
      call deallocate_p(main,n)
      main%p(n)%j2 => datND
      END subroutine 
      
      subroutine assignI3(nn,datND,main,n)
      use gd_mod
      INTEGER nn(3),n
      integer,Target:: datND(nn(1),nn(2),nn(3))
      type(gd) main
      call deallocate_p(main,n)
      main%p(n)%j3 => datND
      END subroutine 
      
      subroutine assignI4(nn,datND,main,n)
      use gd_mod
      INTEGER nn(4),n
      integer,Target:: datND(nn(1),nn(2),nn(3),nn(4))
      type(gd) main
      call deallocate_p(main,n)
      main%p(n)%j4 => datND
      END subroutine 
      
      subroutine assignI5(nn,datND,main,n)
      use gd_mod
      INTEGER nn(5),n
      integer,Target:: datND(nn(1),nn(2),nn(3),nn(4),nn(5))
      type(gd) main
      call deallocate_p(main,n)
      main%p(n)%j5 => datND
      END subroutine 
      
      subroutine assignI6(nn,datND,main,n)
      use gd_mod
      INTEGER nn(6),n
      integer,Target:: datND(nn(1),nn(2),nn(3),nn(4),nn(5),nn(6))
      type(gd) main
      call deallocate_p(main,n)
      main%p(n)%j6 => datND
      END subroutine 
      
      subroutine assignL1(nn,datND,main,n)
      use gd_mod
      INTEGER nn(1),n
      logical*1,Target:: datND(nn(1))
      type(gd) main
      main%p(n)%l1 => datND
      END subroutine 
      
      subroutine assignL2(nn,datND,main,n)
      use gd_mod
      INTEGER nn(2),n
      logical*1,Target:: datND(nn(1),nn(2))
      type(gd) main
      call deallocate_p(main,n)
      main%p(n)%l2 => datND
      END subroutine 
      
      subroutine assignL3(nn,datND,main,n)
      use gd_mod
      INTEGER nn(3),n
      logical*1,Target:: datND(nn(1),nn(2),nn(3))
      type(gd) main
      call deallocate_p(main,n)
      main%p(n)%l3 => datND
      END subroutine 
      
      subroutine assignL4(nn,datND,main,n)
      use gd_mod
      INTEGER nn(4),n
      logical*1,Target:: datND(nn(1),nn(2),nn(3),nn(4))
      type(gd) main
      call deallocate_p(main,n)
      main%p(n)%l4 => datND
      END subroutine 
      
      subroutine assignL5(nn,datND,main,n)
      use gd_mod
      INTEGER nn(5),n
      logical*1,Target:: datND(nn(1),nn(2),nn(3),nn(4),nn(5))
      type(gd) main
      call deallocate_p(main,n)
      main%p(n)%l5 => datND
      END subroutine 
      
      subroutine assignL6(nn,datND,main,n)
      use gd_mod
      INTEGER nn(6),n
      logical*1,Target:: datND(nn(1),nn(2),nn(3),nn(4),nn(5),nn(6))
      type(gd) main
      call deallocate_p(main,n)
      main%p(n)%l6 => datND
      END subroutine 
      
      !#########  end assign pointer routines: dimension and type ############
      !#######################################################################
      
!      
!      interface connectField
!      subroutine connectFieldReal8
!      
!      subroutine connectFieldInteger
!      
!      subroutine connectFieldLogical
!      
!      end subroutine
!      end interface
      
