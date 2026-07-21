      MODULE memoryChunk
      ! CHaopeng
      ! this is a big chunk of storage that uses to hold CLM data
      ! all CLM data are stored in here. These names are not available without a USE statement
      ! only this module can store general data, other modules can only share (PUBLIC) data structure instances!
      ! The advantage is that these pointers can be allocated or directed to other memory (such as passed in from matlab)!
      PRIVATE
      real*8,pointer,public:: R1(:)=>NULL(),R2(:,:)=>NULL(),R3(:,:,:)=>NULL()
      real*8,pointer,public:: R4(:,:,:,:)=>NULL(),R5(:,:,:,:,:)=>NULL()
      integer,pointer,public::L1(:)=>NULL(),L2(:,:)=>NULL(),L3(:,:,:)=>NULL()
      integer,pointer,public::L4(:,:,:,:)=>NULL(),L5(:,:,:,:,:)=>NULL()
      LOGICAL,pointer,public:: TFM(:,:,:)=>NULL()
      
      ! THESE are for subroutine that for whatever reason needs temporary storages (save heap space, or results that needs to be saved from call to call, etc...)
      ! TEMP (TEMPL) are for 2D; TEMPZ(TEMPZL) are for 3D with nz as 3rd dimension, TEMPT(TEMPTL) are for 3D with nRPT as 3rd dimension. Their last dimensions are for layers(items)
      ! They come with an occupancy indicator: IT,ITZ,ITT,ITL,ITZL,ITTL who have the same lengths as last dimensions, which indicate if each layer is occupied (=1) or not (=0)
      ! They are allocated in allocVData
      real*8,pointer,public :: TEMP(:,:,:)=>NULL(),TEMPZ(:,:,:,:)=>NULL(),TEMPT(:,:,:,:)=>NULL()
      integer,pointer,public :: TEMPL(:,:,:)=>NULL(),TEMPZL(:,:,:,:)=>NULL(),TEMPTL(:,:,:,:)=>NULL()
      integer,pointer,public :: IT(:)=>NULL(),ITZ(:)=>NULL(),ITT(:)=>NULL(),ITL(:)=>NULL()
      integer,pointer,public :: ITZL(:)=>NULL(),ITTL(:)=>NULL()
      integer,pointer,public :: cLevelMap(:,:)=> NULL()
      
      integer,target,public :: I1,I2,I3,I4,I5 ! real counter 
      integer,target,public :: J1,J2,J3,J4,J5 ! integer counter
      integer,target,public :: K3
      integer,target,public :: mm(4),mn,nn,iv
      ! mm = /nr,nc,nt,nz/, mn = nr*nc; nn = nr*nc*nz
      
      REAL*8,POINTER,public:: VBB(:,:)=>NULL()  ! A real*8 array (will be used to store pftcon data)
      REAL*8,POINTER,public:: VBCD(:,:)=>NULL() ! This is for VBC (crop.dat)
      
      type, public :: index_counter ! start and ending of indices in the central memory storage
      integer,dimension(2):: cps,pps,pepv
      integer,dimension(2):: ccs,cns,cws,ces
      integer,dimension(2):: ccf,cnf,cwf,cef,cmf
      integer,dimension(2):: g,lu,c,p
      end type index_counter
      
      type(index_counter),public:: R1C,R2C,R3C,R4C,R5C
      type(index_counter),public:: L1C,L2C,L3C,L4C,L5C
      
      
      interface tempStorage      ! call assignPtr(nn,datND,main,n)
      module procedure tempStorage0,tempStorage1,tempStorage2,tempStorage3,tempStorage4
      module procedure tempStorageI2,tempStorageI3,tempStorageI4
      end interface
      
      interface tAVG
      module procedure temporalMovingAverage2d
      end interface 
      
      PUBLIC:: FINDEMPTY, tempStorage, tAVG
      
      CONTAINS
         FUNCTION FINDEMPTY(LIST,N,M)
         use displayMod, Only: display
         ! find a consecutive portion of a list which is = 0
         ! always assume there is enough space!!!!!
         ! this function is used for temporary storages that needs to find a place
         INTEGER LIST(M)
         INTEGER I,N,K,FINDEMPTY
         I =  0; K = 0; FINDEMPTY = 0
         DO WHILE (K < N)
           I = I + 1
           IF (I > M) THEN
             call display( 'not enough temporary storage, allocate more temporary storage in allocVData')
             PAUSE
             call display( 'type to exit program')
             PAUSE
             STOP
           ENDIF
           IF (LIST(I) .eq. 0) THEN
             K = K + 1
           ELSE
             K = 0
           ENDIF
         ENDDO
         FINDEMPTY = I - N + 1
         END FUNCTION
         
         
         !! ######### TEMPORARY STORAGE ROUTINES ########################
         !! interface will support 2D, 3D, 4D REAL*8 and INTEGER
         !! flexible access to temporary storage, not field does not exist, create it
         !! type and dimension will be dependent on input ptr
         !! n is element in last dimension, optional, default to 1 if it is input, then it will be used to 
         ! a. check dimensional correctness if is already there
         ! b. allocate memory if it is not there
         !! HOWEVER, for 2D version, this check will not be done
         !! Notice this subroutine is intended for temporary storage, 
         !and will throw exceptions if you try to link other stuff to it, or have same named accesses, so be careful here
         ! by default, 2D will be (mm(1),mm(2)); 3D will be (mm(1),mm(2),mm)
         
         
         subroutine tempStorage0(f,ptr)
         USE gd_mod,  only: gd_temp, searchStrF, getptr, allocate_gd, addptr
         implicit none
         character*(*) f
         integer i, itype, n
         integer:: a_nz
         
         !! overloaded input
         real*8,pointer :: ptr, tempPtr(:)
         integer :: ndim = 1                     !!overloaded
         integer nd(4)
         i = searchStrF(gd_temp,f); 
         if (i < 0) then
            if (.not. associated(gd_temp%p)) call allocate_gd(gd_temp)
            n = gd_temp%np + 1
            nd(1) = a_nz   !! overloaded
            allocate(tempPtr(1) )          !! overloaded
            call addPtr(nd,f,tempPtr,gd_temp); tempPtr=0
            ptr => tempPtr(1)
         else
            call getPtr(gd_temp,f,tempPtr)
            ptr => tempPtr(1)
         endif
         END subroutine
         
         subroutine tempStorage1(f,ptr,a_nz)
         USE gd_mod,  only: gd_temp, searchStrF, getptr, allocate_gd, addptr
         implicit none
         character*(*) f
         integer i, itype, n
         integer:: a_nz
         
         !! overloaded input
         real*8,dimension(:),pointer :: ptr
         integer :: ndim = 1                     !!overloaded
         integer nd(4)
         i = searchStrF(gd_temp,f); 
         if (i < 0) then
            if (.not. associated(gd_temp%p)) call allocate_gd(gd_temp)
            n = gd_temp%np + 1
            nd(1) = a_nz   !! overloaded
            allocate(ptr(a_nz) )          !! overloaded
            call addPtr(nd,f,ptr,gd_temp); ptr=0
         else
            call getPtr(gd_temp,f,ptr)
         endif
         END subroutine
         
         
         
         
         subroutine tempStorage2(f,ptr,a_nz)
         USE gd_mod,  only: gd_temp, searchStrF, getptr, allocate_gd, addptr
         implicit none
         character*(*) f
         integer i, itype, n
         integer, optional:: a_nz
         
         !! overloaded input
         real*8,dimension(:,:),pointer :: ptr
         integer :: ndim = 2                     !!overloaded
         integer nd(4)
         i = searchStrF(gd_temp,f); 
         if (i < 0) then
            if (.not. associated(gd_temp%p)) call allocate_gd(gd_temp)
            n = gd_temp%np + 1
            nd = (/mm(1), mm(2), mm(4), mm(3)/)   !! overloaded
            allocate(ptr(nd(1),nd(2)) )          !! overloaded
            call addPtr(nd,f,ptr,gd_temp); ptr=0
         else
            call getPtr(gd_temp,f,ptr)
         endif
         END subroutine
      
         subroutine tempStorage3(f,ptr,a_nz)
         USE gd_mod,  only: gd_temp, searchStrF, getptr, allocate_gd, addptr
         implicit none
         character*(*) f
         integer i, itype, n
         integer, optional:: a_nz
         
         !! overloaded input
         real*8,dimension(:,:,:),pointer :: ptr
         integer :: ndim = 3                    !!overloaded
         integer nd(4)
         i = searchStrF(gd_temp,f); 
         if (i < 0) then
            if (.not. associated(gd_temp%p)) call allocate_gd(gd_temp)
            n = gd_temp%np + 1
            nd = (/mm(1), mm(2), a_nz, mm(3)/)
            allocate(ptr(nd(1),nd(2), nd(3)) )          !! overloaded
            call addPtr(nd,f,ptr,gd_temp); ptr=0
         else
            call getPtr(gd_temp,f,ptr)
            if (size(ptr,ndim) .ne. a_nz) call ch_assert('incorrect dimension')
         endif
         END subroutine
      
         subroutine tempStorage4(f,ptr,a_nz)
         USE gd_mod,  only: gd_temp, searchStrF, getptr, allocate_gd, addptr
         implicit none
         character*(*) f
         integer i, itype, n
         integer, optional:: a_nz
         
         !! overloaded input
         real*8,dimension(:,:,:,:),pointer :: ptr
         integer :: ndim = 3                    !!overloaded
         integer nd(4)
         i = searchStrF(gd_temp,f); 
         if (i < 0) then
            if (.not. associated(gd_temp%p)) call allocate_gd(gd_temp)
            n = gd_temp%np + 1
            nd = (/mm(1), mm(2), mm(4), a_nz/)
            allocate(ptr(nd(1),nd(2),nd(3),nd(4)) )          !! overloaded
            call addPtr(nd,f,ptr,gd_temp); ptr=0
         else
            call getPtr(gd_temp,f,ptr)
            if (size(ptr,ndim) .ne. a_nz) call ch_assert('incorrect dimension')
         endif
         END subroutine
      
      
         subroutine tempStorageI2(f,ptr,a_nz)
         USE gd_mod,  only: gd_temp, searchStrF, getptr, allocate_gd, addptr
         implicit none
         character*(*) f
         integer i, itype, n
         integer, optional:: a_nz
         
         !! overloaded input
         integer,dimension(:,:),pointer :: ptr
         integer :: ndim = 2                     !!overloaded
         integer nd(4)
         i = searchStrF(gd_temp,f); 
         if (i < 0) then
            if (.not. associated(gd_temp%p)) call allocate_gd(gd_temp)
            n = gd_temp%np + 1
            nd = (/mm(1), mm(2), mm(4), mm(3)/)   !! overloaded
            allocate(ptr(nd(1),nd(2)) )          !! overloaded
            call addPtr(nd,f,ptr,gd_temp); ptr=0
         else
            call getPtr(gd_temp,f,ptr)
         endif
         END subroutine
      
         subroutine tempStorageI3(f,ptr,a_nz)
         USE gd_mod,  only: gd_temp, searchStrF, getptr, allocate_gd, addptr
         implicit none
         character*(*) f
         integer i, itype, n
         integer, optional:: a_nz
         
         !! overloaded input
         integer,dimension(:,:,:),pointer :: ptr
         integer :: ndim = 3                    !!overloaded
         integer nd(4)
         i = searchStrF(gd_temp,f); 
         if (i < 0) then
            if (.not. associated(gd_temp%p)) call allocate_gd(gd_temp)
            n = gd_temp%np + 1
            nd = (/mm(1), mm(2), a_nz, mm(3)/)
            allocate(ptr(nd(1),nd(2), nd(3)) )          !! overloaded
            call addPtr(nd,f,ptr,gd_temp); ptr=0
         else
            call getPtr(gd_temp,f,ptr)
            if (size(ptr,ndim) .ne. a_nz) call ch_assert('incorrect dimension')
         endif
         END subroutine
      
         subroutine tempStorageI4(f,ptr,a_nz)
         USE gd_mod,  only: gd_temp, searchStrF, getptr, allocate_gd, addptr
         implicit none
         character*(*) f
         integer i, itype, n
         integer, optional:: a_nz
         
         !! overloaded input
         integer,dimension(:,:,:,:),pointer :: ptr
         integer :: ndim = 3                    !!overloaded
         integer nd(4)
         i = searchStrF(gd_temp,f); 
         if (i < 0) then
            if (.not. associated(gd_temp%p)) call allocate_gd(gd_temp)
            n = gd_temp%np + 1
            nd = (/mm(1), mm(2), mm(4), a_nz/)
            allocate(ptr(nd(1),nd(2),nd(3),nd(4)) )          !! overloaded
            call addPtr(nd,f,ptr,gd_temp); ptr=0
         else
            call getPtr(gd_temp,f,ptr)
            if (size(ptr,ndim) .ne. a_nz) call ch_assert('incorrect dimension')
         endif
         END subroutine
         
         subroutine temporalMovingAverage2d(datIn,datOut,name,nstep)
         integer m(2),nstep
         real*8,dimension(:,:,:),pointer:: ptr=>NULL()
!         integer,dimension(:,:),pointer:: Counter
         real*8,dimension(mm(1),mm(2)):: datIn,datOut
         character(*) name
         character(100) nameC,nameP
         integer nd(4)
         nd = (/mm(1), mm(2), 0, 0/)
!         
!         nameC = TRIM(name)//'_Counter'
!         nameP = TRIM(name)//'_Pop'
         call tempStorage(name , ptr, nstep)
!         call tempStorage(nameC, Counter,1) ! but we only use Counter(1,1)
         
!         Counter(1,1) = Counter(1,1) + 1 ! in fact don't need Counter
         
         ptr(:,:,1:nstep-1)=ptr(:,:,2:nstep)
         ptr(:,:,nstep    )=datIn
         datOut            = sum(ptr,3)/nstep
         
         end subroutine
         
      END MODULE
