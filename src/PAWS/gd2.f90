      
      MODULE pp_Mod
      ! A storage-retrieval unit
      PUBLIC:: pp
      PUBLIC:: getR1
      PRIVATE
      
      type :: pp
        real*8, dimension(:), pointer :: r1  => null()
        real*8, dimension(:,:), pointer :: r2 => null()
        real*8, dimension(:,:,:), pointer :: r3 => null()
        real*8, dimension(:,:,:,:), pointer :: r4 => null()
        real*8, dimension(:,:,:,:,:), pointer :: r5 => null()
        real*8, dimension(:,:,:,:,:,:), pointer :: r6 => null()
        integer,dimension(:), pointer :: i1 => null()
        integer,dimension(:,:), pointer :: i2 => null()
        integer,dimension(:,:,:), pointer :: i3 => null()
        integer,dimension(:,:,:,:), pointer :: i4 => null()
        integer,dimension(:,:,:,:,:), pointer :: i5 => null()
        integer,dimension(:,:,:,:,:,:), pointer :: i6 => null()
        logical*1,dimension(:), pointer :: l1 => null()
        logical*1,dimension(:,:), pointer :: l2 => null()
        logical*1,dimension(:,:,:), pointer :: l3 => null()
        logical*1,dimension(:,:,:,:), pointer :: l4 => null()
        logical*1,dimension(:,:,:,:,:), pointer :: l5 => null()
        logical*1,dimension(:,:,:,:,:,:), pointer :: l6 => null()
        CONTAINS
         PROCEDURE:: deallocate => deallocate_pp
         PROCEDURE:: getR1
         PROCEDURE:: getR1s
         PROCEDURE:: setR1s
         PROCEDURE:: returnTypeN
         
      end type pp
      
      Contains
      
      function returnTypeN(obj)
      class(pp) obj
      character*10 returnTypeN
      returnTypeN = 'crap'
      end function
      
      subroutine deallocate_pp(obj,a_option) ! doesn't actually touch the memory pointed to
      class(pp) obj
      integer,optional:: a_option
      integer option
      option = 0
      if (present(a_option)) option = a_option
      if (option .eq. 0) then
        ! shallow deallocate
        obj%r1 => null()
        obj%r2 => null()
        obj%r3 => null()
        obj%r4 => null()
        obj%r5 => null()
        obj%r6 => null()
        obj%i1 => null()
        obj%i2 => null()
        obj%i3 => null()
        obj%i4 => null()
        obj%i5 => null()
        obj%i6 => null()
        obj%l1 => null()
        obj%l2 => null()
        obj%l3 => null()
        obj%l4 => null()
        obj%l5 => null()
        obj%l6 => null()
      elseif (option .eq. 1) then
        ! deep deallocate
        if (associated(obj%r1)) deallocate(obj%r1)
        if (associated(obj%r2)) deallocate(obj%r2)
        if (associated(obj%r3)) deallocate(obj%r3)
        if (associated(obj%r4)) deallocate(obj%r1)
        if (associated(obj%r5)) deallocate(obj%r1)
        if (associated(obj%r6)) deallocate(obj%r1)
        if (associated(obj%i1)) deallocate(obj%i1)
        if (associated(obj%i2)) deallocate(obj%i2)
        if (associated(obj%i3)) deallocate(obj%i3)
        if (associated(obj%i4)) deallocate(obj%i4)
        if (associated(obj%i5)) deallocate(obj%i5)
        if (associated(obj%i6)) deallocate(obj%i6)
        if (associated(obj%l1)) deallocate(obj%l1)
        if (associated(obj%l2)) deallocate(obj%l2)
        if (associated(obj%l3)) deallocate(obj%l3)
        if (associated(obj%l4)) deallocate(obj%l4)
        if (associated(obj%l5)) deallocate(obj%l5)
        if (associated(obj%l6)) deallocate(obj%l6)
      endif
      end subroutine 
      
      function getR1(obj)
      class(pp) obj
      Real*8,DIMENSION(:),pointer::getR1
      if (.not. associated(obj%r1)) then
        call ch_assert('pp::getR1::Requested type R1, but the stored type is '//returnTypeN(obj)) ! reason for not issuing an error is that higher level object can give a better error description
      else
        getR1 => obj%r1
      endif
      end function
      
      subroutine getR1s(obj,a_ptr)
      class(pp) obj
      Real*8,DIMENSION(:),pointer::a_ptr
      a_ptr => getR1(obj)
      end subroutine
      
      subroutine setR1s(obj,a_ptr)
      class(pp) obj
      Real*8,DIMENSION(:),pointer::a_ptr
      call deallocate_pp(obj)
      obj%r1 => a_ptr
      end subroutine
      
      end module
      
      
      MODULE lev_mod
      USE pp_mod, Only: pp
      PUBLIC:: lev
      Private
      type :: lev
        type(pp),dimension(:),pointer:: b
        integer m_maxBoxes
        integer m_datInd
        Contains
        Procedure:: dataIterator
        Procedure:: allocate => deallocate_lev
        Procedure:: get
      endtype lev
      
      CONTAINS
      
      function dataIterator(l)
      class(lev) l
      integer dataIterator
      l%m_datInd = l%m_datInd + 1
      if (l%m_datInd > l%m_maxBoxes) l%m_datInd = 0
      dataIterator = l%m_datInd
      end function
      
      subroutine allocate_lev(l,a_numBoxes)
      class(lev) l
      integer a_numBoxes
      
      call deallocate_lev(l)
      allocate(l%b(a_numBoxes))
      l%m_maxBoxes = a_numBoxes
      end subroutine
      
      subroutine deallocate_lev(l,a_option)
      class(lev) l
      integer I
      integer,optional:: a_option
      integer option
      option = 0
      if (present(a_option)) option = a_option
      
      if (associated(l%b)) then
        DO I=1,size(l%b)
          call l%b(I)%deallocate(option)
        ENDDO
        deallocate(l%b)
      endif
      end subroutine
      
      function get(l,a_ind)
      type(pp),pointer :: get
      class(lev) l
      integer,optional:: a_ind
      integer ind
      ind = 1
      if (present(a_ind)) ind = a_ind
      get => l%b(ind)
      end function
      
      subroutine set(l,b,a_ind)
      type(pp),pointer :: b
      class(lev) l
      integer,optional:: a_ind
      integer ind
      ind = 1
      if (present(a_ind)) ind = a_ind
      l%b(ind) = b
      end subroutine
      
      END MODULE
      
      
      MODULE gd
      USE lev_mod
      integer,parameter :: nTypeStr = 30
      type:: gdLevel
      integer np                                    ! number of associated fields
      type(lev),dimension(:),pointer::l  => null() ! for all pointers
      type(gdLevel),dimension(:),pointer::g  => null() ! for sub-structures
      character*100,dimension(:),pointer::f  => null() ! fieldnames
      
      endtype gdLevel
      END MODULE
      
      
      
