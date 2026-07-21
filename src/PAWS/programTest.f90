#include "fintrf.h"
  
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! THIS IS A TEST PROGRAM !!!!!!
    ! This file should not be included in PAWS build project!
    
    PROGRAM programTest
    use matTable
    use gd_Mod
    use Vdata
    use matRW
    real*8, dimension(:,:), pointer:: ptr
    type(gd),pointer:: gd_tempLocal
    mwpointer, dimension(:), pointer :: arg_pa
    character*(100), dimension(:), pointer :: arg_names
    
    call buildMatTable('Clinton.mat', arg_pa, arg_names) 
    
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
    
    n = totalCounter
    gd_tempLocal => gd_base
    gd_ptr => gd_tempLocal .G. 'Clinton.mat' .G. 'g'
    g%Topo%E => gd_ptr .G. 'Topo' .FF. 'E'
    ptr => gd_ptr .G. 'Topo' .FF. 'E'
    ptr => g%Topo%E
    print*, ptr(20,20)
    gd_ptr => gd_base .G. 'Clinton.mat' .G. 'g' .G. 'Topo'
    print*, size(gd_ptr .F. 'S0p')
    print*, size(gd_ptr .F. 'S0x') ! will give you 0 because type is wrong.
    print*, size(gd_ptr .FF. 'S0x') ! tested 4480
    gd_tempLocal => gd_tempLocal .G. 'Clinton.mat'
    call bridgeGlobalGd(gd_tempLocal)
    !call linkVdataFromGd(gd_tempLocal)
    END
    
    
      
      ! test case to run on different system
      subroutine gd_inHouseTest
      use gd_mod
      type(gd) :: gd_test
      real*8 :: r1(1),r2(2,2),r3(3,3,3),r4(4,4,4,4),r5(5,5,5,5,5)
      real*8 :: r6(6,6,6,6,6,6),j6(6,6,6,6,6,6),l6(6,6,6,6,6,6)
      integer:: j1(1),j2(2,2),j3(3,3,3),j4(4,4,4,4),j5(5,5,5,5,5)
      logical*1:: l1(1),l2(2,2),l3(3,3,3),l4(4,4,4,4),l5(5,5,5,5,5)
      real*8, pointer :: r1p(:)=>NULL(),r2p(:,:)=>NULL(),r3p(:,:,:)=>NULL()
      real*8, pointer :: r4p(:,:,:,:)=>NULL(),r5p(:,:,:,:,:)=>NULL(),        &
     &   r6p(:,:,:,:,:,:)=>NULL()
      integer, pointer:: j1p(:)=>NULL(),j2p(:,:)=>NULL(),j3p(:,:,:)=>NULL()
      integer, pointer:: j4p(:,:,:,:)=>NULL(),j5p(:,:,:,:,:)=>NULL(),        &
     &   j6p(:,:,:,:,:,:)=>NULL()
      logical*1, pointer:: l1p(:)=>NULL(),l2p(:,:)=>NULL(),                  &
     &   l3p(:,:,:)=>NULL(),l4p(:,:,:,:)=>NULL()
      logical*1, pointer:: l5p(:,:,:,:,:)=>NULL(),l6p(:,:,:,:,:,:)=>NULL()
      integer mm(6) 
      
      mm =(/1,2,3,4,5,6/)
      
      r1=1; r2=2; r3=3; r4=4; r5=5; r6=6;
      j1=1; j2=2; j3=3; j4=4; j5=5; j6=6;
      l1=.true.; l2=.false.; l3=.true.; l4=.false.; l5=.true.; l6=.false.;
      
      !interface getPtr         ! call getPtr(main,field,ptr,[num])
      
      !interface assignPtr      ! call assignPtr(nn,datND,main,n)
      
      !interface addPtr         ! call addPtr(nn,field,datND,main)
      print*, associated(r1p)
      if (.not. associated(gd_test%p)) call allocate_gd(gd_test)
      if (.not. associated(gd_test%p)) write(*,*) 'allocate_gd failed1'
      ! add a piece of data into structure, with its name also added
      call addptr(mm,'r1',r1,gd_test)
      call addptr(mm,'r2',r2,gd_test)
      call addptr(mm,'r3',r3,gd_test)
      call addptr(mm,'r4',r4,gd_test)
      call addptr(mm,'r5',r5,gd_test)
      call addptr(mm,'r6',r6,gd_test)
      call addptr(mm,'l1',l1,gd_test)
      call addptr(mm,'l2',l2,gd_test)
      call addptr(mm,'l3',l3,gd_test)
      call addptr(mm,'l4',l4,gd_test)
      call addptr(mm,'l5',l5,gd_test)
      call addptr(mm,'l6',l6,gd_test)
      call addptr(mm,'j1',j1,gd_test)
      call addptr(mm,'j2',j2,gd_test)
      call addptr(mm,'j3',j3,gd_test)
      call addptr(mm,'j4',j4,gd_test)
      call addptr(mm,'j5',j5,gd_test)
      call addptr(mm,'j6',j6,gd_test)
      write(*,*) 'addPtr finished'

      
      call getptr(gd_test,'r1',r1p); if (ANY(r1p .ne. 1D0)) write(*,*) 'r1 failed'
      call getptr(gd_test,'r2',r2p); if (ANY(r2p .ne. 2D0)) write(*,*) 'r2 failed'
      call getptr(gd_test,'r3',r3p); if (ANY(r3p .ne. 3D0)) write(*,*) 'r3 failed'
      call getptr(gd_test,'r4',r4p); if (ANY(r4p .ne. 4D0)) write(*,*) 'r4 failed'
      call getptr(gd_test,'r5',r5p); if (ANY(r5p .ne. 5D0)) write(*,*) 'r5 failed'
      call getptr(gd_test,'r6',r6p); if (ANY(r6p .ne. 6D0)) write(*,*) 'r6 failed'
      call getptr(gd_test,'j1',j1p); if (ANY(j1p .ne. 1D0)) write(*,*) 'j1 failed'
      call getptr(gd_test,'j2',j2p); if (ANY(j2p .ne. 2D0)) write(*,*) 'j2 failed'
      call getptr(gd_test,'j3',j3p); if (ANY(j3p .ne. 3D0)) write(*,*) 'j3 failed'
      call getptr(gd_test,'j4',j4p); if (ANY(j4p .ne. 4D0)) write(*,*) 'j4 failed'
      call getptr(gd_test,'j5',j5p); if (ANY(j5p .ne. 5D0)) write(*,*) 'j5 failed'
      call getptr(gd_test,'j6',j6p); if (ANY(j6p .ne. 6D0)) write(*,*) 'j6 failed'
      call getptr(gd_test,'l1',l1p); if (ANY(.not.l1p))     write(*,*) 'l1 failed'
      call getptr(gd_test,'l2',l2p); if (ANY(l2p))          write(*,*) 'l2 failed'
      call getptr(gd_test,'l3',l3p); if (ANY(.not.l3p))     write(*,*) 'l3 failed'
      call getptr(gd_test,'l4',l4p); if (ANY(l4p))          write(*,*) 'l4 failed'
      call getptr(gd_test,'l5',l5p); if (ANY(.not.l5p))     write(*,*) 'l5 failed'
      call getptr(gd_test,'l6',l6p); if (ANY(l6p))          write(*,*) 'l6 failed'
      write(*,*) 'getPtr verified'

      
      
      call deallocate_gd(gd_test)
      write(*,*) 'deallocate_gd finished'

      if (associated(gd_test%p)) write(*,*) 'deAllocate_gd failed'
      if (.not. associated(gd_test%p)) call allocate_gd(gd_test)
      if (.not. associated(gd_test%p)) write(*,*) 'allocate_gd failed2'
      write(*,*) 'allocate_gd finished'
      
      n =1
      call assignptr(mm,r3,gd_test,n); gd_test%np =1
      k = setFieldName(gd_test,'r3',n); if (k .ne. 1) write(*,*) 'setFieldName failed'
      call getptr(gd_test,'r3',r3p,1); if (ANY(r3p .ne. 3D0)) write(*,*) 'assignPtr failed'
      write(*,*) 'assignPtrR1 verified'

      
      ! WITHOUT THESE TWO STATEMENTS, this test actually fails on Release version, Ifort, 64-bit WIN
      ! My guess is for a characters array, the optimized code does sth fancy and has wrong memory retrieving
      n =4
      call assignptr(mm,r4,gd_test,n); gd_test%np =n
      k = setFieldName(gd_test,'r3',n); if (k .ne. 2) write(*,*) 'setFieldName failed2',k, n ! this test actually failed on Release version 64-bit win
      call getptr(gd_test,'r3',r4p,2); if (ANY(r4p .ne. 4D0)) write(*,*) 'multi-name assignPtr failed'
      write(*,*) 'assignPtrR4 verified'

      
      
      n =1
      call assignptr(mm,j3,gd_test,n); gd_test%np =1
      k = setFieldName(gd_test,'j3',n); if (k .ne. 1) write(*,*) 'setFieldName failed integer'
      call getptr(gd_test,'j3',j3p,1); if (ANY(j3p .ne. 3)) write(*,*) 'assignPtr failed'
      write(*,*) 'assignPtrJ1 verified'

      n =4
      call assignptr(mm,j4,gd_test,n); gd_test%np =n
      k = setFieldName(gd_test,'j3',n); if (k .ne. 2) write(*,*) 'setFieldName failed integer'
      call getptr(gd_test,'j3',j4p,2); if (ANY(j4p .ne. 4)) write(*,*) 'multi-name assignPtr failed integer'
      write(*,*) 'assignPtrJ4 verified'

      
      
      end subroutine
      