!

      subroutine rAvgStages()
      !% (E) nComputing average stage for overland flow- channel exchange calculation
      !% Output:
      !% g.Riv.Stages
      !global r g w
      use Vdata
      implicit none
      real*8,dimension(:,:),POINTER :: M=>NULL()
      integer,dimension(:),POINTER :: Cidx=>NULL()
      integer, dimension(:), pointer :: rL=>NULL()
      real*8, dimension(:), allocatable :: stages
      real*8 :: gid
      integer i,c
      
      rL => w%CR
      gid = w%g
      do i = 1, size(rL)  !% A river may flow through multiple land grids
        c = rL(i)
        M => g%Riv%TM(i)%ptr
        Cidx => g%Riv%Cidx(i)%ptr
        if(allocated(stages)) deallocate(stages)
        allocate(stages(size(Cidx)))
        stages = r(c)%Riv%h(Cidx) + r(c)%Riv%E(Cidx)
        !g(gid).Riv.Stages{i}=(stages*M)';
        g%Riv%Stages(i)%ptr = matmul(stages,M)
        !%g(gid).Riv.W{i}=(r(c).Riv.w(Cidx)*M)'; This is computed only once
      enddo
      
      if(allocated(stages)) deallocate(stages)
      end subroutine rAvgStages
