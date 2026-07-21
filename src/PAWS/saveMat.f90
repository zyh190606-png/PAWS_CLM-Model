#include "fintrf.h"
!

      subroutine saveMat(file_name,ivar)
      use matRW
      use displayMod, Only: display
      implicit none
      
      character(len=200) :: file_name
      mwPointer :: mp
      integer :: status, var
      integer, optional:: ivar
      
      var = 0
      if (present(ivar)) var = ivar


      mp = matOpen(trim(file_name), 'w')
      if (mp .eq. 0) then
         call display( 'Can''t open ''data_out.mat'' for writing.')
         stop
      end if
      if (var .eq. 0 .or. var .eq. 1) then
        status = matPutVariable(mp, 'g', pGlobalG)
        if (status .ne. 0) then
           call display( 'matPutVariable ''g'' failed')
           stop
        end if
      endif
      if (var .eq. 0 .or. var .eq. 2) then
         status = matPutVariable(mp, 'w', pGlobalW)
         if (status .ne. 0) then
           call display( 'matPutVariable ''w'' failed')
           stop
         end if
      endif 
      if (var .eq. 0 .or. var .eq. 3) then
          status = matPutVariable(mp, 'r', pGlobalR)
          if (status .ne. 0) then
             call display( 'matPutVariable ''r'' failed')
             stop
          end if
      endif
      status = matClose(mp)
      if (status .ne. 0) then
         call display( 'Error closing MAT-file')
         stop
      end if
      
      end subroutine saveMat
      

      subroutine saveW(file_name)
      use matRW
      use displayMod, Only: display
      implicit none
      
      character*(*) :: file_name
      !character(len=200),optional :: var_name
      ! This above optional argument does not work because procedures with optional arguments must have explicit interface.
      ! This can be made by include this procedure in a module. However, mex and mat statements cannot be in a module
      ! therefore, we cannot tolerate optional arguments in procedures containing mat- or mx- operations.
      ! http://stackoverflow.com/questions/18300229/fortran-95-optional-statement-does-not-work-using-ftn95-and-plato
      character(len=100) :: var
      mwPointer :: mp
      integer :: status
      mp = matOpen(trim(file_name), 'w')
      if (mp .eq. 0) then
         call ch_assert( 'Can''t open ''data_out.mat'' for writing.')
      end if
      
      var = 'w'
      !if (present(var_name)) var = var_name
      if (trim(var) .eq. 'w') then
         cR3 = pGlobalW
      elseif (trim(var) .eq. 'r') then
         cR3 = pGlobalR
      elseif (trim(var) .eq. 'g') then
         cR3 = pGlobalG
      endif
      !status = matPutVariable(mp, 'w', pGlobalW)
      status = matPutVariable(mp, trim(var), cR3)
      if (status .ne. 0) then
         call ch_assert( 'matPutVariable ''w'' failed')
      end if
      status = matClose(mp)
      if (status .ne. 0) then
         call ch_assert( 'Error closing MAT-file')
      end if
      
      end subroutine saveW
      
      subroutine saveMat2(file_name)
      ! save the PAWS mat data as well as CLM data
      use matRW
      use displayMod, Only: display
      use Vdata, Only: num2str
      use memoryChunk
      implicit none
      integer, parameter :: nField=11
      character*32 :: Fields(nField)
      integer i,j
      
      character(len=*) :: file_name
      character(len=200) :: file_nameCLM
      mwPointer :: mp
      mwSize :: number_of_dim
      mwPointer :: m(6)
      integer :: status, numel


      mp = matOpen(trim(file_name), 'w')
      if (mp .eq. 0) then
         call ch_assert( 'Can''t open ''data_out.mat'' for writing.')
      end if
      status = matPutVariable(mp, 'g', pGlobalG)
      if (status .ne. 0) then
         call ch_assert( 'matPutVariable ''g'' failed')
      end if
      status = matPutVariable(mp, 'w', pGlobalW)
      if (status .ne. 0) then
         call ch_assert( 'matPutVariable ''w'' failed')
      end if
      status = matPutVariable(mp, 'r', pGlobalR)
      if (status .ne. 0) then
         call ch_assert( 'matPutVariable ''r'' failed')
      end if
      status = matClose(mp)
      if (status .ne. 0) then
         call ch_assert( 'Error closing MAT-file')
      end if
      
      !!!! SAVING CLM DATA
      i = index(file_name, '.')
      file_nameCLM =' '
      file_nameCLM(1:i-1) = file_name(1:i-1)
      file_nameCLM(i:i+7)='_CLM.mat'
      
      mp = matOpen(trim(file_nameCLM), 'w')
      if (mp .eq. 0) then
         call ch_assert( 'Can''t open ''data_out_CLM.mat'' for writing.')
      end if
      
     Fields = (/ 'L1','L2','L3','L4','L5','R1','R2','R3','R4','R5','cLevelMap'/)
     do i = 1, nField
         
       numel = mxGetNumberOfElements(cPtr(i))
       
       if (numel*8D0<2D9) then
       status = matPutVariable(mp, Fields(i), cPtr(i))
       if (status .ne. 0) then
          call display( 'savemat2:: matPutVariable ',Fields(i),' failed. Status: ',num2str(status))
          if (i .ne. 9) then
          ! sometimes it failed because it is too large
            stop
          endif
       end if
       else
          call display( 'savemat2:: matPutVariable variable exceeds 2GB, skipped: ',Fields(i))
       endif
     enddo
!      status = matPutVariable(mp, 'R3', cR3)
!      if (status .ne. 0) then
!         call display( 'matPutVariable ''R3'' failed')
!         stop
!      end if
!      status = matPutVariable(mp, 'R4', cR4)
!      if (status .ne. 0) then
!         call display( 'matPutVariable ''R4'' failed')
!         stop
!      end if
!      status = matPutVariable(mp, 'L3', cL3)
!      if (status .ne. 0) then
!         call display( 'matPutVariable ''r'' failed')
!         stop
!      end if
!      status = matPutVariable(mp, 'L4', cL4)
!      if (status .ne. 0) then
!         call display( 'matPutVariable ''L4'' failed')
!         stop
!      end if
      status = matClose(mp)
      if (status .ne. 0) then
         call ch_assert( 'Error closing MAT-file')
      end if
      
      end subroutine saveMat2
      
!====================================================================================================
      SUBROUTINE saveVar(varname,mm,ptr,sig,seq)
      use matRW
      use VData
      implicit none
      integer mm,sig
      mwPointer,SAVE ::mp
      CHARACTER*(32) varname
      CHARACTER*(32) fname,vn
      mwPointer :: var_array
      mwPointer :: px
      mwSize :: m, n
      integer*4 :: status,I
      real*8, dimension(mm) :: ptr
      integer*4 :: classid, ComplexFlag,seq
      ComplexFlag = 0
      m = mm
      n = 1
      
      vn = 'saveVar'
      fname = TRIM(vn)//TRIM(num2str(seq))//'.mat'
      var_array = mxCreateDoubleMatrix(m,n,ComplexFlag)
      px = mxGetPr(var_array)
      call mxCopyReal8ToPtr(ptr, px, m)
      if (sig .eq. 0) THEN
        mp = matOpen(fname, 'w')
        if (mp .eq. 0) then
           call ch_assert( 'Can''t open ''var.mat'' for writing.')
        end if
        sig = 1
      endif
      status = matPutVariable(mp, TRIM(varname), var_array)
      if (status .ne. 0) then
         call ch_assert( 'matPutVariable failed')
      end if
      if (sig .eq. 2 .or. seq .eq. 0) then
        status = matClose(mp)
        if (status .ne. 0) then
          call ch_assert( 'Error closing MAT-file')
        end if
      endif        
      call mxDestroyArray(var_array) ! otherwise causes memory leak
      
      END subroutine saveVar


      SUBROUTINE tempSave(seq)
      use matRW
      use Vdata
      implicit none
      integer*4 :: seq, ComplexFlag, status, I, classID
      mwPointer, save :: mp
      CHARACTER*(32) fname, vn
      mwPointer :: var_array, px
      mwSize :: m, n
      ComplexFlag = 0
      
      m = size(g%VDZ%h)
      n = 1

      vn = 'saveVar'
      fname = TRIM(vn)//TRIM(num2str(seq))//'.mat'

      var_array = mxCreateDoubleMatrix(m,n,ComplexFlag) 
      px = mxGetPr(var_array)

      mp = matOpen(fname, 'w')
      if (mp .eq. 0) then
        call ch_assert( 'Can''t open ''var.mat'' for writing.')
      end if

      if(seq == 0) then
        call mxCopyReal8ToPtr(g%VDZ%N, px, m)
        status = matPutVariable(mp, "N", var_array)
        if (status .ne. 0) then
          call display( 'matPutVariable failed')
          stop
        endif

        call mxCopyReal8ToPtr(g%VDZ%ALPHA, px, m)
        status = matPutVariable(mp, "ALPHA", var_array)
        if (status .ne. 0) then
          call ch_assert( 'matPutVariable failed')
        endif

        call mxCopyReal8ToPtr(g%VDZ%LAMBDA, px, m)
        status = matPutVariable(mp, "LAMBDA", var_array)
        if (status .ne. 0) then
          call ch_assert( 'matPutVariable failed')
        endif
      else
        call mxCopyReal8ToPtr(g%VDZ%h, px, m)
        status = matPutVariable(mp, "h", var_array)
        if (status .ne. 0) then
          call ch_assert( 'matPutVariable failed')
        endif

        call mxCopyReal8ToPtr(g%VDZ%THE, px, m)
        status = matPutVariable(mp, "THE", var_array)
        if (status .ne. 0) then
          call ch_assert( 'matPutVariable failed')
        endif

        call mxCopyReal8ToPtr(g%VDZ%THER, px, m)
        status = matPutVariable(mp, "THER", var_array)
        if (status .ne. 0) then
          call ch_assert( 'matPutVariable failed')
        endif

        call mxCopyReal8ToPtr(g%VDZ%THES, px, m)
        status = matPutVariable(mp, "THES", var_array)
        if (status .ne. 0) then
          call ch_assert( 'matPutVariable failed')
        endif

        call mxCopyReal8ToPtr(g%VDZ%KS, px, m)
        status = matPutVariable(mp, "KS", var_array)
        if (status .ne. 0) then
          call ch_assert( 'matPutVariable failed')
        endif
      endif
      
      status = matClose(mp)
      if (status .ne. 0) then
        call ch_assert( 'Error closing MAT-file')
      end if
        
      call mxDestroyArray(var_array) ! otherwise causes memory leak

      END SUBROUTINE tempSave
