!

      subroutine SaveDebug()
      use matRW
      use Vdata
      implicit none

      real*8, pointer:: wr1(:)=>NULL(),wr2(:,:)=>NULL(),wr3(:,:,:)=>NULL(),writeout(:)=>NULL()
      integer :: sig, nw, I
      integer :: m(3), dsig
      
      nw = 1
      DO I=1,nw
        wr1 => r(25)%Riv%h; wr2 => g%GW(1)%h; wr3 => g%VDZ%h;
        m = w%m(1:3);
        IF (I .eq. 1) dsig = 0
        writeout => oneD_ptr(size(wr1), wr1)
        call saveVar('rh',size(wr1),writeout, dsig)
        writeout => oneD_ptr(size(wr2), wr2)
        call saveVar('gh',size(wr2),writeout, dsig)
        writeout => oneD_ptr(size(wr3), wr3)
        IF (I .eq. nw) dsig = 2
        call saveVar('vh',size(wr3),writeout, dsig)
      ENDDO
      
      end subroutine SaveDebug        
