module exportASCII
    ! a simple ASCII exporter
    use vdata, Only: oneD_ptr
    integer :: matUnit=5871, tsUnit=5872
    integer, parameter :: nM = 4
    public:: writeDimMat

      interface writeDimMat         ! call getPtr(main,field,ptr,[num])
       !writeDimMat(dat,name,time)
       !dat: 1D-4D matrices. 
       !name: a string denoting the name of the variable. could be sth like 'g.OVN.h'; 
       !time: a timestamp for the variable
      module procedure writeDimMat1,writeDimMat2,writeDimMat3,writeDimMat4
      end interface
      
    contains
    
    subroutine writeDimMat1m(dat,name,t,m)
    use vdata, Only: num2str
    real*8 dat(:)
    real*8 :: t
    integer m(nM)
    character*(*) :: name
    integer nr,nc
    logical ok
    character*32 :: fmt
    nr = size(dat)
    inquire(unit = matUnit, OPENED = ok)
    if (.not. ok) open(unit = matUnit, file='matUnit.dat', status='unknown')
    fmt = trim('('// trim(num2str(nr))//'e13.5)')
    write(matUnit,'(A1,A32,e13.5,5(i5))') '#',name,t,m
    !write(matUnit, '(<nr>e23.16)') dat ! (dat(i),i=1,nr)
    write(matUnit, fmt) dat 
    end subroutine
    
    subroutine writeDimMat1(dat,name,t)
    real*8,target :: dat(:)
    real*8,pointer:: dat_p(:)
    real*8, optional :: t
    character*(*) :: name
    integer :: m(nM)=1
    m(1) = size(dat)
    call writeDimMat1m(dat,name,t,m)
    end subroutine
    
    subroutine writeDimMat2(dat,name,t)
    real*8,target :: dat(:,:)
    real*8,pointer:: dat_p(:)
    real*8, optional :: t
    character*(*) :: name
    integer :: m(nM)=1
    integer :: nd=2
    do i=1,nd
        m(i)=size(dat,i)
    enddo    
    dat_p => oneD_ptr(size(dat),dat)
    call writeDimMat1m(dat_p,name,t,m)
    end subroutine
    
    subroutine writeDimMat3(dat,name,t)
    real*8,target :: dat(:,:,:)
    real*8,pointer:: dat_p(:)
    real*8, optional :: t
    character*(*) :: name
    integer :: m(nM)=1
    integer :: nd=3
    do i=1,nd
        m(i)=size(dat,i)
    enddo    
    dat_p => oneD_ptr(size(dat),dat)
    call writeDimMat1m(dat_p,name,t,m)
    end subroutine
    
    subroutine writeDimMat4(dat,name,t)
    real*8,target :: dat(:,:,:,:)
    real*8,pointer:: dat_p(:)
    real*8, optional :: t
    character*(*) :: name
    integer :: m(nM)=1
    integer :: nd=4
    do i=1,nd
        m(i)=size(dat,i)
    enddo    
    dat_p => oneD_ptr(size(dat),dat)
    call writeDimMat1m(dat_p,name,t,m)
    end subroutine
    
end module
      
!====================================================================================
      subroutine find1DL(array,m,filter,n)
      ! this is an emulation of matlab find function
      ! this function finds the 1D indices of 'array' where value is > 0
      ! returns a pointer allocated with indices
      implicit none
      logical array(m)
      integer i,k,m,n
      integer filter(m)
      
      k = 0
      do i = 1,m
        if (array(i)) then
          k = k + 1
          filter(k) = i
        endif
      enddo
      n = k

      END subroutine
      
      
      function digit(n)
      integer n, digit
      if (abs(n) < 10) then
        digit = 0;
      else
        do i=1,100
          if (abs(n) < 10**(I+1) .and. abs(n) >= 10**(I) ) then
            digit = I
            EXIT
          endif
        enddo
      endif
      end function
      
!====================================================================================
    module polyEval
    private
    public:: polyIntVal,polyIntWeights,gauss,polyIntIntpInvMat,polyInterp1,linearInterp1d
    contains
    
    subroutine linearInterp1d(x,y,x1,y1,loc,lastLoc)
    real*8, optional :: lastLoc ! good for a loop of calls with monotonically ascending search position
    real*8, intent(in) :: x(:), y(:), x1
    real*8, intent(out) :: y1, loc
    real*8 d, d1
    integer n, i,j, k, i1
    n = size(x)
    if (present(lastLoc)) then
        i1 = max(int(lastLoc),1)
    else
        i1 = 1
    endif
    if (x1<x(1)) then
        y1=y(1)
        return
    elseif (x1>x(n)) then
        y1=y(n)
        return
    else
      do i=i1, n-1
        if (x(i)<=x1 .and. x(i+1)>x1) then
           d = x(i+1)-x(i)
           d1= x1 - x(i)
           loc = d1/d + i
           y1= y(i)+ (loc-int(loc))*(y(i+1)-y(i))
           exit
        endif
      enddo
    endif  
    end subroutine
    
    subroutine polyIntWeights(n,stencil, mode, A)
    implicit none
      integer n,i,j,k,mode
      real*8 stencil(n+1)
      real*8 A(n,n),coeff(n)
      
      
      do i =1,n
        do j=1,n
            A(i,j) = (stencil(i+1)**(n+1D00-j) - stencil(i)**(n+1D0-j) )/(n+1D0-j)
        enddo    
      enddo
    end subroutine
    
    subroutine polyInterp1(m,t,v,t1,res,mode,n,MX,state)
    ! assuming uniform time step of t
    implicit none
    integer m,mode,n,state
    integer k1,k2,i,j
    real*8, dimension(m):: t,v
    real*8 t1(2) ! interpolate to an interval. if only to a point, the second value can be whatever
    real*8,dimension(4),target:: stencilF3=(/-1D0, 0D0, 1D0, 2D0 /),stencilC3=(/-1.5D0, -0.5D0, 0.5D0, 1.5D0 /)
    real*8,dimension(3),target:: stencilF2=(/0D0, 1D0, 2D0 /),stencilCF2=(/-0.5D0, 0.5D0, 1.5D0 /),stencilC2=(/-1D0, 0D0, 1D0 /)
    real*8,dimension(:),pointer::stencil
    real*8 :: MX(n,n), res
    real*8 :: interval(2),vp(10),P(n)
    real*8 :: PREC=1D-15,dt
    
    if (state .eq. 0) then
      select case(mode)
        case (0) !
            if (mode .eq. 0) then  ! each time stamp means the interval forward. i.e., usgs type time stamp
              if (n .eq. 2) then
                stencil => stencilF2
              elseif (n .eq. 3) then
                stencil => stencilF3
              endif
            else
              if (n .eq. 2) then
                stencil => stencilC2
              elseif (n .eq. 3) then
                stencil => stencilC3
              endif
            endif
            call polyIntIntpInvMat(n,stencil,MX)
            state = 1
      end select
    endif
    
    dt = t(2)-t(1)
    
    if (mode .eq. 0) then  ! each time stamp means the interval forward. i.e., usgs type time stamp
        if (n .eq. 2) then
        stencil => stencilF2
        elseif (n .eq. 3) then
        stencil => stencilF3
        endif
    else
        if (n .eq. 2) then
        stencil => stencilC2
        elseif (n .eq. 3) then
        stencil => stencilC3
        endif
    endif
    if (mode .eq. 0) then ! each value means integral of the forward interval. interpolate to an interval
      k1 = int((t1(1) + PREC - t(1))/dt ) + 1 ! base
      k2 = int((t1(2) + PREC - t(1))/dt ) + 1 !
      if (k2-k1> maxval(stencil)) then
        write(*,*) 'polyInterp1:: this portion not coded yet'
        stop
      else
        interval(1) = max(t1(1) - t(K1),0D0)/dt
        interval(2) = max(t1(2) - t(K1),0D0)/dt
        do i=1,size(stencil)
          vp(i) = v(min(max(1,k1+int(stencil(i))),size(v)))
        enddo
        P = matmul(MX,vp(1:n))
        res = polyIntVal(n,P,interval)
      endif
    end if
    
    end subroutine
    
    function polyIntVal(n,P,interval)
    ! integrate polynomial P of (n-1)-th order (so n coefficients) in interval
    implicit none
    integer n
    real*8 P(n)
    real*8 interval(2)
    real*8 polyIntVal
    
    integer i,j,k
    polyIntVal = 0D0
    do j=1,n
        polyIntVal = polyIntVal+P(j)*(interval(2)**(n+1-j) - interval(1)**(n+1-j) )/(n+1-j)
    enddo    
    end function
    
    subroutine polyIntIntpInvMat(n,stencil,MX)
    ! polynomial integral interpolation inverse matrix
    ! (n-1)-th order polynomal integral interpolation
    ! giving n integrals over stencil intervals
    ! construct n-th polynomial to be evaluated in interval
    ! here the inverse matrix is passed out
    ! to obtain the polynomial, do matmul(MX,integrals)
    ! then for integral evaluation,  res =polyIntVal
    integer n,mode
    real*8:: stencil(n+1),ints(n),interval(2),res
     real*8 BB(N,N),MX(N,N),t1
    integer::INDX(N)
    !real*8, intent(out)::P(3)
    
     call polyIntWeights(n,stencil,1,BB)
     if (n > 2) then
       call MIGS(BB,3,MX,INDX)
     else
       t1 = BB(1,1)*BB(2,2)-BB(1,2)*BB(2,1)
       MX(1,:)= (/BB(2,2)/t1, -BB(1,2)/t1 /) 
       MX(2,:)= (/-BB(2,1)/t1, BB(1,1)/t1 /) 
     endif
     ![  b22/(b11*b22 - b12*b21), -b12/(b11*b22 - b12*b21)]
     ![ -b21/(b11*b22 - b12*b21),  b11/(b11*b22 - b12*b21)]
     !MIGS(A,N,X,INDX)
     !call  gauss(3,BB,b,x)
     !P = matmul(MX,b)
     !res = polyIntVal(3,x,(/0D0,1D0/))
    endsubroutine
    

!########## BELOW ARE SOME FILES OBTAINED AS LIBRARIES ###########################
!########## ORIGINAL COPYRIGHT INFO MAINTAINED ###################################
    subroutine gauss(n,a,bi,x)
    implicit none
      integer n
      real*8 a(n,n),s(n),bi(n),b(n)
      real*8 x(n)
      integer p(n), pk,i,j,k
      real*8 smax,rmax,r,z,sum
      
      b= bi
      do i=1,n
         p(i) = i
         smax = 0.0D0 
         do j=1,n
            smax = max(smax,abs(a(i,j)))
         enddo
         s(i) = smax
     enddo
    
      do k=1,n-1
         rmax = 0.0D0
         do i=k,n
            r = abs(a(p(i),k))/s(p(i))
            if (r .gt. rmax) then
               j = i
               rmax = r
            endif
         enddo
         
         pk = p(j)
         p(j) = p(k)
         p(k) = pk
         
         do i=k+1,n      
            z = a(p(i),k)/a(p(k),k)       
            a(p(i),k) = z
            do j=k+1,n    
               a(p(i),j) = a(p(i),j) - z*a(p(k),j)      
            enddo
         enddo
      enddo   
      
      do k=1,n-1
         do i=k+1,n
            b(p(i)) = b(p(i)) - a(p(i),k)*b(p(k))
         enddo
      enddo
    
      do i=n,1,-1
         sum = b(p(i))
         do j=i+1,n
            sum = sum - a(p(i),j)*x(j)
         enddo
         x(i) = sum/a(p(i),i)
      enddo

    end subroutine
    
    
    ! Two subroutines below were modified from Tao Pang's code
    ! original statements retained below
      SUBROUTINE MIGS(A,N,X,INDX)
      IMPLICIT REAL*8(a-h,o-z)
      IMPLICIT INTEGER(I-N)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c                                                                      c
!c Please Note:                                                         c
!c                                                                      c
!c (1) This computer program is part of the book, "An Introduction to   c
!c     Computational Physics," written by Tao Pang and published and    c
!c     copyrighted by Cambridge University Press in 1997.               c
!c                                                                      c
!c (2) No warranties, express or implied, are made for this program.    c
!c                                                                      c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      REAL*8:: A(N,N),X(N,N),B(N,N)
      integer::INDX(N)
! Subroutine to invert matrix A(N,N) with the inverse stored
! in X(N,N) in the output.
      DO    I = 1, N
        DO     J = 1, N
          B(I,J) = 0.0D0
        ENDDO
      ENDDO
      DO I = 1, N
          B(I,I) = 1.0D0
      ENDDO
      CALL ELGS(A,N,INDX)
      DO I = 1, N-1
        DO  J = I+1, N
          DO   K = 1, N
            B(INDX(J),K) = B(INDX(J),K)-A(INDX(J),I)*B(INDX(I),K)
          ENDDO
        ENDDO
      ENDDO
      DO   I = 1, N
        X(N,I) = B(INDX(N),I)/A(INDX(N),N)
        DO  J = N-1, 1, -1
          X(J,I) = B(INDX(J),I)
          DO K = J+1, N
            X(J,I) = X(J,I)-A(INDX(J),K)*X(K,I)
          ENDDO
          X(J,I) =  X(J,I)/A(INDX(J),J)
        ENDDO
      ENDDO
      RETURN
      ENDsubroutine
      
      SUBROUTINE ELGS(A,N,INDX)
      IMPLICIT REAL*8(a-h,o-z)
      IMPLICIT INTEGER(I-N)
 !Subroutine to perform the partial-pivoting Gaussian elimination.
 !A(N,N) is the original matrix in the input and transformed
 !matrix plus the pivoting element ratios below the diagonal in
 !the output.  INDX(N) records the pivoting order.
      REAL*8:: A(N,N),C(N)
      REAL*8 C1,PI
      integer K,I,j
      INTEGER:: INDX(N)
      DO I = 1, N
        INDX(I) = I
      ENDDO

  ! Find the rescaling factors, one from each row
        DO  I = 1, N
          C1= 0.0D0
          DO J = 1, N
            C1 = max(C1,ABS(A(I,J)))
          ENDDO
          C(I) = C1
        ENDDO

 !C Search the pivoting (largest) element from each column
      DO  J = 1, N-1
        PI1 = 0.0
        DO   I = J, N
          PI = ABS(A(INDX(I),J))/C(INDX(I))
          IF (PI.GT.PI1) THEN
            PI1 = PI
            K   = I
          ELSE
          ENDIF
        ENDDO
 !Interchange the rows via INDX(N) to record pivoting order

        ITMP    = INDX(J)
        INDX(J) = INDX(K)
        INDX(K) = ITMP
        DO   I = J+1, N
          PJ  = A(INDX(I),J)/A(INDX(J),J)

 !Record pivoting ratios below the diagonal

          A(INDX(I),J) = PJ

 !Modify other elements accordingly

          DO  K = J+1, N
            A(INDX(I),K) = A(INDX(I),K)-PJ*A(INDX(J),K)
          ENDDO
        ENDDO
      ENDDO
      RETURN
      ENDSubroutine
    
    end module