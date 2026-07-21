      module displayMod 
      integer :: errorHandle = 0  
#if (defined MEX)
      integer :: printCode = 1 ! cannot print to screen in mex environment 
#else      
      integer :: printCode = 0 
#endif
       ! 0: output to default (screen); 1: write to a file (if printFileName not defined anywhere, fort.272 may be automatically used)
       ! 2: print to both screen and file
      integer :: printID = 272  ! No need to change
      
      private
      public :: errorHandle, printCode, printID
      public :: display
      
      interface num2str
        module procedure num2str
        module procedure num2strI
        module procedure num2strI3
      end interface
      interface display
        module procedure displayStr
        module procedure displaySI,displaySII,displaySS,displaySR
      end interface
      
    contains
    
       FUNCTION num2strI(num)
       ! Converting an integer number into a string variable
       ! Limited to 16 digits, use TRIM outside to trim the trailing blanks
       integer I,N
       parameter (N = 16)
       INTEGER,intent(in):: num
       CHARACTER*(N) num2strI
       CHARACTER*(N) str1
       
       write(str1,'(I16)') num
       ! maybe:
       ! write(str1,'(I<ia+ib>)') num
       num2strI = ADJUSTL(str1)
       
       END FUNCTION

       FUNCTION num2strI3(num)
       ! created specifically for 3-element integer vectors (such as those from datevec2)
       integer I,N
       parameter (N = 16)
       INTEGER,intent(in):: num(3)
       CHARACTER*(N) num2strI3
       CHARACTER*(N) str1
       
       write(str1,'(I16 I16 I16)') num(1), num(2), num(3)
       ! maybe:
       ! write(str1,'(I<ia+ib>)') num
       num2strI3 = ADJUSTL(str1)
       
       END FUNCTION
       
       FUNCTION num2str(num)
       ! Converting an real number into a string variable
       ! string is written in its default format in 16 digits
       ! Limited to F20.10 digits, use TRIM outside to trim the trailing blanks
       integer I,N
       parameter (N = 20)
       real*8,intent(in):: num
       CHARACTER*(N) num2str
       CHARACTER*(N) str1,str2
       
       write(str1,'(F20.10)') num
       ! maybe:
       ! write(str1,'(I<ia+ib>)') num
       num2str = ADJUSTL(str1)
       str2    = num2str
       ! make it more compact:
       DO I = N,1,-1
         IF (str2(I:I) .eq. ' ' .OR. str2(I:I) .eq. '0' .OR. str2(I:I) .eq. '.') THEN
           num2str(I:I) = ' '
         ELSE
           EXIT
         ENDIF
       ENDDO
       
       END FUNCTION
       
      subroutine displayStr(msg)
      character*(*) msg
      if (printCode .eq. 0) then
          write(*,*) trim(msg)
      elseif (printCode .eq. 1) then
          write(printID,*) trim(msg)
      elseif (printCode .eq. 2) then
          write(*,*) trim(msg)
          write(printID,*) trim(msg)
      endif
      end subroutine
    
      subroutine displaySI(msg1,i1,msg2,r2,msg3,r3)
      character*(*) msg1
      character*(*), optional :: msg2,msg3
      integer:: i1
      real*8, optional:: r2,r3
      character*1000 msg
      msg = trim(msg1)//' '//trim(num2str(i1))
      if (present(msg2)) msg = trim(msg)//' '//trim(msg2)
      if (present(r2)) msg = trim(msg)//' '//trim(num2str(r2))
      if (present(msg3)) msg = trim(msg)//' '//trim(msg3)
      if (present(r3)) msg = trim(msg)//' '//trim(num2str(r3))
      call displayStr(msg)
      end subroutine
    
      subroutine displaySII(msg1,i1,msg2,i2,msg3,i3)
      character*(*) msg1,msg2
      integer:: i1
      integer:: i2
      integer, optional:: i3
      character*(*), optional :: msg3
      character*1000 msg
      msg = trim(msg1)//' '//trim(num2str(i1))
      msg = msg//' '//trim(msg2)//' '//trim(num2str(i2))
      if (present(msg3)) msg = trim(msg)//' '//trim(msg3)
      if (present(i3)) msg = trim(msg)//' '//trim(num2str(i3))
      call displayStr(msg)
      end subroutine
      
      subroutine displaySS(msg1,msg2,msg3,msg4)
      character*(*) msg1,msg2
      character*(*), optional :: msg3,msg4
      character*1000 msg
      msg = trim(msg1)//' '//trim(msg2)
      if (present(msg3)) msg = trim(msg)//' '//trim(msg3)
      if (present(msg4)) msg = trim(msg)//' '//trim(msg4)
      call displayStr(msg)
      end subroutine
      
      subroutine displaySR(msg1,r1,msg2,r2,msg3,r3)
      character*(*) msg1
      character*(*), optional :: msg2,msg3
      real*8:: r1
      real*8, optional:: r2,r3
      character*1000 msg
      msg = trim(msg1)//' '//trim(num2str(r1))
      if (present(msg2)) msg = trim(msg)//' '//trim(msg2)
      if (present(r2)) msg = trim(msg)//' '//trim(num2str(r2))
      if (present(msg3)) msg = trim(msg)//' '//trim(msg3)
      if (present(r3)) msg = trim(msg)//' '//trim(num2str(r3))
      call displayStr(msg)
      end subroutine
      
    end module displayMod 
    
    
      subroutine CH_ASSERT(msg)
      use displayMod, Only: errorHandle, display
      character*(*) msg
      character*(10) command
      if (errorHandle .eq. 0) then
          ! Bark the message and stop code
          call display(msg)
          call display('enter C to continue; enter any other key to abort')
          command = ' '
          read (*,*) command !PAUSE is obsolete
          if (TRIM(command) .eq. 'C') then
          else
             STOP
          endif
      elseif (errorHandle .eq. 1) then
          ! do nothing, pass quietly and attempt to continue
      elseif (errorHandle .eq. 2) then
          ! print the error and attempt to continue
          call display(msg)
      else
          call display('ch_assert: Unspecified errorHandle')
          call display(msg)
      endif    
    end subroutine
      
    
    