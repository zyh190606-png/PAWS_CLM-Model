
    !  paws.f90
    !
    !  FUNCTIONS:
    !  paws - Entry point of console application.
    !

    !****************************************************************************
    !
    !  PROGRAM: paws
    !
    !  PURPOSE:  Entry point for the console application.
    !
    !****************************************************************************

    program paws

    use matRW
    use Vdata
    implicit none

    ! Variables
    character(len=length_of_character) ::  file_name
    integer status, nargin,m(3)
    CHARACTER(len=length_of_character)::file,file2,prj

    ! Body of paws
    !file_name = 'data_in.mat'

    call display('version g228 checked out')

    nargin = NARGS( )
    if (nargin > 1) then
        CALL GETARG (1, file2, status)
        file_name = file2
    else
        file_name = 'data_in.mat'
    endif

    file_name = 'Detroit.in'

    call run(file_name)

    !pause

    end program paws
