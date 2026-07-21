!CP This is a dummy module that is created by Chaopeng
MODULE spmdMod
logical, public :: masterproc = .true.
integer :: mpicom = 1
integer :: MPI_REAL8
integer :: MPI_INTEGER
logical :: MPI_LOGICAL
integer :: iam
END MODULE