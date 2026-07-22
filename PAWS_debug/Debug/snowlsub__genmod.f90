        !COMPILER-GENERATED INTERFACE MODULE: Sun Apr 09 23:48:24 2017
        MODULE SNOWLSUB__genmod
          INTERFACE 
            SUBROUTINE SNOWLSUB(DT,NSTEP,INPT,SITEV,STATEV,PARAM,IFLAG, &
     &NSTEPDAY,CUMP,CUME,CUMMR,OUTV,TSBACKUP,TAVEBACKUP,NDEPLETIONPOINTS&
     &,DFC)
              INTEGER(KIND=4) :: NDEPLETIONPOINTS
              REAL(KIND=8) :: DT
              INTEGER(KIND=4) :: NSTEP
              REAL(KIND=8) :: INPT(7,1)
              REAL(KIND=8) :: SITEV(*)
              REAL(KIND=8) :: STATEV(*)
              REAL(KIND=8) :: PARAM(*)
              INTEGER(KIND=4) :: IFLAG(*)
              INTEGER(KIND=4) :: NSTEPDAY
              REAL(KIND=8) :: CUMP
              REAL(KIND=8) :: CUME
              REAL(KIND=8) :: CUMMR
              REAL(KIND=8) :: OUTV(*)
              REAL(KIND=8) :: TSBACKUP(*)
              REAL(KIND=8) :: TAVEBACKUP(*)
              REAL(KIND=8) :: DFC(NDEPLETIONPOINTS,2)
            END SUBROUTINE SNOWLSUB
          END INTERFACE 
        END MODULE SNOWLSUB__genmod
