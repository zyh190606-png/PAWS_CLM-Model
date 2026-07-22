        !COMPILER-GENERATED INTERFACE MODULE: Sun Apr 09 23:48:25 2017
        MODULE COMPTS__genmod
          INTERFACE 
            SUBROUTINE COMPTS(TOBS,TSIM,SD,ED,OPT,STATS)
              USE VDATA
              TYPE (SIM_TYPE) :: TOBS
              TYPE (SIM_TYPE) :: TSIM
              INTEGER(KIND=4) :: SD
              INTEGER(KIND=4) :: ED
              INTEGER(KIND=4) :: OPT
              REAL(KIND=8) :: STATS(4)
            END SUBROUTINE COMPTS
          END INTERFACE 
        END MODULE COMPTS__genmod
