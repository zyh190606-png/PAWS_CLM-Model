        !COMPILER-GENERATED INTERFACE MODULE: Sun Apr 09 23:48:18 2017
        MODULE READPAR__genmod
          INTERFACE 
            SUBROUTINE READPAR(FILE,PAR,NPAR)
              USE VDATA
              CHARACTER(*) :: FILE
              TYPE (PAR_TYPE) :: PAR(100)
              INTEGER(KIND=4) :: NPAR
            END SUBROUTINE READPAR
          END INTERFACE 
        END MODULE READPAR__genmod
