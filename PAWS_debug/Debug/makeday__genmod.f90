        !COMPILER-GENERATED INTERFACE MODULE: Sun Apr 09 23:48:25 2017
        MODULE MAKEDAY__genmod
          INTERFACE 
            SUBROUTINE MAKEDAY(RFRAC,FRAC,YEAR,MM,JD,DAY)
              USE VDATA
              REAL(KIND=8) :: RFRAC
              REAL(KIND=8) :: FRAC
              INTEGER(KIND=4) :: YEAR
              INTEGER(KIND=4) :: MM
              INTEGER(KIND=4) :: JD
              TYPE (STA_TYPE) :: DAY
            END SUBROUTINE MAKEDAY
          END INTERFACE 
        END MODULE MAKEDAY__genmod
