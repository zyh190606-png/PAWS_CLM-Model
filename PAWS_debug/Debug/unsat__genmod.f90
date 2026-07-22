        !COMPILER-GENERATED INTERFACE MODULE: Sun Apr 09 23:48:22 2017
        MODULE UNSAT__genmod
          INTERFACE 
            SUBROUTINE UNSAT(MAPP,V,C,HB,OPT)
              USE VDATA
              TYPE (VDZ_TYPE) ,TARGET :: V
              REAL(KIND=8) :: MAPP(V%M((1)),V%M((2)))
              REAL(KIND=8) :: C(V%M((1)),V%M((2)))
              REAL(KIND=8) :: HB(V%M((1)),V%M((2)))
              INTEGER(KIND=4) ,OPTIONAL :: OPT
            END SUBROUTINE UNSAT
          END INTERFACE 
        END MODULE UNSAT__genmod
