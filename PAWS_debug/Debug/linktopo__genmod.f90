        !COMPILER-GENERATED INTERFACE MODULE: Sun Apr 09 23:47:55 2017
        MODULE LINKTOPO__genmod
          INTERFACE 
            SUBROUTINE LINKTOPO(MS,E,EX,EY,S0,S0X,S0Y)
              INTEGER(KIND=4) :: MS(2)
              REAL(KIND=8) ,TARGET :: E(MS((1)),MS((2)))
              REAL(KIND=8) ,TARGET :: EX(MS((1)),MS((2)))
              REAL(KIND=8) ,TARGET :: EY(MS((1)),MS((2)))
              REAL(KIND=8) ,TARGET :: S0(MS((1)),MS((2)))
              REAL(KIND=8) ,TARGET :: S0X(MS((1)),MS((2)))
              REAL(KIND=8) ,TARGET :: S0Y(MS((1)),MS((2)))
            END SUBROUTINE LINKTOPO
          END INTERFACE 
        END MODULE LINKTOPO__genmod
