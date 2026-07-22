        !COMPILER-GENERATED INTERFACE MODULE: Sun Apr 09 23:47:43 2017
        MODULE ASSIGNR1__genmod
          INTERFACE 
            SUBROUTINE ASSIGNR1(NN,DATND,MAIN,N)
              USE GD_MOD
              INTEGER(KIND=4) :: NN(1)
              REAL(KIND=8) ,TARGET :: DATND(NN((1)))
              TYPE (GD) :: MAIN
              INTEGER(KIND=4) :: N
            END SUBROUTINE ASSIGNR1
          END INTERFACE 
        END MODULE ASSIGNR1__genmod
