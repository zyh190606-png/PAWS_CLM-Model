        !COMPILER-GENERATED INTERFACE MODULE: Sun Apr 09 23:47:43 2017
        MODULE ASSIGNR2__genmod
          INTERFACE 
            SUBROUTINE ASSIGNR2(NN,DATND,MAIN,N)
              USE GD_MOD
              INTEGER(KIND=4) :: NN(2)
              REAL(KIND=8) ,TARGET :: DATND(NN((1)),NN((2)))
              TYPE (GD) :: MAIN
              INTEGER(KIND=4) :: N
            END SUBROUTINE ASSIGNR2
          END INTERFACE 
        END MODULE ASSIGNR2__genmod
