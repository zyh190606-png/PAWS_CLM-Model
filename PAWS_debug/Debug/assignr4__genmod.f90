        !COMPILER-GENERATED INTERFACE MODULE: Sun Apr 09 23:47:43 2017
        MODULE ASSIGNR4__genmod
          INTERFACE 
            SUBROUTINE ASSIGNR4(NN,DATND,MAIN,N)
              USE GD_MOD
              INTEGER(KIND=4) :: NN(4)
              REAL(KIND=8) ,TARGET :: DATND(NN((1)),NN((2)),NN((3)),NN((&
     &4)))
              TYPE (GD) :: MAIN
              INTEGER(KIND=4) :: N
            END SUBROUTINE ASSIGNR4
          END INTERFACE 
        END MODULE ASSIGNR4__genmod
