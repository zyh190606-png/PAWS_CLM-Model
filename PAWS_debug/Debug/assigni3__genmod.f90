        !COMPILER-GENERATED INTERFACE MODULE: Sun Apr 09 23:47:43 2017
        MODULE ASSIGNI3__genmod
          INTERFACE 
            SUBROUTINE ASSIGNI3(NN,DATND,MAIN,N)
              USE GD_MOD
              INTEGER(KIND=4) :: NN(3)
              INTEGER(KIND=4) ,TARGET :: DATND(NN((1)),NN((2)),NN((3)))
              TYPE (GD) :: MAIN
              INTEGER(KIND=4) :: N
            END SUBROUTINE ASSIGNI3
          END INTERFACE 
        END MODULE ASSIGNI3__genmod
