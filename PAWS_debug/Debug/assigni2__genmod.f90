        !COMPILER-GENERATED INTERFACE MODULE: Sun Apr 09 23:47:43 2017
        MODULE ASSIGNI2__genmod
          INTERFACE 
            SUBROUTINE ASSIGNI2(NN,DATND,MAIN,N)
              USE GD_MOD
              INTEGER(KIND=4) :: NN(2)
              INTEGER(KIND=4) ,TARGET :: DATND(NN((1)),NN((2)))
              TYPE (GD) :: MAIN
              INTEGER(KIND=4) :: N
            END SUBROUTINE ASSIGNI2
          END INTERFACE 
        END MODULE ASSIGNI2__genmod
