        !COMPILER-GENERATED INTERFACE MODULE: Sun Apr 09 23:47:43 2017
        MODULE ASSIGNI4__genmod
          INTERFACE 
            SUBROUTINE ASSIGNI4(NN,DATND,MAIN,N)
              USE GD_MOD
              INTEGER(KIND=4) :: NN(4)
              INTEGER(KIND=4) ,TARGET :: DATND(NN((1)),NN((2)),NN((3)), &
     &NN((4)))
              TYPE (GD) :: MAIN
              INTEGER(KIND=4) :: N
            END SUBROUTINE ASSIGNI4
          END INTERFACE 
        END MODULE ASSIGNI4__genmod
