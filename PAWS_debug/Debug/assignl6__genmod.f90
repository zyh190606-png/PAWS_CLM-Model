        !COMPILER-GENERATED INTERFACE MODULE: Sun Apr 09 23:47:44 2017
        MODULE ASSIGNL6__genmod
          INTERFACE 
            SUBROUTINE ASSIGNL6(NN,DATND,MAIN,N)
              USE GD_MOD
              INTEGER(KIND=4) :: NN(6)
              LOGICAL(KIND=1) ,TARGET :: DATND(NN((1)),NN((2)),NN((3)), &
     &NN((4)),NN((5)),NN((6)))
              TYPE (GD) :: MAIN
              INTEGER(KIND=4) :: N
            END SUBROUTINE ASSIGNL6
          END INTERFACE 
        END MODULE ASSIGNL6__genmod
