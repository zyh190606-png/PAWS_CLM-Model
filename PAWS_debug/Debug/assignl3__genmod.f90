        !COMPILER-GENERATED INTERFACE MODULE: Sun Apr 09 23:47:44 2017
        MODULE ASSIGNL3__genmod
          INTERFACE 
            SUBROUTINE ASSIGNL3(NN,DATND,MAIN,N)
              USE GD_MOD
              INTEGER(KIND=4) :: NN(3)
              LOGICAL(KIND=1) ,TARGET :: DATND(NN((1)),NN((2)),NN((3)))
              TYPE (GD) :: MAIN
              INTEGER(KIND=4) :: N
            END SUBROUTINE ASSIGNL3
          END INTERFACE 
        END MODULE ASSIGNL3__genmod
