        !COMPILER-GENERATED INTERFACE MODULE: Sun Apr 09 23:47:44 2017
        MODULE ASSIGNL5__genmod
          INTERFACE 
            SUBROUTINE ASSIGNL5(NN,DATND,MAIN,N)
              USE GD_MOD
              INTEGER(KIND=4) :: NN(5)
              LOGICAL(KIND=1) ,TARGET :: DATND(NN((1)),NN((2)),NN((3)), &
     &NN((4)),NN((5)))
              TYPE (GD) :: MAIN
              INTEGER(KIND=4) :: N
            END SUBROUTINE ASSIGNL5
          END INTERFACE 
        END MODULE ASSIGNL5__genmod
