        !COMPILER-GENERATED INTERFACE MODULE: Sun Apr 09 23:47:44 2017
        MODULE ASSIGNL2__genmod
          INTERFACE 
            SUBROUTINE ASSIGNL2(NN,DATND,MAIN,N)
              USE GD_MOD
              INTEGER(KIND=4) :: NN(2)
              LOGICAL(KIND=1) ,TARGET :: DATND(NN((1)),NN((2)))
              TYPE (GD) :: MAIN
              INTEGER(KIND=4) :: N
            END SUBROUTINE ASSIGNL2
          END INTERFACE 
        END MODULE ASSIGNL2__genmod
