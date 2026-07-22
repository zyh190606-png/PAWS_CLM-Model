        !COMPILER-GENERATED INTERFACE MODULE: Sun Apr 09 23:47:43 2017
        MODULE ASSIGNI5__genmod
          INTERFACE 
            SUBROUTINE ASSIGNI5(NN,DATND,MAIN,N)
              USE GD_MOD
              INTEGER(KIND=4) :: NN(5)
              INTEGER(KIND=4) ,TARGET :: DATND(NN((1)),NN((2)),NN((3)), &
     &NN((4)),NN((5)))
              TYPE (GD) :: MAIN
              INTEGER(KIND=4) :: N
            END SUBROUTINE ASSIGNI5
          END INTERFACE 
        END MODULE ASSIGNI5__genmod
