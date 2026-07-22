        !COMPILER-GENERATED INTERFACE MODULE: Sun Apr 09 23:47:55 2017
        MODULE LINKVEGDAY__genmod
          INTERFACE 
            SUBROUTINE LINKVEGDAY(M,LAI,ROOT,HC,KC)
              INTEGER(KIND=4) :: M(3)
              REAL(KIND=8) ,TARGET :: LAI(M((1)))
              REAL(KIND=8) ,TARGET :: ROOT(M((1)))
              REAL(KIND=8) ,TARGET :: HC(M((1)))
              REAL(KIND=8) ,TARGET :: KC(M((1)))
            END SUBROUTINE LINKVEGDAY
          END INTERFACE 
        END MODULE LINKVEGDAY__genmod
