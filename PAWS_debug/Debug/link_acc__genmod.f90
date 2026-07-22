        !COMPILER-GENERATED INTERFACE MODULE: Sun Apr 09 23:47:54 2017
        MODULE LINK_ACC__genmod
          INTERFACE 
            SUBROUTINE LINK_ACC(M,QOC,QOUT,TRIB)
              INTEGER(KIND=4) :: M(3)
              REAL(KIND=8) ,TARGET :: QOC(M((2)))
              REAL(KIND=8) ,TARGET :: QOUT(M((2)))
              REAL(KIND=8) ,TARGET :: TRIB(M((2)))
            END SUBROUTINE LINK_ACC
          END INTERFACE 
        END MODULE LINK_ACC__genmod
