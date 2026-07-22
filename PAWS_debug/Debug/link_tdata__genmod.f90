        !COMPILER-GENERATED INTERFACE MODULE: Sun Apr 09 23:47:55 2017
        MODULE LINK_TDATA__genmod
          INTERFACE 
            SUBROUTINE LINK_TDATA(DC,H1SAVE,H2SAVE,OBJFUN,M3R,M3E,DIMS, &
     &NFIELD)
              INTEGER(KIND=4) :: NFIELD
              INTEGER(KIND=4) :: DIMS(3,NFIELD)
              REAL(KIND=8) ,TARGET :: DC(DIMS((1,1)),DIMS((2,1)))
              REAL(KIND=8) ,TARGET :: H1SAVE(DIMS((1,2)),DIMS((2,2)))
              REAL(KIND=8) ,TARGET :: H2SAVE(DIMS((1,3)),DIMS((2,3)))
              REAL(KIND=8) ,TARGET :: OBJFUN
              REAL(KIND=8) ,TARGET :: M3R(DIMS((1,5)),DIMS((2,5)))
              REAL(KIND=8) ,TARGET :: M3E(DIMS((1,6)),DIMS((2,6)))
            END SUBROUTINE LINK_TDATA
          END INTERFACE 
        END MODULE LINK_TDATA__genmod
