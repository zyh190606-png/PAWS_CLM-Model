        !COMPILER-GENERATED INTERFACE MODULE: Sun Apr 09 23:47:55 2017
        MODULE LINK_ROBJ__genmod
          INTERFACE 
            SUBROUTINE LINK_ROBJ(TYPE,ID,DIST,DDIST,XY,IX,SIDX,OBJNUM,  &
     &ACC,ACCT,V,CC,WEIGHT,DIMS,NFIELDS,IS)
              INTEGER(KIND=4) :: NFIELDS
              INTEGER(KIND=4) :: DIMS(3,NFIELDS)
              CHARACTER(LEN=1) ,TARGET :: TYPE
              INTEGER(KIND=4) ,TARGET :: ID
              REAL(KIND=8) ,TARGET :: DIST
              REAL(KIND=8) ,TARGET :: DDIST
              REAL(KIND=8) ,TARGET :: XY(DIMS((2,5)))
              INTEGER(KIND=4) ,TARGET :: IX(DIMS((2,6)))
              INTEGER(KIND=4) ,TARGET :: SIDX
              INTEGER(KIND=4) ,TARGET :: OBJNUM
              REAL(KIND=8) ,TARGET :: ACC(DIMS((2,9)))
              REAL(KIND=8) ,TARGET :: ACCT
              REAL(KIND=8) ,TARGET :: V(DIMS((2,11)))
              REAL(KIND=8) ,TARGET :: CC
              REAL(KIND=8) ,TARGET :: WEIGHT
              INTEGER(KIND=4) :: IS
            END SUBROUTINE LINK_ROBJ
          END INTERFACE 
        END MODULE LINK_ROBJ__genmod
