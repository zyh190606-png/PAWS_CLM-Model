        !COMPILER-GENERATED INTERFACE MODULE: Sun Apr 09 23:47:55 2017
        MODULE LINKDM__genmod
          INTERFACE 
            SUBROUTINE LINKDM(D,MSIZE,NGHO,MAP,MAPL,NXY,DIM,ORIGIN,     &
     &ORIGIN_IDX,X,Y,PBOUND,CBOUND,EBOUND,BINDEX,MASK,MASK2,DIMS,NFIELDS&
     &,IS,OBJ)
              INTEGER(KIND=4) :: NFIELDS
              INTEGER(KIND=4) :: DIMS(3,NFIELDS)
              REAL(KIND=8) ,TARGET :: D(DIMS((2,1)))
              REAL(KIND=8) ,TARGET :: MSIZE(DIMS((2,2)))
              REAL(KIND=8) ,TARGET :: NGHO(DIMS((2,3)))
              REAL(KIND=8) ,TARGET :: MAP(DIMS((1,4)),DIMS((2,4)))
              INTEGER(KIND=4) ,TARGET :: MAPL(DIMS((1,5)),DIMS((2,5)))
              REAL(KIND=8) ,TARGET :: NXY(DIMS((2,6)))
              REAL(KIND=8) ,TARGET :: DIM(DIMS((2,7)))
              REAL(KIND=8) ,TARGET :: ORIGIN(DIMS((2,8)))
              REAL(KIND=8) ,TARGET :: ORIGIN_IDX(DIMS((2,9)))
              REAL(KIND=8) ,TARGET :: X(DIMS((2,10)))
              REAL(KIND=8) ,TARGET :: Y(DIMS((2,11)))
              REAL(KIND=8) ,TARGET :: PBOUND(DIMS((2,12)))
              REAL(KIND=8) ,TARGET :: CBOUND(DIMS((2,13)))
              REAL(KIND=8) ,TARGET :: EBOUND(DIMS((2,14)))
              REAL(KIND=8) ,TARGET :: BINDEX(DIMS((2,15)))
              LOGICAL(KIND=1) ,TARGET :: MASK(DIMS((1,16)),DIMS((2,16)))
              REAL(KIND=8) ,TARGET :: MASK2(DIMS((1,17)),DIMS((2,17)))
              INTEGER(KIND=4) :: IS
              CHARACTER(*) :: OBJ
            END SUBROUTINE LINKDM
          END INTERFACE 
        END MODULE LINKDM__genmod
