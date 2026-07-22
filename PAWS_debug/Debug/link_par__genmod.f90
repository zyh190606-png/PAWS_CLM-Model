        !COMPILER-GENERATED INTERFACE MODULE: Sun Apr 09 23:47:55 2017
        MODULE LINK_PAR__genmod
          INTERFACE 
            SUBROUTINE LINK_PAR(ID,IMET,CODE,MAPP,V,UL,LL,DIMS,NFIELDS, &
     &INDEX)
              INTEGER(KIND=4) :: NFIELDS
              INTEGER(KIND=4) :: DIMS(3,NFIELDS)
              INTEGER(KIND=4) ,TARGET :: ID
              INTEGER(KIND=4) ,TARGET :: IMET
              INTEGER(KIND=4) ,TARGET :: CODE
              INTEGER(KIND=4) ,TARGET :: MAPP(DIMS((1,4)))
              REAL(KIND=8) ,TARGET :: V
              REAL(KIND=8) ,TARGET :: UL
              REAL(KIND=8) ,TARGET :: LL
              INTEGER(KIND=4) :: INDEX
            END SUBROUTINE LINK_PAR
          END INTERFACE 
        END MODULE LINK_PAR__genmod
