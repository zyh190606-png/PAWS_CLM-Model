        !COMPILER-GENERATED INTERFACE MODULE: Sun Apr 09 23:47:55 2017
        MODULE LINK_USGS__genmod
          INTERFACE 
            SUBROUTINE LINK_USGS(ID,DATES,V,DATENUMS,T,SD,ED,DIMS,      &
     &NFIELDS,INDEX)
              INTEGER(KIND=4) :: NFIELDS
              INTEGER(KIND=4) :: DIMS(3,NFIELDS)
              REAL(KIND=8) ,TARGET :: ID
              REAL(KIND=8) ,TARGET :: DATES(DIMS((1,3)))
              REAL(KIND=8) ,TARGET :: V(DIMS((1,4)))
              REAL(KIND=8) ,TARGET :: DATENUMS(DIMS((1,5)))
              REAL(KIND=8) ,TARGET :: T(DIMS((1,6)))
              REAL(KIND=8) ,TARGET :: SD
              REAL(KIND=8) ,TARGET :: ED
              INTEGER(KIND=4) :: INDEX
            END SUBROUTINE LINK_USGS
          END INTERFACE 
        END MODULE LINK_USGS__genmod
