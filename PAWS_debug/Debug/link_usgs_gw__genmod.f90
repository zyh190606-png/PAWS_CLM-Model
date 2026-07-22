        !COMPILER-GENERATED INTERFACE MODULE: Sun Apr 09 23:47:55 2017
        MODULE LINK_USGS_GW__genmod
          INTERFACE 
            SUBROUTINE LINK_USGS_GW(DATES,DTWT,DATENUMS,V,T,XY,INDICES, &
     &ELEV,SD,ED,DIMS,NFIELDS,INDEX)
              INTEGER(KIND=4) :: NFIELDS
              INTEGER(KIND=4) :: DIMS(3,NFIELDS)
              REAL(KIND=8) ,TARGET :: DATES(DIMS((1,1)))
              REAL(KIND=8) ,TARGET :: DTWT(DIMS((1,2)))
              REAL(KIND=8) ,TARGET :: DATENUMS(DIMS((1,3)))
              REAL(KIND=8) ,TARGET :: V(DIMS((1,4)))
              REAL(KIND=8) ,TARGET :: T(DIMS((1,5)))
              REAL(KIND=8) ,TARGET :: XY(DIMS((2,6)))
              REAL(KIND=8) ,TARGET :: INDICES(DIMS((2,7)))
              REAL(KIND=8) ,TARGET :: ELEV
              REAL(KIND=8) ,TARGET :: SD
              REAL(KIND=8) ,TARGET :: ED
              INTEGER(KIND=4) :: INDEX
            END SUBROUTINE LINK_USGS_GW
          END INTERFACE 
        END MODULE LINK_USGS_GW__genmod
