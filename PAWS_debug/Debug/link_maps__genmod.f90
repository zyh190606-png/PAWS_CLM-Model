        !COMPILER-GENERATED INTERFACE MODULE: Sun Apr 09 23:47:55 2017
        MODULE LINK_MAPS__genmod
          INTERFACE 
            SUBROUTINE LINK_MAPS(MAPS,PRCP,INF,PET,ET,EVAPG,DPERC,DPERC2&
     &,SATE,INFE,BOTT,QOC,QGC,SUBM,SMELT,GW,OVNH,N,CLMMAP,DIMS,NFIELDS)
              USE VDATA
              INTEGER(KIND=4) :: NFIELDS
              INTEGER(KIND=4) :: DIMS(3,NFIELDS)
              TYPE (MAPS_TYPE) :: MAPS
              REAL(KIND=8) ,TARGET :: PRCP(DIMS((1,1)),DIMS((2,1)))
              REAL(KIND=8) ,TARGET :: INF(DIMS((1,2)),DIMS((2,2)))
              REAL(KIND=8) ,TARGET :: PET(DIMS((1,3)),DIMS((2,3)))
              REAL(KIND=8) ,TARGET :: ET(DIMS((1,4)),DIMS((2,4)))
              REAL(KIND=8) ,TARGET :: EVAPG(DIMS((1,5)),DIMS((2,5)))
              REAL(KIND=8) ,TARGET :: DPERC(DIMS((1,6)),DIMS((2,6)))
              REAL(KIND=8) ,TARGET :: DPERC2(DIMS((1,7)),DIMS((2,7)))
              REAL(KIND=8) ,TARGET :: SATE(DIMS((1,8)),DIMS((2,8)))
              REAL(KIND=8) ,TARGET :: INFE(DIMS((1,9)),DIMS((2,9)))
              REAL(KIND=8) ,TARGET :: BOTT(DIMS((1,10)),DIMS((2,10)))
              REAL(KIND=8) ,TARGET :: QOC(DIMS((1,11)),DIMS((2,11)))
              REAL(KIND=8) ,TARGET :: QGC(DIMS((1,12)),DIMS((2,12)))
              REAL(KIND=8) ,TARGET :: SUBM(DIMS((1,13)),DIMS((2,13)))
              REAL(KIND=8) ,TARGET :: SMELT(DIMS((1,14)),DIMS((2,14)))
              REAL(KIND=8) ,TARGET :: GW(DIMS((1,15)),DIMS((2,15)))
              REAL(KIND=8) ,TARGET :: OVNH(DIMS((1,16)),DIMS((2,16)))
              REAL(KIND=8) ,TARGET :: N
              REAL(KIND=8) ,TARGET :: CLMMAP(DIMS((1,18)),DIMS((2,18)), &
     &DIMS((3,18)))
            END SUBROUTINE LINK_MAPS
          END INTERFACE 
        END MODULE LINK_MAPS__genmod
