        !COMPILER-GENERATED INTERFACE MODULE: Sun Apr 09 23:48:19 2017
        MODULE RASTER2RIV__genmod
          INTERFACE 
            SUBROUTINE RASTER2RIV(M,DAT,J,OUTPUT,OPTION)
              INTEGER(KIND=4) :: J
              INTEGER(KIND=4) :: M(2)
              REAL(KIND=8) :: DAT(M((1)),M((2)))
              REAL(KIND=8) :: OUTPUT(SIZE())
              INTEGER(KIND=4) ,OPTIONAL :: OPTION
            END SUBROUTINE RASTER2RIV
          END INTERFACE 
        END MODULE RASTER2RIV__genmod
