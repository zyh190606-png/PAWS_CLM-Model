        !COMPILER-GENERATED INTERFACE MODULE: Sun Apr 09 23:47:55 2017
        MODULE LINK_NET__genmod
          INTERFACE 
            SUBROUTINE LINK_NET(LEVEL,UPS,UPSC,UPSQ,DS,DSC,INSERT,DSQ,  &
     &PARENT,TRIB,TRIBC,DIMS,NFIELDS,IS)
              INTEGER(KIND=4) :: NFIELDS
              INTEGER(KIND=4) :: DIMS(3,NFIELDS)
              INTEGER(KIND=4) ,TARGET :: LEVEL(DIMS((2,1)))
              REAL(KIND=8) ,TARGET :: UPS
              INTEGER(KIND=4) ,TARGET :: UPSC(DIMS((2,3)))
              REAL(KIND=8) ,TARGET :: UPSQ
              REAL(KIND=8) ,TARGET :: DS
              REAL(KIND=8) ,TARGET :: DSC
              INTEGER(KIND=4) ,TARGET :: INSERT(DIMS((2,7)))
              REAL(KIND=8) ,TARGET :: DSQ
              INTEGER(KIND=4) ,TARGET :: PARENT(DIMS((2,9)))
              INTEGER(KIND=4) ,TARGET :: TRIB(DIMS((2,10)))
              INTEGER(KIND=4) ,TARGET :: TRIBC(DIMS((2,11)))
              INTEGER(KIND=4) :: IS
            END SUBROUTINE LINK_NET
          END INTERFACE 
        END MODULE LINK_NET__genmod
