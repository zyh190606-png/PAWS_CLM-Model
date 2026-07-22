        !COMPILER-GENERATED INTERFACE MODULE: Sun Apr 09 23:47:55 2017
        MODULE LINKGA__genmod
          INTERFACE 
            SUBROUTINE LINKGA(MS,M,NSPAN,KBAR,WFP,FM,CODE,THEW,SW,MM,H,N&
     &,L,DN,DL,HSPAN,GCODE)
              INTEGER(KIND=4) ,TARGET :: NSPAN
              INTEGER(KIND=4) ,TARGET :: M(3)
              INTEGER(KIND=4) :: MS(3)
              REAL(KIND=8) ,TARGET :: KBAR(MS((1)),MS((2)),MS((3)))
              REAL(KIND=8) ,TARGET :: WFP(MS((1)),MS((2)))
              REAL(KIND=8) ,TARGET :: FM(MS((1)),MS((2)),MS((3)))
              REAL(KIND=8) ,TARGET :: CODE(MS((1)),MS((2)))
              REAL(KIND=8) ,TARGET :: THEW(MS((1)),MS((2)))
              REAL(KIND=8) ,TARGET :: SW(M((1)),M((2)),M((3)))
              REAL(KIND=8) ,TARGET :: MM(M((1)),M((2)))
              REAL(KIND=8) ,TARGET :: H(M((3)))
              REAL(KIND=8) ,TARGET :: N(M((1)))
              REAL(KIND=8) ,TARGET :: L(M((2)))
              REAL(KIND=8) ,TARGET :: DN
              REAL(KIND=8) ,TARGET :: DL
              REAL(KIND=8) ,TARGET :: HSPAN(NSPAN,4)
              REAL(KIND=8) ,TARGET :: GCODE(4)
            END SUBROUTINE LINKGA
          END INTERFACE 
        END MODULE LINKGA__genmod
