        !COMPILER-GENERATED INTERFACE MODULE: Sun Apr 09 23:47:55 2017
        MODULE LINKVEG__genmod
          INTERFACE 
            SUBROUTINE LINKVEG(J,M,FRAC,LAI,RPT,EVAP,EVAPCS,BSOIL,ROOT, &
     &RDFN,EDFN,DELTA,HI,T,DT,KC,HC,CPC,CS,TP,ET,SMAP,RPET,PET,ALBEDO,  &
     &ASP,SWE,GPRCP,EVAPG,SALB,SNOCOV,SNOM,SNOCOEF,SNOW,NEWS,H,IMP,LUIDX&
     &)
              INTEGER(KIND=4) :: M(3)
              INTEGER(KIND=4) :: J
              REAL(KIND=8) ,TARGET :: FRAC(M((1)),M((2)),M((3)))
              REAL(KIND=8) ,TARGET :: LAI(M((1)),M((2)),M((3)))
              INTEGER(KIND=4) ,TARGET :: RPT(M((1)),M((2)),M((3)))
              REAL(KIND=8) ,TARGET :: EVAP(M((1)),M((2)))
              REAL(KIND=8) ,TARGET :: EVAPCS(M((1)),M((2)))
              REAL(KIND=8) ,TARGET :: BSOIL(M((1)),M((2)),M((3)))
              REAL(KIND=8) ,TARGET :: ROOT(M((1)),M((2)),M((3)))
              INTEGER(KIND=4) ,TARGET :: RDFN(M((1)),M((2)),M((3)))
              INTEGER(KIND=4) ,TARGET :: EDFN(M((1)),M((2)))
              REAL(KIND=8) ,TARGET :: DELTA(M((1)),M((2)),M((3)))
              REAL(KIND=8) ,TARGET :: HI(M((1)),M((2)))
              REAL(KIND=8) ,TARGET :: T
              REAL(KIND=8) ,TARGET :: DT
              REAL(KIND=8) ,TARGET :: KC(M((1)),M((2)),M((3)))
              REAL(KIND=8) ,TARGET :: HC(M((1)),M((2)),M((3)))
              REAL(KIND=8) ,TARGET :: CPC(M((1)),M((2)))
              REAL(KIND=8) ,TARGET :: CS(M((1)),M((2)))
              REAL(KIND=8) ,TARGET :: TP(M((1)),M((2)))
              REAL(KIND=8) ,TARGET :: ET(M((1)),M((2)))
              INTEGER(KIND=4) ,TARGET :: SMAP(M((1)),M((2)))
              REAL(KIND=8) ,TARGET :: RPET(M((1)),M((2)))
              REAL(KIND=8) ,TARGET :: PET(M((1)),M((2)))
              REAL(KIND=8) ,TARGET :: ALBEDO(M((1)),M((2)))
              REAL(KIND=8) ,TARGET :: ASP(M((1)),M((2)))
              REAL(KIND=8) ,TARGET :: SWE(M((1)),M((2)))
              REAL(KIND=8) ,TARGET :: GPRCP(M((1)),M((2)))
              REAL(KIND=8) ,TARGET :: EVAPG(M((1)),M((2)))
              REAL(KIND=8) ,TARGET :: SALB(M((1)),M((2)))
              REAL(KIND=8) ,TARGET :: SNOCOV(M((1)),M((2)))
              REAL(KIND=8) ,TARGET :: SNOM(M((1)),M((2)))
              REAL(KIND=8) ,TARGET :: SNOCOEF(2)
              REAL(KIND=8) ,TARGET :: SNOW(M((1)),M((2)))
              REAL(KIND=8) ,TARGET :: NEWS(M((1)),M((2)))
              REAL(KIND=8) ,TARGET :: H(M((1)),M((2)))
              REAL(KIND=8) ,TARGET :: IMP(M((1)),M((2)))
              REAL(KIND=8) ,TARGET :: LUIDX(5)
            END SUBROUTINE LINKVEG
          END INTERFACE 
        END MODULE LINKVEG__genmod
