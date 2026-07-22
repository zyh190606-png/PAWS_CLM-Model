        !COMPILER-GENERATED INTERFACE MODULE: Sun Apr 09 23:47:55 2017
        MODULE LINKSTATIONS__genmod
          INTERFACE 
            SUBROUTINE LINKSTATIONS(ID,XYELEV,LATLONG,GINDICES,WPCT,    &
     &DATES,DATENUMS,VRATIO,PRCP,TMIN,TMAX,SOL,WND,HMD,DPTP,AWND,EVAP,  &
     &RAD,RD,RI,TAU,RRAD,S1,S0,SNOW,MNRH,MNTP,MXRH,NEIGHBOR,REF,PTR,    &
     &DDATA,PCT,DAT_TIME,COSZ,DECL,DIMS,NFIELDS,IS)
              INTEGER(KIND=4) :: NFIELDS
              INTEGER(KIND=4) :: DIMS(3,NFIELDS)
              INTEGER(KIND=4) ,TARGET :: ID(DIMS((2,1)))
              REAL(KIND=8) ,TARGET :: XYELEV(DIMS((2,2)))
              REAL(KIND=8) ,TARGET :: LATLONG(DIMS((2,3)))
              INTEGER(KIND=4) ,TARGET :: GINDICES(DIMS((1,4)))
              REAL(KIND=8) ,TARGET :: WPCT(DIMS((1,5)))
              REAL(KIND=8) ,TARGET :: DATES(DIMS((1,6)))
              REAL(KIND=8) ,TARGET :: DATENUMS(DIMS((1,7)))
              REAL(KIND=8) ,TARGET :: VRATIO(DIMS((1,8)))
              REAL(KIND=8) ,TARGET :: PRCP(DIMS((1,9)))
              REAL(KIND=8) ,TARGET :: TMIN(DIMS((1,10)))
              REAL(KIND=8) ,TARGET :: TMAX(DIMS((1,11)))
              REAL(KIND=8) ,TARGET :: SOL(DIMS((1,12)))
              REAL(KIND=8) ,TARGET :: WND(DIMS((1,13)))
              REAL(KIND=8) ,TARGET :: HMD(DIMS((1,14)))
              REAL(KIND=8) ,TARGET :: DPTP(DIMS((1,15)))
              REAL(KIND=8) ,TARGET :: AWND(DIMS((1,16)))
              REAL(KIND=8) ,TARGET :: EVAP(DIMS((1,17)))
              REAL(KIND=8) ,TARGET :: RAD(DIMS((1,18)))
              REAL(KIND=8) ,TARGET :: RD(DIMS((1,19)))
              REAL(KIND=8) ,TARGET :: RI(DIMS((1,20)))
              REAL(KIND=8) ,TARGET :: TAU(DIMS((1,21)))
              REAL(KIND=8) ,TARGET :: RRAD(DIMS((1,22)))
              REAL(KIND=8) ,TARGET :: S1(DIMS((1,23)))
              REAL(KIND=8) ,TARGET :: S0(DIMS((1,24)))
              REAL(KIND=8) ,TARGET :: SNOW(DIMS((1,25)))
              REAL(KIND=8) ,TARGET :: MNRH(DIMS((1,26)))
              REAL(KIND=8) ,TARGET :: MNTP(DIMS((1,27)))
              REAL(KIND=8) ,TARGET :: MXRH(DIMS((1,28)))
              INTEGER(KIND=4) ,TARGET :: NEIGHBOR(DIMS((1,29)),DIMS((2, &
     &29)))
              REAL(KIND=8) ,TARGET :: REF(DIMS((1,30)))
              REAL(KIND=8) ,TARGET :: PTR(DIMS((1,31)))
              REAL(KIND=8) ,TARGET :: DDATA(DIMS((1,32)))
              REAL(KIND=8) ,TARGET :: PCT(DIMS((1,33)))
              REAL(KIND=8) ,TARGET :: DAT_TIME(DIMS((1,34)))
              REAL(KIND=8) ,TARGET :: COSZ(DIMS((1,35)))
              REAL(KIND=8) ,TARGET :: DECL(DIMS((1,35)))
              INTEGER(KIND=4) :: IS
            END SUBROUTINE LINKSTATIONS
          END INTERFACE 
        END MODULE LINKSTATIONS__genmod
