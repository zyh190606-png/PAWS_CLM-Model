        !COMPILER-GENERATED INTERFACE MODULE: Sun Apr 09 23:47:55 2017
        MODULE LINKGW__genmod
          INTERFACE 
            SUBROUTINE LINKGW(J,NY,NX,NN,DD,H,EB,ES,E,D,K,W1,W2,ST,T,   &
     &HNEW,DZ,SEEP,DLDR,DR,DT,QPERC,A1B,A2,A3,A4,A5,TT,MAPP,STO,NLC,BOTK&
     &,BOTD,STATE,TOPBOUND)
              INTEGER(KIND=4) :: NN(3)
              INTEGER(KIND=4) :: NX
              INTEGER(KIND=4) :: NY
              INTEGER(KIND=4) :: J
              REAL(KIND=8) ,TARGET :: DD(2)
              REAL(KIND=8) ,TARGET :: H(NY,NX)
              REAL(KIND=8) ,TARGET :: EB(NY,NX)
              REAL(KIND=8) ,TARGET :: ES(NY,NX)
              REAL(KIND=8) ,TARGET :: E(NY,NX)
              REAL(KIND=8) ,TARGET :: D(NY,NX)
              REAL(KIND=8) ,TARGET :: K(NY,NX)
              REAL(KIND=8) ,TARGET :: W1(NY,NX)
              REAL(KIND=8) ,TARGET :: W2(NY,NX)
              REAL(KIND=8) ,TARGET :: ST(NY,NX)
              REAL(KIND=8) ,TARGET :: T(NY,NX)
              REAL(KIND=8) ,TARGET :: HNEW(NY,NX)
              REAL(KIND=8) ,TARGET :: DZ(NY,NX)
              REAL(KIND=8) ,TARGET :: SEEP(NY,NX)
              REAL(KIND=8) ,TARGET :: DLDR(NY,NX)
              REAL(KIND=8) ,TARGET :: DR(NY,NX)
              REAL(KIND=8) ,TARGET :: DT
              REAL(KIND=8) ,TARGET :: QPERC(NY,NX)
              REAL(KIND=8) ,TARGET :: A1B(NY,NX)
              REAL(KIND=8) ,TARGET :: A2(NY,NX)
              REAL(KIND=8) ,TARGET :: A3(NY,NX)
              REAL(KIND=8) ,TARGET :: A4(NY,NX)
              REAL(KIND=8) ,TARGET :: A5(NY,NX)
              REAL(KIND=8) ,TARGET :: TT
              REAL(KIND=8) ,TARGET :: MAPP(NY,NX)
              REAL(KIND=8) ,TARGET :: STO(NY,NX)
              REAL(KIND=8) ,TARGET :: NLC(NN((1)),NN((2)),NN((3)))
              REAL(KIND=8) ,TARGET :: BOTK(NY,NX)
              REAL(KIND=8) ,TARGET :: BOTD(NY,NX)
              REAL(KIND=8) ,TARGET :: STATE
              REAL(KIND=8) ,TARGET :: TOPBOUND(NY,NX)
            END SUBROUTINE LINKGW
          END INTERFACE 
        END MODULE LINKGW__genmod
