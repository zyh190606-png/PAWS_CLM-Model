        !COMPILER-GENERATED INTERFACE MODULE: Sun Apr 09 23:47:55 2017
        MODULE LINKVDZ__genmod
          INTERFACE 
            SUBROUTINE LINKVDZ(IS,M,N,NZ,NP,NZC,DZ,DDZ,H,SRC,JTYPE,WTL,K&
     &,CI,THE,DF,DPERC,T,DT,DTP,EB,E,RM,THER,THES,KS,KO,NN,LAMBDA,ALPHA,&
     &ET,FC,WP,MAPP,COND,PRATIO,VPARM,ND,GWL,STATE,SMR,KSSAVE,NR,REP,   &
     &NREP,SAND,CLAY,OM,ICE,SOIC,GTI)
              INTEGER(KIND=4) :: NREP
              INTEGER(KIND=4) :: NR
              INTEGER(KIND=4) :: ND
              INTEGER(KIND=4) :: NP
              INTEGER(KIND=4) :: NZ
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: M
              INTEGER(KIND=4) :: IS
              INTEGER(KIND=4) ,TARGET :: NZC(M,N)
              REAL(KIND=8) ,TARGET :: DZ(M,N,NZ)
              REAL(KIND=8) ,TARGET :: DDZ(M,N,NZ)
              REAL(KIND=8) ,TARGET :: H(M,N,NZ)
              REAL(KIND=8) ,TARGET :: SRC(M,N,NZ)
              INTEGER(KIND=4) ,TARGET :: JTYPE(M,N)
              INTEGER(KIND=4) ,TARGET :: WTL(M,N)
              REAL(KIND=8) ,TARGET :: K(M,N,NZ)
              REAL(KIND=8) ,TARGET :: CI(M,N,NZ)
              REAL(KIND=8) ,TARGET :: THE(M,N,NZ)
              REAL(KIND=8) ,TARGET :: DF(M,N)
              REAL(KIND=8) ,TARGET :: DPERC(M,N)
              REAL(KIND=8) ,TARGET :: T
              REAL(KIND=8) ,TARGET :: DT(M,N)
              REAL(KIND=8) ,TARGET :: DTP(3)
              REAL(KIND=8) ,TARGET :: EB(M,N,NZ)
              REAL(KIND=8) ,TARGET :: E(M,N,NZ)
              REAL(KIND=8) ,TARGET :: RM(M,N,NR)
              REAL(KIND=8) ,TARGET :: THER(M,N,NZ)
              REAL(KIND=8) ,TARGET :: THES(M,N,NZ)
              REAL(KIND=8) ,TARGET :: KS(M,N,NZ)
              REAL(KIND=8) ,TARGET :: KO(M,N,NZ)
              REAL(KIND=8) ,TARGET :: NN(M,N,NZ)
              REAL(KIND=8) ,TARGET :: LAMBDA(M,N,NZ)
              REAL(KIND=8) ,TARGET :: ALPHA(M,N,NZ)
              REAL(KIND=8) ,TARGET :: ET(M,N,NZ)
              REAL(KIND=8) ,TARGET :: FC(M,N,NZ)
              REAL(KIND=8) ,TARGET :: WP(M,N,NZ)
              REAL(KIND=8) ,TARGET :: MAPP(M,N)
              REAL(KIND=8) ,TARGET :: COND(M,N,ND)
              REAL(KIND=8) ,TARGET :: PRATIO(M,N)
              REAL(KIND=8) ,TARGET :: VPARM(NP)
              REAL(KIND=8) ,TARGET :: GWL(M,N)
              REAL(KIND=8) ,TARGET :: STATE
              REAL(KIND=8) ,TARGET :: SMR(M,N)
              REAL(KIND=8) ,TARGET :: KSSAVE(M,N,NZ)
              REAL(KIND=8) ,TARGET :: REP(M,N,NREP)
              REAL(KIND=8) ,TARGET :: SAND(M,N,NZ)
              REAL(KIND=8) ,TARGET :: CLAY(M,N,NZ)
              REAL(KIND=8) ,TARGET :: OM(M,N,NZ)
              REAL(KIND=8) ,TARGET :: ICE(M,N,NZ)
              INTEGER(KIND=4) ,TARGET :: SOIC(M,N)
              REAL(KIND=8) ,TARGET :: GTI(M,N)
            END SUBROUTINE LINKVDZ
          END INTERFACE 
        END MODULE LINKVDZ__genmod
