        !COMPILER-GENERATED INTERFACE MODULE: Sun Apr 09 23:47:55 2017
        MODULE LINKSOILS__genmod
          INTERFACE 
            SUBROUTINE LINKSOILS(M,N,ND,KS,THER,FC,THES,ALPHA,LAMBDA,NN,&
     &DP,DPN,DPUL,DPUU,STATE)
              INTEGER(KIND=4) :: ND
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: M
              REAL(KIND=8) ,TARGET :: KS(M,N)
              REAL(KIND=8) ,TARGET :: THER(M,N)
              REAL(KIND=8) ,TARGET :: FC(M,N)
              REAL(KIND=8) ,TARGET :: THES(M,N)
              REAL(KIND=8) ,TARGET :: ALPHA(M,N)
              REAL(KIND=8) ,TARGET :: LAMBDA(M,N)
              REAL(KIND=8) ,TARGET :: NN(M,N)
              REAL(KIND=8) ,TARGET :: DP(ND,N)
              INTEGER(KIND=4) ,TARGET :: DPN(N)
              REAL(KIND=8) ,TARGET :: DPUL(N)
              REAL(KIND=8) ,TARGET :: DPUU(N)
              REAL(KIND=8) ,TARGET :: STATE
            END SUBROUTINE LINKSOILS
          END INTERFACE 
        END MODULE LINKSOILS__genmod
