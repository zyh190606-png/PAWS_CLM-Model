        !COMPILER-GENERATED INTERFACE MODULE: Sun Apr 09 23:47:55 2017
        MODULE LINKTH__genmod
          INTERFACE 
            SUBROUTINE LINKTH(J,M,A,C,BBASE,TEMP,KT,CT,SSTATEV,SSITEV,  &
     &SOUTV,TSS,TSA,SPARM,DFC,SMR,SNOP,SNOE,SALB,BSUB,ND,NP,STATE)
              INTEGER(KIND=4) :: NP
              INTEGER(KIND=4) :: ND(2)
              INTEGER(KIND=4) :: M(3)
              INTEGER(KIND=4) :: J
              REAL(KIND=8) ,TARGET :: A(M((1)),M((2)),M((3)))
              REAL(KIND=8) ,TARGET :: C(M((1)),M((2)),M((3)))
              REAL(KIND=8) ,TARGET :: BBASE(M((1)),M((2)),M((3)))
              REAL(KIND=8) ,TARGET :: TEMP(M((1)),M((2)),M((3)))
              REAL(KIND=8) ,TARGET :: KT(M((1)),M((2)),M((3)))
              REAL(KIND=8) ,TARGET :: CT(M((1)),M((2)),M((3)))
              REAL(KIND=8) ,TARGET :: SSTATEV(M((1)),M((2)),11)
              REAL(KIND=8) ,TARGET :: SSITEV(M((1)),M((2)),6)
              REAL(KIND=8) ,TARGET :: SOUTV(M((1)),M((2)),23)
              REAL(KIND=8) ,TARGET :: TSS(M((1)),M((2)),24)
              REAL(KIND=8) ,TARGET :: TSA(M((1)),M((2)),24)
              REAL(KIND=8) ,TARGET :: SPARM(NP)
              REAL(KIND=8) ,TARGET :: DFC(ND((1)),ND((2)))
              REAL(KIND=8) ,TARGET :: SMR(M((1)),M((2)))
              REAL(KIND=8) ,TARGET :: SNOP(M((1)),M((2)))
              REAL(KIND=8) ,TARGET :: SNOE(M((1)),M((2)))
              REAL(KIND=8) ,TARGET :: SALB(M((1)),M((2)))
              REAL(KIND=8) ,TARGET :: BSUB(M((1)),M((2)))
              REAL(KIND=8) ,TARGET :: STATE
            END SUBROUTINE LINKTH
          END INTERFACE 
        END MODULE LINKTH__genmod
