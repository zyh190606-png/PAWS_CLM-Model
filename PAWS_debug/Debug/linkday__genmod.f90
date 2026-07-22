        !COMPILER-GENERATED INTERFACE MODULE: Sun Apr 09 23:47:55 2017
        MODULE LINKDAY__genmod
          INTERFACE 
            SUBROUTINE LINKDAY(J,M,PRCP,TEMP,RAD,RET,RETS,WND,HMD,RN,   &
     &DPTP,HB,HLI,COSZ)
              INTEGER(KIND=4) :: M(2)
              INTEGER(KIND=4) :: J
              REAL(KIND=8) ,TARGET :: PRCP(M((1)),M((2)))
              REAL(KIND=8) ,TARGET :: TEMP(M((1)),M((2)))
              REAL(KIND=8) ,TARGET :: RAD(M((1)),M((2)))
              REAL(KIND=8) ,TARGET :: RET(M((1)),M((2)))
              REAL(KIND=8) ,TARGET :: RETS(M((1)),M((2)))
              REAL(KIND=8) ,TARGET :: WND(M((1)),M((2)))
              REAL(KIND=8) ,TARGET :: HMD(M((1)),M((2)))
              REAL(KIND=8) ,TARGET :: RN(M((1)),M((2)))
              REAL(KIND=8) ,TARGET :: DPTP(M((1)),M((2)))
              REAL(KIND=8) ,TARGET :: HB(M((1)),M((2)))
              REAL(KIND=8) ,TARGET :: HLI(M((1)),M((2)))
              REAL(KIND=8) ,TARGET :: COSZ(M((1)),M((2)))
            END SUBROUTINE LINKDAY
          END INTERFACE 
        END MODULE LINKDAY__genmod
