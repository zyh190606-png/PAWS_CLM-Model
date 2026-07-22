        !COMPILER-GENERATED INTERFACE MODULE: Sun Apr 09 23:47:43 2017
        MODULE PENMANMONTEITH__genmod
          INTERFACE 
            SUBROUTINE PENMANMONTEITH(NT,RAD,HMD,E0,WND,T,ELEV,TAU,REF, &
     &RET,RETS,RN,HB,HLI,RSS)
              INTEGER(KIND=4) :: NT
              REAL(KIND=8) :: RAD(NT)
              REAL(KIND=8) :: HMD(NT)
              REAL(KIND=8) :: E0(NT)
              REAL(KIND=8) :: WND(NT)
              REAL(KIND=8) :: T(NT)
              REAL(KIND=8) :: ELEV
              REAL(KIND=8) :: TAU
              INTEGER(KIND=4) :: REF
              REAL(KIND=8) :: RET(NT)
              REAL(KIND=8) :: RETS(NT)
              REAL(KIND=8) :: RN(NT)
              REAL(KIND=8) :: HB(NT)
              REAL(KIND=8) :: HLI(NT)
              REAL(KIND=8) :: RSS(NT,2)
            END SUBROUTINE PENMANMONTEITH
          END INTERFACE 
        END MODULE PENMANMONTEITH__genmod
