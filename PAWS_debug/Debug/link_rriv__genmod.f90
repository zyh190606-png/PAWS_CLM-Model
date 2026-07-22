        !COMPILER-GENERATED INTERFACE MODULE: Sun Apr 09 23:47:55 2017
        MODULE LINK_RRIV__genmod
          INTERFACE 
            SUBROUTINE LINK_RRIV(H,HX,AC,U,WM,WF,RM,E,ETA,EX,ZBANK,QX,  &
     &ESTE,DX,MANN,S,SE,S0,SM,CODE,L,CTYPE,UPST,UPSV,DST,DSV,K,HG,MR,FG,&
     &QOC,QGC,QL,PRCP,EVAP,TRIB,COND,AB,RTR,T,DT,STATE,UMAX,X,Y,DDIST,  &
     &DIMS,NFIELDS,IS)
              INTEGER(KIND=4) :: NFIELDS
              INTEGER(KIND=4) :: DIMS(3,NFIELDS)
              REAL(KIND=8) ,TARGET :: H(DIMS((2,1)))
              REAL(KIND=8) ,TARGET :: HX(DIMS((2,2)))
              REAL(KIND=8) ,TARGET :: AC(DIMS((2,3)))
              REAL(KIND=8) ,TARGET :: U(DIMS((2,4)))
              REAL(KIND=8) ,TARGET :: WM(DIMS((2,5)))
              REAL(KIND=8) ,TARGET :: WF(DIMS((2,6)))
              REAL(KIND=8) ,TARGET :: RM(DIMS((2,7)))
              REAL(KIND=8) ,TARGET :: E(DIMS((2,8)))
              REAL(KIND=8) ,TARGET :: ETA(DIMS((2,9)))
              REAL(KIND=8) ,TARGET :: EX(DIMS((2,10)))
              REAL(KIND=8) ,TARGET :: ZBANK(DIMS((2,11)))
              REAL(KIND=8) ,TARGET :: QX(DIMS((2,12)))
              REAL(KIND=8) ,TARGET :: ESTE(DIMS((2,13)))
              REAL(KIND=8) ,TARGET :: DX
              REAL(KIND=8) ,TARGET :: MANN(DIMS((2,15)))
              REAL(KIND=8) ,TARGET :: S(DIMS((2,16)))
              REAL(KIND=8) ,TARGET :: SE(DIMS((2,17)))
              REAL(KIND=8) ,TARGET :: S0(DIMS((2,18)))
              REAL(KIND=8) ,TARGET :: SM(DIMS((2,19)))
              REAL(KIND=8) ,TARGET :: CODE
              REAL(KIND=8) ,TARGET :: L
              INTEGER(KIND=4) ,TARGET :: CTYPE(DIMS((2,22)))
              REAL(KIND=8) ,TARGET :: UPST
              REAL(KIND=8) ,TARGET :: UPSV
              REAL(KIND=8) ,TARGET :: DST
              REAL(KIND=8) ,TARGET :: DSV(DIMS((2,26)))
              REAL(KIND=8) ,TARGET :: K(DIMS((2,27)))
              REAL(KIND=8) ,TARGET :: HG(DIMS((2,28)))
              REAL(KIND=8) ,TARGET :: MR(DIMS((2,29)))
              REAL(KIND=8) ,TARGET :: FG(DIMS((2,30)))
              REAL(KIND=8) ,TARGET :: QOC(DIMS((2,31)))
              REAL(KIND=8) ,TARGET :: QGC(DIMS((2,32)))
              REAL(KIND=8) ,TARGET :: QL(DIMS((2,33)))
              REAL(KIND=8) ,TARGET :: PRCP(DIMS((2,34)))
              REAL(KIND=8) ,TARGET :: EVAP(DIMS((2,35)))
              REAL(KIND=8) ,TARGET :: TRIB(DIMS((2,36)))
              REAL(KIND=8) ,TARGET :: COND(DIMS((2,37)))
              REAL(KIND=8) ,TARGET :: AB(DIMS((2,38)))
              REAL(KIND=8) ,TARGET :: RTR(DIMS((2,39)))
              REAL(KIND=8) ,TARGET :: T
              REAL(KIND=8) ,TARGET :: DT
              INTEGER(KIND=4) ,TARGET :: STATE
              REAL(KIND=8) ,TARGET :: UMAX
              REAL(KIND=8) ,TARGET :: X(DIMS((2,44)))
              REAL(KIND=8) ,TARGET :: Y(DIMS((2,45)))
              REAL(KIND=8) ,TARGET :: DDIST(DIMS((2,46)))
              INTEGER(KIND=4) :: IS
            END SUBROUTINE LINK_RRIV
          END INTERFACE 
        END MODULE LINK_RRIV__genmod
