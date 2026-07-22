        !COMPILER-GENERATED INTERFACE MODULE: Sun Apr 09 23:47:54 2017
        MODULE LINK_REC__genmod
          INTERFACE 
            SUBROUTINE LINK_REC(HG,MEM,MEMR,MEMQGC,HSAVE,TSS,M,NFIELD)
              INTEGER(KIND=4) :: NFIELD
              INTEGER(KIND=4) :: M(3,NFIELD)
              REAL(KIND=8) ,TARGET :: HG(M((1,1)),M((2,1)))
              REAL(KIND=8) ,TARGET :: MEM(M((2,2)))
              REAL(KIND=8) ,TARGET :: MEMR(M((2,3)))
              REAL(KIND=8) ,TARGET :: MEMQGC(M((2,4)))
              REAL(KIND=8) ,TARGET :: HSAVE(M((1,5)),M((2,5)))
              REAL(KIND=8) ,TARGET :: TSS(M((2,6)))
            END SUBROUTINE LINK_REC
          END INTERFACE 
        END MODULE LINK_REC__genmod
