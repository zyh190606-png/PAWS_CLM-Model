        !COMPILER-GENERATED INTERFACE MODULE: Sun Apr 09 23:47:54 2017
        MODULE LINK_MCHUNK__genmod
          INTERFACE 
            SUBROUTINE LINK_MCHUNK(L1M,L2M,L3M,L4M,L5M,R1M,R2M,R3M,R4M, &
     &R5M,CLEVELMAPM,M,NFIELD)
              INTEGER(KIND=4) :: NFIELD
              INTEGER(KIND=4) :: M(5,NFIELD)
              INTEGER(KIND=4) ,TARGET :: L1M(M((1,1)))
              INTEGER(KIND=4) ,TARGET :: L2M(M((1,2)),M((2,2)))
              INTEGER(KIND=4) ,TARGET :: L3M(M((1,3)),M((2,3)),M((3,3)))
              INTEGER(KIND=4) ,TARGET :: L4M(M((1,4)),M((2,4)),M((3,4)),&
     &M((4,4)))
              INTEGER(KIND=4) ,TARGET :: L5M(M((1,5)),M((2,5)),M((3,5)),&
     &M((4,5)),M((5,5)))
              REAL(KIND=8) ,TARGET :: R1M(M((1,6)))
              REAL(KIND=8) ,TARGET :: R2M(M((1,7)),M((2,7)))
              REAL(KIND=8) ,TARGET :: R3M(M((1,8)),M((2,8)),M((3,8)))
              REAL(KIND=8) ,TARGET :: R4M(M((1,9)),M((2,9)),M((3,9)),M((&
     &4,9)))
              REAL(KIND=8) ,TARGET :: R5M(M((1,10)),M((2,10)),M((3,10)),&
     &M((4,10)),M((5,10)))
              INTEGER(KIND=4) ,TARGET :: CLEVELMAPM(M((1,11)),M((2,11)))
            END SUBROUTINE LINK_MCHUNK
          END INTERFACE 
        END MODULE LINK_MCHUNK__genmod
