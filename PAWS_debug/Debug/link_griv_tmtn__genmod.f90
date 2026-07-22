        !COMPILER-GENERATED INTERFACE MODULE: Sun Apr 09 23:47:55 2017
        MODULE LINK_GRIV_TMTN__genmod
          INTERFACE 
            SUBROUTINE LINK_GRIV_TMTN(VALUES,IR,NZMAX,JC,COLUMNS,M,DIMS,&
     &INDEX,FIELD)
              INTEGER(KIND=4) :: COLUMNS
              INTEGER(KIND=4) :: NZMAX
              REAL(KIND=8) :: VALUES(NZMAX)
              INTEGER(KIND=4) :: IR(NZMAX)
              INTEGER(KIND=4) :: JC(COLUMNS)
              INTEGER(KIND=4) :: M(3)
              INTEGER(KIND=4) :: DIMS(3)
              INTEGER(KIND=4) :: INDEX
              CHARACTER(*) :: FIELD
            END SUBROUTINE LINK_GRIV_TMTN
          END INTERFACE 
        END MODULE LINK_GRIV_TMTN__genmod
