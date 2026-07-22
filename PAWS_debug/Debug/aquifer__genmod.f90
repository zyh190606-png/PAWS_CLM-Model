        !COMPILER-GENERATED INTERFACE MODULE: Sun Apr 09 23:48:22 2017
        MODULE AQUIFER__genmod
          INTERFACE 
            SUBROUTINE AQUIFER(GW,MAPP,BBC,NBC,OPT)
              USE VDATA
              INTEGER(KIND=4) :: NBC
              TYPE (GW2D_TYPE) :: GW
              REAL(KIND=8) :: MAPP(GW%NY,GW%NX)
              TYPE (BC_TYPE) :: BBC(NBC)
              INTEGER(KIND=4) :: OPT
            END SUBROUTINE AQUIFER
          END INTERFACE 
        END MODULE AQUIFER__genmod
