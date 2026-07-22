        !COMPILER-GENERATED INTERFACE MODULE: Sun Apr 09 23:48:22 2017
        MODULE UPDATEVEG__genmod
          INTERFACE 
            SUBROUTINE UPDATEVEG(VEG,MAPP,V,DD)
              USE VDATA
              TYPE (VEGE_TYPE) :: VEG
              REAL(KIND=8) :: MAPP(VEG%NY,VEG%NX)
              TYPE (VDZ_TYPE) :: V
              TYPE (STA_TYPE) :: DD
            END SUBROUTINE UPDATEVEG
          END INTERFACE 
        END MODULE UPDATEVEG__genmod
