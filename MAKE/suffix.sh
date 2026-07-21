#!/bin/bash
# Change suffix of fortran code from F90 to f90

cd ../src/PAWS/
rename F90 f90 *.F90
cd ../CLM/
cd biogeochem/
rename F90 f90 *.F90
cd ../biogeophys/
rename F90 f90 *.F90
cd ../main/
rename F90 f90 *.F90
cd ../../obsolete/
rename F90 f90 *.F90
cd ../../MAKE