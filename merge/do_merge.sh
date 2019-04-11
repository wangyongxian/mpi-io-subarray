#!/usr/bin/bash

# mpif90 -fpp -D REAL8 merge-files-into-one.f90 -o merge-real8.exe
mpif90 -fpp -D REAL4 merge-files-into-one.f90 -o merge-real4.exe

yhrun -n 432 -N 16 -p nmdis_meit1 ./merge-real4.exe

