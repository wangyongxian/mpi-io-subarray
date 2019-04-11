#!/usr/bin/bash

# case 1: OK!

echo "------ first run: "
rm test.*.dat
mpif90 test-ser-write-files.f90                     -o ser-write.exe
./ser-write.exe
mpif90 test-par-read-file.f90                       -o par-read.exe
yhrun -n 12 -N 2 -p nmdis_meit1 ./par-read.exe

# case 2: failed!!!

echo "------ second run: "
rm test.*.dat
mpif90 test-ser-write-files.f90 -convert big_endian -o ser-write.exe
./ser-write.exe
mpif90 test-par-read-file.f90   -convert big_endian -o par-read.exe
yhrun -n 12 -N 2 -p nmdis_meit1 ./par-read.exe

