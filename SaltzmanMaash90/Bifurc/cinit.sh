#!/bin/sh
gfortran cinit.f90  sal3_forced.f90
./a.out > vdp.dat
head -n 1 vdp.dat && tail -n 1 vdp.dat

