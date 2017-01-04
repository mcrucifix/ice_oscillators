#!/bin/sh
gfortran cinit.f90  fort_pp04_forced.f90
./a.out > vdp.dat
head -n 1 vdp.dat && tail -n 1 vdp.dat

