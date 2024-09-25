#!/bin/bash
gfortran -g general.mod.f90 genetic.mod.f90 io.mod.f90 geompack3.f90 make_ics_dynamic.f90 -o mkics_dyn.e; rm *.mod
./mkics_dyn.e
chmod 777 mkics_dyn.e

