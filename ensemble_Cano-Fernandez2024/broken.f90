program copy_network


! gfortran -g general.mod.f90 genetic.mod.f90 io.mod.f90 geompack3.f90 neighboring.mod.f90 analysis.mod.f90 broken.f90 -o broken.e; rm *.mod

  use general
  use genetic
  use io
  use neighboring
  use analysis

Implicit none
integer, parameter :: read_unit = 99
character*314 :: copyit
integer ios
real micha, michi, diffi, diffa, mua, mui
real goldenratio, myrand
character*400 :: input,output, rangfile, when, used
character*500 :: inpu1
character*152 :: name_dat2
integer :: ihc,jhc,lhc, khc, phc, shc, ord,ord2, nhc, mhc  !!>> HC 24-11-2020 counters (we use i for the gene that is going to mutate and j for the other gene in the interaction)
real*8, dimension(1:nga) :: max_elim, min_elim  !!>> HC 27-11-2020 These vectors store the limits of e matrix 
integer, dimension(1:nga) :: rembeh !>> HC 28-11-2020 Unused cellular behaviors/properties CHECK DIMENSION MANUALLY
integer, dimension(1:nga+1) :: ens_rembeh
real*8, dimension(1:5) :: max_glim, min_glim
real*8 :: intensa, intensb, intensc, intensd
character( len= 4 ) :: ithc
real prob_enact, doesEnact, regulmax, regulij
real mumin,mumax,effmin,effmax
real*8 :: myEffect, mymax
character*140 originNet, mygradient
integer ich, jch, kch, lch
integer :: mygen, myNWas, myWa
real*8 :: myval, mySign
integer :: totcounter, icounter, bcounter

nseed=33

call getarg(1,input)

!print*,trim(input)
call iniread
call readsnap(input)

call neighbor_build
!call fill_node_arrays
!call OPCval

bcounter=0; totcounter=0
do ich=1,nd
  if(node(ich)%tipus .gt. 2)cycle
  totcounter=totcounter+1
  icounter=0
  do jch=1,nneigh(ich)
    if(node(neigh(ich,jch))%tipus .gt. 2)cycle
    icounter=icounter+1
  end do
  if(icounter .lt. 3)then
    bcounter=bcounter+1
  end if
end do

if(bcounter .gt. 1)then
  print*,"THIS STUFF IS BROKEN"
  print*,bcounter
else
  print*,"EVERYTHING IS FINE"
end if

call do_analysis(1)

end program
