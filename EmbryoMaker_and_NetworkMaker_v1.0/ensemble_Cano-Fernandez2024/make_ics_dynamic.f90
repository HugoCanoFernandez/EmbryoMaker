program copy_network


! gfortran -g general.mod.f90 genetic.mod.f90 io.mod.f90 geompack3.f90 make_ics_dynamic.f90 -o mkics_dyn.e; rm *.mod

  use general
  use genetic
  use io

Implicit none
integer, parameter :: read_unit = 99
character*314 :: copyit, patternfile, pinput, prefix
integer ios
real micha, michi, diffi, diffa, mua, mui
real goldenratio, myrand
character*400 :: input,output, rangfile, when, used
character*500 :: inpu1
character*152 :: name_dat2
character(len=2)::jthc
character(len=1)::kthc
character(len=4)::jithc
character(len=3)::kithc
character(len=1)::ggg1,ggg2
character(len=2)::gggg1,gggg2
integer :: ihc,jhc,lhc, khc, phc, shc, ord,ord2, nhc, mhc,jich  !!>> HC 24-11-2020 counters (we use i for the gene that is going to mutate and j for the other gene in the interaction)
real*8, dimension(1:nga) :: max_elim, min_elim  !!>> HC 27-11-2020 These vectors store the limits of e matrix 
integer, dimension(1:nga) :: rembeh !>> HC 28-11-2020 Unused cellular behaviors/properties CHECK DIMENSION MANUALLY
integer, dimension(1:nga+1) :: ens_rembeh
real*8, dimension(1:5) :: max_glim, min_glim
real*8 :: intensa, intensb, intensc, intensd
character( len= 5 ) :: ithc
character( len= 5 ) :: ach
character( len= 5 ) :: bch
character( len= 5 ) :: dch
character( len= 5 ) :: ech
character( len= 20 ) :: idtt
real prob_enact, doesEnact, regulmax, regulij
real mumin,mumax,effmin,effmax
real*8 :: myEffect, mymax
character*140 originNet, mygradient
integer ich, jch, kch, lch
integer :: mygen, myNWas, myWa
real*8 :: myval, mySign
integer :: listParam(9)
character*500 :: listPatterns(19)
real*8 :: maxgex
real*8,allocatable :: tabGex(:,:,:)
integer, allocatable :: listBeh(:)
integer :: tmpNGen, g1, g2
integer :: mygenes(3), mywas(3), isgood
integer :: nact
integer, allocatable :: actstf(:), isneg(:), networks(:), netgens(:,:)
integer :: gen1, gen2, ntot, npattern, maxgen
real*8 :: myval1, myval1bis, myval2, myval2bis
integer, allocatable :: listAct(:,:), myngs(:)
real*8, allocatable :: listEnact(:,:), maxexp(:,:)
real*8, allocatable :: tO(:,:,:), kindofO(:,:), muO(:,:), diffO(:,:), eO(:,:,:), michO(:,:)
integer, allocatable :: prugenes(:,:)
integer, allocatable :: newics(:)
!integer :: gnet(1035)

!open(unit=read_unit, file="list_good_indices_dyna.dat", iostat=ios)
!do ich=1, 1035
!read(read_unit,*)gnet(ich)
!end do
!close(read_unit)

nseed=33



! INITIAL FILES
input="icoepi_nuclei.dat"
rangfile="ranges_ENSEMBLE_epi.dat"
patternfile="list_pruned_networks.txt"


! READ RANGES
call read_rang(rangfile, max_elim, min_elim, max_glim, min_glim, rembeh)  

!call fdate(when)
!used="cp "//trim(rangfile)//" used_ranges"//trim(when)//".dat"
!call system(used)

if(allocated(actstf))deallocate(actstf)
nact=52-sum(rembeh)
allocate(actstf(1:nact)); actstf=0
ord=0
do ich=1,52
  if(rembeh(ich) .eq. 0)then
    ord=ord+1
    actstf(ord)=ich
  end if
end do
if(allocated(isneg))deallocate(isneg)
allocate(isneg(1:nact)); isneg=0
do ich=1,nact
  if(min_elim(actstf(ich)) < 0d0)isneg(ich)=1
end do

! READ LIST OF PATTERNS
open(unit=read_unit, file=patternfile, iostat=ios)
read(read_unit,*) npattern
read(read_unit,*) maxgen
if(allocated(networks))deallocate(networks)
if(allocated(netgens))deallocate(netgens)
if(allocated(prugenes))deallocate(prugenes)
if(allocated(newics))deallocate(newics)
allocate(networks(1:npattern)); networks=0
allocate(netgens(1:npattern,1:maxgen)); netgens=0
allocate(prugenes(1:npattern,2)); prugenes=0
allocate(newics(1:npattern)); newics=0
do ich=1, npattern
  read(read_unit,*) networks(ich)
  read(read_unit,*) prugenes(ich,1:2)
  read(read_unit,*) netgens(ich,1:2)
  read(read_unit,*) newics(ich)
end do
close(read_unit)

print*,sum(newics)


! IMPORT NETWORKS AND INITIAL GENE PATTERNS
if(allocated(tO))deallocate(tO)
allocate(tO(1:npattern,1:30,1:30)); tO=0d0
if(allocated(kindofO))deallocate(kindofO)
allocate(kindofO(1:npattern,1:30)); kindofO=0d0
if(allocated(muO))deallocate(muO)
allocate(muO(1:npattern,1:30)); muO=0d0
if(allocated(michO))deallocate(michO)
allocate(michO(1:npattern,1:30)); michO=0d0
if(allocated(diffO))deallocate(diffO)
allocate(diffO(1:npattern,1:30)); diffO=0d0
if(allocated(maxexp))deallocate(maxexp)
allocate(maxexp(1:npattern,1:30)); maxexp=0d0
if(allocated(tabGex))deallocate(tabGex)
allocate(tabGex(1:npattern,1:744,1:2)); tabGex=0d0

!if(allocated(eO))deallocate(eO)
!allocate(eO(1:npattern,1:30,1:nparam_per_node+ngcb)); eO=0d0
if(allocated(myngs))deallocate(myngs)
allocate(myngs(1:npattern)); myngs=0
do ich=1, npattern
  write( ithc, "(I5.5)" ) networks(ich)
  if (newics(ich) .eq. 1) then
    prefix="new_ics/"
  else
    prefix="ics/"
  end if
  if(prugenes(ich,1) < 10 .and. prugenes(ich,2) < 10)then
    write(ggg1, "(I1)")prugenes(ich,1)
    write(ggg2, "(I1)")prugenes(ich,2)
    pinput=trim(prefix)//"icoepi_nuclei.datmod"//ithc//".dat_prune_"//trim(ggg1)//"_"//trim(ggg2)//".dat"
  else if(prugenes(ich,1) .ge. 10 .and. prugenes(ich,2) < 10)then
    write(gggg1, "(I2)")prugenes(ich,1)
    write(ggg2, "(I1)")prugenes(ich,2)
    pinput=trim(prefix)//"icoepi_nuclei.datmod"//ithc//".dat_prune_"//trim(gggg1)//"_"//trim(ggg2)//".dat"
  else if(prugenes(ich,1) < 10 .and. prugenes(ich,2) .ge. 10)then
    write(ggg1, "(I1)")prugenes(ich,1)
    write(gggg2, "(I2)")prugenes(ich,2)
    pinput=trim(prefix)//"icoepi_nuclei.datmod"//ithc//".dat_prune_"//trim(ggg1)//"_"//trim(gggg2)//".dat"
  else
    write(gggg1, "(I2)")prugenes(ich,1)
    write(gggg2, "(I2)")prugenes(ich,2)
    pinput=trim(prefix)//"icoepi_nuclei.datmod"//ithc//".dat_prune_"//trim(gggg1)//"_"//trim(gggg2)//".dat"
  end if
  print*,trim(pinput)
  call iniread
  call readsnap(pinput)
  myngs(ich)=ng
  do jch=1,ng
    maxexp(ich,jch) = maxval(gex(:,jch))
    kindofO(ich,jch) = gen(jch)%kindof
    muO(ich,jch) = gen(jch)%mu
    diffO(ich,jch) = gen(jch)%diffu
    michO(ich,jch) = gen(jch)%mich
    do kch=1,ng
      tO(ich,jch,kch)=gen(jch)%t(kch)
    end do
  end do
  do jch=1,nd
    tabGex(ich,jch,1)=gex(jch,1)
    tabGex(ich,jch,2)=gex(jch,2)
  end do
end do


! GET MEAN MAX EXPRESSION
do ich=1, npattern
  write( ithc, "(I5.5)" ) networks(ich)
  do jch=1,10
    if (newics(ich) .eq. 0) then
      jich=2*jch
      print*,"checking network",ich,"at time",jich
      ! Making ICs name
      if(jich .ge. 10)then
        write(jthc, "(I2)") jich
        if(prugenes(ich,1) < 10 .and. prugenes(ich,2) < 10)then
          write(ggg1, "(I1)")prugenes(ich,1)
          write(ggg2, "(I1)")prugenes(ich,2)
          pinput="ics/icoepi_nuclei.datmod"//ithc//".dat_prune_"//trim(ggg1)//"_"//trim(ggg2)//".dat"//trim(jthc)//"000.dat"
        else if(prugenes(ich,1) .ge. 10 .and. prugenes(ich,2) < 10)then
          write(gggg1, "(I2)")prugenes(ich,1)
          write(ggg2, "(I1)")prugenes(ich,2)
          pinput="ics/icoepi_nuclei.datmod"//ithc//".dat_prune_"//trim(gggg1)//"_"//trim(ggg2)//".dat"//trim(jthc)//"000.dat"
        else if(prugenes(ich,1) < 10 .and. prugenes(ich,2) .ge. 10)then
          write(ggg1, "(I1)")prugenes(ich,1)
          write(gggg2, "(I2)")prugenes(ich,2)
          pinput="ics/icoepi_nuclei.datmod"//ithc//".dat_prune_"//trim(ggg1)//"_"//trim(gggg2)//".dat"//trim(jthc)//"000.dat"
        else
          write(gggg1, "(I2)")prugenes(ich,1)
          write(gggg2, "(I2)")prugenes(ich,2)
          pinput="ics/icoepi_nuclei.datmod"//ithc//".dat_prune_"//trim(gggg1)//"_"//trim(gggg2)//".dat"//trim(jthc)//"000.dat"
        end if
      else
        write(kthc,"(I1)") jich
        if(prugenes(ich,1) < 10 .and. prugenes(ich,2) < 10)then
          write(ggg1, "(I1)")prugenes(ich,1)
          write(ggg2, "(I1)")prugenes(ich,2)
          pinput="ics/icoepi_nuclei.datmod"//ithc//".dat_prune_"//ggg1//"_"//ggg2//".dat"//trim(kthc)//"000.dat"
        else if(prugenes(ich,1) .ge. 10 .and. prugenes(ich,2) < 10)then
          write(gggg1, "(I2)")prugenes(ich,1)
          write(ggg2, "(I1)")prugenes(ich,2)
          pinput="ics/icoepi_nuclei.datmod"//ithc//".dat_prune_"//gggg1//"_"//ggg2//".dat"//trim(kthc)//"000.dat"
        else if(prugenes(ich,1) < 10 .and. prugenes(ich,2) .ge. 10)then
          write(ggg1, "(I1)")prugenes(ich,1)
          write(gggg2, "(I2)")prugenes(ich,2)
          pinput="ics/icoepi_nuclei.datmod"//ithc//".dat_prune_"//ggg1//"_"//gggg2//".dat"//trim(kthc)//"000.dat"
        else
          write(gggg1, "(I2)")prugenes(ich,1)
          write(gggg2, "(I2)")prugenes(ich,2)
          pinput="ics/icoepi_nuclei.datmod"//ithc//".dat_prune_"//gggg1//"_"//gggg2//".dat"//trim(kthc)//"000.dat"
        end if
      end if
    else
      jich=332*jch
      print*,"checking network",ich,"at time",jich
      ! Making ICs name
      if(jich .ge. 1000)then
        write(jithc, "(I4)") jich
        if(prugenes(ich,1) < 10 .and. prugenes(ich,2) < 10)then
          write(ggg1, "(I1)")prugenes(ich,1)
          write(ggg2, "(I1)")prugenes(ich,2)
          pinput="new_ics/icoepi_nuclei.datmod"//ithc//".dat_prune_"//ggg1//"_"//ggg2//".dat"//trim(jithc)//".dat"
        else if(prugenes(ich,1) .ge. 10 .and. prugenes(ich,2) < 10)then
          write(gggg1, "(I2)")prugenes(ich,1)
          write(ggg2, "(I1)")prugenes(ich,2)
          pinput="new_ics/icoepi_nuclei.datmod"//ithc//".dat_prune_"//gggg1//"_"//ggg2//".dat"//trim(jithc)//".dat"
        else if(prugenes(ich,1) < 10 .and. prugenes(ich,2) .ge. 10)then
          write(ggg1, "(I1)")prugenes(ich,1)
          write(gggg2, "(I2)")prugenes(ich,2)
          pinput="new_ics/icoepi_nuclei.datmod"//ithc//".dat_prune_"//ggg1//"_"//gggg2//".dat"//trim(jithc)//".dat"
        else
          write(gggg1, "(I2)")prugenes(ich,1)
          write(gggg2, "(I2)")prugenes(ich,2)
          pinput="new_ics/icoepi_nuclei.datmod"//ithc//".dat_prune_"//gggg1//"_"//gggg2//".dat"//trim(jithc)//".dat"
        end if
      else
        write(kithc,"(I3)") jich
        if(prugenes(ich,1) < 10 .and. prugenes(ich,2) < 10)then
          write(ggg1, "(I1)")prugenes(ich,1)
          write(ggg2, "(I1)")prugenes(ich,2)
          pinput="new_ics/icoepi_nuclei.datmod"//ithc//".dat_prune_"//ggg1//"_"//ggg2//".dat"//trim(kithc)//".dat"
        else if(prugenes(ich,1) .ge. 10 .and. prugenes(ich,2) < 10)then
          write(gggg1, "(I2)")prugenes(ich,1)
          write(ggg2, "(I1)")prugenes(ich,2)
          pinput="new_ics/icoepi_nuclei.datmod"//ithc//".dat_prune_"//gggg1//"_"//ggg2//".dat"//trim(kithc)//".dat"
        else if(prugenes(ich,1) < 10 .and. prugenes(ich,2) .ge. 10)then
          write(ggg1, "(I1)")prugenes(ich,1)
          write(gggg2, "(I2)")prugenes(ich,2)
          pinput="new_ics/icoepi_nuclei.datmod"//ithc//".dat_prune_"//ggg1//"_"//gggg2//".dat"//trim(kithc)//".dat"
        else
          write(gggg1, "(I2)")prugenes(ich,1)
          write(gggg2, "(I2)")prugenes(ich,2)
          pinput="new_ics/icoepi_nuclei.datmod"//ithc//".dat_prune_"//gggg1//"_"//gggg2//".dat"//trim(kithc)//".dat"
        end if
      end if
    end if
    ! Importing network and reading max gene expression
    call iniread
    call readsnap(pinput)
    do kch=1,myngs(ich)
      !maxexp(ich,kch) = maxexp(ich,kch) + maxval(gex(:,kch))
      if(maxval(gex(:,kch)) > maxexp(ich,kch))then
        maxexp(ich,kch) = maxval(gex(:,kch))
      end if
    end do
  end do
  !do kch=1,myngs(ich)
  !  maxexp(ich,kch) = maxexp(ich,kch)/11
  !end do
end do


! LIST NUMBER OF ICS
ntot=0
do ich=1, npattern
  do jch=1, nact
    myval1bis=0d0
    myval1=max_elim(actstf(jch))
    if(isneg(jch) .eq. 1)myval1bis=min_elim(actstf(jch))
    do kch=1, nact
      if(kch .eq. jch)cycle
      myval2bis=0d0
      myval2=max_elim(actstf(kch))
      if(isneg(kch) .eq. 1)myval2bis=min_elim(actstf(kch))
      if(myval1bis > -1d-8)then
        if(myval2bis > -1d-8)then
          ntot=ntot+1
        else
          ntot=ntot+2
        end if
      else
        if(myval2bis > -1d-8)then
          ntot=ntot+2
        else
          ntot=ntot+4
        end if
      end if
    end do
  end do
end do

print*,"There will be",ntot,"IC files"

if(allocated(listAct))deallocate(listAct)
if(allocated(listEnact))deallocate(listEnact)
allocate(listAct(1:ntot,6)); listAct=0
allocate(listEnact(1:ntot,2)); listEnact=0d0

! REGISTER ALL POSSIBLE ICs COMBINATIONS
ord=0
do ich=1, npattern
  do jch=1, nact
    myval1bis=0d0
    myval1=max_elim(actstf(jch))
    if(isneg(jch) .eq. 1)myval1bis=min_elim(actstf(jch))
    do kch=1, nact
      if(kch .eq. jch)cycle
      myval2bis=0d0
      myval2=max_elim(actstf(kch))
      if(isneg(kch) .eq. 1)myval2bis=min_elim(actstf(kch))
      if(myval1bis > -1d-8)then
        if(myval2bis > -1d-8)then
          ord=ord+1
          listAct(ord,1)=ich; listAct(ord,2)=1; listAct(ord,3)=2
          listAct(ord,4)=actstf(jch); listAct(ord,5)=actstf(kch)
          listEnact(ord,1)=myval1; listEnact(ord,2)=myval2
          listAct(ord,6)=networks(ich)
        else
          ord=ord+1
          listAct(ord,1)=ich; listAct(ord,2)=1; listAct(ord,3)=2
          listAct(ord,4)=actstf(jch); listAct(ord,5)=actstf(kch)
          listEnact(ord,1)=myval1; listEnact(ord,2)=myval2
          listAct(ord,6)=networks(ich)
          ord=ord+1
          listAct(ord,1)=ich; listAct(ord,2)=1; listAct(ord,3)=2
          listAct(ord,4)=actstf(jch); listAct(ord,5)=actstf(kch)
          listEnact(ord,1)=myval1; listEnact(ord,2)=myval2bis
          listAct(ord,6)=networks(ich)
        end if
      else
        if(myval2bis > -1d-8)then
          ord=ord+1
          listAct(ord,1)=ich; listAct(ord,2)=1; listAct(ord,3)=2
          listAct(ord,4)=actstf(jch); listAct(ord,5)=actstf(kch)
          listEnact(ord,1)=myval1; listEnact(ord,2)=myval2
          listAct(ord,6)=networks(ich)
          ord=ord+1
          listAct(ord,1)=ich; listAct(ord,2)=1; listAct(ord,3)=2
          listAct(ord,4)=actstf(jch); listAct(ord,5)=actstf(kch)
          listEnact(ord,1)=myval1bis; listEnact(ord,2)=myval2
          listAct(ord,6)=networks(ich)
        else
          ord=ord+1
          listAct(ord,1)=ich; listAct(ord,2)=1; listAct(ord,3)=2
          listAct(ord,4)=actstf(jch); listAct(ord,5)=actstf(kch)
          listEnact(ord,1)=myval1; listEnact(ord,2)=myval2
          listAct(ord,6)=networks(ich)
          ord=ord+1
          listAct(ord,1)=ich; listAct(ord,2)=1; listAct(ord,3)=2
          listAct(ord,4)=actstf(jch); listAct(ord,5)=actstf(kch)
          listEnact(ord,1)=myval1bis; listEnact(ord,2)=myval2
          listAct(ord,6)=networks(ich)
          ord=ord+1
          listAct(ord,1)=ich; listAct(ord,2)=1; listAct(ord,3)=2
          listAct(ord,4)=actstf(jch); listAct(ord,5)=actstf(kch)
          listEnact(ord,1)=myval1; listEnact(ord,2)=myval2bis
          listAct(ord,6)=networks(ich)
          ord=ord+1
          listAct(ord,1)=ich; listAct(ord,2)=1; listAct(ord,3)=2
          listAct(ord,4)=actstf(jch); listAct(ord,5)=actstf(kch)
          listEnact(ord,1)=myval1bis; listEnact(ord,2)=myval2bis
          listAct(ord,6)=networks(ich)
        end if
      end if
    end do
  end do
end do


! MAKE ICs AND EXPORT ICs FILES

do ich=1, ntot
  print*, listAct(ich,:)
enddo

do ich=1, ntot
  !if(all(gnet .ne. ich))cycle
  !if(ich .ne. 1)cycle
  call iniread
  call readsnap(input)
  
  
  ng=31
  call initiate_gene
  if(any(listAct(ich,:) .eq. 0))cycle
  
  !if(allocated(kadh))deallocate(kadh) !!!!!!!!
  !allocate(kadh(2,2)) !!!!!!!!!
  !kadh=0d0 !!!!!!!!!!
  !kadh(1,1)=4d0 !!!!!!!!
  !kadh(2,2)=4d0 !!!!!!!!!!
  !ntipusadh=2 !!!!!!!!!!!
  ffu(14)=0
  ffu(22)=0
  
  
  
  gex=0d0
  do jch=1,nd
    gex(jch,1) = tabGex(listAct(ich,1),jch,1)
    gex(jch,2) = tabGex(listAct(ich,1),jch,2)
    if(newics(listAct(ich,1)) .eq. 1)then
      gex(jch,31) = 1.0d0
    else
      gex(jch,31) = 1.0d0
    end if
  end do
  do jch=1, myngs(listAct(ich,1))
    gen(jch)%diffu = diffO(listAct(ich,1),jch)
    gen(jch)%mu = muO(listAct(ich,1),jch)
    gen(jch)%kindof = kindofO(listAct(ich,1),jch)
    gen(jch)%mich = michO(listAct(ich,1),jch)
    do kch=1, myngs(listAct(ich,1))
      gen(jch)%t(kch) = tO(listAct(ich,1),jch,kch)
    end do
  end do
  g1 = netgens(listAct(ich,1),listAct(ich,2))
  g2 = netgens(listAct(ich,1),listAct(ich,3))
  
  print*,"maxexp for g1",maxexp(listAct(ich,1),g1)
  print*,"and for g2",maxexp(listAct(ich,1),g2)
  gen(g1)%e(listAct(ich,4)) = (listEnact(ich,1)/maxexp(listAct(ich,1),g1))*0.750d0
  gen(g2)%e(listAct(ich,5)) = (listEnact(ich,2)/maxexp(listAct(ich,1),g2))*0.750d0
  !gen(g1)%e(1) = 1d0/maxexp(listAct(ich,1),g1) !!!!!!!!!!!!!!!
  !gen(g2)%e(2) = 1d0/maxexp(listAct(ich,1),g2) !!!!!!!!!!!!!!!
  if(gen(g1)%e(nparam_per_node+8) > 0d0)then
    gen(31)%e(nparam_per_node+17)=(0.1d0/maxexp(listAct(ich,1),g1))*0.750d0
  end if
  if(gen(g2)%e(nparam_per_node+8) > 0d0)then
    gen(31)%e(nparam_per_node+17)=(0.1d0/maxexp(listAct(ich,1),g2))*0.750d0
  end if
  
  gen(31)%e(nparam_per_node+2) = 0.50d0
  
  
  call iniio
  call writesnap
  write( idtt, "(4I5.5)" ) ich, listAct(ich,6), listAct(ich,4), listAct(ich,5)  !This saves the number of the file to use it then in the name of the file
  write( ach, "(I5.5)" ) ich
  write( bch, "(I5.5)" ) listAct(ich,6)
  write( ech, "(I5.5)" ) listAct(ich,4)
  write( dch, "(I5.5)" ) listAct(ich,5)
  output=trim(input)//"dyn"//ach//"_"//bch//"_"//ech//"_"//dch//".dat"
  print*,output
	!This actually makes the new file
  name_dat2=trim("cp fort.1 "//adjustl(output)) 
  call system(name_dat2)
end do
open(unit=read_unit, file="epi_patterns_dynamic.dat", iostat=ios)
do ich=1, ntot
  write(read_unit,'(I5,A,I5,A,I5,A,I5,A,I5,A,I5,A,E20.8,A,E20.8)') ich,",",networks(listAct(ich,1)),",",&
  & netgens(listAct(ich,1),listAct(ich,2)),",",&
  & netgens(listAct(ich,1),listAct(ich,3)),",",&
  & listAct(ich,4),",", listAct(ich,5),",",&
  & listEnact(ich,1)/maxexp(listAct(ich,1),netgens(listAct(ich,1),listAct(ich,2))),&
  & ",",listEnact(ich,2)/maxexp(listAct(ich,1),netgens(listAct(ich,1),listAct(ich,3)))
end do
close(read_unit)
if(allocated(actstf))deallocate(actstf)
if(allocated(networks))deallocate(networks)
if(allocated(netgens))deallocate(netgens)

end program


