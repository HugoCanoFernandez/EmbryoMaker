program rewrite


! gfortran -g general.mod.f90 neighboring.mod.f90 genetic.mod.f90 io.mod.f90 geompack3.f90 check_genes.f90 -o chkgenes.e; rm *.mod

use io
use general
use genetic
use neighboring

Implicit none
character*314 :: copyit
character*500 :: input,output, input1
character*500 :: inpu1
character*152 :: name_dat2
integer iji, ihc, jhc, khc, lhc, nhc, mhc, phc, ord, jch, ich, kch
real micha, michi, diffi, diffa, mua, mui, thedist
real goldenratio
character( len= 80 ) :: ithc
integer :: kevincounter,iqpas
real*8 :: guardaene, kevinaux, kevintol
integer, allocatable :: list_change(:), list_change1(:)
integer :: sumlistchange, oldng
integer, allocatable :: tabkev(:,:), tabkev1(:,:)
real*8, allocatable :: tabgex(:,:), tabgex1(:,:), distkev(:,:)
integer :: sumbigexpr, sumepi

!******************************************************************************************************!
nseed=33
call getarg(1,input)

  call iniread
  call readsnap(trim(input))
!This program is used to change the mesenchyme adhesion, noise and deltamin
!of selected initial conditions 
call neighbor_build

if(allocated(list_change))deallocate(list_change)
allocate(list_change(1:ng)); list_change=0

do iqpas=1,ng

!print*,"which gene?"
!	read(*,*)iqpas   ! iqpas which gene you want to visualize
    sumbigexpr=0; sumepi=0
    do ich=1, nd
      if(node(ich)%tipus .gt. 2)cycle
      sumepi=sumepi+1
      if(gex(ich,iqpas) > 1d-8)sumbigexpr=sumbigexpr+1
    end do
    if(sumbigexpr .lt. 5.0d-2 * sumepi)cycle
    kevintol=5.0d0*abs((maxval(gex(:,iqpas))-minval(gex(:,iqpas)))/100d0)
    node(:)%orix=0.0
	  do i=1,nd
	    node(i)%orix=0d0
	  	
	    if (node(i)%tipus==1) then !apical epithelium	
	 	    kevinaux=0d0
	 	    kevincounter=0
	 	    do j=1,nneigh(i)	
	  		  ii=neigh(i,j)
	  		  if (node(ii)%tipus==1) then
	  		  	kevincounter=kevincounter+1
				    if (gex(i,iqpas)-gex(ii,iqpas)<-kevintol) then			
					    kevinaux=kevinaux-1d0
				    elseif (gex(i,iqpas)-gex(ii,iqpas)>kevintol) then 
					    kevinaux=kevinaux+1d0
		 		    endif
		 	    endif
		    enddo
		
		    if (kevincounter==0) then
			    kevinaux=0d0
		    else
	 		    kevinaux=2d0*kevinaux/real(kevincounter)
	 	    endif
	 		
	 	    if (kevinaux==2d0) then
	  		  !node(i)%bor=2
	  		  node(i)%orix=2d0
	  	  elseif (Kevinaux<2d0 .and. kevinaux>0d0) then
	  		  !node(i)%bor=1
	  		  node(i)%orix=1d0
	  	  elseif (kevinaux==0d0) then
	  		  !node(i)%bor=0
	  		  node(i)%orix=0d0
	  	  elseif (kevinaux<0d0 .and. kevinaux>-2d0) then
	  		  !node(i)%bor=-1
	  		  node(i)%orix=-1d0
	  	  elseif (kevinaux==-2d0) then
	  		  !node(i)%bor=-2
	  		  node(i)%orix=-2d0
	  	  endif		  		 	
	    endif
	  enddo 
	
    if(any(node(1:nd)%orix .ne. 0)) then
      list_change(iqpas)=1
      !tabkev(iqpas,1:nd)=node(1:nd)%orix
    end if
end do


sumlistchange=sum(list_change)
oldng=ng
if(sumlistchange>0)then
  print*,"yes"
else
  print*,"no"
end if


if(allocated(tabkev))deallocate(tabkev)
if(allocated(tabgex))deallocate(tabgex)
if(allocated(list_change))deallocate(list_change)
if(allocated(list_change1))deallocate(list_change1)
if(allocated(distkev))deallocate(distkev)

end program


