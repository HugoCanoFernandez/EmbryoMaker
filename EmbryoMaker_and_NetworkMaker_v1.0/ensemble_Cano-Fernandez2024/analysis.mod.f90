!    EmbryoMaker software (General Node Model)
!    Computational model to simulate morphogenetic processes in living organs and tissues.
!    Copyright (C) 2014 Miquel Marin-Riera, Miguel Brun-Usan, Roland Zimm, Tommi Välikangas & Isaac Salazar-Ciudad

!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    any later version.

!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.

!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <http://www.gnu.org/licenses/>.

! Created by Renske 24-04-2018
! This module contains functions for analysing embryonic development. They could be called by auto in automaticon, 
! which would however require one to divide an elli run into subruns (otherwise there is not much to analyse)
! 

module analysis

  use general
  use genetic
  use io
  use neighboring
  !parameters
  implicit none 
  real*8, allocatable :: storprev(:,:),stornow(:,:) !to generate a centered tissue
    
  real*8, allocatable ::maxgex(:,:) !store max expression per gene (and the time of maximum
  integer ::epi_prev, epi_now
  
  
  integer :: nrnods !the previous nr of nodes
  integer :: grofile=942
  integer :: OPCdistfile=337
  integer :: OPCfile=338
  integer :: gefile1=863  !file handle
  integer :: gefile2=864  !file handle
  integer :: gefile3=865 !file handle
  
  
  !************************************

  
  contains
  
  ! gathers the several subroutines in one handy bundle
  subroutine do_analysis(which)
    integer :: which
    
    call fill_node_arrays
    
    if (which==1 .or. which==3) then 
      call analyse_growth
      call OPCval
    end if
    if (which==2 .or. which==3) then 
      call analyse_gene_expression(1)
      call analyse_gene_expression(2)
      call maxexpr
    else if (which/=1) then
      print *, "Error in analysis.mod: invalid option. Please use 1, 2 or 3. Exiting..."
      stop
    end if
     
    !copy the situation now for later 
    !if (which==1 .or. which==3) call copy_nodes
    
  end subroutine

  !to save the initial conditions
  subroutine initialise_analysis(filestart,which)
    integer which
    character*300 filestart
    
    print *, "starting analysis"  
    call fill_node_arrays  
      
    if (which==1 .or. which==3) then !analyse tissue growth and shape
      open(grofile, file=trim(filestart)//"growthprofile.dat")
      write (grofile, *) 0, nd, 0.0, 0, 0, 0, 0, 0, 0 
      open(OPCdistfile, file=trim(filestart)//"patchsizedist.dat")
      open(OPCfile, file=trim(filestart)//"patchcount.dat")
      call OPCval
      call copy_nodes
    end if
    
    if (which==2 .or. which==3) then !analyse gene expression
      !call copy_gex
      open(gefile1, file=trim(filestart)//"genexpression1.dat")
      call analyse_gene_expression(1)
      open(gefile2, file=trim(filestart)//"genexpression2.dat")
      call analyse_gene_expression(2)
      
      !store max expression
      open(gefile3,file=trim(filestart)//"maxgex.dat")
      if (allocated(maxgex)) deallocate(maxgex)
      allocate(maxgex(ng,2))
      maxgex=0.0d0
      call maxexpr
      
    else if (which/=1) then
      print *, "Error in analysis.mod: invalid option. Please use 1, 2 or 3. Exiting..."
      stop
    end if
    

    
  end subroutine 
  
  subroutine copy_nodes
  
    if(allocated(storprev)) deallocate(storprev)
    allocate(storprev(epi_now,6))
    
    storprev=stornow
    
  end subroutine
  
  
  !close the files that were written into, delete the used memory
  subroutine close_analysis(which)
    integer which
    print *, "ending analysis"
    if (which==1 .or. which==3) then
      close(grofile)
      close(OPCfile)
      close(OPCdistfile)
      !if (allocated(prevnods)) deallocate(prevnods)
    end if
    if (which==2 .or. which==3) then 
      close(gefile1)
      close(gefile2)
      !if (allocated(prevgex)) deallocate(prevgex)
      call writemax
      close(gefile3)
    end if
    
    if (allocated(storprev)) deallocate(storprev)
    if (allocated(stornow)) deallocate(stornow)
    if (allocated(maxgex)) deallocate(maxgex)
    
  end subroutine
  

  subroutine maxexpr
    
    do i=1,nd
      do j=1, ng
        if (gex(i,j)>maxgex(j,1)) then
          maxgex(j,1)=gex(i,j)
          maxgex(j,2)=rtime
        endif
      end do
    end do
  
  end subroutine
  
  subroutine writemax
    print *, "here"
    write(gefile3, *) "gene"//char(9)//"maxex"//char(9)//"timept"
    do j=1, ng
      write(gefile3, '(A,es10.3,es10.1)')gen(j)%label, maxgex(j,1), maxgex(j,2)
    end do
    
  end subroutine
  
  !new version: gotta deal with the fact that nodes are in different layers
  subroutine analyse_gene_expression(layer)
    integer::layer !are we looking at apical or basal layer?
    !real*8, allocatable linegex(:,:,:) !to store the expression on a line in a nr of directions
    integer :: i,j, gefile
    
    !get the right file handle
    if (layer==1) gefile=gefile1
    if (layer==2) gefile=gefile2
    
    
    !x is 0, y and z vary
    do i=1, epi_now
    !!.and. stornow(i,2)>(miny_x+cell-stornow(i,4)*0.5) & .and. stornow(i,2)<(miny_x+cell + stornow(i,4)*0.5)
      if (stornow(i,5)==layer .and. stornow(i,1)<stornow(i,4)*1.5 .and. stornow(i,1)>-stornow(i,4)*1.5  ) then 
        write(gefile, *) rtime, stornow(i,1), stornow(i,2), stornow(i,3), gex(int(stornow(i,6)),:)
        !cell=cell+stornow(i,4)
      end if
    end do
    write(gefile, *)
    
    !y is 0, x and z vary
    do i=1, epi_now
    !.and. stornow(i,1)>minx_y+cell-stornow(i,4)*0.5 & .and. stornow(i,1)<minx_y+cell+stornow(i,4)*0.5
      if(stornow(i,5)==layer .and. stornow(i,2)<stornow(i,4)*1.5 .and. stornow(i,2)>-stornow(i,4)*1.5 ) then 
        write(gefile, *) rtime, stornow(i,1), stornow(i,2), stornow(i,3), gex(int(stornow(i,6)),:)
        !cell=cell+stornow(i,4)
      end if
    end do
    write(gefile, *)
    
    !z is 0, x and y vary
    do i=1, epi_now
    !.and. stornow(i,1)>minx_z+cell-stornow(i,4)*0.5 & .and. stornow(i,1)<minx_z+cell+stornow(i,4)*0.5
      if(stornow(i,5)==layer .and.  stornow(i,3)<stornow(i,4)*1.5 .and. stornow(i,3)>-stornow(i,4)*1.5 ) then 
        write(gefile, *) rtime, stornow(i,1), stornow(i,2), stornow(i,3), gex(int(stornow(i,6)),:)
        !cell=cell+stornow(i,4)
      end if
    end do
    write(gefile, *)
    write(gefile, *)
    
  end subroutine
  
  
  subroutine analyse_growth
    real*8 :: d, dd, aa, a, b, c, req1, req2 !for emd
    integer :: rostral1, caudal1, dorsal1, ventral1, left1, right1 !to assess growth or displacement 
    integer :: rostral2, caudal2, dorsal2, ventral2, left2, right2 !in different areas
    integer :: i,j
    integer :: patchcount
      !now EMD (this does not seem to be working too well: always prints 1000000)
      
      d=0.0d0
      do i=1,epi_prev
        ! if(node(i)%tipus>2)cycle
        a=storprev(i,1) ; b=storprev(i,2) ; c=storprev(i,3) ; req1=storprev(i,4)
        aa=1000000
        do j=1,epi_now
          dd=((sqrt((stornow(j,1)-a)**2+(stornow(j,2)-b)**2+(stornow(j,3)-c)**2)))
          if (dd<aa) aa=dd
        end do          
        d=d+aa
      end do
      do i=1,epi_now
        a=stornow(i,1) ; b=stornow(i,2) ; c=stornow(i,3) ; req2=stornow(i,4)
        aa=1000000
        do j=1,epi_prev
          dd=((sqrt((storprev(j,1)-a)**2+(storprev(j,2)-b)**2+(storprev(j,3)-c)**2)))
          if (dd<aa) aa=dd
        end do          
        d=d+aa
      enddo
      d=d/real(epi_prev+epi_now)
    
      !use the centered tissues to identify where more growth has occurred: use the fact that new nodes are added at the end of the array!
      rostral1=0
      caudal1=0
      dorsal1=0
      ventral1=0
      left1=0
      right1=0
      do i=epi_prev+1, epi_now
        if (stornow(i,1)>0.0) then ; right1=right1+1; else; left1=left1+1; end if
        if (stornow(i,2)>0.0) then ; dorsal1=dorsal1+1; else; ventral1=ventral1+1; end if
        if (stornow(i,3)>0.0) then ; rostral1=rostral1+1; else; caudal1=caudal1+1; end if  
      end do
      
!       !use the centered tissues to identify where more growth has occurred
!       rostral2=0
!       caudal2=0
!       dorsal2=0
!       ventral2=0
!       left2=0
!       right2=0
!       do i=1, epi_now
!         if (stornow(i,1)>0.0) then ; right2=right2+1; else; left2=left2+1; end if
!         if (stornow(i,2)>0.0) then ; dorsal2=dorsal2+1; else; ventral2=ventral2+1; end if
!         if (stornow(i,3)>0.0) then ; rostral2=rostral2+1; else; caudal2=caudal2+1; end if  
!       end do
    
      
      write (grofile, *) rtime, nd, d, rostral1,caudal1,dorsal1,ventral1,left1,right1
    
     
  end subroutine

  
  !The arrays with the orientation of each node, and their neighbours, have already been filled in fill_node_arrays
  subroutine OPCval
    integer :: patch_count
    integer, allocatable :: checks(:), store1(:),store2(:)
    integer, allocatable :: orientation(:), epi_node_list(:),translate(:) !to store the epithelium orientation
    integer::patch_size,nstore1,nstore2,min_patch_size,max_patch, maxmaxpatch 
    integer,allocatable::epi_neigh_list(:,:), patchsize_dist(:)
    integer :: nepi, maxneigh !nr of epithelial cells (by counting only apical nodes), and nr of neighbouring cells
    real*8 :: a, b, c, d
    integer :: i, j, ii, jj, jjj, k, kkk
    
    maxmaxpatch=500000000 !!>>HC 3-9-2020 with large polyps we can have a large OPC sizes
    
    max_patch=0
    min_patch_size=3
    !print *, "starting OPC"
    
    !allocating memory
    if (allocated(orientation))deallocate(orientation)
    allocate(orientation(ncels))

    if (allocated(epi_node_list))deallocate(epi_node_list)
    allocate(epi_node_list(ncels))
    
    if (allocated(store1))deallocate(store1)
    allocate(store1(ncels))

    if (allocated(store2))deallocate(store2)
    allocate(store2(ncels))
    
    if (allocated(checks))deallocate(checks)
    allocate(checks(ncels))    

    if (allocated(epi_neigh_list)) deallocate(epi_neigh_list)
    maxneigh=maxval(nneigh) !; print*,"maxneigh",maxneigh
    allocate(epi_neigh_list(ncels,maxneigh))
    
    if (allocated(patchsize_dist))deallocate(patchsize_dist)
    allocate(patchsize_dist(ncels))    
    
    if (allocated(translate))deallocate(translate)
    allocate(translate(nda))
      
    nepi=0
    do i=1,nd
      if (node(i)%tipus==2) then
        nepi=nepi+1
        epi_node_list(nepi)=i
        translate(i)=nepi
        ii=0
        ii=node(i)%altre
        a=node(i)%x-node(ii)%x  !subtraction gives the vector spanning the epithelial cylinder
        b=node(i)%y-node(ii)%y
        c=node(i)%z-node(ii)%z
        d=1/sqrt(a**2+b**2+c**2)
        a=a*d ; b=b*d ; c=c*d
           
       !we determine the orientation of the vector  
       if(c>=-0.95d0.and.c<=-epsilod)then  
         if(b<=-epsilod)then
           if(a<=-epsilod)then
             orientation(nepi)=1
           else
             orientation(nepi)=2
           end if
         else
           if(a<=-epsilod)then
             orientation(nepi)=3
           else
             orientation(nepi)=4
           end if
         end if
       elseif(c<=0.95d0.and.c>=epsilod)then
         if(b<=-epsilod)then
           if(a<=-epsilod)then
             orientation(nepi)=5
           else
             orientation(nepi)=6
           end if
         else
           if(a<=-epsilod)then
             orientation(nepi)=7
           else
             orientation(nepi)=8
           end if
         end if
       else
         orientation(nepi)=9
       end if !end of orientation determination
         
     end if !end of apical nodes conditional
        
    enddo
    
    epi_neigh_list=0
    !fill list of neighbouring nodes
    do ii=1,nepi
      jj=0
      i=epi_node_list(ii)
      do j=1,nneigh(i)
        if(node(neigh(i,j))%tipus==2)then
          jj=jj+1
          epi_neigh_list(ii,jj)=translate(neigh(i,j))
        end if 
      end do
    end do
   
    !count patches and record distribution
    checks=0
    patch_count=0
    patchsize_dist=0
    do i=1,nepi
    ii=epi_node_list(i)
    if(checks(i)==1) cycle
    checks(i)=1
    if(orientation(i)==0) cycle !it has not proper orientation
    !we start a patch
    store1=0 ; store2=0 ; nstore1=0 ; nstore2=0
    patch_size=1
    do j=1,maxneigh
      k=epi_neigh_list(i,j) 
      if(k==0) exit
    !  print*,"i",i,"j",j,"k",k
      if(checks(k)==1) cycle
      if(orientation(i)/=orientation(k))cycle
      patch_size=patch_size+1
      nstore2=nstore2+1
      store2(nstore2)=k
      checks(k)=1
    end do
    do while (nstore2/=0)
      store1=store2 ; store2=0
      nstore1=nstore2 ; nstore2=0
      do j=1,nstore1
        jjj=store1(j)
        do k=1,maxneigh
          kkk=epi_neigh_list(jjj,k)
          if(kkk==0) exit
          if(checks(kkk)==1) cycle
          if(orientation(jjj)/=orientation(kkk))cycle
          nstore2=nstore2+1
          patch_size=patch_size+1
          store2(nstore2)=kkk
          checks(kkk)=1
        end do
      end do
    end do
 
    if(patch_size>=min_patch_size)then 
      patch_count=patch_count+1
      patchsize_dist(patch_size)=patchsize_dist(patch_size)+1
      if(patch_size>max_patch) max_patch=patch_size
    endif
    
  end do
  
  if (max_patch>maxmaxpatch) print *, "OPC warning: exceeding maximum defined patch size"
  
  
  !write the size distribution to a file
  do i=min_patch_size,max_patch
    if (patchsize_dist(i)>0) then
      do j=1,patchsize_dist(i)
        write(OPCdistfile, fmt="(es10.2,I4)") rtime, i
      end do
    end if
  end do
!     do i=max_patch+1,maxmaxpatch
!       write(OPCdistfile, fmt="(I8,I4,I4)") timepoint, i, 0 
!     end do
!     write(OPCdistfile, *) ""

!    print*, "OPC", patch_count  !!>> HC 30-11-2020

      open(202,file="individual.datfitness")                              !!>> HC 11-11-2021 WRITING OUTPUT
      write(202,*) real(patch_count)                                      !!>> HC 11-11-2021
      close(202)                                                          !!>> HC 11-11-2021
      open(203,file="OPC.val")                                            !!>> HC 11-11-2021
      write(203,*) real(patch_count)                                      !!>> HC 11-11-2021

    
    write(OPCfile, fmt="(es10.2,I4)") rtime, patch_count 
    
    !deallocating memory, this important, otherwise mem leaks may occur!
    if (allocated(orientation))deallocate(orientation)

    if (allocated(epi_node_list))deallocate(epi_node_list)
     
    if (allocated(store1))deallocate(store1)
   
    if (allocated(store2))deallocate(store2)
      
    if (allocated(checks))deallocate(checks)
   
    if (allocated(epi_neigh_list)) deallocate(epi_neigh_list)
       
    if (allocated(patchsize_dist))deallocate(patchsize_dist)
        
    if (allocated(translate))deallocate(translate)
   
    
  end subroutine
  
    
  subroutine fill_node_arrays
    real*8 :: centx, centy, centz
    integer ::epnd
    integer :: i
     
    
    !store the node positions of the current time point
    epi_now=0
    do i=1,nd
      if(node(i)%tipus>2)cycle  
      epi_now=epi_now+1
    enddo
    
    if (epi_now>0) then 
    !now we transfer the coordinates of the individual to matrix storprev
      if (allocated(stornow)) deallocate(stornow)
      allocate(stornow(epi_now,6))
      epnd=0
      
      do i=1,nd
        if(node(i)%tipus>2)cycle  
        epnd=epnd+1
        stornow(epnd,1)=node(i)%x
        stornow(epnd,2)=node(i)%y
        stornow(epnd,3)=node(i)%z
        stornow(epnd,4)=node(i)%eqd  !!>>HC 30-6-2020
        stornow(epnd,5)=node(i)%tipus
        stornow(epnd,6)=i
      
      end do
      
    endif
  
    centx=sum(stornow(:,1))/epnd     !centroid
    centy=sum(stornow(:,2))/epnd
    centz=sum(stornow(:,3))/epnd
   
    do i=1,epnd
      stornow(i,1)=stornow(i,1)-centx   !centering 
      stornow(i,2)=stornow(i,2)-centy
      stornow(i,3)=stornow(i,3)-centz
    end do
  
 
  end subroutine
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
  subroutine joint_entropy_fourth_original_neighbors(fitelli)  !!>> HC 28-9-2021 This subroutine was mady by Hugo Cano and it calculates complexity
                                                      !!>> HC 28-9-2021 as the joint entropy of the curvature classes 
                                                      !!>> HC 28-9-2021 of the original nodes of the embryo with their 4th original neighbors
                                                      !!>> HC 28-9-2021 We assume that there is not apoptosis
                                                      !!>> HC 28-9-2021 We need to know the original number of apical nodes and
                                                      !!>> HC 28-9-2021 the original number of  relationships node-neighbor
                                                      !!>> HC 28-9-2021 These relationships must also be written in two files in the working directory
                                                      !!>> HC 28-9-2021 fourth_original_neighs.dat for the tipe 1 nodes
                                                      !!>> HC 28-9-2021 fourth_original_neighs2.dat for the tipe 2 nodes
   integer, parameter :: finitnn = 8280               !!>> HC 28-9-2021 Initial number of association node-4th neighbour
   integer, parameter :: napical = 372                !!>> HC 28-9-2021 Number of original nodes
   real*8 :: angle, upv, modu, modv, sumangle, numn, medianglee, entropy, entropy2, entropy1 !!>> HC 28-9-2021
   real*8 ::  pij,  classtemsum, rnapical, totalneightshc                                    !!>> HC 28-9-2021
   integer ::  ihc, jhc, khc, tiponodhc, order, order2, node2_1                            !!>> HC 28-9-2021
   integer:: tempclass, lowlimit, highlimit, steps, scaling, nclasseshc                    !!>> HC 28-9-2021
   integer, dimension(1:napical) :: classangles, classtem, oldorder                        !!>> HC 28-9-2021
   real*8, dimension(1:3) :: u, v                                                            !!>> HC 28-9-2021
   real*8,  dimension(1:napical) :: avangles, probclass                                      !!>> HC 28-9-2021
   real*8, allocatable, dimension(:,:) :: cooccurrences, cooprobs, classprod                 !!>> HC 28-9-2021
   integer, dimension(1:finitnn) :: foldnodes, foldnodesneight, angles1, angles2           !!>> HC 28-9-2021
   real*8 :: fitelli                                                                       !!>> HC 11-11-2021
   
   call neighbor_build                                                                     !!>> HC 28-9-2021



   lowlimit = -10                                            !!>> HC 28-9-2021 Lower limit of the interval in classes
   highlimit= 9                                              !!>> HC 28-9-2021 Higher limit of the interval in classes
   steps = 1                                                 !!>> HC 28-9-2021 width of the interval in classes
   nclasseshc = 20                                           !!>> HC 28-9-2021 Number of classes
   scaling = 10                                              !!>> HC 28-9-2021 Scaling of the avangle vector to fit into the classes

   if(allocated(cooccurrences))deallocate(cooccurrences)     !!>> HC 28-9-2021
   allocate( cooccurrences( nclasseshc, nclasseshc ) )       !!>> HC 28-9-2021
   if(allocated(cooprobs))deallocate(cooprobs)               !!>> HC 28-9-2021
   allocate( cooprobs( nclasseshc, nclasseshc ) )            !!>> HC 28-9-2021
   if(allocated(classprod))deallocate(classprod)             !!>> HC 28-9-2021
   allocate( classprod( nclasseshc, nclasseshc ) )           !!>> HC 28-9-2021
   rnapical = napical !!>> HC 28-9-2021 In case we have to operate with napical as a real number 

   entropy1=0; entropy2=0                                    !!>> HC 28-9-2021

   do tiponodhc=1, 2                                         !!>> HC 28-9-2021
      if (tiponodhc==1) then                                 !!>> HC 28-9-2021
      !!>> HC 28-9-2021 This reads the file were the initial 4th neightbour associations are written down
	 open (777,file="fourth_original_neighs.dat")             !!>> HC 28-9-2021 
	      do ihc=1, finitnn                                   !!>> HC 28-9-2021
	         read(777,*) foldnodes(ihc), foldnodesneight(ihc) !!>> HC 28-9-2021
	      enddo                                               !!>> HC 28-9-2021
	 close(777)                                               !!>> HC 28-9-2021
      else                                                        !!>> HC 28-9-2021
        !!>> HC 28-9-2021 This reads the file were the initial 4th neightbour associations are written down
	 open ( 777, file = "fourth_original_neighs2.dat" )       !!>> HC 28-9-2021
	      do ihc = 1, finitnn                                 !!>> HC 28-9-2021
		 read(777,*) foldnodes(ihc), foldnodesneight(ihc) !!>> HC 28-9-2021
	      enddo                                               !!>> HC 28-9-2021
	 close(777)                                               !!>> HC 28-9-2021
      endif                                                       !!>> HC 28-9-2021

 	!!>> HC 28-9-2021 this sets the dimension of the vector of average angles to the number of apical nodes

	!initializing 
      avangles=0; classangles=0; classtem=0; probclass=0           !!>> HC 28-9-2021
      cooccurrences=0                                              !!>> HC 28-9-2021
      oldorder=0; cooccurrences=0; cooprobs=0; tempclass=0         !!>> HC 28-9-2021
      totalneightshc=0; medianglee=0;  sumangle=0 ; numn=0         !!>> HC 28-9-2021
      pij=0; u=0; v=0; order=0; order2=0 ; angles1=0; angles2=0    !!>> HC 28-9-2021
	             

	!!!!!>> HC 28-9-2021 NEIGHBOURS !!!
	!!!!!>> HC 28-9-2021 This is the neightbour finding/angle calculating loop!!!
 
      do ihc=1, finitnn                                                   !!>> HC 28-9-2021 Look in the list of initial 2nd neightbourd associations
         u(1)=(node(node(foldnodes(ihc))%altre)%x-node(foldnodes(ihc))%x) !!>> HC 28-9-2021 Apical-basal vector i cell  
         u(2)=(node(node(foldnodes(ihc))%altre)%y-node(foldnodes(ihc))%y) !!>> HC 28-9-2021
         u(3)=(node(node(foldnodes(ihc))%altre)%z-node(foldnodes(ihc))%z) !!>> HC 28-9-2021
         v(1)=(node(foldnodesneight(ihc))%x-node(foldnodes(ihc))%x)       !!>> HC 28-9-2021 Vector between apical i cell node and the apical j cell node	
         v(2)=(node(foldnodesneight(ihc))%y-node(foldnodes(ihc))%y)       !!>> HC 28-9-2021
         v(3)=(node(foldnodesneight(ihc))%z-node(foldnodes(ihc))%z)       !!>> HC 28-9-2021

         upv=u(1)*v(1)+u(2)*v(2)+u(3)*v(3)                                !!>> HC 28-9-2021 vectorial product u * v
         modu=sqrt((u(1))**2+(u(2))**2+(u(3))**2)                         !!>> HC 28-9-2021 modulus u
         modv=sqrt((v(1))**2+(v(2))**2+(v(3))**2)                         !!>> HC 28-9-2021 modulus v
         angle=upv/(modu*modv)                                            !!>> HC 28-9-2021 getting the dot product between u and v ranges from -1 to 1 and is directly proportional to the angle
		
         if (ihc==finitnn)then                     !!>> HC 28-9-2021
            sumangle=sumangle+angle                !!>> HC 28-9-2021 sum of angles of the neighbouts in a cell i (last node of the class must be taken into account)
            numn=numn +1                           !!>> HC 28-9-2021 Counter of the real number of neights (last node of the class must be taken into account)		
            order=order +1                         !!>> HC 28-9-2021 This is actually the order of the nodes, fom 1st to 271st
            medianglee=sumangle/numn               !!>> HC 28-9-2021 This obtains the average dot product of the cell i with its neightbours
            avangles(order)=medianglee             !!>> HC 28-9-2021 This stores this average in the avangles vector
            numn=0                                 !!>> HC 28-9-2021 This resets the number of neights
            sumangle=0                             !!>> HC 28-9-2021 This resets the summation of dot products
            oldorder(order)=foldnodes(ihc)         !!>> HC 28-9-2021 Saves the old node number to use it later when calculating the cooccurrence                          
         else
            if (foldnodes(ihc)==foldnodes(ihc+1)) then    !!>> HC 28-9-2021 If the next node of the list is still the same
               sumangle=sumangle+angle                    !!>> HC 28-9-2021 sum of angles of the neighbouts in a cell i 
               numn=numn+1                                !!>> HC 28-9-2021 Counter of the real number of neights 
            else                                          !!>> HC 28-9-2021 If the next node is different
               sumangle=sumangle+angle                    !!>> HC 28-9-2021 sum of angles of the neighbouts in a cell i (last node of the class must be taken into account)
               numn=numn+1                                !!>> HC 28-9-2021 Counter of the real number of neights (last node of the class must be taken into account)		
               order=order+1                              !!>> HC 28-9-2021 This is actually the order of the nodes, fom 1st to 271st
               medianglee=sumangle/numn                   !!>> HC 28-9-2021 This obtains the average dot product of the cell i with its neightbours
               avangles(order)=medianglee                 !!>> HC 28-9-2021 This stores this average in the avangles vector
               numn=0                                     !!>> HC 28-9-2021 This resets the number of neights
               sumangle=0                                 !!>> HC 28-9-2021 This resets the summation of dot products
               oldorder(order)=foldnodes(ihc)             !!>> HC 28-9-2021 Saves the old node number to use it later when calculating the cooccurrence
            endif                                         !!>> HC 28-9-2021
         endif                                            !!>> HC 28-9-2021
      enddo                                               !!>> HC 28-9-2021


	!!!!>> HC 28-9-2021 ANGLE CLASSES!!
	!!!>> HC 28-9-2021 This divides the averange angles in classes each 0.1
      do ihc=1,napical                                      !!>> HC 28-9-2021
         do jhc=lowlimit, highlimit, steps                  !!>> HC 28-9-2021 The class -10 is form -1 to -0.9 and so on... The class 9 is form 0.9 to 1
            if (avangles(ihc)*scaling .ge. jhc ) then       !!>> HC 28-9-2021
               if (avangles(ihc)*scaling < jhc +1 ) then    !!>> HC 28-9-2021
                  classangles(ihc)=jhc                      !!>> HC 28-9-2021
               endif                                        !!>> HC 28-9-2021
            endif                                           !!>> HC 28-9-2021
         enddo                                              !!>> HC 28-9-2021
      enddo                                                 !!>> HC 28-9-2021


	!!! COOCURRENCE !!!!!>> HC 28-9-2021
	!!! This is the neightbour finding/ cooccurrence calculating loop!!!!!>> HC 28-9-2021

      order=0                                               !!>> HC 28-9-2021
      !!!>> HC 28-9-2021 This loop replaces the number of node by its curvature class
      !!!>> HC 28-9-2021 In the original lists of nodes and neigbours
      do ihc=1, finitnn                                     !!>> HC 28-9-2021
         do jhc=1, napical                                  !!>> HC 28-9-2021
	     if (foldnodes(ihc)==oldorder(jhc))then         !!>> HC 28-9-2021
                order=order+1                               !!>> HC 28-9-2021
                angles1(order)=classangles(jhc)             !!>> HC 28-9-2021
             endif                                          !!>> HC 28-9-2021
             if (foldnodesneight(ihc)==oldorder(jhc))then   !!>> HC 28-9-2021
                order2=order2+1                             !!>> HC 28-9-2021
                angles2(order2)=classangles(jhc)            !!>> HC 28-9-2021
             endif                                          !!>> HC 28-9-2021
         enddo                                              !!>> HC 28-9-2021
      enddo                                                 !!>> HC 28-9-2021

      ! We have to add 11 to transform the curvature class into the   !!>> HC 28-9-2021
      ! number of row/column in the matrix cooccurrence matrix        !!>> HC 28-9-2021
      angles1=angles1+11                                              !!>> HC 28-9-2021
      angles2=angles2+11                                              !!>> HC 28-9-2021

      !This loop assign each combination of curvatures to a cell in the coocurrence matrix   !!>> HC 28-9-2021
      do ihc=1, finitnn                                                                      !!>> HC 28-9-2021
         cooccurrences(angles1(ihc),angles2(ihc))=cooccurrences(angles1(ihc),angles2(ihc))+1 !!>> HC 28-9-2021
      enddo                                                                                  !!>> HC 28-9-2021

      !This loop gets the total number of node-neights in the embryo !!>> HC 28-9-2021
      do ihc=1, nclasseshc                                           !!>> HC 28-9-2021
         do jhc=1, nclasseshc                                        !!>> HC 28-9-2021
            totalneightshc=totalneightshc+cooccurrences(ihc,jhc)     !!>> HC 28-9-2021
         enddo                                                       !!>> HC 28-9-2021
      enddo                                                          !!>> HC 28-9-2021

      !The probabilities are the number of neights with a given combination of curvatures between the total num of neights !!>> HC 28-9-2021
      cooprobs=cooccurrences/totalneightshc                                                                                !!>> HC 28-9-2021

      if (tiponodhc==1)then                                                       !!>> HC 28-9-2021 CALCULATING JOINT ENTROPY
         !This calculates the joint entropy as -sum( pij*log(pij)) (Sole 2004)    !!>> HC 28-9-2021
         do ihc=1, nclasseshc                                                     !!>> HC 28-9-2021
            do jhc=1, nclasseshc                                                  !!>> HC 28-9-2021
	       if (jhc>ihc)cycle                                                  !!>> HC 11-11-2021
	       if (cooprobs(ihc,jhc)==0.0d0) cycle                                    !!>> HC 11-11-2021
	       if (jhc==ihc)then                                                  !!>> HC 11-11-2021
	          pij=cooprobs(ihc,jhc)                                           !!>> HC 11-11-2021
	       else                                                               !!>> HC 11-11-2021
	          pij=cooprobs(ihc,jhc)+cooprobs(jhc,ihc)                         !!>> HC 11-11-2021
               endif                                                              !!>> HC 28-9-2021
               entropy1=entropy1+(pij*log(pij)/log(2d0))                          !!>> HC 28-9-2021
            enddo                                                                 !!>> HC 28-9-2021
         enddo                                                                    !!>> HC 28-9-2021
         entropy1=-entropy1                                                       !!>> HC 28-9-2021
      else                                                                        !!>> HC 28-9-2021
         !This calculates the joint entropy as -sum( pij*log(pij)) (Sole 2004)    !!>> HC 28-9-2021
         do ihc=1, nclasseshc                                                     !!>> HC 28-9-2021
            do jhc=1, nclasseshc                                                  !!>> HC 28-9-2021
	       if (jhc>ihc)cycle                                                  !!>> HC 11-11-2021
	       if (cooprobs(ihc,jhc)==0.0d0) cycle                                !!>> HC 11-11-2021
	       if (jhc==ihc)then                                                  !!>> HC 11-11-2021
	          pij=cooprobs(ihc,jhc)                                           !!>> HC 11-11-2021
	       else                                                               !!>> HC 11-11-2021
	          pij=cooprobs(ihc,jhc)+cooprobs(jhc,ihc)                         !!>> HC 11-11-2021
               endif                                                              !!>> HC 28-9-2021
               pij=cooprobs(ihc,jhc)                                              !!>> HC 28-9-2021
               entropy2=entropy2+(pij*log(pij)/log(2d0))                          !!>> HC 28-9-2021
            enddo                                                                 !!>> HC 28-9-2021
         enddo                                                                    !!>> HC 28-9-2021
         entropy2=-entropy2                                                       !!>> HC 28-9-2021
      endif                                                                       !!>> HC 28-9-2021
  enddo                                                                           !!>> HC 28-9-2021

    entropy=(entropy2+entropy1)/2          !!>> HC 20-10-2021
    fitelli=entropy                        !!>> HC 11-11-2021
    
  end subroutine                           !!>> HC 28-9-2021
 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  !!>> HC 11-11-2021
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  !!>> HC 11-11-2021
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  !!>> HC 11-11-2021
 
  subroutine flux_resistance_decrease(fitelli)                         !!>> HC 11-11-2021 This subroutine was written by HUgo Cano and it calculates the flux resistance of a morphology
    implicit none                                                      !!>> HC 11-11-2021 calculated as the average dot product between an arbritary flux vector and the spring vectors.
    integer ::  ihc, jhc, khc, nnodes, naltech                         !!>> HC 11-11-2021
    real*8, dimension(1:3) :: flux, u                                  !!>> HC 11-11-2021
    real*8 :: upv, sumd, average_flux_resistance, modu, fitelli        !!>> HC 11-11-2021
 
    call neighbor_build                                                !!>> HC 11-11-2021

    flux(1)=0.0d0; flux(2)=0.0d0; flux(3)=1.0d0                        !!>> HC 11-11-2021 This is the FLUX VECTOR in the environment (arbitrary)
    nnodes=0; sumd=0.0d0; average_flux_resistance=0.0d0; modu=0.0d0    !!>> HC 11-11-2021 it has to be a unit vector

    do ihc=1, nd                                                       !!>> HC 11-11-2021
       if(node(ihc)%tipus.ne.1)cycle                                   !!>> HC 11-11-2021
       nnodes=nnodes+1                                                 !!>> HC 11-11-2021
       naltech=node(ihc)%altre                                         !!>> HC 11-11-2021
       u(1)=node(ihc)%x-node(naltech)%x                                !!>> HC 11-11-2021
       u(2)=node(ihc)%y-node(naltech)%y                                !!>> HC 11-11-2021
       u(3)=node(ihc)%z-node(naltech)%z                                !!>> HC 11-11-2021
       modu=sqrt(u(1)**2+u(2)**2+u(3)**2)                              !!>> HC 11-11-2021
       u(1)=u(1)/modu                                                  !!>> HC 11-11-2021 u is the unit spring vector of all the epithelial cells
       u(2)=u(2)/modu                                                  !!>> HC 11-11-2021
       u(3)=u(3)/modu                                                  !!>> HC 11-11-2021
       upv=flux(1)*u(1)+flux(2)*u(2)+flux(3)*u(3)                      !!>> HC 11-11-2021 this is the dot producti between the spring vector and the flux vector
       sumd=sumd+abs(upv)                                              !!>> HC 11-11-2021
    enddo                                                              !!>> HC 11-11-2021

    average_flux_resistance=(1-(sumd/nnodes))                          !!>> HC 11-11-2021 this is the average fux resistance
    fitelli=average_flux_resistance
      
  end subroutine                                                       !!>> HC 11-11-2021
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  !!>> HC 11-11-2021
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  !!>> HC 11-11-2021
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  !!>> HC 11-11-2021
 
  subroutine flux_resistance_increase(fitelli)                         !!>> HC 11-11-2021 This subroutine was written by HUgo Cano and it calculates the flux resistance of a morphology
    implicit none                                                      !!>> HC 11-11-2021 calculated as the average dot product between an arbritary flux vector and the spring vectors.
    integer ::  ihc, jhc, khc, nnodes, naltech                         !!>> HC 11-11-2021
    real*8, dimension(1:3) :: flux, u                                  !!>> HC 11-11-2021
    real*8 :: upv, sumd, average_flux_resistance, modu, fitelli        !!>> HC 11-11-2021
 
    call neighbor_build                                                !!>> HC 11-11-2021

    flux(1)=0.0d0; flux(2)=0.0d0; flux(3)=1.0d0                        !!>> HC 11-11-2021 This is the FLUX VECTOR in the environment (arbitrary)
    nnodes=0; sumd=0.0d0; average_flux_resistance=0.0d0; modu=0.0d0    !!>> HC 11-11-2021 it has to be a unit vector

    do ihc=1, nd                                                       !!>> HC 11-11-2021
       if(node(ihc)%tipus.ne.1)cycle                                   !!>> HC 11-11-2021
       nnodes=nnodes+1                                                 !!>> HC 11-11-2021
       naltech=node(ihc)%altre                                         !!>> HC 11-11-2021
       u(1)=node(ihc)%x-node(naltech)%x                                !!>> HC 11-11-2021
       u(2)=node(ihc)%y-node(naltech)%y                                !!>> HC 11-11-2021
       u(3)=node(ihc)%z-node(naltech)%z                                !!>> HC 11-11-2021
       modu=sqrt(u(1)**2+u(2)**2+u(3)**2)                              !!>> HC 11-11-2021
       u(1)=u(1)/modu                                                  !!>> HC 11-11-2021 u is the unit spring vector of all the epithelial cells
       u(2)=u(2)/modu                                                  !!>> HC 11-11-2021
       u(3)=u(3)/modu                                                  !!>> HC 11-11-2021
       upv=flux(1)*u(1)+flux(2)*u(2)+flux(3)*u(3)                      !!>> HC 11-11-2021 this is the dot producti between the spring vector and the flux vector
       sumd=sumd+abs(upv)                                              !!>> HC 11-11-2021
    enddo                                                              !!>> HC 11-11-2021

    average_flux_resistance=(sumd/nnodes)                              !!>> HC 11-11-2021 this is the average fux resistance
    fitelli=average_flux_resistance                                    !!>> HC 11-11-2021
  
  end subroutine                                                       !!>> HC 11-11-2021
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  !!>> HC 11-11-2021
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  !!>> HC 11-11-2021
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  !!>> HC 11-11-2021
  
  subroutine local_joint_entropy(fitelli)                                                   !!>> HC 11-11-2021 This subroutine was created by Hugo Cano
    Implicit none                                                                           !!>> HC 11-11-2021 and it calculates the local joint entropy
    real*8 :: angle, upv, modu, modv, sumangle, medianglee, entropy, entropy1                 !!>> HC 11-11-2021 of a morphology using the first and second order neighbors.
    real*8 :: entropy2, pij, entropy3, totalneightshc,fitelli                                 !!>> HC 11-11-2021
    integer ::  napicale, napicali, ihc, jhc, khc, lhc, mhc, order, order2                  !!>> HC 11-11-2021
    integer :: node2_1, numn, ndouble, neichi, ord                                          !!>> HC 11-11-2021
    integer:: tempclass, lowlimit, highlimit, steps, tiponodhc, scaling, nclasseshc, neich  !!>> HC 11-11-2021
    integer, allocatable, dimension(:) :: classangles, oldorder                             !!>> HC 11-11-2021
    integer, allocatable, dimension(:,:) :: double_neighs, trans_double_neighs              !!>> HC 11-11-2021
    integer, allocatable, dimension(:) :: double_nneighs                                    !!>> HC 11-11-2021
    real*8, dimension(1:3) :: u, v                                                            !!>> HC 11-11-2021
    real*8, allocatable, dimension(:) :: avangles                                             !!>> HC 11-11-2021
    real*8, allocatable, dimension(:,:) :: cooccurrences, cooprobs                            !!>> HC 11-11-2021
    
 !   call neighbor_build

    tiponodhc=2                       !!>> HC 11-11-2021 Type of node we are aiming 2= apical 1 = basal
    lowlimit=-10                      !!>> HC 11-11-2021 Lower limit of the interval in classes
    highlimit=9                       !!>> HC 11-11-2021 Higher limit of the interval in classes
    steps=1                           !!>> HC 11-11-2021 width of the interval in classes
    nclasseshc=20                     !!>> HC 11-11-2021 Number of classes
    scaling=10                        !!>> HC 11-11-2021 Scaling of the avangle vector to fit into the classes
    ndouble=maxval(nneigh(1:nd))*2    !!>> HC 11-11-2021

    !!>> HC 11-11-2021 This calculates the number of apical cells
    napicale=0                             !!>> HC 11-11-2021
    do ihc=1,nd                            !!>> HC 11-11-2021 This loop calculates the number of apical nodes that are going to be inspected
       if (node(ihc)%tipus==tiponodhc)then !!>> HC 11-11-2021 Select only the apical nodes
         napicale=napicale+1               !!>> HC 11-11-2021 Counter of the number of apical nodes inspected
       endif
    enddo


     !!>> HC 11-11-2021 this sets the dimension of the vector of average angles to the number of apical nodes
 
 
      if(allocated(avangles))deallocate(avangles)                           !!>> HC 11-11-2021
      allocate(avangles(napicale))                                          !!>> HC 11-11-2021

      if(allocated(classangles))deallocate(classangles)                     !!>> HC 11-11-2021
      allocate(classangles(napicale))                                       !!>> HC 11-11-2021
 
      if(allocated(oldorder))deallocate(oldorder)                           !!>> HC 11-11-2021
      allocate(oldorder(napicale))                                          !!>> HC 11-11-2021

      if(allocated(cooccurrences))deallocate(cooccurrences)                 !!>> HC 11-11-2021
      allocate(cooccurrences(nclasseshc,nclasseshc))                        !!>> HC 11-11-2021
  
      if(allocated(cooprobs))deallocate(cooprobs)                           !!>> HC 11-11-2021
      allocate(cooprobs(nclasseshc,nclasseshc))                             !!>> HC 11-11-2021

      if(allocated(double_neighs))deallocate(double_neighs)                 !!>> HC 11-11-2021
      allocate(double_neighs(napicale,ndouble))                             !!>> HC 11-11-2021

      if(allocated(trans_double_neighs))deallocate(trans_double_neighs)     !!>> HC 11-11-2021
      allocate(trans_double_neighs(napicale,ndouble))                       !!>> HC 11-11-2021

      if(allocated(double_nneighs))deallocate(double_nneighs)               !!>> HC 11-11-2021
      allocate(double_nneighs(napicale))                                    !!>> HC 11-11-2021
  
      entropy=0.0d0; entropy1=0.0d0; entropy2=0.0d0                         !!>> HC 11-11-2021
  
      do tiponodhc=1,2                                                      !!>> HC 11-11-2021
         avangles=666.0d0; classangles=666                                  !!>> HC 11-11-2021
         napicali=0                                                         !!>> HC 11-11-2021
         oldorder=0; cooccurrences=0.0d0; cooprobs=0.0d0; tempclass=0       !!>> HC 11-11-2021
         totalneightshc=0.0d0;                                              !!>> HC 11-11-2021
         pij=0.0d0; order=0                                                 !!>> HC 11-11-2021
         double_neighs=0; double_nneighs=0; trans_double_neighs=0           !!>> HC 11-11-2021
         

         !!>> HC 11-11-2021 Saving first order neighbors
         order=0                                                            !!>> HC 11-11-2021
         do ihc=1,nd                                                        !!>> HC 11-11-2021
            if(node(ihc)%tipus.ne.tiponodhc)cycle                           !!>> HC 11-11-2021
            order=order+1                                                   !!>> HC 11-11-2021
            double_neighs(order,1:nneigh(ihc))=neigh(ihc,1:nneigh(ihc))     !!>> HC 11-11-2021
            double_nneighs(order)=nneigh(ihc)                               !!>> HC 11-11-2021
         enddo                                                              !!>> HC 11-11-2021

         !!>> HC 11-11-2021 Finding second order neighbors (can be optimized...)
         order2=0                                                                                      !!>> HC 11-11-2021
         do ihc=1,nd                                                                                   !!>> HC 11-11-2021
            if(node(ihc)%tipus.ne.tiponodhc)cycle                                                      !!>> HC 11-11-2021
            order2=order2+1                                                                            !!>> HC 11-11-2021
            order=nneigh(ihc)                                                                          !!>> HC 11-11-2021
            do jhc=1, nneigh(ihc)                                                                      !!>> HC 11-11-2021
               neich=neigh(ihc,jhc)                                                                    !!>> HC 11-11-2021
               if(neich==ihc)cycle                                                                     !!>> HC 11-11-2021
               if(node(neich)%tipus.ne.tiponodhc)cycle                                                 !!>> HC 11-11-2021
               do khc=1,nneigh(neich)                                                                  !!>> HC 11-11-2021
                  neichi=neigh(neich,khc)                                                              !!>> HC 11-11-2021
                  if(node(neichi)%tipus.ne.tiponodhc)cycle                                             !!>> HC 11-11-2021
                  if(neichi==neich)cycle                                                               !!>> HC 11-11-2021
                  if(neichi==ihc)cycle                                                                 !!>> HC 11-11-2021
                  if(any(double_neighs(order2,1:double_nneighs(order2))==neichi))cycle                 !!>> HC 11-11-2021
                  order=order+1                                                                        !!>> HC 11-11-2021
                  if (order>ndouble)then                                                               !!>> HC 11-11-2021
                     if(allocated(trans_double_neighs))deallocate(trans_double_neighs)                 !!>> HC 11-11-2021
                     allocate(trans_double_neighs(napicale,ndouble))                                   !!>> HC 11-11-2021
                     trans_double_neighs=0                                                             !!>> HC 11-11-2021
                     trans_double_neighs(1:napicale,1:ndouble)=double_neighs(1:napicale,1:ndouble)     !!>> HC 11-11-2021
                     ndouble=ndouble+1                                                                 !!>> HC 11-11-2021
                     if(allocated(double_neighs))deallocate(double_neighs)                             !!>> HC 11-11-2021
                     allocate(double_neighs(napicale,ndouble))                                         !!>> HC 11-11-2021
                     double_neighs=0                                                                   !!>> HC 11-11-2021
                     double_neighs(1:napicale,1:ndouble-1)=trans_double_neighs(1:napicale,1:ndouble-1) !!>> HC 11-11-2021   
                  endif                                                                                !!>> HC 11-11-2021
                  double_neighs(order2,order)=neichi                                                   !!>> HC 11-11-2021
                  double_nneighs(order2)=double_nneighs(order2)+1                                      !!>> HC 11-11-2021
               enddo                                                                                   !!>> HC 11-11-2021
            enddo                                                                                      !!>> HC 11-11-2021
         enddo                                                                                         !!>> HC 11-11-2021


         !!>> HC 11-11-2021 This is th neightbour finding/angle calculating loop!!!
         ord=0; order=0                                    !!>> HC 11-11-2021
         do ihc=1,nd                                       !!>> HC 11-11-2021 Inspecting all the nodes i of the embryo
            sumangle=0.0d0 ; numn=0                        !!>> HC 11-11-2021 initializing      
            if (node(ihc)%tipus.ne.tiponodhc )cycle        !!>> HC 11-11-2021 We are just looking for apical cells 
            ord=ord+1                                      !!>> HC 11-11-2021
            do jhc=1,double_nneighs(ord)                   !!>> HC 11-11-2021 Comaring them with its neighbors
               neich=double_neighs(ord,jhc)                !!>> HC 11-11-2021
               if(ihc==neich)cycle                         !!>> HC 11-11-2021 Not themselves
               if(node(neich)%tipus.ne.tiponodhc)cycle     !!>> HC 11-11-2021 We are just looking for apical cells 
               v(1)=(node(neich)%x-node(ihc)%x)            !!>> HC 11-11-2021 Vector between apical i cell node and the apical j cell node
               v(2)=(node(neich)%y-node(ihc)%y)            !!>> HC 11-11-2021
               v(3)=(node(neich)%z-node(ihc)%z)	     !!>> HC 11-11-2021
               modv=sqrt( (v(1))**2+(v(2))**2+(v(3))**2)   !!>> HC 11-11-2021 modulus v
               u(1)=(node(node(ihc)%altre)%x-node(ihc)%x)  !!>> HC 11-11-2021 Apical-basal vector i cell  
               u(2)=(node(node(ihc)%altre)%y-node(ihc)%y)  !!>> HC 11-11-2021
               u(3)=(node(node(ihc)%altre)%z-node(ihc)%z)  !!>> HC 11-11-2021
               modu=sqrt((u(1))**2+(u(2))**2+(u(3))**2)    !!>> HC 11-11-2021 modulus u 
               upv=u(1)*v(1)+u(2)*v(2)+u(3)*v(3)           !!>> HC 11-11-2021 vectorial product u * v
               angle=upv/(modu*modv)                       !!>> HC 11-11-2021 getting the dot product between u and v ranges from -1 to 1 and is directly proportional to the angle
               sumangle=sumangle+angle                     !!>> HC 11-11-2021 sum of angles of the neighbouts in a cell i 
               numn=numn+1                                 !!>> HC 11-11-2021 Counter of the real number of neights (just apical nodes) of the cell i
            enddo                                          !!>> HC 11-11-2021
            if (numn.ne.0) then                            !!>> HC 11-11-2021 Just in case a cell has no neights
               napicali=napicali+1                         !!>> HC 11-11-2021 This is a counter of the number of cells i with neights
                                                           !!>> HC 11-11-2021 it will be used after to get the probability of classes
               medianglee=sumangle/real(numn)              !!>> HC 11-11-2021 This obtains the average angle of the cell i with its neightbours
               avangles(napicali)=medianglee               !!>> HC 11-11-2021
            endif                                          !!>> HC 11-11-2021
            order=order+1                                  !!>> HC 11-11-2021
            oldorder(order)=ihc                            !!>> HC 11-11-2021 Saves the old node number to use it later when calculating the cooccurrence    
         enddo                                             !!>> HC 11-11-2021
         
        
         order=0                                           !!>> HC 11-11-2021 
         !!>> HC 11-11-2021 This divides the average angles in classes each 0.1
         do ihc=1,napicale                                 !!>> HC 11-11-2021
            do jhc=lowlimit,highlimit, steps               !!>> HC 11-11-2021 The class -10 is form -1 to -0.9 and so on... The class 9 is form 0.9 to 1
               if (avangles(ihc)*scaling.ge.jhc)then       !!>> HC 11-11-2021
                  if (avangles(ihc)*scaling<jhc+1)then     !!>> HC 11-11-2021
                     classangles(ihc)=jhc                  !!>> HC 11-11-2021
                  endif                                    !!>> HC 11-11-2021
               endif                                       !!>> HC 11-11-2021
            enddo                                          !!>> HC 11-11-2021
         enddo                                             !!>> HC 11-11-2021


         !!! COOCURRENCE !!!!!>> HC 11-11-2021

         !!! This is the neightbour finding/ cooccurrence calculating loop!!! !!>> HC 11-11-2021
         v=0.0d0                                                 !!>> HC 11-11-2021
         order=0; order2=0                                       !!>> HC 11-11-2021
         do khc=lowlimit, highlimit, steps                       !!>> HC 11-11-2021 For all the classes
            order=order+1                                        !!>> HC 11-11-2021 this is the order in one of the axes of the coocurrence matrix
            do ihc=1,napicale                                    !!>> HC 11-11-2021 Inspecting all the nodes i with neightbours     
               if (classangles(ihc).ne.khc)cycle                 !!>> HC 11-11-2021 Whose classangle is equal to khc      
               do jhc=1,double_nneighs(ihc)                      !!>> HC 11-11-2021 look all the rest of the cells in the embryo
                  neich=double_neighs(ihc,jhc)                   !!>> HC 11-11-2021
                  if (oldorder(ihc)==neich)cycle                 !!>> HC 11-11-2021 except themselves (you cannot be your own neight)
                  if (node(neich)%tipus.ne.tiponodhc)cycle       !!>> HC 11-11-2021 We are just looking for apical cells 
                  do mhc=1,napicale                              !!>> HC 11-11-2021 For all the nodes whose class has been calculated         
                     if(oldorder(mhc).ne.neich)cycle             !!>> HC 11-11-2021 This correlates the old node number with the new possition in classanges vector
                     order2=0                                    !!>> HC 11-11-2021 This resets the order in the other axis of the correlation matrix
                     tempclass=classangles(mhc)                  !!>> HC 11-11-2021 This stores temporally the curvature value of the neight
                     do lhc=lowlimit, highlimit, steps           !!>> HC 11-11-2021 This loop assigns the neightbour to a position in the array of cooccurences
                        order2 = order2 +1                       !!>> HC 11-11-2021 This is the order the other axe of the coocurrence matrix
                        if (tempclass.ne.lhc)cycle               !!>> HC 11-11-2021 "for all the classes, when the class is equal to the neight class 
                                                                 !!>> HC 11-11-2021 Then it assigns it to its the corresponding place in the array"
                        cooccurrences(order,order2)=cooccurrences(order,order2)+1.0d0   !!>> HC 11-11-2021             	               
                     enddo                                       !!>> HC 11-11-2021
                  enddo	                                  !!>> HC 11-11-2021
               enddo                                             !!>> HC 11-11-2021
            enddo                                                !!>> HC 11-11-2021
         enddo                                                   !!>> HC 11-11-2021
     
         !!>> HC 11-11-2021 This loop gets the total number of node-neights in the embryo
         do ihc=1, nclasseshc                                          !!>> HC 11-11-2021
            do jhc=1, nclasseshc                                       !!>> HC 11-11-2021
               totalneightshc=totalneightshc+cooccurrences(ihc,jhc)    !!>> HC 11-11-2021
            enddo                                                      !!>> HC 11-11-2021
         enddo                                                         !!>> HC 11-11-2021

         !The probabilities are the number of neights with a given combination of curvatures between the total num of neights 
         do ihc=1, nclasseshc                                               !!>> HC 11-11-2021
            do jhc=1, nclasseshc                                            !!>> HC 11-11-2021
               cooprobs(ihc,jhc)=cooccurrences(ihc,jhc)/totalneightshc      !!>> HC 11-11-2021
            enddo                                                           !!>> HC 11-11-2021
         enddo                                                              !!>> HC 11-11-2021

         do ihc=1,nclasseshc                                                !!>> HC 11-11-2021
            do jhc=1,nclasseshc                                             !!>> HC 11-11-2021
               if(ihc==jhc)cycle                                            !!>> HC 11-11-2021
               if(ihc>jhc)cycle                                             !!>> HC 11-11-2021
               cooprobs(ihc,jhc)=cooprobs(ihc,jhc)+cooprobs(jhc,ihc)        !!>> HC 11-11-2021
               cooprobs(jhc,ihc)=0.0d0                                      !!>> HC 11-11-2021
            enddo                                                           !!>> HC 11-11-2021
         enddo                                                              !!>> HC 11-11-2021


         !!>> HC 11-11-2021 This calculates the joint entropy as -sum( pij*log(pij)) (Sole 2004)
         if (tiponodhc==1)then                                          !!>> HC 11-11-2021
            do ihc=1, nclasseshc                                        !!>> HC 11-11-2021
               do jhc=1, nclasseshc                                     !!>> HC 11-11-2021
                  if (cooprobs(ihc,jhc)==0.0d0)cycle                    !!>> HC 11-11-2021
                  pij=cooprobs(ihc,jhc)                                 !!>> HC 11-11-2021
                  entropy1=entropy1+(pij*log(pij)/log(2.0000000000))    !!>> HC 11-11-2021
               enddo                                                    !!>> HC 11-11-2021
            enddo                                                       !!>> HC 11-11-2021
            entropy1=-entropy1                                          !!>> HC 11-11-2021
         else                                                           !!>> HC 11-11-2021
            do ihc=1, nclasseshc                                        !!>> HC 11-11-2021
               do jhc=1, nclasseshc                                     !!>> HC 11-11-2021
                  if (cooprobs(ihc,jhc)==0.0d0)cycle                    !!>> HC 11-11-2021
                  pij=cooprobs(ihc,jhc)                                 !!>> HC 11-11-2021
                  entropy2=entropy2+(pij*log(pij)/log(2.0000000000))    !!>> HC 11-11-2021
               enddo                                                    !!>> HC 11-11-2021
            enddo                                                       !!>> HC 11-11-2021
            entropy2=-entropy2                                          !!>> HC 11-11-2021
         endif                                                          !!>> HC 11-11-2021

      enddo                                                             !!>> HC 11-11-2021

      entropy=(entropy2+entropy1)/2.0d0                                 !!>> HC 11-11-2021
      fitelli=entropy                                                   !!>> HC 11-11-2021

  end subroutine                                                        !!>> HC 11-11-2021
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  !!>> HC 11-11-2021
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  !!>> HC 11-11-2021
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  !!>> HC 11-11-2021
  
  subroutine surface_volume_ratio(fitelli)                                                       !!>> HC 11-11-2021 This subroutine was created by Hugo Cano
    Implicit none                                                                                !!>> HC 11-11-2021 and if calculates the surface/volume ratio 
    integer ::  ihc, jhc, khc, nnodes, naltech, surface, volume, outside, total, inside          !!>> HC 11-11-2021 of a morphology using a volume integration method
    integer ::  iihc,jjhc,kkhc, ord, lhc                                                         !!>> HC 11-11-2021
    real*8  :: suvol, fitelli                                                                    !!>> HC 11-11-2021
    real*8 :: upv, sumd, average_flux_resistance, modu                                           !!>> HC 11-11-2021
    integer, allocatable, dimension(:,:,:) :: bfillz, bfillz2, bfillx, bfillx2, bfilly, bfilly2  !!>> HC 11-11-2021
    call extrem                                                                                  !!>> HC 11-11-2021

    a=2*maxval(node(:nd)%add)                                           !!>> HC 11-11-2021 maximal interaction distance between mesenchymal nodes
    rv=a+1d-3                                                           !!>> HC 11-11-2021
    urv=1.0d0/a                                                         !!>> HC 11-11-2021

    nboxes=nint(extre*urv)+int(maxval(node(:nd)%dmo)+dmax)+1            !!>> HC 11-11-2021
    if (allocated(boxes)) deallocate(boxes)                             !!>> HC 11-11-2021
    allocate(boxes(-nboxes:nboxes,-nboxes:nboxes,-nboxes:nboxes))       !!>> HC 11-11-2021
    boxes=0                                                             !!>> HC 11-11-2021
    do ihc=1,nd                                                         !!>> HC 11-11-2021
      if(node(ihc)%tipus.ne.1)cycle                                     !!>> HC 11-11-2021
      iihc=nint(node(ihc)%x*urv);jjhc=nint(node(ihc)%y*urv);kkhc=nint(node(ihc)%z*urv)  !!>> HC 11-11-2021
      boxes(iihc,jjhc,kkhc)=ihc                                         !!>> HC 11-11-2021
    end do                                                              !!>> HC 11-11-2021

    if(allocated(bfillz))deallocate(bfillz)                             !!>> HC 11-11-2021 This saves whether the box is
    allocate(bfillz(-nboxes:nboxes,-nboxes:nboxes,-nboxes:nboxes))      !!>> HC 11-11-2021 1=surface; 2=inside (volume); 0=outside
    bfillz=666                                                          !!>> HC 11-11-2021

    if(allocated(bfillz2))deallocate(bfillz2)                           !!>> HC 11-11-2021
    allocate(bfillz2(-nboxes:nboxes,-nboxes:nboxes,-nboxes:nboxes))     !!>> HC 11-11-2021
    bfillz2=666                                                         !!>> HC 11-11-2021
    surface=0; volume=0; total=0; outside=0                             !!>> HC 11-11-2021

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                      !!>> HC 11-11-2021
    !!!!!!!!!!!!!!!!!!  Z AXIS  !!!!!!!!!!!!!!!!!!                      !!>> HC 11-11-2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                      !!>> HC 11-11-2021
    do ihc=-nboxes,nboxes                                               !!>> HC 11-11-2021 This goes over all the boxes in
       do jhc=-nboxes,nboxes                                            !!>> HC 11-11-2021 one of the directions of the z axis
          inside=0                                                      !!>> HC 11-11-2021
          do khc=-nboxes,nboxes                                         !!>> HC 11-11-2021
             if (boxes(ihc,jhc,khc)>0)then                              !!>> HC 11-11-2021 if the box is full of nodes, then it is surface
                bfillz(ihc,jhc,khc)=1                                   !!>> HC 11-11-2021
                if (inside==0)then; inside=1; else; inside=0; endif     !!>> HC 11-11-2021 and if it is the firs one we find, then we start counting
             else                                                       !!>> HC 11-11-2021 the volume (inner boxes)
                if (inside==1)then                                      !!>> HC 11-11-2021 If the box is empty, it can be
                   bfillz(ihc,jhc,khc)=2                                !!>> HC 11-11-2021 inside the embryo = volume
                else                                                    !!>> HC 11-11-2021
                   bfillz(ihc,jhc,khc)=0                                !!>> HC 11-11-2021 outside the embryo = outer space
                endif                                                   !!>> HC 11-11-2021
             endif                                                      !!>> HC 11-11-2021
          enddo                                                         !!>> HC 11-11-2021
       enddo                                                            !!>> HC 11-11-2021
    enddo                                                               !!>> HC 11-11-2021


    do ihc=-nboxes,nboxes                                               !!>> HC 11-11-2021  This goes over all the boxes in
       do jhc=-nboxes,nboxes                                            !!>> HC 11-11-2021  the other direction of the z axis
          inside=0; ord=0                                               !!>> HC 11-11-2021  an it does the same as above
          do lhc=-nboxes,nboxes                                         !!>> HC 11-11-2021
             khc=nboxes-ord                                             !!>> HC 11-11-2021
             ord=ord+1                                                  !!>> HC 11-11-2021
             if (boxes(ihc,jhc,khc)>0)then                              !!>> HC 11-11-2021
                bfillz2(ihc,jhc,khc)=1                                  !!>> HC 11-11-2021
                if (inside==0)then; inside=1; else; inside=0; endif     !!>> HC 11-11-2021
             else                                                       !!>> HC 11-11-2021
                if (inside==1)then                                      !!>> HC 11-11-2021
                   bfillz2(ihc,jhc,khc)=2                               !!>> HC 11-11-2021
                else                                                    !!>> HC 11-11-2021
                   bfillz2(ihc,jhc,khc)=0                               !!>> HC 11-11-2021
                endif                                                   !!>> HC 11-11-2021
             endif                                                      !!>> HC 11-11-2021
          enddo                                                         !!>> HC 11-11-2021
       enddo                                                            !!>> HC 11-11-2021
    enddo                                                               !!>> HC 11-11-2021


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                !!>> HC 11-11-2021
    !!!!!!!!!!!!!!!!!!  CALCULATIONS  !!!!!!!!!!!!!!!!!!                !!>> HC 11-11-2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                !!>> HC 11-11-2021

    do ihc=-nboxes,nboxes                                               !!>> HC 11-11-2021
       do jhc=-nboxes,nboxes                                            !!>> HC 11-11-2021
          do khc=-nboxes,nboxes                                         !!>> HC 11-11-2021
             total=total+1                                              !!>> HC 11-11-2021
             if(bfillz2(ihc,jhc,khc)==0) bfillz(ihc,jhc,khc)=0          !!>> HC 11-11-2021 the outer space has to be the same in both directions
             if(bfillz(ihc,jhc,khc)==0) outside=outside+1               !!>> HC 11-11-2021 we count the boxes in each category
             if(bfillz(ihc,jhc,khc)==2) volume=volume+1                 !!>> HC 11-11-2021
             if(bfillz(ihc,jhc,khc)==1) surface=surface+1               !!>> HC 11-11-2021
          enddo                                                         !!>> HC 11-11-2021
       enddo                                                            !!>> HC 11-11-2021
    enddo                                                               !!>> HC 11-11-2021
    
    if (volume>100)then                                                 !!>> HC 11-11-2021
       suvol=real(surface)/real(volume)                                 !!>> HC 11-11-2021 SURFACE/VOLUME
       fitelli=suvol                                                    !!>> HC 11-11-2021
    else                                                                !!>> HC 11-11-2021
       fitelli=0.0d0                                                    !!>> HC 11-11-2021 This is completely broken
    endif                                                               !!>> HC 11-11-2021

  
  end subroutine                                                        !!>> HC 11-11-2021
  
  
  
  
end module  
  
