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




!***************************************************************************
!***************  MODUL VEINATGE ********************************************
!***************************************************************************

module neighboring
use general

integer,public, allocatable  :: boxes(:,:,:),list(:)
integer,public               :: nboxes,iextre
integer,public, allocatable  :: borderface(:,:)  !>>Miquel3-2-14
integer,public               :: nborder
real*8 ,private              :: maxlen    ! >>> Is 7-2-15

contains

!**************************************************************************
subroutine extrem 
    ! troba quina es la cel que esta mes allunyada del centre del embryo
    extre=0.0d0 ; iextre=1 
    do i=1,nd
      a=sqrt((node(i)%x)**2+(node(i)%y)**2+(node(i)%z)**2)
      if (a>extre) then ; extre=a ; iextre=i ; end if ;
      !if(extre>10000) print*,"es passa",extre,iextre,"coord",node(iextre)%x,node(iextre)%y,node(iextre)%z
    end do
end subroutine extrem

!**************************************************************************
subroutine iniboxes
integer ic,ii,jj,kk
integer,allocatable:: cboxes(:,:,:)

    oextre=0	
    call extrem
 !print*,"extre",extre
    !>>> Is 7-2-15
    maxlen=sqrt((2*maxval(node(:nd)%reqs))**2+(2*maxval(node(:nd)%da)**2)) !this is the maximal interaction distance between epithelial nodes
    a=2*maxval(node(:nd)%da) !maximal interaction distance between mesenchymal nodes
    if(a>maxlen)then ; rv=a ; urv=1.0d0/a ; else ; rv=maxlen ; urv=1.0d0/maxlen ;end if
    !>>> Is 7-2-15

    nboxes=nint(extre*urv)+maxval(node(:nd)%dmo)+dmax+1
    if (allocated(list)) deallocate(list)
    allocate(list(nda))	!>>Miquel 14-10-12
    list=0
    if (allocated(boxes)) deallocate(boxes)
    if (allocated(cboxes)) deallocate(cboxes)
    allocate(boxes(-nboxes:nboxes,-nboxes:nboxes,-nboxes:nboxes))
    allocate(cboxes(-nboxes:nboxes,-nboxes:nboxes,-nboxes:nboxes))
    boxes=0
    cboxes=0
    do i=1,nd
      ii=nint(node(i)%x*urv);jj=nint(node(i)%y*urv);kk=nint(node(i)%z*urv)
      list(i)=boxes(ii,jj,kk)
      boxes(ii,jj,kk)=i
      cboxes(ii,jj,kk)=cboxes(ii,jj,kk)+1
    end do
    mnn_dyn=maxval(cboxes)
    !print*,"dynamic mnn pre",mnn_dyn
end subroutine iniboxes

!**************************************************************************
subroutine iniboxesll
    oextre=0
    call extrem

    nboxes=nint(extre*urv)+maxval(node(:nd)%dmo)+dmax+1 !;if(nboxes<7) print*,"NBOXESll",nboxes,"extre",extre
    if (allocated(list)) deallocate(list)
    allocate(list(nda))	!>>Miquel 14-10-12
    list=0
    if (allocated(boxes)) deallocate(boxes)
    allocate(boxes(-nboxes:nboxes,-nboxes:nboxes,-nboxes:nboxes))
    boxes=0

    do i=1,nd 
      ii=nint(node(i)%x*urv);jj=nint(node(i)%y*urv);kk=nint(node(i)%z*urv) !;print*,"ii",ii,"jj",jj,"kk",kk
      list(i)=boxes(ii,jj,kk)
      boxes(ii,jj,kk)=i                 
    end do
end subroutine iniboxesll

!************************************************

subroutine neighbor_build !this subroutine is called at each iteration to assess the neighbors of all the nodes !>>>Miquel24-2-14
integer:: ivv,ii1,ii2,ii3,nbo,iii1,iii2,iii3,ie,ierr,ij,ik,jk,iv
real*8:: ix,iy,iz,dist,udist,nbh

integer::tipi  !>>Miquel31-12-14
real*8::dai,maxlen    !>>Miquel31-12-14

!!!triangulation variables
integer:: npt !number of points
integer:: sizht !size of hash table
integer:: maxbf,maxfc,nbf,nfc,nface,ntetra   !size of arrays
integer,allocatable :: vecinod(:,:),vecic(:)!miguel
real*8,allocatable :: dvecinod(:,:)
real*8,allocatable :: vcl(:,:) !point coordinates
integer,allocatable :: vm(:),ht(:),bf(:,:),fc(:,:) !point indices (the algorithm reorders)
integer,dimension(:) :: border(nd)
real*8::jx,jy,jz,cx,cy,cz  !>>Miquel6-3-14
integer::o                 !>>Miquel6-3-14
!integer,dimension(:)::sneigh(mnn)
!real*8,dimension(:)::sdneigh(mnn)
integer,allocatable::osneigh(:),sneigh(:),trans_neigh(:,:)
real*8,allocatable::osdneigh(:),sdneigh(:),trans_dneigh(:,:)
integer,allocatable::oosneigh(:)
real*8,allocatable::oosdneigh(:)
integer::snneigh,error
integer::L,RRA2,IR
real*8::RRA
!integer::mnn_dynam

  !>>Miquel31-12-14
  !maxlen=maxval(node(:)%da)*2
  !maxlen=sqrt((2*maxval(node(:nd)%reqs))**2+(2*maxval(node(:nd)%da)**2)) !this is the maximal interaction distance between epithelial nodes
  !a=2*maxval(node(:nd)%da) !maximal interaction distance between mesenchymal nodes
  !if(a>maxlen)then ; rv=a ; urv=1.0d0/a ; else ; rv=maxlen ; urv=1.0d0/maxlen ;end if
  !if(a>maxlen)then ; rv=a ; else ; rv=maxlen ;end if
  !>>Miquel31-12-14

  !rv=2*maxval(node(:nd)%da);urv=1d0/rv !;print*,"RV",rv

!  a=maxval(node(:nd)%reqs)
!  if (rv<a) then ; rv=a;urv=1d0/rv ;end if !this is mostly to allow diffusion between the two faces of the epithelium
  rdiffmax=2*maxval(node(:nd)%da)*dmax

  !omnn=0
  call iniboxes
    
  mnn_dynam=mnn_dyn*(2*nint(rdiffmax*urv)+1)**3
  !print*,"mnn_dynam global",mnn_dynam,"(mnn_dyn",mnn_dyn,")"
  allocate(trans_neigh(nd,mnn_dynam),trans_dneigh(nd,mnn_dynam))
  !trans_neigh=0 ; trans_dneigh=0d0  

  if (ffu(13)==0)then !normal neighboring, extensive search of the boxes
    omnn=0
    do i=1,nd
      if (rdiffmax<2*node(i)%da)then
        nbh=2*node(i)%da  !the neighbor search range for diffusion is the same as for node interactions
      else
        nbh=rdiffmax
      end if
      
      ix=node(i)%x     ; iy=node(i)%y     ; iz=node(i)%z   
      ii1=nint(iz*urv) ; ii2=nint(iy*urv) ; ii3=nint(ix*urv)
      ivv=node(i)%altre

      nbo=nint(nbh*urv) !;print*,"nbo",nbo,"rdiffmax",rdiffmax
      mnn_dynam=mnn_dyn*(2*nbo+1)**3 !;print*,"mnn_dyn def",mnn_dynam !calculating the alleged maximal width of the neigh matrix

      allocate(sdneigh(mnn_dynam),sneigh(mnn_dynam))
      snneigh=0
      sneigh=0 ; sdneigh=0d0

      if(node(i)%tipus<3)then
        snneigh=1
        sneigh(1)=ivv
        sdneigh(1)=sqrt((node(ivv)%x-ix)**2+(node(ivv)%y-iy)**2+(node(ivv)%z-iz)**2)
      end if


      !>>> Is 18-4-15
      if (ffu(24)==0) then
        tipi=node(i)%tipus
        dai=node(i)%da
        do i1=-nbo,nbo 
          iii1=ii1+i1
          do i2=-nbo,nbo
            iii2=ii2+i2
            do i3=-nbo,nbo
              iii3=ii3+i3
              ie=boxes(iii3,iii2,iii1)
              do while(ie.ne.0)
                if (ie==i) then ; ie=list(ie) ; cycle ; end if
                if (ie==ivv) then ; ie=list(ie) ; cycle ; end if
                if (node(ie)%tipus==1.and.tipi==2) then ; ie=list(ie) ; cycle ; end if !>>> Is 18-4-15
                if (node(ie)%tipus==2.and.tipi==1) then ; ie=list(ie) ; cycle ; end if !>>> Is 18-4-15
                dist=sqrt((node(ie)%x-ix)**2+(node(ie)%y-iy)**2+(node(ie)%z-iz)**2)
                !>>Miquel31-12-14
                  a=dai+node(ie)%da
                if(tipi>=3)then
                  if(node(ie)%tipus>=3)then !mesench/ECM vs mesench/ECM
                    if(dist>a)then ; ie=list(ie) ; cycle ; end if
                  else  !mesench/ECM vs epithelial
                    b=sqrt(2*(a**2))
                    if(dist>b)then ; ie=list(ie) ; cycle ; end if
                  end if
                else !epithelial vs epithelial
                  b=sqrt(2*(a**2))
                  if(tipi==node(ie)%tipus)then !same face epithelials
                    if(b<maxlen) b=maxlen
                    if(dist>b)then ; ie=list(ie) ; cycle ; end if
                  else
                    if(dist>b)then ; ie=list(ie) ; cycle ; end if    
                  end if
                end if
                !>>Miquel31-12-14
                snneigh=snneigh+1
                sneigh(snneigh)=ie
                sdneigh(snneigh)=dist
                ie=list(ie)
              end do
            end do
          end do
        end do
      else
        tipi=node(i)%tipus
        dai=node(i)%da
        do i1=-nbo,nbo 
          iii1=ii1+i1
          do i2=-nbo,nbo
            iii2=ii2+i2
            do i3=-nbo,nbo
              iii3=ii3+i3
              ie=boxes(iii3,iii2,iii1)
              do while(ie.ne.0)
                if (ie==i) then ; ie=list(ie) ; cycle ; end if
                if (ie==ivv) then ; ie=list(ie) ; cycle ; end if
                !if (node(ie)%tipus==1.and.tipi==2) then ; ie=list(ie) ; cycle ; end if !>>> Is 18-4-15
                !if (node(ie)%tipus==2.and.tipi==1) then ; ie=list(ie) ; cycle ; end if !>>> Is 18-4-15
                dist=sqrt((node(ie)%x-ix)**2+(node(ie)%y-iy)**2+(node(ie)%z-iz)**2)
                !>>Miquel31-12-14
                  a=dai+node(ie)%da
                if(tipi>=3)then
                  if(node(ie)%tipus>=3)then !mesench/ECM vs mesench/ECM
                    if(dist>a)then ; ie=list(ie) ; cycle ; end if
                  else  !mesench/ECM vs epithelial
                    b=sqrt(2*(a**2))
                    if(dist>b)then ; ie=list(ie) ; cycle ; end if
                  end if
                else !epithelial vs epithelial
                  b=sqrt(2*(a**2))
                  if(tipi==node(ie)%tipus)then !same face epithelials
                    if(b<maxlen) b=maxlen
                    if(dist>b)then ; ie=list(ie) ; cycle ; end if
                  else
                    if(dist>b)then ; ie=list(ie) ; cycle ; end if    
                  end if
                end if
                !>>Miquel31-12-14
                snneigh=snneigh+1
                sneigh(snneigh)=ie
                sdneigh(snneigh)=dist
                ie=list(ie)
              end do
            end do
          end do
        end do
      end if

      allocate(osdneigh(snneigh),osneigh(snneigh))
      osdneigh(1:snneigh)=sdneigh(1:snneigh)
      osneigh(1:snneigh)=sneigh(1:snneigh)

   !print*,"sneigh",sneigh(1:nneigh(i))
      if(ffu(3)==1)then
        !screening by Gabriel graph !>>Miquel6-3-14
        !a neighbor connection is deleted if the sphere which diameter is the vector connecting the two nodes contains any other node
        !ordering the neighors with increasing distance

        allocate(oosdneigh(snneigh),oosneigh(snneigh))

        !sorting algorithm by selection, it's ok
        do j=1,snneigh-1  
          b=osdneigh(j)
          ii=0
          do k=j+1,snneigh
            c=osdneigh(k)
             if(b>c)then
              ii=k ; b=osdneigh(k)
            end if
          end do
          if(ii/=0)then
            !jj=osneigh(j)
            kk=osneigh(ii) ; c=osdneigh(ii) !the swap
            osneigh(ii)=osneigh(j) ; osdneigh(ii)=osdneigh(j)
            osneigh(j)=kk ; osdneigh(j)=c
          end if
        end do

        !sorting algorithm by insertion, more efficient
!        do j=2,snneigh
!          a=osdneigh(j) ; jj=osneigh(j)
!          do ii=j-1,1,-1
!            if(osdneigh(ii)<=a) goto 157
!            osdneigh(ii+1)=osdneigh(ii)
!            osneigh(ii+1)=osneigh(ii)
!          end do
!          ii=0
!157       osdneigh(ii+1)=a ; osneigh(ii+1)=jj
!        end do


        !heapsort algorithm, more efficient for large datasets
        !The index L will be decremented from its initial value during the
        !"hiring" (heap creation) phase. Once it reaches 1, the index IR 
        !will be decremented from its initial value down to 1 during the
        !"retirement-and-promotion" (heap selection) phase.
!        L=snneigh/2+1
!        IR=snneigh
!        888 continue
!        if(L > 1)then
!          L=L-1
!          RRA=osdneigh(L)
!          RRA2=osneigh(L)
!        else
!          RRA=osdneigh(IR)
!          RRA2=osneigh(IR)
!          osdneigh(IR)=osdneigh(1)
!          osneigh(IR)=osneigh(1)
!          IR=IR-1
!          if(IR.eq.1)then
!            osdneigh(1)=RRA
!            osneigh(1)=RRA2
!            goto 157
!          end if
!        end if
!        II=L
!        JJ=L+L
!899     if(JJ.le.IR)then
!          if(JJ < IR)then
!            if(osdneigh(JJ) < osdneigh(JJ+1))  JJ=JJ+1
!          end if
!          if(RRA < osdneigh(JJ))then
!            osdneigh(II)=osdneigh(JJ)
!            osneigh(II)=osneigh(JJ)
!            II=JJ; JJ=JJ+JJ
!          else
!            JJ=IR+1
!          end if
!          goto 899
!        end if
!        osdneigh(II)=RRA
!        osneigh(II)=RRA2
!        goto 888
!157     continue

        !the screening
        !ii=0 !the number of eliminated connections
        do j=snneigh,1,-1
          jj=osneigh(j)
          if(jj==ivv) cycle
          a=osdneigh(j)*0.5*screen_radius !the radius of the sphere !>>Miquel28-7-14
          jx=node(jj)%x ; jy=node(jj)%y ; jz=node(jj)%z
          cx=(ix+jx)*0.5 ; cy=(iy+jy)*0.5 ; cz=(iz+jz)*0.5 !the midpoint
          do k=j-1,1,-1
            kk=osneigh(k)
            d=sqrt((cx-node(kk)%x)**2+(cy-node(kk)%y)**2+(cz-node(kk)%z)**2)
            if(d-a<epsilod)then !there is one node within the sphere, we must delete this connection
              osneigh(j)=0
              !do l=j,snneigh-ii-1
              !  osneigh(l)=osneigh(l+1)
              !  osdneigh(l)=osdneigh(l+1)
              !end do
              !osneigh(snneigh-ii)=0
              !osdneigh(snneigh-ii)=0
              !ii=ii+1
              exit
            end if
          end do
        end do
        ii=0
        do j=1,snneigh
          jj=osneigh(j)
          if(jj/=0)then
            ii=ii+1
            oosneigh(ii)=jj ; oosdneigh(ii)=osdneigh(j)
          end if
        end do
        snneigh=ii
        !snneigh=snneigh-ii
        !if(i==1) print*,nneigh(1),"neigh1",neigh(1,1:nneigh(1))

        trans_neigh(i,1:snneigh)=oosneigh(1:snneigh)
        trans_dneigh(i,1:snneigh)=oosdneigh(1:snneigh)
      else                                                 !>>> Is 23-4-15
        trans_neigh(i,1:snneigh)=osneigh(1:snneigh)        !>>> Is 23-4-15
        trans_dneigh(i,1:snneigh)=osdneigh(1:snneigh)      !>>> Is 23-4-15
      end if                                            
      
      if(snneigh>omnn) omnn=snneigh
      nneigh(i)=snneigh
      deallocate(osdneigh,osneigh,sdneigh,sneigh)
      if (ffu(3)==1) deallocate(oosneigh,oosdneigh)

    end do
!    print*,"omnn",omnn
    !omnn=maxval(nneigh,dim=1,mask=nneigh<=nd) ; print*,"omnn",omnn
    if(allocated(neigh)) deallocate(neigh)
    if(allocated(dneigh)) deallocate(dneigh)
    allocate(neigh(nda,omnn),dneigh(nda,omnn))
    neigh(1:nd,1:omnn)=trans_neigh(1:nd,1:omnn)
    dneigh(1:nd,1:omnn)=trans_dneigh(1:nd,1:omnn)
   
  else   !*****3D triangulation neighooring***************************************************************
    
    npt=nd
!    sizht=3/2*npt
    sizht=3*npt
    maxfc=npt**2
    maxbf=npt**2
    allocate(vcl(3,npt),vm(npt))
    allocate(bf(1:3,maxbf),fc(1:7,maxfc))
    allocate(ht(0:sizht-1))
    
    
    !call extrem !not necessary, but sets some variables that later may be used in pinta and creix !>>Miquel21-3-14
    
    
    do i=1,npt  !building the input arrays for the triangualtion function
      vcl(1,i)=node(i)%x ; vcl(2,i)=node(i)%y ; vcl(3,i)=node(i)%z
      vm(i)=i
    end do !;print*,"pillat pre delau"
    !calling the triangulation subroutine
    call dtriw3(npt, sizht, maxbf, maxfc, vcl, vm, nbf, nfc, nface, ntetra, bf, fc, ht, ierr)
    !if(ierr/=0)then; print*,"error in the triangulation, avort or something",ierr,getot ; endif !call exit(24); end if
    !translating the ouptut into a neighbor matrix

    do i=1,nd !initialize neighbor matrix
!      neigh(i,1:nneigh(i))=0
      neigh(i,1:)=0  !IS 29-4-14
      nneigh(i)=0
    end do
    do j=1,nfc               ! it passes through all triangles 
      if((fc(1,j)>0).and.(fc(2,j)>0).and.(fc(3,j)>0))then ! it is a "valid" triangle
        ii=fc(1,j) ; jj=fc(2,j) ; kk=fc(3,j)
        iii=vm(ii);jjj=vm(jj);kkk=vm(kk) !; print*,"iii jjj kkk",iii,jjj,kkk
        tiii=node(iii)%tipus ; tjjj=node(jjj)%tipus ; tkkk=node(kkk)%tipus

        ivv=nneigh(iii)
        ij=0 ; ik=0 ; jk=0
        !connection iii-jjj*****
        if(ivv==0)then
          d=sqrt((vcl(1,iii)-vcl(1,jjj))**2+(vcl(2,iii)-vcl(2,jjj))**2+(vcl(3,iii)-vcl(3,jjj))**2)
          ivv=1
          trans_neigh(iii,ivv)=jjj ; trans_dneigh(iii,ivv)=d
        else
          do i=1,nneigh(iii)
            if(trans_neigh(iii,i)==jjj)then; ij=1;exit;end if
          end do
          if(ij==0)then
            d=sqrt((vcl(1,iii)-vcl(1,jjj))**2+(vcl(2,iii)-vcl(2,jjj))**2+(vcl(3,iii)-vcl(3,jjj))**2)
            ivv=ivv+1
            trans_neigh(iii,ivv)=jjj ; trans_dneigh(iii,ivv)=d
          end if
        end if
        if(ij==0)then
          !connection jjj-iii*******
          if(nneigh(jjj)==0)then ; nneigh(jjj)=1 ; trans_neigh(jjj,1)=iii ; trans_dneigh(jjj,1)=d
          else; iv=nneigh(jjj)+1 ; trans_neigh(jjj,iv)=iii ; trans_dneigh(jjj,iv)=d ; nneigh(jjj)=iv ; end if
        end if
        !connection iii-kkk*****
        do i=1,nneigh(iii)
          if(trans_neigh(iii,i)==kkk)then; ik=1;exit;end if
        end do
        if(ik==0)then
          d=sqrt((vcl(1,iii)-vcl(1,kkk))**2+(vcl(2,iii)-vcl(2,kkk))**2+(vcl(3,iii)-vcl(3,kkk))**2)
          ivv=ivv+1
          trans_neigh(iii,ivv)=kkk ; trans_dneigh(iii,ivv)=d
          !connection kkk-iii*******
          if(nneigh(kkk)==0)then ; nneigh(kkk)=1 ; trans_neigh(kkk,1)=iii ; trans_dneigh(kkk,1)=d
          else; iv=nneigh(kkk)+1 ; trans_neigh(kkk,iv)=iii ; trans_dneigh(kkk,iv)=d ; nneigh(kkk)=iv ; end if
        end if
        nneigh(iii)=ivv
        !connection jjj-kkk*****
        do i=1,nneigh(jjj)
          if(trans_neigh(jjj,i)==kkk)then; jk=1;exit;end if
        end do
        if(jk==0)then
          d=sqrt((vcl(1,jjj)-vcl(1,kkk))**2+(vcl(2,jjj)-vcl(2,kkk))**2+(vcl(3,jjj)-vcl(3,kkk))**2)
          iv=nneigh(jjj)+1
          trans_neigh(jjj,iv)=kkk ; trans_dneigh(jjj,iv)=d ;nneigh(jjj)=iv
          !connection kkk-jjj*******
          iv=nneigh(kkk)+1 ; trans_neigh(kkk,iv)=jjj ; trans_dneigh(kkk,iv)=d ; nneigh(kkk)=iv
        end if
      end if
    end do


    if(ffu(3)==1)then
    
      do i=1,nd
        !Screening by Gabriel graph !>>Miquel6-3-14
        !A neighbor connection is deleted if the sphere which diameter is the vector connecting the two nodes contains any other node
        !Ordering the neighors with increasing distance

        !sorting algorithm by selection, it's ok
        do j=1,nneigh(i)-1  
          b=trans_dneigh(i,j)
          ii=0
          do k=j+1,nneigh(i)
            c=trans_dneigh(i,k)
             if(b>c)then
              ii=k ; b=trans_dneigh(i,k)
            end if
          end do
          if(ii/=0)then
            !jj=osneigh(j)
            kk=trans_neigh(i,ii) ; c=trans_dneigh(i,ii) !the swap
            trans_neigh(i,ii)=trans_neigh(i,j) ; trans_dneigh(i,ii)=trans_dneigh(i,j)
            trans_neigh(i,j)=kk ; trans_dneigh(i,j)=c
          end if
        end do
        
        !the screening
        ii=0 !the number of eliminated connections
        do j=nneigh(i),1,-1
          jj=trans_neigh(i,j)
          if(jj==ivv) cycle
          a=trans_dneigh(i,j)*0.5*screen_radius !the radius of the sphere
          jx=node(jj)%x ; jy=node(jj)%y ; jz=node(jj)%z
          cx=(ix+jx)*0.5 ; cy=(iy+jy)*0.5 ; cz=(iz+jz)*0.5 !the midpoint
          do k=j-1,1,-1
            kk=trans_neigh(i,k)
            d=sqrt((cx-node(kk)%x)**2+(cy-node(kk)%y)**2+(cz-node(kk)%z)**2)
            if(d<a)then !there is one node within the sphere, we must delete this connection
              do l=j,nneigh(i)-ii-1
                trans_neigh(i,l)=trans_neigh(i,l+1)
                trans_dneigh(i,l)=trans_dneigh(i,l+1)
              end do
              trans_neigh(i,nneigh(i)-ii)=0
              trans_dneigh(i,nneigh(i)-ii)=0
              ii=ii+1
              exit
            end if
          end do
        end do
        nneigh(i)=nneigh(i)-ii
      end do
    !else                                                 !>>> Is 23-4-15
    !  trans_neigh(i,1:snneigh)=osneigh(1:snneigh)        !>>> Is 23-4-15
    !  trans_dneigh(i,1:snneigh)=osdneigh(1:snneigh)      !>>> Is 23-4-15      
    end if
    
    omnn=0
    do i=1,nd
      if(nneigh(i)>omnn) omnn=nneigh(i)
    end do
    !print*,"omnn",omnn
    if(allocated(neigh)) deallocate(neigh)
    if(allocated(dneigh)) deallocate(dneigh)
    allocate(neigh(nda,omnn),dneigh(nda,omnn))
    neigh(1:nd,1:omnn)=trans_neigh(1:nd,1:omnn)
    dneigh(1:nd,1:omnn)=trans_dneigh(1:nd,1:omnn)
  end if
!!do i=1,nd
!print *,i,neigh(i,:nneigh(i)),"ui"
!end do  
end subroutine neighbor_build
!************************************************

subroutine neighbor_build_node(i) !this calculates the neighbors for one node onlly, used for random noise !>>Miquel16-12-14
integer:: ivv,ii1,ii2,ii3,nbo,iii1,iii2,iii3,ie,ierr,ij,ik,jk,iv
real*8:: ix,iy,iz,dist,udist,nbh
integer::i

!!!triangulation variables
integer:: npt !number of points  
integer:: sizht !size of hash table
integer:: maxbf,maxfc,nbf,nfc,nface,ntetra   !size of arrays
integer,allocatable :: vecinod(:,:),vecic(:)!miguel
real*8,allocatable :: dvecinod(:,:)
real*8,allocatable :: vcl(:,:) !point coordinates
integer,allocatable :: vm(:),ht(:),bf(:,:),fc(:,:) !point indices (the algorithm reorders)
integer,dimension(:) :: border(nd)
real*8::jx,jy,jz,cx,cy,cz  !>>Miquel6-3-14
integer::o                 !>>Miquel6-3-14
real*8::dai,maxlen
integer,allocatable::osneigh(:),sneigh(:),trans_neigh(:,:)
real*8,allocatable::osdneigh(:),sdneigh(:),trans_dneigh(:,:)
integer::snneigh,error
integer::L,RRA2,IR
real*8::RRA
integer::mnn_dynam


  maxlen=sqrt((2*maxval(node(:nd)%reqs))**2+(2*maxval(node(:nd)%da)**2)) !this is the maximal interaction distance between epithelial nodes
  a=2*maxval(node(:nd)%da) !maximal interaction distance between mesenchymal nodes
  if(a>maxlen)then ; rv=a ; urv=1.0d0/a ; else ; rv=maxlen ; urv=1.0d0/maxlen ;end if


!  rv=2*maxval(node(:nd)%da);urv=1d0/rv !;print*,"RV",rv
!  a=maxval(node(:nd)%reqs)
!  if (rv<a) then ; rv=a;urv=1d0/rv ;end if !this is mostly to allow diffusion between the two faces of the epithelium
  rdiffmax=2*maxval(node(:nd)%da)*dmax
  
  call iniboxes
 
    if (rdiffmax<2*node(i)%da)then
      nbh=2*node(i)%da  !the neighbor search range for diffusion is the same as for node interactions
    else
      nbh=rdiffmax
    end if
    ix=node(i)%x     ; iy=node(i)%y     ; iz=node(i)%z   
    ii1=nint(iz*urv) ; ii2=nint(iy*urv) ; ii3=nint(ix*urv)
    ivv=node(i)%altre
    
    nbo=nint(nbh*urv) !;print*,"nbo",nbo,"rdiffmax",rdiffmax
    mnn_dynam=mnn_dyn*(2*nbo+1)**3 !;print*,"mnn_dyn def",mnn_dynam !calculating the alleged maximal width of the neigh matrix
    
    do j=1,nneigh(i)                            !!This to erase node i from the neighbor matrix, as it will be ovewritten later  !! >>>Miguel17-12-14
      k=neigh(i,j)                              ! "k" neighbor
      do ii=1,nneigh(k)                         ! it searches in the neighborhood of "k"         
        if(i.eq.neigh(k,ii))then                ! if it (i) appears    
          neigh(k,ii:nneigh(k)-1)=neigh(k,ii+1:nneigh(k))   ; neigh(k,nneigh(k))=0   !  neigh matrix is displaced
          dneigh(k,ii:nneigh(k)-1)=dneigh(k,ii+1:nneigh(k)) ; dneigh(k,nneigh(k))=0  !  dneigh matrix is diaplced
          nneigh(k)=nneigh(k)-1 ; exit                                               !  neighbor counter
        end if 
      end do        
    end do
    neigh(i,:)=0 ; nneigh(i)=0 ; dneigh(i,:)=0d0  !! >>>Miguel17-12-14
    
    allocate(sdneigh(mnn_dynam),sneigh(mnn_dynam))
    snneigh=0
    sneigh=0 ; sdneigh=0d0
    
    
    tipi=node(i)%tipus
    if(node(i)%tipus<3)then
      snneigh=1
      sneigh(1)=ivv 
      sdneigh(1)=sqrt((node(ivv)%x-ix)**2+(node(ivv)%y-iy)**2+(node(ivv)%z-iz)**2)
    end if

    !>>> Is 18-4-15
    if (ffu(24)==0) then
      dai=node(i)%da
      do i1=-nbo,nbo 
        iii1=ii1+i1
        do i2=-nbo,nbo
          iii2=ii2+i2
          do i3=-nbo,nbo
            iii3=ii3+i3
            ie=boxes(iii3,iii2,iii1)
            do while(ie.ne.0)
              if (ie==i) then ; ie=list(ie) ; cycle ; end if
              if (ie==ivv) then ; ie=list(ie) ; cycle ; end if   
              if (node(ie)%tipus==1.and.tipi==2) then ; ie=list(ie) ; cycle ; end if !>>> Is 18-4-15
              if (node(ie)%tipus==2.and.tipi==1) then ; ie=list(ie) ; cycle ; end if !>>> Is 18-4-15
              dist=sqrt((node(ie)%x-ix)**2+(node(ie)%y-iy)**2+(node(ie)%z-iz)**2)
              !>>Miquel31-12-14
              a=dai+node(ie)%da
              if(tipi>=3)then
                if(node(ie)%tipus>=3)then !mesench/ECM vs mesench/ECM
                  if(dist>a)then ; ie=list(ie) ; cycle ; end if
                else  !mesench/ECM vs epithelial
                  b=sqrt(2*(a**2))
                  if(dist>b)then ; ie=list(ie) ; cycle ; end if
                end if
              else !epithelial vs epithelial
                b=sqrt(2*(a**2))
                if(tipi==node(ie)%tipus)then !same face epithelials
                  if(b<maxlen) b=maxlen
                  if(dist>b)then ; ie=list(ie) ; cycle ; end if
                else
                  ie=list(ie) ; cycle !>>> Is 18-4-15
                end if
              end if
              !>>Miquel31-12-14
              snneigh=snneigh+1
              sneigh(snneigh)=ie
              !dneigh(i,nneigh(i))=sqrt((node(ie)%x-ix)**2+(node(ie)%y-iy)**2+(node(ie)%z-iz)**2)
              sdneigh(snneigh)=dist
              ie=list(ie)
            end do
          end do
        end do
      end do
      !<<< Is 18-4-15
    else
      dai=node(i)%da
      do i1=-nbo,nbo 
        iii1=ii1+i1
        do i2=-nbo,nbo
          iii2=ii2+i2
          do i3=-nbo,nbo
            iii3=ii3+i3
            ie=boxes(iii3,iii2,iii1)
            do while(ie.ne.0)
              if (ie==i) then ; ie=list(ie) ; cycle ; end if
              if (ie==ivv) then ; ie=list(ie) ; cycle ; end if   
              dist=sqrt((node(ie)%x-ix)**2+(node(ie)%y-iy)**2+(node(ie)%z-iz)**2)
              !>>Miquel31-12-14
              a=dai+node(ie)%da
              if(tipi>=3)then
                if(node(ie)%tipus>=3)then !mesench/ECM vs mesench/ECM
                  if(dist>a)then ; ie=list(ie) ; cycle ; end if
                else  !mesench/ECM vs epithelial
                  b=sqrt(2*(a**2))
                  if(dist>b)then ; ie=list(ie) ; cycle ; end if
                end if
              else !epithelial vs epithelial
                b=sqrt(2*(a**2))
                if(tipi==node(ie)%tipus)then !same face epithelials
                  if(b<maxlen) b=maxlen
                  if(dist>b)then ; ie=list(ie) ; cycle ; end if
                else
                  if(dist>b)then ; ie=list(ie) ; cycle ; end if
                end if
              end if
              !>>Miquel31-12-14
              snneigh=snneigh+1
              sneigh(snneigh)=ie
              !dneigh(i,nneigh(i))=sqrt((node(ie)%x-ix)**2+(node(ie)%y-iy)**2+(node(ie)%z-iz)**2)
              sdneigh(snneigh)=dist
              ie=list(ie)
            end do
          end do
        end do
      end do
    end if

    allocate(osdneigh(snneigh),osneigh(snneigh))
    osdneigh(1:snneigh)=sdneigh(1:snneigh)
    osneigh(1:snneigh)=sneigh(1:snneigh)
    
    if((ffu(3)==1).or.(ffu(13)==1))then !! >>>Miguel17-12-14
      !screening by Gabriel graph !>>Miquel6-3-14
      !a neighbor connection is deleted if the sphere which diameter is the vector connecting the two nodes contains any other node
      do j=1,snneigh-1  !ordering the neighors with increasing distance
        b=osdneigh(j)
        ii=0
        do k=j+1,snneigh
          c=osdneigh(k)
           if(b>c)then
            ii=k ; b=osdneigh(k)
          end if
        end do
        if(ii/=0)then
          !jj=neigh(i,j)
          kk=osneigh(ii) ; c=osdneigh(ii) !the swap
          osneigh(ii)=osneigh(j) ; osdneigh(ii)=osdneigh(j)
          osneigh(j)=kk ; osdneigh(j)=c
        end if
      end do
      if(ffu(13).eq.1)then;b=0.75d0;else;b=screen_radius;endif ! >>>Miguel17-12-14
      !the screening
      ii=0 !the number of eliminated connections
      do j=snneigh,1,-1
        jj=osneigh(j)
        if(jj==ivv) cycle
        a=osdneigh(j)*0.5*b !the radius of the sphere !>>Miquel28-7-14   ! >>>Miguel17-12-14       
        jx=node(jj)%x ; jy=node(jj)%y ; jz=node(jj)%z
        cx=(ix+jx)*0.5 ; cy=(iy+jy)*0.5 ; cz=(iz+jz)*0.5 !the midpoint
        do k=j-1,1,-1
          kk=osneigh(k)
          d=sqrt((cx-node(kk)%x)**2+(cy-node(kk)%y)**2+(cz-node(kk)%z)**2)
          if(d-a<epsilod)then !there is one node within the sphere, we must delete this connection
            do l=j,snneigh-ii-1
              osneigh(l)=osneigh(l+1)
              osdneigh(l)=osdneigh(l+1)
            end do
            osneigh(snneigh-ii)=0
            osdneigh(snneigh-ii)=0
            ii=ii+1
            exit
          end if
        end do
      end do
      snneigh=snneigh-ii       
    else                                                 !>>> Is 23-4-15
      trans_neigh(i,1:snneigh)=osneigh(1:snneigh)        !>>> Is 23-4-15
      trans_dneigh(i,1:snneigh)=osdneigh(1:snneigh)      !>>> Is 23-4-15 
    end if
    
    if(snneigh>omnn)then !the whole matrix has to be reallocated
      ii=omnn
      allocate(trans_neigh(nd,omnn),trans_dneigh(nd,omnn))
      trans_neigh(1:nd,1:omnn)=neigh(1:nd,1:omnn)
      trans_dneigh(1:nd,1:omnn)=dneigh(1:nd,1:omnn)
      deallocate(neigh,dneigh)
      omnn=snneigh
      allocate(neigh(nda,omnn),dneigh(nda,omnn))
      neigh(1:nd,1:ii)=trans_neigh(1:nd,1:ii)
      dneigh(1:nd,1:ii)=trans_dneigh(1:nd,1:ii)
      deallocate(trans_neigh,trans_dneigh)
    end if
    
    nneigh(i)=snneigh
    neigh(i,1:snneigh)=osneigh(1:snneigh)
    dneigh(i,1:snneigh)=osdneigh(1:snneigh)
    
    do j=1,snneigh !here we put i on the reciprocal neighborhood of its neighbors
      k=osneigh(j)
      nneigh(k)=nneigh(k)+1
      if(nneigh(k)>omnn)then  !the whole matrix has to be reallocated
        ii=omnn
        allocate(trans_neigh(nd,omnn),trans_dneigh(nd,omnn))
        trans_neigh(1:nd,1:omnn)=neigh(1:nd,1:omnn)
        trans_dneigh(1:nd,1:omnn)=dneigh(1:nd,1:omnn)
        deallocate(neigh,dneigh)
        omnn=nneigh(k)
        allocate(neigh(nda,omnn),dneigh(nda,omnn))
        neigh(1:nd,1:ii)=trans_neigh(1:nd,1:ii)
        dneigh(1:nd,1:ii)=trans_dneigh(1:nd,1:ii)
        deallocate(trans_neigh,trans_dneigh)
      end if
      neigh(k,nneigh(k))=i
      dneigh(k,nneigh(k))=osdneigh(j)
    end do

end subroutine neighbor_build_node

end module neighboring
