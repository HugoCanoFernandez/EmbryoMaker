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




module biomechanic	!>>>>>>>> by Miquel 17-5-13

use general
use genetic
use neighboring
use shell       ! miguel4-11-13
!use nexus       ! Is 3-1-14

contains

subroutine iterdiferencial
integer::nodmo,i,j,k,ii
real*8::a,b,c
real*8::ox,oy,oz !miguel4-1-13

  !CALCULATING FORCES AND MOVEMENT VECTORS

  call forces
  !print*,"start iteration",getot
  !print*,"check*** p1",px(1),py(1),pz(1)
  !print*,"check*** p2",px(2),py(2),pz(2)
  !print*,"check*** p3",px(3),py(3),pz(3)
  !print*,getot,"maxval fmeanl",maxval(fmeanl(1:nd))


  a=0.0d0
  do i=1,nd
    if (dex(i)>a) then ; a=dex(i) ; ii=i ; end if
  end do
  
  if (a<epsilod) then ; rtime=rtime+delta ; delta=deltamin; return ; end if  !>>> Is 29-8-13
  if(ffu(12)==0)then
    delta=resmax/a      !delta is adjusted so the node that has the greatest force
                      !will always move the same absolute distance, and the rest will
                      !move proportionally to that force
  else
    delta=deltamin
!    delta=deltamax 
  end if
                      
  if (delta>deltamax) delta=deltamax
  if (delta<deltamin) delta=deltamin
end subroutine iterdiferencial

!***************************************************************************************************

subroutine rungekutta4(d)  ! Runge-Kutta fourth order integration
real*8 d,halfd,sixthd
real*8 ox(nd),oy(nd),oz(nd)
real*8 kux(nd),kuy(nd),kuz(nd)
real*8 kdx(nd),kdy(nd),kdz(nd)
real*8 ktx(nd),kty(nd),ktz(nd)
real*8 kqx(nd),kqy(nd),kqz(nd)

halfd=d*0.5d0
sixthd=d/6.0d0

ox=node(:nd)%x ; oy=node(:nd)%y ; oz=node(:nd)%z 

!k1
kux=px(:nd) ; kuy=py(:nd) ; kuz=pz(:nd)

!k2
node(:nd)%x=node(:nd)%x+halfd*px(:nd)
node(:nd)%y=node(:nd)%y+halfd*py(:nd)
node(:nd)%z=node(:nd)%z+halfd*pz(:nd)

if(nd>1) call neighbor_build
call forces

kdx=px(:nd) ; kdy=py(:nd) ; kdz=pz(:nd)

!k3
node(:nd)%x=node(:nd)%x+halfd*px(:nd)
node(:nd)%y=node(:nd)%y+halfd*py(:nd)
node(:nd)%z=node(:nd)%z+halfd*pz(:nd)

if(nd>1) call neighbor_build
call forces

ktx=px(:nd) ; kty=py(:nd) ; ktz=pz(:nd)

!k4 
node(:nd)%x=node(:nd)%x+d*px(:nd)
node(:nd)%y=node(:nd)%y+d*py(:nd)
node(:nd)%z=node(:nd)%z+d*pz(:nd)

if(nd>1) call neighbor_build
call forces

kqx=px(:nd) ; kqy=py(:nd) ; kqz=pz(:nd)

!final
node(:nd)%x=ox+sixthd*(kux+2*kdx+2*ktx+kqx)
node(:nd)%y=oy+sixthd*(kuy+2*kdy+2*kty+kqy)
node(:nd)%z=oz+sixthd*(kuz+2*kdz+2*ktz+kqz)

end subroutine

!***************************************************************************************************

subroutine adaptive_rungekutta
real*8 ox(nd),oy(nd),oz(nd)
real*8 aux(nd),auy(nd),auz(nd)
real*8 adx(nd),ady(nd),adz(nd)
real*8 r,halfdelta,invdelta,mdx,mdy,mdz,suggesteddelta

37 continue

halfdelta=0.5d0*delta
invdelta=1d0/delta

ox=node(:nd)%x ; oy=node(:nd)%y ; oz=node(:nd)%z

call rungekutta4(delta)

aux=node(:nd)%x ; auy=node(:nd)%y ; auz=node(:nd)%z
node(:nd)%x=ox  ; node(:nd)%y=oy  ; node(:nd)%z=oz

call rungekutta4(halfdelta)
call rungekutta4(halfdelta)

adx=node(:nd)%x ; ady=node(:nd)%y ; adz=node(:nd)%z

mdx=maxval(abs(aux-adx)*invdelta)
mdy=maxval(abs(auy-ady)*invdelta)
mdz=maxval(abs(auz-adz)*invdelta)

!mdx=maxval(abs(aux-adx)/ox)
!mdy=maxval(abs(auy-ady)/oy)
!mdz=maxval(abs(auz-adz)/oz)

if (mdx>=mdy.and.mdx>=mdz) r=mdx
if (mdy>=mdz.and.mdy>=mdz) r=mdy
if (mdz>=mdy.and.mdz>=mdx) r=mdz

suggesteddelta=0.9d0*prec*delta/r

if (r>prec) then !the step is too large
  delta=suggesteddelta
!print *,"NO",delta,r,prec
  goto 37
else
!print *,"SI",delta,r,prec
  ! in theory we should make delta=suggested delta but we prefer delta to be decided based on resmax in each step
  node(:nd)%x=adx
  node(:nd)%y=ady
  node(:nd)%z=adz
end if

end subroutine

!***************************************************************************************************

subroutine forces
real*8   ::ix,iy,iz,dd
real*8   ::a,b,c,d,e,f,g
integer  ::celi,celj,nod
real*8   ::youe,repe,adhe,adho,repcele,deqe,ideqe
real*8   ::younod,repnod,adhnod,repcelnod,reqnod,tornod,stornod    !>>>> Miquel 16-8-13
real*8   ::ax,ay,az,bx,by,bz
real*8   ::ud,udd,uddd
real*8   ::cx,cy,cz,ccx,ccy,ccz,dotp,pesco
real*8   ::icx,icy,icz,idd,iudd,id !>>>>>>>> MIQUEL 4-3-13
real*8   ::mcx,mcy,mcz,md,umd					!>>>>>>>> MIQUEL 30-4-13
real*8   ::nodda,posca
integer  ::ivv		!>>>>>>>> MIQUEL 4-3-13
integer  ::i,j,ii,jj,kk,ic,iii,jjj,kkk,iiii,jjjj,kkkk,iv,kjjj,jkkk
integer  ::nuve,inuve
integer  ::tipi,tipic																!>>>>>>>>>>>>>>>>>>>>>Miquel 23-4-13
integer  ::switch,twoep

integer  ::lateral,vertical !flags that tell if there is a lateral or vertical component to take into account !>>Miquel28-1-14

real*8   ::rvx,rvy,rvz   !the resulting force vector
real*8   ::uvx,uvy,uvz   !unit vector
real*8   ::pox,poy,poz   !polarisation vector (from the cell)

real*8   ::ad,fd  !>>Miquel28-1-14

real*8   ::upr(nd)
integer  ::iupr(nd)
real*8   ::r(nd),er(nd)

real*8   ::rcilx(nd),rcily(nd),rcilz(nd)
real*8   ::rtorx(nd),rtory(nd),rtorz(nd)
real*8   ::rstorx(nd),rstory(nd),rstorz(nd)
real*8   ::rsprx,rspry,rsprz


!integer, parameter :: mnn=200  !declared in general module !>>Miquel24-2-14

                              ! ACHTUNG, WE ASSUME THAT THE MAXIMAL NUMBER OF NODES INTERACTING WITH A NODE IS mnn, otherwise these 
                              ! matrices become unbareably large for large systems
!integer  ::vei(omnn,nd)        ! force storing arrays for optimization    >>>>>>Miquel 7-8-16
real*8   ::vforce(omnn,nd)     ! so we don't calculate each interaction twice
real*8   ::vuvx(omnn,nd)       ! (for now we do it only for cylinder and normal interactions,
real*8   ::vuvy(omnn,nd)       ! not torsion nor springs)
real*8   ::vuvz(omnn,nd)       !
real*8   ::vtorforce(omnn,nd)  ! NOTE THE NON-CONVENTIONAL FORM OF THE ARRAYS:
real*8   ::vtoruvx(omnn,nd)    ! DIM 1 IS THE LISTS THE NEIGHBORS,
real*8   ::vtoruvy(omnn,nd)    ! DIM 2 LISTS THE NODE ID
real*8   ::vtoruvz(omnn,nd)    !
real*8   ::vstorforce(omnn,nd) !
real*8   ::vstoruvx(omnn,nd)   !
real*8   ::vstoruvy(omnn,nd)   !
real*8   ::vstoruvz(omnn,nd)   !
!integer ::nveins(nd)         !
integer  ::epinveins(nd)      ! we store how many same-side epithelial neighbors an epithelial node has !>>>Miquel4-4-14
integer  ::alone              ! 0 if the node is really alone

  vcilx=0 ; vcily=0 ; vcilz=0 ; vsprx=0     !force vectors storage matrices, for different components
  vtorx=0 ; vtory=0 ; vtorz=0 ; vspry=0     !of the resulting force
  vstorx=0 ; vstory=0 ; vstorz=0 ; vsprz=0
  vei=0 ; vforce=0 ; vuvx=0 ; vuvy=0 ; vuvz=0 ; nveins=0  !>>>>>>Miquel 7-8-16
  vtorforce=0 ; vtoruvx=0 ; vtoruvy=0 ; vtoruvz=0         !
  vstorforce=0 ; vstoruvx=0 ; vstoruvy=0 ; vstoruvz=0     !

  fmeanl=0 ; fmeanv=0  !storage vector that makes the balance between compressive and tensile forces within a node    !>>>Miquel23-1-14

  rcilx=0d0  ; rcily=0d0  ; rcilz=0d0  !they store the force components for all the nodes !>>Miquel4-4-14
  rtorx=0d0  ; rtory=0d0  ; rtorz=0d0
  rstorx=0d0 ; rstory=0d0 ; rstorz=0d0
  epinveins=0
  
  do nod=1,nd
    !lonely=0 !fossile? >>Miquel27-2-14

    if (node(nod)%hold==2) then                ! >>> Is 30-6-14
      dex(nod)=0  !module of the force vector ! >>> Is 30-6-14
      px(nod)=0 ; py(nod)=0 ; pz(nod)=0       ! >>> Is 30-6-14
      cycle                                   ! >>> Is 30-6-14
    end if                                    ! >>> Is 30-6-14

    rsprx=0.0d0  ; rspry=0.0d0; rsprz=0.0d0
    ix=node(nod)%x ; iy=node(nod)%y ; iz=node(nod)%z
    tipi=node(nod)%tipus ; celi=node(nod)%icel
    iii=nint(ix*urv)    ; jjj=nint(iy*urv)   ; kkk=nint(iz*urv)	
    rvx=0d0    ; rvy=0d0    ; rvz=0d0
    nuve=0!nuve=nveins(nod) !>>>>>Miquel 7-8-13 : this is not always zero since we fill this matrix from its neighbours
    switch=0
    alone=0

    !SPRINGS
    if(tipi<3)then
      iv=node(nod)%altre
      ax=node(iv)%x   ; ay=node(iv)%y    ; az=node(iv)%z
      cx=ax-ix        ; cy=ay-iy         ; cz=az-iz
      dd=sqrt(cx**2+cy**2+cz**2)
      udd=1d0/dd
      !if(ffu(1)==1)then
        uvx=cx*udd ; uvy=cy*udd ; uvz=cz*udd  !the unit vector
        ddd=dd-node(nod)%reqs-node(iv)%reqs  !>>Miquel5-2-14
        f=2*node(nod)%ke*ddd    !the force
        rsprx=f*uvx ; rspry=f*uvy ; rsprz=f*uvz
        fmeanv(nod)=ddd !>>>Miquel5-2-14
      !end if
    else
      iv=0    !>>>>Miquel17-1-14
    end if
    younod=node(nod)%you                                            !>>>>>Miquel 16-8-13
    repnod=node(nod)%rep                                            !
    adhnod=node(nod)%adh !default inespecific adhesion of node nodmo!
    repcelnod=node(nod)%repcel                                      !
    tornod=node(nod)%tor                                            !
    stornod=node(nod)%stor                                          !
    reqnod=node(nod)%req                                            !
    nodda=node(nod)%da

    !NODE's REPULSIONS AND ADHESIONS

    do i=1,nneigh(nod)
      ic=neigh(nod,i)
      if(ic<nod.and.node(ic)%hold/=2)then !we have calculated that interaction already !>>Miquel6-5-15
        alone=1
        cycle
      end if
      !do iiii=1,nuve                                     !>>>>>>>>>>>Miquel 7-8-13
      !  if(ic==vei(iiii,nod))then                        !This means that nod and ic have already interacted 
      !    switch=0                                       !so we do not need to make all the calculations again
      !    f=vforce(iiii,nod)                             !
      !    rcilx(iiii)=rcilx(iiii)+f*vuvx(iiii,nod)       !IS THIS RIGHT, do we need rcilx in its own sum??????????????
      !    rcily(iiii)=rcily(iiii)+f*vuvy(iiii,nod)       !
      !    rcilz(iiii)=rcilz(iiii)+f*vuvz(iiii,nod)       !
      !    f=vtorforce(iiii,nod)                          !
      !    rtorx(iiii)=rtorx(iiii)+f*vtoruvx(iiii,nod)    !
      !    rtory(iiii)=rtory(iiii)+f*vtoruvy(iiii,nod)    !
      !    rtorz(iiii)=rtorz(iiii)+f*vtoruvz(iiii,nod)    !
      !    f=vstorforce(iiii,nod)                         !
      !    rstorx(iiii)=rstorx(iiii)+f*vstoruvx(iiii,nod) !
      !    rstory(iiii)=rstory(iiii)+f*vstoruvy(iiii,nod) !
      !    rstorz(iiii)=rstorz(iiii)+f*vstoruvz(iiii,nod) !
      !    alone=1
      !    goto 324 !this sends you at the end of the loop!
      !  end if                                           !
      !end do                                             !

      !so it turns out that the nod-ic interactions has not been calculated before

      bx=node(ic)%x   ; by=node(ic)%y    ; bz=node(ic)%z
      ccx=bx-ix       ; ccy=by-iy        ; ccz=bz-iz		
      !d=sqrt(ccx**2+ccy**2+ccz**2)
      d=dneigh(nod,i)
      ud=1d0/d
      tipic=node(ic)%tipus
      twoep=0
      ad=0 ; fd=0 ; ddd=0     !>>Miquel28-1-14

      alone=1      !this is crappy but fast, it makes that lonely nodes are eliminated in squares
      if (tipi<3)then
        if(tipic<3)then
          ivv=node(ic)%altre
          icx=node(ivv)%x-bx ; icy=node(ivv)%y-by ; icz=node(ivv)%z-bz 
          idd=sqrt(icx**2+icy**2+icz**2)      ; iudd=1d0/idd	!ic's spring vector	
          posca=icx*cx+icy*cy+icz*cz
          if (tipi==tipic) then           ! BOTH NODES ARE EPITHELIAL IN THE SAME SIDE we are equal so we must have lateral adhesion
            if (posca>epsilod) then    !SO WE THE NEIGHBOUR IS IN A CONTIGUOUS PART OF THE EPITHELIUM
              !switch=1      !lateral interaction
              !icx=node(ivv)%x-bx ; icy=node(ivv)%y-by ; icz=node(ivv)%z-bz 
              !idd=sqrt(icx**2+icy**2+icz**2)      ; iudd=1d0/idd	!ic's spring vector	
              mcx=icx+cx; mcy=icy+cy; mcz=icz+cz; md=1d0/sqrt(mcx**2+mcy**2+mcz**2)  !we take as spring vector the sum of ic's and nod's spring vectors
              dotp=(mcx*ccx+mcy*ccy+mcz*ccz)
              ddd=abs(dotp)*md  !vertical component
              !a=nodda+node(ic)%da
              ad=d**2-ddd**2 ; if(ad<epsilod) ad=epsilod
              ad=sqrt(ad)              !lateral component
              if (ad-nodda-node(ic)%da>epsilod) cycle
              pesco=dotp*md**2
              uvx=(ccx-mcx*pesco) ; uvy=(ccy-mcy*pesco) ; uvz=(ccz-mcz*pesco)  !unit vector of the within cilinder force !el modul és el mateix que el vector c
              !a=1d0/sqrt(uvx**2+uvy**2+uvz**2)
              a=1/ad
              uvx=uvx*a ; uvy=uvy*a ; uvz=uvz*a
              fd=ad  !fd as the distance we will use to assert the force modulus  !>>Miquel28-1-14
              twoep=1
              epinveins(nod)=epinveins(nod)+1
              epinveins(ic)=epinveins(ic)+1
            else
              ! apical/basal contact, two epithelia from the same side
              if(ffu(16)==0)then
                mcx=icx+cx; mcy=icy+cy; mcz=icz+cz; md=sqrt(mcx**2+mcy**2+mcz**2)  !we take as spring vector the sum of ic's and nod's spring vectors
                if(md<epsilod)then !the two cylinders are parallel, the vector used is the spring !>>Miquel20-8-14
                  dotp=(cx*ccx+cy*ccy+cz*ccz)
                  if (dotp<epsilod) then         ! projection of the vector from nod to ic into the vector from nod to iv
                    ddd=-dotp*udd                 ! that is the distance UP in the direction of altre
                    a=nodda+node(ic)%da
                    if (ddd-a<epsilod) then
                      ddd=d**2-ad**2 ;if(ddd<epsilod) ddd=epsilod
                      ddd=sqrt(ddd)              !lateral component
                      if (ad-a>epsilod) cycle
                      uvx=-cx*udd ; uvy=-cy*udd ; uvz=-cz*udd  !unit vector vertical
                      fd=ddd     !distance used, vertical component
                      twoep=2
                      goto 300
                    end if
                    cycle 
                  end if  
                else  !proper apical cylindric interface applied !>>Miquel20-8-14
                  md=1d0/md
                  dotp=(mcx*ccx+mcy*ccy+mcz*ccz)
                  ad=abs(dotp)*md  !vertical component
                  a=nodda+node(ic)%da
                  if(ad-a<epsilod)then
                    ddd=d**2-ad**2 ;if(ddd<epsilod) ddd=epsilod
                    ddd=sqrt(ddd)              !lateral component
                    if (ddd-a>epsilod) cycle
                    pesco=dotp*md**2
                    uvx=(ccx-mcx*pesco) ; uvy=(ccy-mcy*pesco) ; uvz=(ccz-mcz*pesco)  !unit vector of the within cilinder force !el modul és el mateix que el vector c
                    a=1/ddd
                    uvx=uvx*a ; uvy=uvy*a ; uvz=uvz*a
                    fd=ddd     !distance used, vertical component
                    twoep=2
                    goto 300
                  end if
                  cycle 
                end if  
              else  !apical-apical contact from the same face, sphere-sphere interface used !>>Miquel20-8-14
                if(d-nodda-node(ic)%da<epsilod)then
                  uvx=ccx*ud ; uvy=ccy*ud ; uvz=ccz*ud
                  fd=d
                  twoep=2
                  goto 300
                end if
                cycle 
              end if  

            end if
          else
            if (posca<0.0d0) then           
              ! we are from contrary sides so we must have face adhesion
              ! this is adhesion and repulsion between epithelial sides
              dotp=(cx*ccx+cy*ccy+cz*ccz)                     
              if (dotp<epsilod) then         ! projection of the vector from nod to ic into the vector from nod to iv
                ddd=-dotp*udd                 ! that is the distance UP in the direction of altre
                a=nodda+node(ic)%da
                if (ddd-a<epsilod) then
                  ad=d**2-ddd**2 ;if(ad<epsilod) ad=epsilod
                  ad=sqrt(ad)              !lateral component
                  if (ad-a>epsilod) cycle
                  uvx=cx*udd ; uvy=cy*udd ; uvz=cz*udd  !unit vector vertical
                  fd=ddd     !distance used, vertical component
                  twoep=2
                  goto 300
                end if
              end if
            end if
            cycle
          end if
        else
          ! nod IS EPITHELIAL and ic is not
          ! we check the distance to radial distance in the plane of the ic cylinder 
          dotp=(cx*ccx+cy*ccy+cz*ccz)!;print*,"dotp epi-mes",dotp,nod,ic !hi ha un vector del revés, per tant això està al revés també
          if (dotp<0.0) then
            ddd=abs(dotp*udd)    !vertical component
            a=nodda+node(ic)%da
            if (ddd-a<epsilod) then
              ad=d**2-ddd**2 ;if(ad<epsilod) ad=epsilod
              ad=sqrt(ad)              !lateral component
              if (ad-a>epsilod) cycle
              !uvx=ccx*ud ; uvy=ccy*ud ; uvz=ccz*ud  
              !fd=d  !now interactions between epith. and mesench. work in a sphere-sphere manner !>>Miquel21-2-14
              
              uvx=-cx*udd ; uvy=-cy*udd ; uvz=-cz*udd  !unit vector vertical
              fd=ddd     !distance used, vertical component
            else
              cycle
            end if
          else
            cycle
          end if
        end if
      else
        if(tipic<3)then          ! IC IS EPITHELIAL and nod is not
          ivv=node(ic)%altre     ! we check the distance to radial distance in the plane of the ic cylinder 
          icx=node(ivv)%x-bx ; icy=node(ivv)%y-by ; icz=node(ivv)%z-bz   
          idd=1d0/sqrt(icx**2+icy**2+icz**2)
          dotp=(icx*ccx+icy*ccy+icz*ccz)  !hi ha un vector del revés, per tant això està al revés també
          if (dotp>0.0) then
            ddd=dotp*idd        !vertical component
            a=nodda+node(ic)%da
            if (ddd-a<0.0) then
              ad=d**2-ddd**2 ;if(ad<epsilod) ad=epsilod
              ad=sqrt(ad)              !lateral component
              if (ad-a>epsilod) cycle
              !uvx=ccx*ud ; uvy=ccy*ud ; uvz=ccz*ud
              !fd=d  !now interactions between epith. and mesench. work in a sphere-sphere manner !>>Miquel21-2-14
              
              uvx=icx*idd ; uvy=icy*idd ; uvz=icz*idd  !unit vector vertical
              fd=ddd     !distance used, vertical component
            else
              cycle
            end if
          else
            cycle
          end if
        else
          fd=d !BOTH NODES ARE NON-EPITHELIAL: we just consider the interactions between nodes as such
          if (fd-nodda-node(ic)%da>epsilod) cycle
          uvx=ccx*ud ; uvy=ccy*ud ; uvz=ccz*ud  !unit vector of the within cilinder force !el modul és el mateix que el vector c
    !if(nod==260.and.ic==267) print*,"mesdist",d,ccx,ccy,ccz
    !if(nod==260.and.ic==267) print*,"pos260",node(nod)%x,node(nod)%y,node(nod)%z,"pos267",node(ic)%x,node(ic)%y,node(ic)%z
        end if
      end if

300   nuve=nuve+1

      if (nuve>mnn) then; 
        print *,"PANIC!!! PANIC!!!PANIC!!!; this is going to crash because there is too many " ; 
        print *,"neighbors per a node to interact, and vuvx matrix and those are too small"
      end if 
  
      !nveins(nod)=nuve ; 
      !vuvx(nuve,nod)=uvx ; vuvy(nuve,nod)=uvy ; vuvz(nuve,nod)=uvz   !>>>>>Miquel 7-8-13

      !if (twoep/=2) then 
      !  vei(nuve,nod)=ic         
      !  jjjj=nveins(ic)+1   !>>>>>Miquel 7-8-13 We save that nod and ic are neighbours
      !  nveins(ic)=jjjj ; 
      !  vuvx(jjjj,ic)=-uvx ; vuvy(jjjj,ic)=-uvy ; vuvz(jjjj,ic)=-uvz   !>>>>>Miquel 7-8-13
      !  vei(jjjj,ic)=nod                            !>>>>>Miquel 7-8-13 
      !end if
      !r(nuve)=fd

      !ALL THAT WAS JUST TO CALCULATE THE RIGHT DISTANCE BETWEEN NODES, fd, NOW WE CALCULATE THE ACTUAL ENERGIES
      !if (ffu(2)==1) then IS 25-12-13
        if(node(nod)%icel==node(ic)%icel)then
          deqe=(reqnod+node(ic)%req)             !>>>> Miquel 16-8-13
          if(fd-deqe<-epsilod)then 				
            repe=0.5d0*(repnod+node(ic)%rep)     !>>>> Miquel 16-8-13
            f=2*repe*(fd-deqe)
          else
            youe=(younod+node(ic)%you)          !>>>> Miquel 16-8-13
            f=youe*(fd-deqe)
          end if
        else
          deqe=(reqnod+node(ic)%req)          !>>>> Miquel 16-8-13
          if(fd-deqe<-epsilod)then 
            repcele=(repcelnod+node(ic)%repcel) !>>>> Miquel 16-8-13
            f=repcele*(fd-deqe)       !in fact that is the derivative of the energy
          else
            adhe=0.5d0*(adhnod+node(ic)%adh)          !>>>> Miquel 16-8-13
            if (npag(1)>0) then ! we have adhesion molecules
              do j=1,npag(1)
                k=whonpag(1,j)
                do kjjj=1,npag(1)
                  jkkk=whonpag(1,kjjj)
                  if (gex(nod,k)>0.0d0.and.gex(ic,jkkk)>0.0d0) then     
                    adhe=adhe+gex(nod,k)*gex(ic,jkkk)*kadh(int(gen(k)%wa(1)),int(gen(jkkk)%wa(1))) ! >>> Is 7-6-14    !this is specific adhesion
                  end if
                end do
              end do
            end if
            f=2*adhe*(fd-deqe)          !in fact that is the derivative of the energy
          end if
        end if !;print*,"f",f,"fd",fd
        !vforce(nuve,nod)=f     !>>>>>Miquel 7-8-13
        !if (twoep/=2) vforce(jjjj,ic)=f      !>>>>>Miquel 7-8-13
        rcilx(nod)=rcilx(nod)+f*uvx ; rcily(nod)=rcily(nod)+f*uvy ; rcilz(nod)=rcilz(nod)+f*uvz  !>>>Miquel 7-8-13
        rcilx(ic)=rcilx(ic)-f*uvx ; rcily(ic)=rcily(ic)-f*uvy ; rcilz(ic)=rcilz(ic)-f*uvz  !>>>Miquel 4-4-14
        !if(tipi<3 .and. tipic<3)then 
           
        fmeanl(nod)=fmeanl(nod)+(fd-deqe) ; fmeanl(ic)=fmeanl(ic)+(fd-deqe) ! >>> Is 21-6-14

        !endif  !tensile forces positive !>>Miquel23-1-14
        !er(nuve)=f
      !end if

     !if(nod==260.or.ic==260) print*,"dzero",nod,ic,fd

      !TORSION
      if(ffu(4)==0 .and.twoep==1 ) then !it is only between epithelial nodes
        !if(switch==0)then	!this is something that we need here but is calculated in the internal repulsion section, if that is off, you have to calculate it here
        !  ivv=node(ic)%altre
        !  icx=node(ivv)%x-bx ; icy=node(ivv)%y-by ; icz=node(ivv)%z-bz; 
        !  idd=sqrt(icx**2+icy**2+icz**2)      ; iudd=1d0/idd	!ic's spring vector	
        !else
        !  switch=0
        !end if

        !surface tension-like torsion (original)
        mcx=icx+cx; mcy=icy+cy; mcz=icz+cz; md=1d0/sqrt(mcx**2+mcy**2+mcz**2)  !we take as spring vector the sum of ic's and nod's spring vectors
        dotp=(mcx*ccx+mcy*ccy+mcz*ccz)*md !vertical projection, more stable than the angle
        !cal fer la mitja de les 2 molles, la d'ic i la de nod, perquè despres no hi hagi assimetries
        if(abs(dotp)-angletor*d>epsilod)then
          uvx=mcx*md ; uvy=mcy*md ; uvz=mcz*md
          f=(tornod+node(ic)%tor)*dotp         !>>>>>Miquel 7-8-13
          rstorx(nod)=rstorx(nod)+f*uvx ; rstory(nod)=rstory(nod)+f*uvy ; rstorz(nod)=rstorz(nod)+f*uvz
          rstorx(ic)=rstorx(ic)-f*uvx ; rstory(ic)=rstory(ic)-f*uvy ; rstorz(ic)=rstorz(ic)-f*uvz
          !vstorforce(nuve,nod)=f    !>>>>>Miquel 7-8-13
          !vstorforce(jjjj,ic)=f     !>>>>>Miquel 7-8-13
          !vstoruvx(nuve,nod)=uvx ; vstoruvy(nuve,nod)=uvy ; vstoruvz(nuve,nod)=uvz   !>>>>>Miquel 7-8-13
          !vstoruvx(jjjj,ic)=-uvx ; vstoruvy(jjjj,ic)=-uvy ; vstoruvz(jjjj,ic)=-uvz   !>>>>>Miquel 7-8-13
          !fmeanv(nod)=fmeanv(nod)-f ; fmeanv(ic)=fmeanv(ic)-f  !stor also contributes to apical-basal compression !>>Miquel28-1-14
!        end if


        !parallel springs torsion   !new  !>>Miquel14-3-14
        dotp=(cx*ccx+cy*ccy+cz*ccz)*udd !vertical projection, more stable than the angle
        !if(abs(dotp)-angltormax>epsilod)then
 !       if(abs(dotp)-0.3>epsilod)then
          f=(tornod+node(ic)%tor)*dotp
          uvx=ccx*ud ; uvy=ccy*ud ; uvz=ccz*ud
          rtorx(nod)=rtorx(nod)+f*uvx ; rtory(nod)=rtory(nod)+f*uvy ; rtorz(nod)=rtorz(nod)+f*uvz
          !vtorforce(nuve,nod)=f
          !vtoruvx(nuve,nod)=uvx ; vtoruvy(nuve,nod)=uvy ; vtoruvz(nuve,nod)=uvz
  !  end if
        !this is the neighbor, this force is not simmetrical, but it's no problem
        dotp=-(icx*ccx+icy*ccy+icz*ccz)*iudd !vertical projection, more stable than the angle
   !     if(abs(dotp)-0.3>epsilod)then
          f=(tornod+node(ic)%tor)*dotp
          rtorx(ic)=rtorx(ic)-f*uvx ; rtory(ic)=rtory(ic)-f*uvy ; rtorz(ic)=rtorz(ic)-f*uvz
          !vtorforce(jjjj,ic)=f
          !vtoruvx(jjjj,ic)=-uvx ; vtoruvy(jjjj,ic)=-uvy ; vtoruvz(jjjj,ic)=-uvz
   end if
        !end if
        !parallel springs torsion   !old 
        !dotp=abs(cx*icx+cy*icy+cz*icz)*udd*iudd
        !f=(tornod+node(ic)%tor)*dotp         !>>>>>Miquel 7-8-13
        !uvx=cx-icx ; uvy=cy-icy ; uvz=cz-icz ; ddd=sqrt(uvx**2+uvy**2+uvz**2)
        !if(ddd>0)then
        !  ddd=1d0/ddd
        !  uvx=uvx*ddd ; uvy=uvy*ddd ; uvz=uvz*ddd     !>>>>>Miquel 7-8-13
        !  rtorx(nuve)=rtorx(nuve)+f*uvx ; rtory(nuve)=rtory(nuve)+f*uvy ; rtorz(nuve)=rtorz(nuve)+f*uvz
        !  vtorforce(nuve,nod)=f    !>>>>>Miquel 7-8-13
        !  vtorforce(jjjj,ic)=f     !>>>>>Miquel 7-8-13
        !  vtoruvx(nuve,nod)=uvx ; vtoruvy(nuve,nod)=uvy ; vtoruvz(nuve,nod)=uvz   !>>>>>Miquel 7-8-13
        !  vtoruvx(jjjj,ic)=-uvx ; vtoruvy(jjjj,ic)=-uvy ; vtoruvz(jjjj,ic)=-uvz   !>>>>>Miquel 7-8-13
        !end if


      end if
324   continue
    end do

    if(epinveins(nod)>0)then                  !>>>Miquel24-3-14
      !fmeanl(nod)=0
    !else
      fmeanl(nod)=fmeanl(nod)/epinveins(nod)
    end if

789 if (ffu(8)==1) then   !this is kind of crappy in the sense that 
      if (alone==0) then
        node(nod)%talone=node(nod)%talone+1
      else
        node(nod)%talone=0
      end if
    end if

    ! IF THERE IS NOT SCREENING THINGS ARE SIMPLER
    !rvx=rsprx ; rvy=rspry ; rvz=rsprz
    !do i=1,nuve
    !  rvx=rvx+rcilx(i)+rtorx(i)+rstorx(i) ; rvy=rvy+rcily(i)+rtory(i)+rstory(i) ; rvz=rvz+rcilz(i)+rtorz(i)+rstorz(i)
    !  vcilx(nod)=vcilx(nod)+rcilx(i)    ; vcily(nod)=vcily(nod)+rcily(i)    ; vcilz(nod)=vcilz(nod)+rcilz(i)         !this part should be turned down
    !  vtorx(nod)=vtorx(nod)+rtorx(i)    ; vtory(nod)=vtory(nod)+rtory(i)    ; vtorz(nod)=vtorz(nod)+rtorz(i)         !when there is no display
    !  vstorx(nod)=vstorx(nod)+rstorx(i) ; vstory(nod)=vstory(nod)+rstory(i) ; vstorz(nod)=vstorz(nod)+rstorz(i)!
    !end do

    if (aut==0) then ! if aut==0 there is not display   
      vsprx(nod)=rsprx ; vspry(nod)=rspry ; vsprz(nod)=rsprz                     !putting the force components into storage vectors for display
      vcilx(nod)=rcilx(nod)    ; vcily(nod)=rcily(nod)    ; vcilz(nod)=rcilz(nod)         !this part should be turned down
      vtorx(nod)=rtorx(nod)    ; vtory(nod)=rtory(nod)    ; vtorz(nod)=rtorz(nod)         !when there is no display
      vstorx(nod)=rstorx(nod)  ; vstory(nod)=rstory(nod)  ; vstorz(nod)=rstorz(nod)!
    end if

    if (epinveins(nod)>0) then  ! >>> Is 7-6-14
      a=epinveins(nod) ! >>> Is 7-6-14
      a=1.0d0/a        ! >>> Is 7-6-14
      rvx=rsprx+rcilx(nod)+(rtorx(nod)+rstorx(nod))*a  ! >>> Is 7-6-14 !summing the force components into a net force vector
      rvy=rspry+rcily(nod)+(rtory(nod)+rstory(nod))*a  ! >>> Is 7-6-14
      rvz=rsprz+rcilz(nod)+(rtorz(nod)+rstorz(nod))*a  ! >>> Is 7-6-14
    else
      rvx=rsprx+rcilx(nod)  ! >>> Is 7-6-14 !summing the force components into a net force vector
      rvy=rspry+rcily(nod)  ! >>> Is 7-6-14
      rvz=rsprz+rcilz(nod)  ! >>> Is 7-6-14
    end if

    !Physical boundaries force  !>>>>>Miquel9-1-14
    if(node(nod)%hold==1)then
      uvx=node(nod)%orix-ix
      uvy=node(nod)%oriy-iy
      uvz=node(nod)%oriz-iz
      d=uvx**2+uvy**2+uvz**2 ; if(d>epsilod)then ; ud=1d0/sqrt(d);else;d=0;end if
      a=node(nod)%khold
      rvx=rvx+uvx*a*ud    !now the force is constant        !no need to calculate the unit vector, because here you have the product of the unit vector and the distance,
      rvy=rvy+uvy*a*ud    !to make it a spring remove ud    !wich is the original vector
      rvz=rvz+uvz*a*ud
    end if

    dex(nod)=sqrt(rvx**2+rvy**2+rvz**2)  !module of the force vector
    px(nod)=rvx ; py(nod)=rvy ; pz(nod)=rvz
  end do
  
end subroutine forces

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ordenarepeq(ma,mt,rang)
  integer rang
  real*8 ma(rang)
  integer mt(rang)
  integer i,j,k
  real*8 a
    mt=0
el: do i=1,rang
      a=ma(i) ; k=1
      do j=1,rang ; if (a>ma(j)) k=k+1 ; end do 
      do j=k,rang ; if (mt(j)==0) then ; mt(j)=i ; cycle el ; end if ; end do
    end do el 
end subroutine ordenarepeq


end module biomechanic
