!    GNOMO software (General Node Model)
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




module single_node
use general
use neighboring
use growth
use mitosis

contains

subroutine should_I_divide_single
integer ick,j,k,ii,jj,kk
real*8 a,b,c,s,sx,sy,sz,d 

do ick=1,ncels
  s=0.0d0
  ii=cels(ick)%node(1)
  if(node(ii)%hold==1) cycle !hold nodes doesn't perform cell behaviors
  if(cels(ick)%ctipus<3)then
    jj=node(ii)%altre
    c=1-0.5*(node(ii)%diffe+node(jj)%diffe)
    do k=1,npag(nparam_per_node+2)
      kk=whonpag(nparam_per_node+2,k)
      if (gex(ii,kk)>0.0d0.or.gex(jj,kk)>0.0d0) then
        s=s+(gex(ii,kk)+gex(jj,kk))*gen(kk)%wa(nparam_per_node+2)*c
      end if
    end do
    s=s*0.5 !this way %fase is independent of cell size !>>>>Miquel2-12-13
  else
    c=1-node(ii)%diffe
    do k=1,npag(nparam_per_node+2)
      kk=whonpag(nparam_per_node+2,k)
      if (gex(ii,kk)>0.0d0) then
        s=s+gex(ii,kk)*gen(kk)%wa(nparam_per_node+2)*c
      end if
    end do
  end if
  cels(ick)%fase=cels(ick)%fase+s*delta
  if (cels(ick)%fase>=1.0d0) then
    if (ffu(25)==1) then  !>>> Is 18-4-15
      if (node(cels(ick)%node(1))%req>=mmae) then !>>> Is 18-4-15 !>> Miquel20-4-15
        call division_single_node(ick)!>>> Is 18-4-15
        cels(ick)%fase=0d0  !>>> Is 18-4-15
        cels(ncels)%fase=0d0  !>>> Is 18-4-15
      end if
    else
      call division_single_node(ick)    !>>> Is 18-4-15
      cels(ick)%fase=0d0    !>>> Is 18-4-15
      cels(ncels)%fase=0d0  !>>> Is 18-4-15
    end if   !>>> Is 18-4-15
  end if  !>>> Is 18-4-15
end do

end subroutine should_I_divide_single

!**************************************************************************************************

subroutine division_single_node(celd)   ! IMPORTANT, THIS IMPLIES AN INELUDIBLE BIAS IF THE INITIAL CONDITIONS ARE REGULAR SINCE THE NEW CELL HAS TO ARISE IN A SPECIFIC DIRECTION
integer::celd,nnod,tipi,nnoda,nnodb,i,j,k,ii,jj,kk
real*8::a,b,c,d,ax,ay,az,bx,by,bz,cx,cy,cz,dotp,ix,iy,iz,nodda
integer,dimension(:)::nodea(cels(celd)%nunodes),nodeb(cels(celd)%nunodes)

  ii=cels(celd)%node(1) !the node
  tipi=node(ii)%tipus
  nodda=node(ii)%da
  !determining the physicial division plane (in relation with the neighboring cells)
  !just taking the longest vector connecting two of its neighbors
  if(tipi<3)then !epithelial : only looks for neighbors of the same %tipus
    aa=0
    bx=0;by=0;bz=0
    if(nneigh(ii)>2)then !>>Miquel20-3-15
      do i=1,nneigh(ii)-1
        jj=neigh(ii,i)
        if(node(jj)%tipus==tipi.and.(dneigh(ii,i)<nodda+node(jj)%da))then
          do j=i+1,nneigh(ii)
            kk=neigh(ii,j)
            if(node(kk)%tipus==tipi.and.dneigh(ii,j)<nodda+node(kk)%da)then
              a=node(kk)%x-node(jj)%x ; b=node(kk)%y-node(jj)%y ; c=node(kk)%z-node(jj)%z
              d=sqrt(a**2+b**2+c**2)
             if(d>aa)then
                aa=d
                bx=a ; by=b ; bz=c
              end if
            end if
          end do
        end if
      end do
    else
      call random_number(bx)
      call random_number(by)
      call random_number(bz)
      aa=sqrt(bx**2+by**2+bz**2)
    end if

    d=1/aa
    bx=bx*d ; by=by*d ; bz=bz*d
    jj=node(ii)%altre
    ax=node(jj)%x-node(ii)%x ; ay=node(jj)%y-node(ii)%y ; az=node(jj)%z-node(ii)%z
    d=1/sqrt(ax**2+ay**2+az**2)
    ax=ax*d;ay=ay*d;az=az*d
    dotp=ax*bx+ay*by+az*bz
    cels(celd)%hpolx=bx-ax*dotp ; cels(celd)%hpoly=by-ay*dotp ; cels(celd)%hpolz=bz-az*dotp    !physical vector
  else    !mesenchymal : looks for any neighbor
    aa=0
    bx=0;by=0;bz=0
    if(nneigh(ii)>1)then !>>Miquel20-3-15
      do i=1,nneigh(ii)-1
        jj=neigh(ii,i)
        if(dneigh(ii,i)<nodda+node(jj)%da)then
          do j=i+1,nneigh(ii)
            kk=neigh(ii,j)
            if(dneigh(ii,j)<nodda+node(kk)%da)then
              a=node(kk)%x-node(jj)%x ; b=node(kk)%y-node(jj)%y ; c=node(kk)%z-node(jj)%z
              d=sqrt(a**2+b**2+c**2)
              if(d>aa)then
                aa=d
                bx=a ; by=b ; bz=c
              end if
            end if
          end do
        end if
      end do
    else
      call random_number(bx)
      call random_number(by)
      call random_number(bz)
      aa=sqrt(bx**2+by**2+bz**2)
    end if
    cels(celd)%hpolx=bx ; cels(celd)%hpoly=by ; cels(celd)%hpolz=bz    !physical vector
  end if
  
  d=0.0d0
  do j=1,npag(nparam_per_node+11)    !number of genes affecting growth
    k=whonpag(nparam_per_node+11,j)  !which are those genes
    if(tipi<3)then
      d=d+(gex(ii,k)+gex(jj,k))*gen(k)%wa(nparam_per_node+11) !this is the differential of growth for the node
    else
      d=d+gex(ii,k)*gen(k)%wa(nparam_per_node+11) !this is the differential of growth for the node
    end if
  end do
  d=1d0/(1d0+d)  !ponderacion entre vector fisico i quimico
 !print*,"d ponderation",d,node(nd)%tipus
  !d=0d0                                 
  !d    ; dependence of the gradient vector (how many it affects to polarization vector)
  !d=0  ; polarization vector comes only from its shape  (default mode) 
  !d=1  ; polarization vector comes only from the gradient       
  !d=0.5; polarization vector comes equally from both the gradient and the shape    
  cx=((1-d)*cels(celd)%polx)+(d*cels(celd)%hpolx)   !vector resultante para division
  cy=((1-d)*cels(celd)%poly)+(d*cels(celd)%hpoly)
  cz=((1-d)*cels(celd)%polz)+(d*cels(celd)%hpolz)
  d=1/sqrt(cx**2+cy**2+cz**2) 
  cx=cx*d ; cy=cy*d ; cz=cz*d
  ncels=ncels+1
  if(tipi<3)then
    ndepi=ndepi+2
    nd=nd+2
    nnodb=2
    nnoda=2
    !i=jj
    !jj=ii
    !ii=i
  else
    ndmes=ndmes+1
    nd=nd+1
    nnodb=1
    nnoda=1
  end if

  call addanode(ii)
  node(nd)%icel=ncels
!  node(nd)%icel=ncels!pfh17-3-15

  node(nd)%x=node(ii)%x+desmax*cx ; node(nd)%y=node(ii)%y+desmax*cy ; node(nd)%z=node(ii)%z+desmax*cz
 !print *,cx,cy,cz,cels(celd)%polx,cels(celd)%poly,cels(celd)%polz,"koop",desmax
  node(nd)%marge=node(ii)%marge !>>> Is 1-15
  nodeo(nd)%x=node(nd)%x ; nodeo(nd)%y=node(nd)%y ; nodeo(nd)%z=node(nd)%z !pfh18-3-15

  if (ffu(25)==0) then !>>> Is 18-4-15
    node(nd)%reqcr=node(ii)%reqcr !>>> Is 18-4-15
  else                            !>>> Is 18-4-15
    node(nd)%reqcr=reqmin         !>>> Is 18-4-15
  end if                          !>>> Is 18-4-15

  if (prop_noise==0d0) node(nd)%e=0d0 !pfh 20-3-15
  if(tipi<3)then
    node(nd-1)%icel=ncels
    node(nd-1)%x=node(jj)%x+desmax*cx ; node(nd-1)%y=node(jj)%y+desmax*cy ; node(nd-1)%z=node(jj)%z+desmax*cz
    if (ffu(25)==0) then !>>> Is 18-4-15
      node(nd-1)%reqcr=node(ii)%reqcr !>>> Is 18-4-15
    else                              !>>> Is 18-4-15
      node(nd-1)%reqcr=reqmin         !>>> Is 18-4-15
    end if                            !>>> Is 18-4-15
    ! node(nd)%marge=0 >>> Is 1-15
    node(nd-1)%marge=1-node(nd)%marge !>>> Is 1-15
    nodeo(nd-1)%x=node(nd-1)%x ; nodeo(nd-1)%y=node(nd-1)%y ; nodeo(nd-1)%z=node(nd-1)%z !pfh18-3-15
    if (prop_noise==0d0) node(nd-1)%e=0d0  !pfh 20-3-15
  else
    node(nd)%marge=0
  end if

  !cels(celd)%nunodes=nnoda
  cels(ncels)%nunodes=nnodb
  !deallocate(cels(celd)%node)
  !cels(celd)%nodela=nnoda+10
  cels(ncels)%nodela=nnodb
  allocate(cels(ncels)%node(nnodb))
  cels(celd)%nunodes=nnoda
  cels(celd)%node(:)=0
  cels(celd)%node(1)=ii

  cels(ncels)%node(:)=0
  cels(ncels)%node(1)=nd

  cels(ncels)%hpolx=0 ; cels(ncels)%hpoly=0 ; cels(ncels)%hpolz=0  !>>Miquel12-9-14
  cels(ncels)%polx=0 ; cels(ncels)%poly=0 ; cels(ncels)%polz=0     !>>Miquel12-9-14

  if(tipi<3)then
    cels(celd)%node(2)=node(ii)%altre
    cels(ncels)%node(2)=nd-1
    cels(ncels)%cex=node(nd-1)%x ; cels(ncels)%cey=node(nd-1)%y ; cels(ncels)%cez=node(nd-1)%z
  else
    cels(ncels)%cex=node(nd)%x ; cels(ncels)%cey=node(nd)%y ; cels(ncels)%cez=node(nd)%z
  end if
  cels(ncels)%ctipus=cels(celd)%ctipus
  cels(ncels)%maxsize_for_div=cels(celd)%maxsize_for_div!pfh-17-3-15
  cels(ncels)%minsize_for_div=cels(celd)%minsize_for_div!pfh-17-3-15
  cels(ncels)%temt=cels(celd)%temt			!pfh-17-3-15

  if (ncels>=ncals) call recels
  
end subroutine division_single_node

!***************************************************************************************

subroutine polarization_single  !it only works for diffusible molecules, should be implemented for membrane-tethered molecules
real*8::vx,vy,vz,ix,iy,iz,ivx,ivy,ivz,ax,ay,az,d,alfa
integer::iv,celi,i,j,k,jj,kk

  !print*,"POLARISATION_SINGLE"
  do i=1,nd
    if(node(i)%tipus==2)then !epithelial
      vx=0 ; vy=0 ; vz=0
      iv=node(i)%altre
      ix=node(i)%x ; iy=node(i)%y ; iz=node(i)%z
      ivx=node(iv)%x ; ivy=node(iv)%y ; ivz=node(iv)%z
      celi=node(i)%icel
      do j=1,nneigh(i)
        jj=neigh(i,j)
        if(node(jj)%tipus/=node(i)%tipus) cycle !only compare nodes of the same layer
        ax=node(jj)%x-ix ; ay=node(jj)%y-iy ; az=node(jj)%z-iz
        d=1/dneigh(i,j)
        do k=1,npag(nparam_per_node+8)
          kk=whonpag(nparam_per_node+8,k)
          if (gex(jj,kk)>0.0d0) then
            alfa=gex(jj,kk)*gen(kk)%wa(nparam_per_node+8)*d!*gen(kk)%diffu !proportional to the diffusion differential
            vx=vx+ax*alfa ; vy=vy+ay*alfa ; vz=vz+az*alfa
          else
            cycle
          end if                                                                         !just a bit different (for simplification)
        end do
      end do
      do j=1,nneigh(iv)
        jj=neigh(iv,j)
        if(node(jj)%tipus/=node(iv)%tipus) cycle !only compare nodes of the same layer
        ax=node(jj)%x-ivx ; ay=node(jj)%y-ivy ; az=node(jj)%z-ivz
        d=1/dneigh(iv,j)
        do k=1,npag(nparam_per_node+8)
          kk=whonpag(nparam_per_node+8,k)
          if (gex(jj,kk)>0.0d0) then
            alfa=gex(jj,kk)*gen(kk)%wa(nparam_per_node+8)*d!*gen(kk)%diffu !proportional to the diffusion differential
            vx=vx+ax*alfa ; vy=vy+ay*alfa ; vz=vz+az*alfa
          else
            cycle
          end if                                                                         !just a bit different (for simplification)
        end do
      end do
      if(vx==0.and.vy==0.and.vz==0)then
        d=0
      else
        d=1/sqrt(vx**2+vy**2+vz**2)
      end if
      cels(celi)%polx=vx*d ; cels(celi)%poly=vy*d ; cels(celi)%polz=vz*d    
!      cels(celi)%polx=0.0d0 ; cels(celi)%poly=0.0d0 ; cels(celi)%polz=1.0d0 !ACHTUNG  ACHTUNG  ACHTUNG  ACHTUNG  ACHTUNG  ACHTUNG  ACHTUNG  ACHTUNG  
    elseif(node(i)%tipus==3)then  !mesenchymal
      vx=0 ; vy=0 ; vz=0
      ix=node(i)%x ; iy=node(i)%y ; iz=node(i)%z
      celi=node(i)%icel
      do j=1,nneigh(i)
        jj=neigh(i,j)
        ax=node(jj)%x-ix ; ay=node(jj)%y-iy ; az=node(jj)%z-iz
        d=1/dneigh(i,j)
        do k=1,npag(nparam_per_node+8)
          kk=whonpag(nparam_per_node+8,k)
          if (gex(jj,kk)>0.0d0) then
            alfa=gex(jj,kk)*gen(kk)%wa(nparam_per_node+8)*d!*gen(kk)%diffu !proportional to the diffusion differential
            vx=vx+ax*alfa ; vy=vy+ay*alfa ; vz=vz+az*alfa
          else
            cycle
          end if                                                                         !just a bit different (for simplification)
        end do
      end do
      if(vx==0.and.vy==0.and.vz==0)then
        d=0
      else
        d=1/sqrt(vx**2+vy**2+vz**2)
      end if
      cels(celi)%polx=vx*d ; cels(celi)%poly=vy*d ; cels(celi)%polz=vz*d
      !print*,i,"pol",cels(celi)%polx,cels(celi)%poly,cels(celi)%polz
    end if
  end do


end subroutine polarization_single

!****************************************************************************************************************

subroutine emt_singleold  !AIXO NO VA MIQUEL >>> Is 4-10-14
integer:: i,j,k,ii,jj,kk,tipi
real*8:: a,b,c,aa,bb,cc,aaa,bbb,ccc,dotp,d,dd,mean,an

  an=-pi/6
  do i=1,nd
    tipi=node(i)%tipus
    if(tipi==2)then
      ii=0 !tipus 1 counter
      jj=0 !tipus 2 counter
      mean=0
      k=node(i)%altre
      a=node(i)%x ; b=node(i)%y ; c=node(i)%z
      aaa=node(k)%x ; bbb=node(k)%y ; ccc=node(k)%z
      do j=1,nneigh(i)
        kk=neigh(i,j)
        !if(node(kk)%tipus==1)then
        !  ii=ii+1
        !  a=a+dneigh(i,j)
        !elseif(node(kk)%tipus==2)then
        !  jj=jj+1
        !  b=b+dneigh(i,j)
        !end if
        if(node(kk)%tipus==2)then
          jj=jj+1
          aa=node(kk)%x ; bb=node(kk)%y ; cc=node(kk)%z
          d=sqrt((aa-a)**2+(bb-b)**2+(cc-c)**2)
          dd=sqrt((aaa-a)**2+(bbb-b)**2+(ccc-c)**2)
          dotp=((aa-a)*(aaa-a)+(bb-b)*(bbb-b)+(cc-c)*(ccc-c))/(dd*d)
          mean=mean+dotp
        end if
      end do
      !a=a/real(ii)
      !b=b/real(jj)
      mean=mean/real(jj) ; print*,i,"TRANS mean",mean
      
      if(mean<sin(an))then
      !print*,i,"TRANS","a",a,"b",b
      !if(a<b)then !the basal node is getting out of the basal epithelial lining
      node(i)%tipus=3 ; node(k)%tipus=3
      node(i)%altre=0 ; node(k)%altre=0
      cels(node(i)%icel)%ctipus=3
        !gex(i,3)=0  !THIS IS PROVISIONAL!!!!!!!!!!!!!!!!!!!!!!!!!!
        print*,"TRANSITIOOOOOOOOOOON",i,k,"a",a,"b",b
      end if
    end if
  end do

end subroutine emt_singleold

subroutine emt_single  !epithelial-mensenchymal transitions

               ! this simply makes the tipus equal to 3 but it needs to approach both sides of the originally epithelial cell so that 
               ! the two parts do not drift away since they are more far away than da (they are at reqs normally)
               ! this transition has to be sudden otherwise other existing forces may impede the nodes from getting close enough to each other
               ! that wouldnt be realistic since it would be as if the cell would explode.

integer ick,j,k,ii,kk,iii,jjj,kkk,iv
real*8 a,b,c,d,aa,bb,cc,dd,s

e: do ick=1,ncels
  if (cels(ick)%ctipus==3) cycle
  s=0.0
  do j=1,cels(ick)%nunodes
    ii=cels(ick)%node(j)
    do k=1,npag(nparam_per_node+13)
      kk=whonpag(nparam_per_node+13,k)
      s=s+gex(ii,kk)*gen(kk)%wa(nparam_per_node+13) 
    end do
  end do
  cels(ick)%temt=cels(ick)%temt+s*delta/real(cels(ick)%nunodes)
  
  if (cels(ick)%temt>=1.0d0) then       ! wa in units of probability: this is an arbitrary value but it is fine since the rest can be re-scaled 
    !so this cell goes to emt
    do jjj=1,cels(ick)%nunodes
      iii=cels(ick)%node(jjj)          
      if (node(iii)%tipus==1) then
        iv=node(iii)%altre
        a=node(iv)%x-node(iii)%x ; b=node(iv)%y-node(iii)%y ; c=node(iv)%z-node(iii)%z 
        d=1d0/sqrt(a**2+b**2+c**2)
        a=a*d ; b=b*d ; c=c*d
        aa=0.5d0*(node(iv)%x+node(iii)%x) ; bb=0.5d0*(node(iv)%y+node(iii)%y) ; cc=0.5d0*(node(iv)%z+node(iii)%z) 
        dd=(node(iii)%da+node(iv)%da)*0.5d0
        node(iii)%x=aa-a*dd 
        node(iii)%y=bb-b*dd 
        node(iii)%z=cc-c*dd
        node(iv)%x=aa+a*dd 
        node(iv)%y=bb+b*dd 
        node(iv)%z=cc+c*dd
        node(iii)%altre=0
        node(iv)%altre=0
        node(iii)%tipus=3
        node(iv)%tipus=3
      end if
    end do
  end if
end do e

end subroutine

end module single_node
