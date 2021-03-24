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




module nexus

  use general
  use genetic
  use growth  !>>>> Miquel 17-6-13
  use death    !>>>> Miquel 18-6-13
  use ecm !>>>> Is 31-8-13 
  use single_node
  use mitosis


contains

!******************************************************************

! Declare the interface for POSIX fsync function

subroutine nexe

!interface
!  function fsync (fd) bind(c,name="fsync")
!  use iso_c_binding, only: c_int
!     integer(c_int), value :: fd
!     integer(c_int) :: fsync
!   end function fsync
!end interface


  integer i,j,k,ii,jj,kk,ik,ikk,ick,ret
  real*8 a
  real*8 differ
  real*8 kplast,kvol,dvol
  character*300 cx

! GENETIC REGULATION OF CELL BEHAVIOURS

  ! extracellular matrix secretion

  if (npag(nparam_per_node+4)>0) then ; call should_I_secrete   ; end if

  if(ffu(1)==0)then

    ! cell polarization  THIS ONE SHOULD BE THE FIRST IN HERE
    if (npag(nparam_per_node+8)>0) then ; call polarization ; end if
    !nparam_per_node+9 is to tell that growth is polarized 
    ! cell growth
    if (npag(nparam_per_node+1)>0) then ; call should_I_grow        ; end if  !It considers also polar growth in it, with nparam_per_node+9

    !cell division
    !if (npag(nparam_per_node+2)>0) then ; 
    call should_I_divide ; 
    !end if
    !if (ffu(9)==1) then                                        !>>> Is 5-2-14                                        !>>> Is 5-2-14
    !  do ick=1,ncels                                           !>>> Is 5-2-14
    !    if (cels(ick)%nunodes>cels(ick)%maxsize_for_div) then  !>>> Is 5-2-14 If the cell has too many nodes it divides si o si
    !       call division(ick)                                   !>>> Is 5-2-14
    !      cels(ick)%fase=cels(ick)%fase-1d0                    !>>> Is 5-2-14
    !    end if                                                 !>>> Is 5-2-14
    !  end do                                                   !>>> Is 5-2-14
    !end if                                                     !>>> Is 5-2-14

    ! cell apoptosis
    if (npag(nparam_per_node+3)>0.or.npag(nparam_per_node+14)>0) then ; call should_I_die  ; end if

    ! change the size of the cell required for division
    if (npag(nparam_per_node+10)>0) then; call change_minsize_for_div ; end if

    !nparam_per_node+11 is to orient division according to the chemical polarization and not according to the physical (hertwig) one

    !nparam_per_node+12 the larger the most asymetric (in mass) is the plane of division

    ! change the maximal number of nodes per cell before the cell divides
    if (npag(nparam_per_node+15)>0) then; call change_maxsize_for_div ; end if  !>>> Is 23-3-14

  else
    if (ffu(25)==1) then ; if (npag(nparam_per_node+1)>0) then ; call should_I_grow ; end if  ; end if  !>>> Is 18-4-15 !It considers also polar growth in it, with nparam_per_node+9 ACHTUNG
    if (npag(nparam_per_node+8)>0) then ; call polarization_single ; end if
    if (npag(nparam_per_node+2)>0) then ; call should_I_divide_single ; end if
    ! cell apoptosis
    if (npag(nparam_per_node+3)>0.or.npag(nparam_per_node+14)>0) then ; call should_I_die  ; end if  !>>Miquel20-3-14
    if (npag(nparam_per_node+13)>0) then; call emt ; end if   ! >>> Is 4-10-14
   !call emt_single
  end if

  ! epithelial-mesenchymal transition
  if (npag(nparam_per_node+13)>0) then; call emt ; end if

  if (ffu(5)==1) then ; agex(:nd,1)=abs(node(:nd)%x) ;end if; !external source trick >>> Is 29-6-14
  if (ffu(5)==2) then ; agex(:nd,1)=node(:nd)%z**4    ;end if!((1.1d0+maxval(node(:nd)%z)+node(:nd)%z)**4) !external source trick miguel >>> Is 29-6-14

  if(ffu(5)==3)then !horizontal gradient (x-y plane) centered on the origin of coordinates
    do i=1,nd
      agex(i,6)=1/(1+(node(i)%x**2+node(i)%y**2)) !; print*,"i",i,"gex",agex(i,5)
    end do
  end if  
  
  if (ffu(15)==1)then; !especial conditions, gene 1 is always expressing in the borders !>>Miquel12-5-14
    do i=1,nd
      if(node(i)%hold==1)then
        agex(i,1)=1d0
        !if(node(i)%tipus<3)then
        !  agex(i,1)=1d0
        !else
        !  agex(i,4)=1d0
        !end if
      end if
    end do
  end if

  if (ffu(8)==1.and.nd>2) then   ! IS 23-4-13 this eliminates the nodes that get alone for too long
    ik=1
    do while(ik<=nd) !;print*,"node(",ik,")%talone=",node(ik)%talone
      if (node(ik)%talone>ttalone) then !;print*,ik,"entra mort",node(ik)%tipus
        if (node(ik)%tipus>2) then  !we only delete an epithelial node if its altre is also lonely
          ikk=node(ik)%icel !;print*,"entra mesenq"
          call apoptosis(ik)
          if (ikk>0) then   !notice that if ikk is negative it means that is not a celular node
            if(cels(ikk)%nunodes==0) call eliminate_cell(ikk)  !we eliminate the cell if it has no nodes left
          end if
        else
          if (node(node(ik)%altre)%talone>ttalone) then
            ikk=node(ik)%icel
            call apoptosis(ik)
            if (ikk>0) then   !notice that if ikk is negative it means that is not a celular node
              if(cels(ikk)%nunodes==0) call eliminate_cell(ikk)  !we eliminate the cell if it has no nodes left
            end if
          end if
        end if
        ik=ik+1 ! Is it? >>> Is 16-1-14
      else
        ik=ik+1   ! I know, it is kind of funky but it should be this way, a loop with nd wont do because nd decreases because of apoptosis
      end if
    end do
!    call cellbreak ! to see if the cell is split in two 
  end if

! GENETIC REGULATION OF NODE PROPERTIES
  do i=1,nd         ! we update that parameter in each cell that expresses the gene
                    ! WE ONLY UPDATES DE NODES IN WHICH THE GENE IS EXPRESSED, OTHERWISE WE LEAVE IT IS AS IT WAS 
    if(node(i)%hold==1) cycle  !we don't want the border cells to perform behaviours because that would alter and possibly break the border  !>>>>Miquel9-1-14
    differ=1-node(i)%diffe !;print*,"differ",differ

    ! DIFFERENTIATION
    if (npag(25)>0) then
      ii=25;a=0.0; 
      do k=1,npag(ii) 
        kk=whonpag(ii,k) ; 
        if (gex(i,kk)>0.0d0) then ; 
          a=a+gex(i,kk)*gen(kk)%wa(ii)  
        endif
      end do
      node(i)%diffe=node(i)%diffe+a*delta
      if (node(i)%diffe>1.0) node(i)%diffe=1.0 
      if (node(i)%diffe<0.0d0) node(i)%diffe=0.0

      nodeo(i)%diffe=node(i)%diffe !>>Miquel17-9-14
    end if

    if(node(i)%tipus<2)then  !this for epithelial nodes  !>>Miquel12-5-14
      j=node(i)%altre
      if(ffu(11)==1) then        !plastic deformation   !>>Miquel5-2-14
        !kplast=node(i)%kplast
        !if(fmeanl(i)<epsilod)then;ki=1/node(i)%rep;else;ki=1/node(i)%you;endif
        !if(fmeanl(j)<epsilod)then;kj=1/node(j)%rep;else;kj=1/node(j)%you;endif
        node(i)%reqp=node(i)%reqp+(node(i)%kplast*fmeanl(i))*delta
        node(j)%reqp=node(j)%reqp+(node(j)%kplast*fmeanl(j))*delta !; print*,"fmeanl",fmeanl(i),fmeanl(j)

        nodeo(i)%reqp=node(i)%reqp  !>>Miquel17-9-14
        nodeo(j)%reqp=node(j)%reqp  !>>Miquel17-9-14

      else
        node(i)%reqp=0
        node(j)%reqp=0
      end if
    
      if (npag(21)>0) then  ! Contraction by genes 
        ii=21 ; a=0 ; b=0
        do k=1,npag(ii) ; 
          kk=whonpag(ii,k) !;print*,"ij",i,j,"gex",gex(i,kk),gex(j,kk),"kk",kk
          if (gex(i,kk)>0.0d0.or.gex(j,kk)>0.0d0) then ; 
            a=a+gex(i,kk)*gen(kk)%wa(ii)
            b=b+gex(j,kk)*gen(kk)%wa(ii)
          endif 
        enddo
        node(i)%reqc=nodeo(i)%reqc+a*differ  !;print*,"a",a,"b",b
        node(j)%reqc=nodeo(j)%reqc+b*differ
        !;print*,"reqc",node(i)%reqc,node(j)%reqc
      else
        node(i)%reqc=0 ;node(j)%reqc=0
      end if

      if(ffu(17)==1) then        !volume conservation   !>>Miquel6-5-14 ! >>> Is 24-5-14
        !kvol=node(i)%kvol
        dvol=0.5*(node(i)%reqcr+node(j)%reqcr-(node(i)%req+node(j)%req))
        node(i)%reqv=node(i)%reqv+node(i)%kvol*dvol*delta
        node(j)%reqv=node(j)%reqv+node(j)%kvol*dvol*delta
        nodeo(i)%reqv=node(i)%reqv  !>>Miquel17-9-14
        nodeo(j)%reqv=node(j)%reqv  !>>Miquel17-9-14

      else
        node(i)%reqv=0 ; node(j)%reqv=0
      end if

      if (ffu(18)==1) call diffusion_of_reqcr       !diffusion of reqcr ! >>> Is 25-5-14

      a=node(i)%da-node(i)%req
      node(i)%req=node(i)%reqcr+node(i)%reqc+node(i)%reqp+node(i)%reqv  !now req is the sum of the req components: growth/apoptosis and contraction/deformation
      if(node(i)%req>df_reqmax) node(i)%req=df_reqmax !put an upper an lower boundary on how much  !>>Miquel28-7-14
      if(node(i)%req<reqmin) node(i)%req=reqmin !the req can be deformed
      node(i)%da=node(i)%req+a

      b=node(j)%da-node(j)%req
      node(j)%req=node(j)%reqcr+node(j)%reqc+node(j)%reqp+node(j)%reqv  !now req is the sum of the req components: growth/apoptosis and contraction/deformation
      if(node(j)%req>df_reqmax) node(j)%req=df_reqmax !put an upper an lower boundary on how much  !>>Miquel28-7-14
      if(node(j)%req<reqmin) node(j)%req=reqmin !the req can be deformed
      node(j)%da=node(j)%req+b
      if (npag(6)>0) then;ii=6;a=0 ; b=0 ; do k=1,npag(ii) ; kk=whonpag(ii,k) ; !wa in req-space units  
      if (gex(i,kk)>0.0d0) then ; a=a+gex(i,kk)*gen(kk)%wa(ii) ;endif;
      if (gex(j,kk)>0.0d0) then ; b=b+gex(j,kk)*gen(kk)%wa(ii) ;endif;enddo;
      node(i)%da=nodeo(i)%da+a*differ;node(j)%da=nodeo(j)%da+b*differ;
      if (node(i)%da<0) node(i)%da=0.0;if (node(j)%da<0) node(j)%da=0.0;end if
      if (npag(7)>0) then;ii=7;a=0 ; b=0 ; do k=1,npag(ii) ; kk=whonpag(ii,k) ; !wa in req-space units  
      if (gex(i,kk)>0.0d0) then ; a=a+gex(i,kk)*gen(kk)%wa(ii) ;endif;
      if (gex(j,kk)>0.0d0) then ; b=b+gex(j,kk)*gen(kk)%wa(ii) ;endif;enddo;
      node(i)%you=nodeo(i)%you+a*differ;node(j)%you=nodeo(j)%you+b*differ;
      if (node(i)%you<0) node(i)%you=0.0;if (node(j)%you<0) node(j)%you=0.0;end if
      if (npag(8)>0) then;ii=8;a=0 ; b=0 ; do k=1,npag(ii) ; kk=whonpag(ii,k) ; !wa in force units 
      if (gex(i,kk)>0.0d0) then ; a=a+gex(i,kk)*gen(kk)%wa(ii) ;endif;
      if (gex(j,kk)>0.0d0) then ; b=b+gex(j,kk)*gen(kk)%wa(ii) ;endif;enddo;
      node(i)%adh=nodeo(i)%adh+a*differ;node(j)%adh=nodeo(j)%adh+b*differ;
      if (node(i)%adh<0) node(i)%adh=0.0;if (node(j)%adh<0) node(j)%adh=0.0;end if
      if (npag(9)>0) then;ii=9;a=0 ; b=0 ; do k=1,npag(ii) ; kk=whonpag(ii,k) ; !wa in force units 
      if (gex(i,kk)>0.0d0) then ; a=a+gex(i,kk)*gen(kk)%wa(ii) ;endif;
      if (gex(j,kk)>0.0d0) then ; b=b+gex(j,kk)*gen(kk)%wa(ii) ;endif;enddo;
      node(i)%rep=nodeo(i)%rep+a*differ;node(j)%rep=nodeo(j)%rep+b*differ;
      if (node(i)%rep<0) node(i)%rep=0.0;if (node(j)%rep<0) node(j)%rep=0.0;end if
      if (npag(10)>0) then;ii=10;a=0 ; b=0 ; do k=1,npag(ii) ; kk=whonpag(ii,k) ; !wa in force units 
      if (gex(i,kk)>0.0d0) then ; a=a+gex(i,kk)*gen(kk)%wa(ii) ;endif;
      if (gex(j,kk)>0.0d0) then ; b=b+gex(j,kk)*gen(kk)%wa(ii) ;endif;enddo;
      node(i)%repcel=nodeo(i)%repcel+a*differ;node(j)%repcel=nodeo(j)%repcel+b*differ;
      if (node(i)%repcel<0) node(i)%repcel=0.0;if (node(j)%repcel<0) node(j)%repcel=0.0;end if
      if (npag(11)>0) then;ii=11;a=0 ; b=0 ; do k=1,npag(ii) ; kk=whonpag(ii,k) ; !wa in force units 
      if (gex(i,kk)>0.0d0) then ; a=a+gex(i,kk)*gen(kk)%wa(ii) ;endif;
      if (gex(j,kk)>0.0d0) then ; b=b+gex(j,kk)*gen(kk)%wa(ii) ;endif;enddo;
      node(i)%tor=nodeo(i)%tor+a*differ;node(j)%tor=nodeo(j)%tor+b*differ;
      if (node(i)%tor<0) node(i)%tor=0.0;if (node(j)%tor<0) node(j)%tor=0.0;end if
      if (npag(12)>0) then;ii=12;a=0 ; b=0 ; do k=1,npag(ii) ; kk=whonpag(ii,k) ; !wa in force units 
      if (gex(i,kk)>0.0d0) then ; a=a+gex(i,kk)*gen(kk)%wa(ii) ;endif;
      if (gex(j,kk)>0.0d0) then ; b=b+gex(j,kk)*gen(kk)%wa(ii) ;endif;enddo;
      node(i)%stor=nodeo(i)%stor+a*differ;node(j)%stor=nodeo(j)%stor+b*differ;
      if (node(i)%stor<0) node(i)%stor=0.0;if (node(j)%stor<0) node(j)%stor=0.0;end if
      if (npag(13)>0) then;ii=13;a=0 ; b=0 ; do k=1,npag(ii) ; kk=whonpag(ii,k) ; !wa in force units 
      if (gex(i,kk)>0.0d0) then ; a=a+gex(i,kk)*gen(kk)%wa(ii) ;endif;
      if (gex(j,kk)>0.0d0) then ; b=b+gex(j,kk)*gen(kk)%wa(ii) ;endif;enddo;
      node(i)%reqs=nodeo(i)%reqs+a*differ;node(j)%reqs=nodeo(j)%reqs+b*differ;
      if (node(i)%reqs<0) node(i)%reqs=0.0;if (node(j)%reqs<0) node(j)%reqs=0.0;end if
      if (npag(14)>0) then;ii=14;a=0 ; b=0 ; do k=1,npag(ii) ; kk=whonpag(ii,k) ; !wa in force units 
      if (gex(i,kk)>0.0d0) then ; a=a+gex(i,kk)*gen(kk)%wa(ii) ;endif;
      if (gex(j,kk)>0.0d0) then ; b=b+gex(j,kk)*gen(kk)%wa(ii) ;endif;enddo;
      node(i)%ke=nodeo(i)%ke+a*differ;node(j)%ke=nodeo(j)%ke+b*differ;
      if (node(i)%ke<0) node(i)%ke=0.0;if (node(j)%ke<0) node(j)%ke=0.0;end if
      if (npag(15)>0) then;ii=15;a=0 ; b=0 ; do k=1,npag(ii) ; kk=whonpag(ii,k) ; !wa in force units 
      if (gex(i,kk)>0.0d0) then ; a=a+gex(i,kk)*gen(kk)%wa(ii) ;endif;
      if (gex(j,kk)>0.0d0) then ; b=b+gex(j,kk)*gen(kk)%wa(ii) ;endif;enddo;
      node(i)%mo=nodeo(i)%mo+a*differ;node(j)%mo=nodeo(j)%mo+b*differ;
      if (node(i)%mo<0) node(i)%mo=0.0;if (node(j)%mo<0) node(j)%mo=0.0;end if
      if (npag(16)>0) then;ii=16;a=0 ; b=0 ; do k=1,npag(ii) ; kk=whonpag(ii,k) ; !wa in req-space units  
      if (gex(i,kk)>0.0d0) then ; a=a+gex(i,kk)*gen(kk)%wa(ii) ;endif;
      if (gex(j,kk)>0.0d0) then ; b=b+gex(j,kk)*gen(kk)%wa(ii) ;endif;enddo;
      node(i)%dmo=nodeo(i)%dmo+a*differ;node(j)%dmo=nodeo(j)%dmo+b*differ;
      if (node(i)%dmo<0) node(i)%dmo=0.0;if (node(j)%dmo<0) node(j)%dmo=0.0;end if
      if (npag(27)>0) then; ii=27;a=0 ; b=0 ; do k=1,npag(ii) ; kk=whonpag(ii,k) ; !wa in req-space units  
      if (gex(i,kk)>0.0d0) then ; a=a+gex(i,kk)*gen(kk)%wa(ii) ;endif;
      if (gex(j,kk)>0.0d0) then ; b=b+gex(j,kk)*gen(kk)%wa(ii) ;endif;enddo;
      node(i)%kplast=nodeo(i)%kplast+a*differ;node(j)%kplast=nodeo(j)%kplast+b*differ;
      if (node(i)%kplast<0) node(i)%kplast=0.0;if (node(j)%kplast<0) node(j)%kplast=0.0;end if
      if (npag(28)>0) then;ii=28;a=0 ; b=0 ; do k=1,npag(ii) ; kk=whonpag(ii,k) ; !wa in req-space units 
      if (gex(i,kk)>0.0d0) then ; a=a+gex(i,kk)*gen(kk)%wa(ii) ;endif;
      if (gex(j,kk)>0.0d0) then ; b=b+gex(j,kk)*gen(kk)%wa(ii) ;endif;enddo;
      node(i)%kvol=nodeo(i)%kvol+a*differ;node(j)%kvol=nodeo(j)%kvol+b*differ;
      if (node(i)%kvol<0) node(i)%kvol=0.0;if (node(j)%kvol<0) node(j)%kvol=0.0;end if


    else if (node(i)%tipus>2)then !this for mesenchyme and ECM

      if (npag(21)>0) then 
        ii=21 ; a=0
        do k=1,npag(ii) ; 
          kk=whonpag(ii,k) 
          if (gex(i,kk)>0.0d0) then ; 
            a=a+gex(i,kk)*gen(kk)%wa(ii)
          endif 
        enddo
        node(i)%reqc=nodeo(i)%reqc+a*differ*delta !wa in req-space units
      else
        node(i)%reqc=0
      end if

      a=node(i)%da-node(i)%req
      node(i)%req=node(i)%reqcr+node(i)%reqc  !now req is the sum of the req components: growth/apoptosis and contraction/deformation
      if(node(i)%req>df_reqmax) node(i)%req=df_reqmax !put an upper an lower boundary on how much  !>>Miquel28-7-14
      if(node(i)%req<reqmin) node(i)%req=reqmin !the req can be deformed
      node(i)%da=node(i)%req+a

      if (npag(6)>0) then;ii=6;a=0 ; do k=1,npag(ii) ; kk=whonpag(ii,k) ; if (gex(i,kk)>0.0d0) then ; 
      a=a+gex(i,kk)*gen(kk)%wa(ii) ;endif;enddo;node(i)%da=nodeo(i)%da+a*differ;if (node(i)%da<0) node(i)%da=0.0;end if
      if (npag(7)>0) then;ii=7;a=0 ; do k=1,npag(ii) ; kk=whonpag(ii,k) ; if (gex(i,kk)>0.0d0) then ; 
      a=a+gex(i,kk)*gen(kk)%wa(ii) ;endif;enddo;node(i)%you=nodeo(i)%you+a*differ;if (node(i)%you<0) node(i)%you=0.0;end if
      if (npag(8)>0) then;ii=8;a=0 ; do k=1,npag(ii) ; kk=whonpag(ii,k) ; if (gex(i,kk)>0.0d0) then ; 
      a=a+gex(i,kk)*gen(kk)%wa(ii) ;endif;enddo;node(i)%adh=nodeo(i)%adh+a*differ;if (node(i)%adh<0) node(i)%adh=0.0;end if
      if (npag(9)>0) then;ii=9;a=0 ; do k=1,npag(ii) ; kk=whonpag(ii,k) ; if (gex(i,kk)>0.0d0) then ; 
      a=a+gex(i,kk)*gen(kk)%wa(ii) ;endif;enddo;node(i)%rep=nodeo(i)%rep+a*differ;if (node(i)%rep<0) node(i)%rep=0.0;end if
      if (npag(10)>0) then;ii=10;a=0 ; do k=1,npag(ii) ; kk=whonpag(ii,k) ; if (gex(i,kk)>0.0d0) then ; 
      a=a+gex(i,kk)*gen(kk)%wa(ii) ;endif;enddo;node(i)%repcel=nodeo(i)%repcel+a;if (node(i)%repcel<0) node(i)%repcel=0.0;end if
      if (npag(11)>0) then;ii=11;a=0 ; do k=1,npag(ii) ; kk=whonpag(ii,k) ; if (gex(i,kk)>0.0d0) then ; 
      a=a+gex(i,kk)*gen(kk)%wa(ii) ;endif;enddo;node(i)%tor=nodeo(i)%tor+a*differ;if (node(i)%tor<0) node(i)%tor=0.0;end if
      if (npag(12)>0) then;ii=12;a=0 ; do k=1,npag(ii) ; kk=whonpag(ii,k) ; if (gex(i,kk)>0.0d0) then ; 
      a=a+gex(i,kk)*gen(kk)%wa(ii) ;endif;enddo;node(i)%stor=nodeo(i)%stor+a*differ;if (node(i)%stor<0) node(i)%stor=0.0;end if
      if (npag(13)>0) then;ii=13;a=0 ; do k=1,npag(ii) ; kk=whonpag(ii,k) ; if (gex(i,kk)>0.0d0) then ; 
      a=a+gex(i,kk)*gen(kk)%wa(ii) ;endif;enddo;node(i)%reqs=nodeo(i)%reqs+a*differ;if (node(i)%reqs<0) node(i)%reqs=0.0;end if
      if (npag(14)>0) then;ii=14;a=0 ; do k=1,npag(ii) ; kk=whonpag(ii,k) ; if (gex(i,kk)>0.0d0) then ; 
      a=a+gex(i,kk)*gen(kk)%wa(ii) ;endif;enddo;node(i)%ke=nodeo(i)%ke+a*differ;if (node(i)%ke<0) node(i)%ke=0.0;end if
      if (npag(15)>0) then;ii=15;a=0 ; do k=1,npag(ii) ; kk=whonpag(ii,k) ; if (gex(i,kk)>0.0d0) then ; 
      a=a+gex(i,kk)*gen(kk)%wa(ii) ;endif;enddo;node(i)%mo=nodeo(i)%mo+a*differ;if (node(i)%mo<0) node(i)%mo=0.0;end if
      if (npag(16)>0) then;ii=16;a=0 ; do k=1,npag(ii) ; kk=whonpag(ii,k) ; if (gex(i,kk)>0.0d0) then ; 
      a=a+gex(i,kk)*gen(kk)%wa(ii);endif;enddo;node(i)%dmo=nodeo(i)%dmo+a*differ;if (node(i)%dmo<0) node(i)%dmo=0.0;end if


    end if
  end do

  ! if a node gets 3 times the original size of node 1 we fucking kill the program
  if (mod(getot,100).eq.0) then !>>> Is 4-2-14
    if (ffu(7)==1) then  ! if a node gets twice its original size we fucking kill the program ! >>> Is 4-2-14
      do i=1,nd                                                                               ! >>> Is 4-2-14 
        if (node(i)%da>ramax) then                        ! >>> Is 4-2-14
          node(i)%da=ramax
        end if                                                                                ! >>> Is 4-2-14
        !if (node(i)%req>ramax) then   !maximal req now is controlled above in the same sub. !>>Miquel7-8-14
        !  node(i)%req=ramax
        !end if
      end do                                                                                  ! >>> Is 4-2-14
    end if                                                                                    ! >>> Is 4-2-14
  end if                                                                                      ! >>> Is 4-2-14

! We check in here if all cells are already differentiated, if they are we stop the simulations

!>>> Is 4-4-14
if (ffu(10)==1) then
  do i=1,nd
    if (node(i)%tipus==4) cycle !pfh, ignore ecm
    if (node(i)%diffe<1.0d0) return
  end do
  print *,""
  print *," THIS IS THE END all cells are differentiated and then the simulations stop",trim(carg)//trim(noff),trim(nofi)
  print *,""
  call writesnap                                                                            !>>> Is 25-2-14
  if (len_trim(carg)/=0) then
    print *,trim(carg)//trim(noff),"ki"
    print *,trim(carg)//trim(nofi),"kii"
    open(23,file=trim(carg)//"t",iostat=i)
    print *,"making...",trim(carg)//"t"
    write(23,*,ERR=46) trim(carg)//trim(nofi)
    flush(23)
    !ret = fsync(fnum(23))
          
    ! Handle possible error
    !if (ret /= 0) stop "Error calling FSYNC"
    open(23,file=trim(carg)//"t",iostat=ii)
    print *,trim(carg)//"t",ii,"iostat"
    read(23,*,END=45,ERR=48) cx
print *,""
print *,cx,"cx HERE"
print *,""
    !close(23)
    call flush(23)
    print *,"done with iostat",ii
    cx="ls -alrt "//trim(carg)//"t" 
    call system(cx) 
  end if
  stop
45 print *,"end of file error"
  stop
46 print *,"error in writing",trim(carg)//"t"
  stop 
48 print *,"other end"
  stop
end if
!>>> Is 4-4-14

end subroutine

!*******************************************************************

subroutine cellbreak   !>>> 17-1-14
  real*8  ax,ay,az,ida
  integer doapo(nd)
  integer nocon(nd,nd)
  integer i,j,k,ii,jj,kk,iii,jjj,kkk,ik

  ! now we check that the cell is not split in two or more parts
  doapo=0
  nocon=0

  do i=1,ncels

    kkk=0
    do j=1,cels(i)%nunodes
      ii=cels(i)%node(j)
      if (node(ii)%marge==0) then ; kkk=ii ; ida=node(kkk)%da ; exit ; end if
    end do

    if (kkk==0) then !it means that the cell has no nucleus and then IT MUST DIE!!!!   !>>> Is 11-6-14
      do j=1,cels(i)%nunodes  !>>> Is 11-6-14
        ii=cels(i)%node(j)    !>>> Is 29-6-14
!print *,ii,i,j,cels(i)%node(:cels(i)%nunodes),ii,"ii",nd
        doapo(ii)=1           !>>> Is 11-6-14
      end do                  !>>> Is 11-6-14
      kkk=1  ! >>> Is 11-6-14
    else
      do j=1,cels(i)%nunodes
        ii=cels(i)%node(j)
        ax=node(ii)%x ; ay=node(ii)%y ; az=node(ii)%z
        if (sqrt((ax-node(kkk)%x)**2+(ay-node(kkk)%y)**2+(az-node(ii)%z)**2)<node(ii)%da+ida) then
          doapo(ii)=1
        end if
      end do
    end if   ! >>> Is 11-6-14

    doapo(kkk)=1
    do j=1,cels(i)%nunodes
      ii=cels(i)%node(j)
      ax=node(ii)%x ; ay=node(ii)%y ; az=node(ii)%z
      ida=node(ii)%da
      do jj=1,cels(i)%nunodes
        if (j==jj) cycle
        iii=cels(i)%node(jj)
        if (sqrt((ax-node(iii)%x)**2+(ay-node(iii)%y)**2+(az-node(iii)%z)**2)<node(iii)%da+ida) then
          nocon(ii,iii)=1
          nocon(iii,ii)=1
          if (doapo(ii)==1) then
            doapo(iii)=1
          else
            if (doapo(iii)==1) doapo(ii)=1
          end if
        end if
      end do
    end do

    do k=1,cels(i)%nunodes/2
      do j=1,cels(i)%nunodes
        ii=cels(i)%node(j)
        do jj=1,cels(i)%nunodes
          if (j==jj) cycle
          iii=cels(i)%node(jj)
          if (nocon(ii,iii)==1) then
            if (doapo(ii)==1) then
              doapo(iii)=1
            else
              if (doapo(iii)==1) doapo(ii)=1
            end if
          end if
        end do
      end do
    end do
  end do

  ik=1
  do while(ik<=nd)
    if (doapo(ik)==0) then
      if (node(ik)%tipus>2) then  !we only delete an epithelial node if its altre is also lonely
        call apoptosis(ik)
      else
        if (doapo(node(ik)%altre)==0) then
          call apoptosis(ik)
        end if
      end if
      ik=ik+1
    else
      ik=ik+1   ! I know, it is kind of funky but it should be this way, a loop with nd wont do because nd decreases because of apoptosis
    end if
  end do


end subroutine

!********************************************************************

subroutine polarization
integer:: celd,nnod,tipi,ggr,ccen
real*8::a,b,c,d,e,ax,ay,az,bx,by,bz,cx,cy,cz,ix,iy,iz,alfa,s
      do celd=1,ncels
        tipi=cels(celd)%ctipus
        nnod=cels(celd)%nunodes        
        if (nnod==0) cycle      ! >>> Is 10-5-14
        iy=1d10 ; cx=0d0 ; cy=0d0 ; cz=0d0
	a=cels(celd)%cex ; b=cels(celd)%cey ; c=cels(celd)%cez   

        do i=1,nnod                                                     ! (gen) in the centroid (in the closest node)
          j=cels(celd)%node(i)
          if(node(j)%tipus==1.or.node(j)%tipus==3)then !in epithelial cells, polarity is planar, so we only take one layer of nodes >>>Miquel22-10-13
            d=sqrt((node(j)%x-a)**2+(node(j)%y-b)**2+(node(j)%z-c)**2)
          end if
          if(d.le.iy)then;iy=d;ccen=j;endif             
        end do   

        alfa=0.0d0                                                 ! concentration in the central node
        do k=1,npag(nparam_per_node+8)    
          kk=whonpag(nparam_per_node+8,k)
          if (gex(ccen,kk)>0.0d0) then
            alfa=alfa+gex(ccen,kk)*gen(kk)%wa(nparam_per_node+8)   ! wa in units of probability such that it makes things to go from 0 to 1
          end if
        end do  

        ix=0d0 ; iy=0d0 ; iz=0d0                                        ! vector of the gradient within a cell
        do i=1,nnod                                                     
            j=cels(celd)%node(i)
            if(node(j)%tipus==1.or.node(j)%tipus==3)then          
              d=sqrt((node(j)%x-a)**2+(node(j)%y-b)**2+(node(j)%z-c)**2)
              if (d<epsilod) cycle
              d=1d0/d                                                   ! module of radial vectors to get unitary vectors     
              s=0.0d0
              do k=1,npag(nparam_per_node+8)
                kk=whonpag(nparam_per_node+8,k)
                if (gex(j,kk)>0.0d0) then
                  s=s+gex(j,kk)*gen(kk)%wa(nparam_per_node+8)
                end if
              end do
              ix=ix+((node(j)%x-a)*d)*(s-alfa)                   ! and ignore shape/size effects
              iy=iy+((node(j)%y-b)*d)*(s-alfa)
              iz=iz+((node(j)%z-c)*d)*(s-alfa)
            end if
        end do

        if((ix.eq.0).and.(iy.eq.0).and.(iz.eq.0))then            ! if the gene has uniform expresion, the vector is random ! >>>Miguel1-7-14
          call random_number(a)                                  ! >>>Miguel1-7-14
          k=int(a*nvaloq)+1                                      ! >>>Miguel1-7-14
          cels(celd)%polx=particions_esfera(k,1)                 ! >>>Miguel1-7-14
          cels(celd)%poly=particions_esfera(k,2)                 ! >>>Miguel1-7-14
          cels(celd)%polz=particions_esfera(k,3)                 ! >>>Miguel1-7-14
        else                                                     ! >>>Miguel1-7-14
          a=ix**2+iy**2+iz**2 
          if(a==0)then
            cels(celd)%polx=0d0 ; cels(celd)%poly=0d0 ; cels(celd)%polz=0d0	! unitary resultant vector (gradient polarization)
          else
            d=1d0/sqrt(a)
            cels(celd)%polx=ix*d ; cels(celd)%poly=iy*d ; cels(celd)%polz=iz*d	! unitary resultant vector (gradient polarization)
          end if
          if((ix.eq.0d0).and.(iy.eq.0d0).and.(iz.eq.0d0))then                     ! miguel27-11-13
            cels(celd)%polx=0d0 ; cels(celd)%poly=0d0 ; cels(celd)%polz=0d0
          endif   ! miguel27-11-13
        endif                                                    ! >>>Miguel1-7-14
      end do
end subroutine

!*******************************************************************

subroutine polarizationisaac
integer i,j,k,ii,kk
real*8 a,b,c,s,sx,sy,sz,d,aa,bb,cc 

do i=1,ncels
  a=cels(i)%cex ; b=cels(i)%cey ; c=cels(i)%cez
  sx=0.0d0      ; sy=0.0d0      ; sz=0.0d0
  if (node(cels(i)%node(1))%tipus<3) then
    do j=1,cels(i)%nunodes
      ii=cels(i)%node(j)
      if (node(ii)%tipus==1) then
        s=0.0d0
        do k=1,npag(nparam_per_node+8)
          kk=whonpag(nparam_per_node+8,k)
          if (gex(ii,kk)>0.0d0) then
            s=s+gex(ii,kk)*gen(kk)%wa(nparam_per_node+8)*delta
          end if
        end do
        if (s/=0.0d0) then
          aa=node(ii)%x-a ; bb=node(ii)%y-b ; cc=node(ii)%z-c
          d=s/sqrt(aa**2+bb**2+cc**2)
          sx=sx+d*aa ; sy=sy+d*bb ; sz=sz+d*cc
        end if
      end if
    end do
    aa=sqrt(sx**2+sy**2+sz**2)
    if (aa>0.0d0) then
      d=1d0/aa
      cels(i)%polx=sx*d ; cels(i)%poly=sy*d ; cels(i)%polz=sz*d  
    else
      cels(i)%polx=0.0d0 ; cels(i)%poly=0.0d0 ; cels(i)%polz=0.0d0   
    end if
  else
    if (node(cels(i)%node(1))%tipus==3) then
      do j=1,cels(i)%nunodes
        ii=cels(i)%node(j)
        s=0.0d0
        do k=1,npag(nparam_per_node+8)
          kk=whonpag(nparam_per_node+8,k)
          if (gex(ii,kk)>0.0d0) then
            s=s+gex(ii,kk)*gen(kk)%wa(nparam_per_node+8)
          end if
        end do
        if (s/=0.0d0) then
          aa=node(ii)%x-a ; bb=node(ii)%y-b ; cc=node(ii)%z-c
          sx=sx+s*aa ; sy=sy+s*bb ; sz=sz+s*cc
        end if
      end do
      aa=sqrt(sx**2+sy**2+sz**2)
      if (aa>0.0d0) then  !epsilod) then
        d=1d0/aa
        cels(i)%polx=sx*d ; cels(i)%poly=sy*d ; cels(i)%polz=sz*d  
      else
        cels(i)%polx=0.0d0 ; cels(i)%poly=0.0d0 ; cels(i)%polz=0.0d0   
      end if
    end if
  end if
end do
end subroutine

!*************************************************************************************

subroutine change_minsize_for_div  ! updates the size required for dividint according to gene expression nparam_per_node+10
integer ick,j,k,ii,kk
real*8 a,b,c,s,sx,sy,sz,d 

do ick=1,ncels
  s=0.0d0
  do j=1,cels(ick)%nunodes
    ii=cels(ick)%node(j)
    do k=1,npag(nparam_per_node+10)
      kk=whonpag(nparam_per_node+10,k)
      if (gex(ii,kk)>0.0d0) then
        s=s+gex(ii,kk)*gen(kk)%wa(nparam_per_node+10)  !wa in units of number of nodes but it can be roughly understood as space-req units
      end if                                           ! THIS WA CAN BE NEGATIVE
    end do
  end do
  s=s*delta
!print *,delta,"ACHTUNG"
!  if (s>0.0d0) then
    s=s/cels(ick)%nunodes !this way %minsize is independent of cell size !>>>>Miquel2-12-13
    cels(ick)%minsize_for_div=cels(ick)%minsize_for_div+s  !this can make a SUDDEN CHANGE
    if (cels(ick)%minsize_for_div<1) cels(ick)%minsize_for_div=1
!  end if
end do

end subroutine

!*************************************************************************************

!>>> Is 5-2-14
subroutine change_maxsize_for_div  ! updates the size required for dividint according to gene expression nparam_per_node+10
integer ick,j,k,ii,kk
real*8 a,b,c,s,sx,sy,sz,d 

do ick=1,ncels
  s=0.0d0
  do j=1,cels(ick)%nunodes
    ii=cels(ick)%node(j)
    do k=1,npag(nparam_per_node+15)
      kk=whonpag(nparam_per_node+15,k)
      if (gex(ii,kk)>0.0d0) then
        s=s+gex(ii,kk)*gen(kk)%wa(nparam_per_node+15)  ! THIS WA MAY BE NEGATIVE
      end if
    end do
  end do
  s=s*delta
!  if (s>0.0d0) then
    s=s/cels(ick)%nunodes !this way %minsize is independent of cell size !>>>>Miquel2-12-13
    cels(ick)%maxsize_for_div=cels(ick)%maxsize_for_div+s  !wa in units of space-req roughly this can make a SUDDEN CHANGE
    if (cels(ick)%maxsize_for_div<1) cels(ick)%maxsize_for_div=1
!  end if
end do

!>>> Is 5-2-14

end subroutine


!**************************************************************************************

subroutine emt !epithelial-mensenchymal transitions

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
      a=gex(ii,kk) !don't change this, weird stuff will happen if you do !>>Miquel18-12-14
      s=s+a*gen(kk)%wa(nparam_per_node+13)
   !if(gex(ii,kk)>epsilod) print*,"emt gene",kk,"node",ii,"cell",ick,"conc",gex(ii,kk)
    end do
  end do
  cels(ick)%temt=cels(ick)%temt+s*delta/real(cels(ick)%nunodes)
  
  if (cels(ick)%temt>=1.0d0) then       ! wa in units of probability: this is an arbitrary value but it is fine since the rest can be re-scaled 
    !so this cell goes to emt
 !print*,"EMT HAPPENS CELL",ick
    cels(ick)%ctipus=3  ! >>> Is 4-1-14
    do jjj=1,cels(ick)%nunodes
      iii=cels(ick)%node(jjj)          
      if (node(iii)%tipus==1) then
        iv=node(iii)%altre
        a=node(iv)%x-node(iii)%x ; b=node(iv)%y-node(iii)%y ; c=node(iv)%z-node(iii)%z 
        d=1d0/sqrt(a**2+b**2+c**2)
        a=a*d ; b=b*d ; c=c*d
        aa=0.5d0*(node(iv)%x+node(iii)%x) ; bb=0.5d0*(node(iv)%y+node(iii)%y) ; cc=0.5d0*(node(iv)%z+node(iii)%z) 
        dd=(node(iii)%da+node(iv)%da)*0.1d0
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
        !if only one node per cell both nodes of the epithelial cell get a nucleus
        if (ffu(1)==1) then ; node(iii)%marge=0 ; node(iv)%marge=0 ; end if ! >>> Is 10-10-14
      end if
    end do
  end if
end do e

end subroutine

!*******************************************************************************************************

subroutine diffusion_of_reqcr
  real*8 hreqcr(nd),hreqp(nd),hreqc(nd)
  integer i,j,k,ii,jj,kk
  real*8 a,b,c,d

  do i=1,nd
    a=0.0d0 ; b=0.0 ; c=0.0d0
    do ii=1,nneigh(i)
      k=neigh(i,ii)
      if (node(i)%icel/=node(k)%icel) cycle    ! only within the same cell
      if (node(i)%tipus/=node(k)%tipus) cycle  ! only within the same side of the cell
      if (node(i)%tipus>2) cycle               ! only for epithelial cells
      d=dneigh(i,ii)
      a=a+(node(k)%reqcr-node(i)%reqcr) !/(d+1d0)
      !b=b+(node(k)%reqc-node(i)%reqc) !/(d+1d0)
      !c=c+(node(k)%reqp-node(i)%reqp) !/(d+1d0)
    end do
    hreqcr(i)=a*dif_req ! >>> 11-6-14
    !hreqc(i)=b*dif_req  ! >>> 11-6-14
    !hreqp(i)=c*dif_req  ! >>> 11-6-14
  end do  
  do i=1,nd
    node(i)%reqcr=node(i)%reqcr+delta*hreqcr(i)
    !node(i)%reqc=node(i)%reqc+delta*hreqc(i)
    !node(i)%reqp=node(i)%reqp+delta*hreqp(i)
  end do
end subroutine

end module nexus
