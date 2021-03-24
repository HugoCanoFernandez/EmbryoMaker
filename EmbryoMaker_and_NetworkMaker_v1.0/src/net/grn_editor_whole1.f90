!    NetworkMaker software (General Node Model)
!    Computational model to simulate morphogenetic processes in living organs and tissues.
!    Copyright (C) 2014 Roland Zimm & Isaac Salazar-Ciudad

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


module some_widgets

  use iso_c_binding, only: c_ptr

type(c_ptr) :: notebook
integer :: w1, w2, w3, w4, w0, aaa0, bbb0, ccc0, ddd0, showup
integer :: choiceboxthere

end module some_widgets

module global_widgets

  use iso_c_binding, only: c_ptr, c_char, c_int
  use some_widgets

  type(c_ptr) :: table, button1, button2, button3, button4, box1, box2, box3, boxL, table2
  type(c_ptr) :: button01,button02,button03,button04,button05,button06,button07,button08,button09,&
  & button10,button11,button12,button13
  type(c_ptr) :: label1, label2, label3, label4, label01, label02, label03
  type(c_ptr) :: label00, label001, label002, label003, label004, label005, label006
  type(c_ptr) :: label04, label05, label06, label07
  type(c_ptr) :: expander,expander2,my_choice_window2,my_choice_windown,my_choice_windows,my_choice_windowp
  type(c_ptr) :: notebookLabel1, notebookLabel2, notebookLabel3, notebookLabel4, notebookLabel5, notebookLabel6
  type(c_ptr) :: linkButton, eventabbox1,checkf,entrysi, entrypi, add_window, my_choice_window
  type(c_ptr) :: labelt(29),spinbutt(29)
  type(c_ptr),allocatable :: entryna(:),labelna(:),ffut(:)

  type(c_ptr) :: my_pixbuf, spinButton1, spinButton2, spinButton3
  type(c_ptr) :: spinButton01, spinButton02, spinButton03
  type(c_ptr) :: spinButton04, spinButton05, spinButton06, spinButtona, spinButtona2, spinbutt9, spinbutt8
  type(c_ptr) :: textView, buffer, scrolled_window, statusBar, combo1
  character(kind=c_char), dimension(:), pointer :: pixel
  integer(kind=c_int) :: nch, rowstride, pixwidth_e, pixheight_e
  logical :: computing = .false.
  character(LEN=80) :: string
  integer :: cual, ffuval, saveposn
  character(LEN=380) :: str, str2, str4
  character*12,allocatable :: genename(:)
  integer,allocatable :: savepos(:,:) 
  integer :: inii=0, yetthere=0, yetthere1=0, yetthere2=0, yetthere3=0, yetthere5=0, yetthere6=0
  integer :: style1, deleteif=0, rompetout=0, namefieldopen=0, modify=1
  real*8 :: quot=0.0d0
  
end module global_widgets

module read0

use general
use genetic
use io
use global_widgets !R!

contains
subroutine readsnap0 ! this starts the network editor with default conditions

    nvarglobal_out=5
    winame=".NOTHING-NOTHING-NOTHING-NOTHING-NOTHING."
    if (allocated(ffu)) deallocate(ffu)
    allocate(ffu(nfu))
    if (allocated(param)) deallocate(param)
    allocate(param(nparam))
    if (allocated(varglobal_out)) deallocate(varglobal_out)
    allocate(varglobal_out(nvarglobal_out))
    if(allocated(genename)) deallocate(genename)
    allocate(genename(1))
    genename="1"
    
    ffu=0; ffu(3)=1; ffu(5)=1
 
    param=0.0000000000000000E+00
    param(2)=1.0000000000000001E-05
    param(3)=1.0000000000000000E+03
    param(4)=9.9999997764825821E-03
    param(5)=2.5000000000000000E-01
    param(6)=1.0000000000000000E+00
    param(7)=6.9999998807907104E-01
    param(8)=-1.0000000000000001E-01
    param(9)=0.0000000000000000E+00
    param(10)=1.0000000000000000E+00
    param(11)=1.0000000000000000E-03
    param(12)=1.0000000000000000E+05
    param(14)=1.0000000000000000E+01
    param(15)=5.0000000000000003E-0
    param(16)=2.5000000000000000E-01
    param(17)=1.0000000000000000E-02
    param(18)=1.0000000000000000E+00
    param(19)=1.0000000000000000E+00
    param(20)=1.0000000000000000E+00
    param(22)=5.0000000000000000E+02
    param(23)=1.0499999821186066E+00
    param(24)=1.0000000149011612E-01
    param(25)=1.0000000000000001E-03
    param(26)=1.0000000000000001E+00
    param(27)=1.0000000000000001E+00
    param(28)=5.0000000000000000E-01
    param(30)=2.2204460492503131E-15

    idumoriginal(1)= -11111
    idum= -11111

    call random_seed(put=idum)
    call get_param_from_matrix_read(param)

    varglobal_out=0d0
    nd=1; ng=1; ncels=1

      call initiate_gene
      if (allocated(gex)) deallocate(gex)
      allocate(gex(nda,ng))    
      if (allocated(gen)) deallocate(gen)
      allocate(gen(ng))
      if(ntipusadh>0)then                    !>>>>>Miquel14-11-13
        if (allocated(kadh)) deallocate(kadh)!
        allocate(kadh(ntipusadh,ntipusadh))  !
      end if                                 !
      do i=1,ng
        if (allocated(gen(i)%w)) deallocate(gen(i)%w)
        allocate(gen(i)%w(ng))    
        if (allocated(gen(i)%ww)) deallocate(gen(i)%ww)
        allocate(gen(i)%ww(ng*ng,3))    
        if (allocated(gen(i)%wa)) deallocate(gen(i)%wa)
        allocate(gen(i)%wa(nga))    
      end do
      gen(1)%w=0.0d0
      gen(1)%nww=0
      gen(1)%wa=0.0d0    
      gen(1)%diffu=0.0d0
      gen(1)%mu=0.0d0
      gen(1)%kindof=1
      gen(1)%npre=0
      gen(1)%npost=0
      gex(1,1)=0.1d0  

      call update_npag

      if (allocated(node)) deallocate(node)
      if (allocated(nodeo)) deallocate(nodeo)
      if (allocated(cels)) deallocate(cels)

      ncals=ncels+10
      allocate(node(nd))
      allocate(nodeo(nd))
      allocate(cels(ncels))

    !call iniarrays

     node(:)%x=0d0 ; node(:)%y=0d0 ; node(:)%z=0d0
     node(:)%e=0d0
     node(:)%req=0.25d0
     node(:)%da=0.35d0
     node(:)%you=10d0
     node(:)%adh=0d0
     node(:)%rep=10d0
     node(:)%repcel=10d0
     node(:)%tor=1d0
     node(:)%stor=0d0
     node(:)%reqs=0d0
     node(:)%ke=0d0
     node(:)%mo=0.00001d0
     node(:)%dmo=0.01d0
     node(:)%orix=0d0 ; node(:)%oriy=0d0 ; node(:)%oriz=0d0
     node(:)%acecm=0d0
     node(:)%reqc=0d0
     node(:)%reqcr=0.25d0
     node(:)%reqp=0d0
     node(:)%reqv=0d0
     node(:)%diffe=0d0
     node(:)%khold=0d0
     node(:)%kplast=0d0
     node(:)%kvol=0d0
     node(:)%tipus=3
     node(:)%icel=1
     node(:)%altre=0
     node(:)%marge=1
     node(:)%talone=0
     node(:)%hold=0
     node(:)%border=0

     cels(:)%minsize_for_div=8
     cels(:)%maxsize_for_div=16
     cels(:)%cex=0d0 ; cels(:)%cey=0d0 ; cels(:)%cez=0d0 ;
     cels(:)%polx=0d0 ; cels(:)%poly=0d0 ; cels(:)%polz=0d0 ;
     cels(:)%hpolx=0d0 ; cels(:)%hpoly=0d0 ; cels(:)%hpolz=0d0 ;
     cels(:)%fase=0d0
     cels(:)%temt=0d0
     cels(:)%nunodes=1
     cels(:)%nodela=17
     cels(:)%ctipus=3

if (allocated(cels(1)%node)) deallocate(cels(1)%node)
allocate(cels(1)%node(1))

     cels(1)%node(1)=1

    nodeo=node !>>Miquel17-9-14

  end subroutine readsnap0

end module read0

module basic_handlers

use iso_c_binding, only: c_int
  
use gtk, only: gtk_container_add, gtk_drawing_area_new, gtk_events_pending, gtk&
  &_main, gtk_main_iteration, gtk_main_iteration_do, gtk_widget_get_window, gtk_w&
  &idget_show, gtk_window_new, gtk_window_set_default, gtk_window_set_default_siz&
  &e, gtk_window_set_title, gtk_widget_queue_draw, gtk_window_set_title, gtk_label_new,&
  & gtk_spin_button_new, gtk_table_new, gtk_button_new, gtk_adjustment_new, gtk_entry_get_text, &
  &gtk_entry_set_text, gtk_table_attach_defaults, gtk_widget_show_all, gtk_button_new_with_mnemonic, &
  &gtk_spin_button_get_value, gtk_widget_destroy, gtk_widget_destroyed, gtk_entry_new, &
  &TRUE, FALSE, c_null_char, GTK_window_TOPLEVEL, gtk_init, g_signal_connect, &
  &gtk_widget_get_toplevel, gtk_widget_is_toplevel, gtk_notebook_get_n_pages,&
  &gtk_notebook_get_current_page, gtk_notebook_prev_page, gtk_notebook_next_page,&
  &CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL, gtk_disable_setlocale
     
use cairo, only: cairo_arc, cairo_create, cairo_curve_to, cairo_destroy, cairo_&
  &get_target, cairo_line_to, cairo_move_to, cairo_new_sub_path, cairo_select_fon&
  &t_face, cairo_set_font_size, cairo_set_line_width, cairo_set_source, cairo_set&
  &_source_rgb, cairo_show_text, cairo_stroke, cairo_surface_write_to_png, cairo_rotate
  
use gdk, only: gdk_cairo_create
use general
use io
use genetic
use read0

use iso_c_binding, only: c_ptr, c_int, c_char
type(c_ptr) :: my_drawing_area, my_drawing_area2, my_window, my_choice_window3

integer(c_int) :: run_status = TRUE
integer :: whichsave, which2, whatsave, what2, whichii, wcl ! just save a value we need later
integer :: test0, counts, back, cmax, nmax, undone1t, togglen, faster, ngc, ng_o
integer :: saveindex2(9)
character*55 :: cellparams(16)
character*120 :: name_dat,name_png

integer, allocatable :: saveindex(:,:), mi(:,:) ! 2.1 is number of matrix called, 2.2 the nth element
real*8, allocatable :: memor1(:,:) ! 2 is number of global params, :.1 the old, :.2 the new value
integer, allocatable :: memor2(:), m0(:), ngback(:), ngback2(:) ! tells simply which gene got copied
integer, allocatable :: undone1(:), undone2(:) ! references which is the next step to go back
real*8, allocatable :: memor3(:,:) ! gene params. 2.1 is the gene affected
real*8, allocatable :: memor4(:,:) ! pre/post matrix changes
real*8, allocatable :: memor5(:,:) ! w matrix changes
real*8, allocatable :: memor6(:,:), mr(:,:) ! ww matrix changes
real*8, allocatable :: memor7(:,:) ! wa matrix changes
real*8, allocatable :: npos(:,:),kpos(:,:)
  integer :: recf
  character*140 :: record(9)
  character(len=8) :: fmt
  character*1 :: eddp

  integer :: rb_arrange_tog

contains

subroutine del(delgen)

use genetic
use general

  type(genes) :: cgen(ng-1)
  integer, intent(in) :: delgen
  integer pres(ng),posts(ng)
  real*8 cgex(nda,ng-1)
  real*8 cw(ng-1),ckadh(ntipusadh,ntipusadh),cww(ng-1,3),ccww(ng*ng,3)
  integer moldadh

  pres=0 ; posts=0

  !for the contents of WW
  do i=1,ng
    k=gen(i)%nww
    
    ccww(:ng,:)=gen(i)%ww(:ng,:)

    kk=0
    do j=1,k
      if (gen(i)%ww(j,1)==delgen.or.gen(i)%ww(j,2)==delgen) then
        kk=kk+1
        do jj=j,k-1
          ccww(jj,:)=ccww(jj+1,:)
        end do
        ccww(k,:)=0.0d0
      else
        if (gen(i)%ww(j,1)>delgen) then
          ccww(j,1)=ccww(j,1)-1
        end if
        if (gen(i)%ww(j,2)>delgen) then
          ccww(j,2)=ccww(j,2)-1
        end if
      end if
    end do
    gen(i)%nww=gen(i)%nww-kk
    gen(i)%ww(:ng,:)=ccww(:ng,:)
  end do

  ! if delgen is an adhesion molecule we have to re-size kadh
  ! now if the form is an adhesion molecule we have to re-size kadh and renumber wa(1)
  if (gen(delgen)%wa(1)/=0) then
    moldadh=gen(delgen)%wa(1)
    do j=1,ng
      if (gen(j)%wa(1)>moldadh) gen(j)%wa(1)=gen(j)%wa(1)-1 
    end do    
    ntipusadh=ntipusadh-1
    if (ntipusadh>0) then
      ckadh(:ntipusadh+1,:ntipusadh+1)=kadh(:ntipusadh+1,:ntipusadh+1)
      deallocate(kadh)
      allocate(kadh(ntipusadh,ntipusadh))
      kadh=0.0d0
      kadh(:moldadh-1,:ntipusadh)=ckadh(:moldadh-1,:ntipusadh)
      kadh(moldadh:ntipusadh,:ntipusadh)=ckadh(moldadh+1:ntipusadh+1,:ntipusadh)
      kadh(:ntipusadh,:moldadh-1)=ckadh(:ntipusadh,:moldadh-1)
      kadh(:ntipusadh,moldadh:ntipusadh)=ckadh(:ntipusadh,moldadh+1:ntipusadh+1)
    else
      deallocate(kadh)
    end if
  end if

  ! we take away from pre and post the references to the deleted gene
  do i=1,ng-1
    do j=1,gen(i)%npre
      if (gen(i)%pre(j)==delgen) then
        pres(:j-1)=gen(i)%pre(:j-1)
        pres(j:gen(i)%npre-1)=gen(i)%pre(j+1:gen(i)%npre)
        deallocate(gen(i)%pre)
        gen(i)%npre=gen(i)%npre-1
        allocate(gen(i)%pre(gen(i)%npre))
        gen(i)%pre(:gen(i)%npre)=pres(:gen(i)%npre)
        exit
      end if
    end do
    do j=1,gen(i)%npost
      if (gen(i)%post(j)==delgen) then
        posts(:j-1)=gen(i)%post(:j-1)
        posts(j:gen(i)%npost-1)=gen(i)%post(j+1:gen(i)%npost)
        deallocate(gen(i)%post)
        gen(i)%npost=gen(i)%npost-1
        allocate(gen(i)%post(gen(i)%npost))
        gen(i)%post(:gen(i)%npost)=posts(:gen(i)%npost)
        exit
      end if
    end do
  end do

  ! the indices in pre and post need to be -1 for the genes after delgen
  do i=1,ng
    do j=1,gen(i)%npre
      if (gen(i)%pre(j)>=delgen) gen(i)%pre(j)=gen(i)%pre(j)-1
    end do
    do j=1,gen(i)%npost
      if (gen(i)%post(j)>=delgen) gen(i)%post(j)=gen(i)%post(j)-1
    end do
  end do

  cgex(:,:delgen-1)=gex(:,:delgen-1)
  cgex(:,delgen:ng-1)=gex(:,delgen+1:ng)

  deallocate(gex)
  allocate(gex(nda,ng-1))

  gex=cgex

  do i=1,delgen-1
    cgen(i)=gen(i)

    deallocate(cgen(i)%w)
    allocate(cgen(i)%w(ng-1))
    cgen(i)%w(:delgen-1)=gen(i)%w(:delgen-1)
    if (delgen<ng) then
      cgen(i)%w(delgen:ng-1)=gen(i)%w(delgen+1:ng)
    end if

    deallocate(cgen(i)%ww)
    allocate(cgen(i)%ww(ng-1,3))
    cgen(i)%ww(:ng-1,:)=gen(i)%ww(:ng-1,:)

    deallocate(cgen(i)%wa)
    allocate(cgen(i)%wa(nga))
    cgen(i)%wa=gen(i)%wa

    if (gen(i)%npre>0) then
      if (allocated(cgen(i)%pre)) deallocate(cgen(i)%pre)
      allocate(cgen(i)%pre(gen(i)%npre))
      cgen(i)%pre=gen(i)%pre
    end if
    if (gen(i)%npost>0) then
      if (allocated(cgen(i)%post)) deallocate(cgen(i)%post)
      allocate(cgen(i)%post(gen(i)%npost))
      cgen(i)%post=gen(i)%post
    end if
  end do

  do i=delgen,ng-1

    cgen(i)=gen(i+1)
    deallocate(cgen(i)%w)
    allocate(cgen(i)%w(ng-1))
    cgen(i)%w(:delgen-1)=gen(i+1)%w(:delgen-1)
    if (delgen<ng) then
      cgen(i)%w(delgen:ng-1)=gen(i+1)%w(delgen+1:ng)
    end if

    deallocate(cgen(i)%ww)
    allocate(cgen(i)%ww(ng-1,3))
    cgen(i)%ww(:ng-1,:)=gen(i+1)%ww(:ng-1,:)

    deallocate(cgen(i)%wa)
    allocate(cgen(i)%wa(nga))
    cgen(i)%wa=gen(i+1)%wa

    if (gen(i+1)%npre>0) then
      if (allocated(cgen(i)%pre)) deallocate(cgen(i)%pre)
      allocate(cgen(i)%pre(gen(i+1)%npre))
      cgen(i)%pre=gen(i+1)%pre
    end if
    if (gen(i+1)%npost>0) then
      if (allocated(cgen(i)%post)) deallocate(cgen(i)%post)
      allocate(cgen(i)%post(gen(i+1)%npost))
      cgen(i)%post=gen(i+1)%post
    end if
  end do
  deallocate(gen)
  allocate(gen(ng-1))

  do i=1,ng-1
    gen(i)=cgen(i)
    if (allocated(gen(i)%w)) deallocate(gen(i)%w)
    allocate(gen(i)%w(ng-1))
    gen(i)%w=cgen(i)%w
    if (allocated(gen(i)%ww)) deallocate(gen(i)%ww)
    allocate(gen(i)%ww(ng-1,3))
    gen(i)%ww=cgen(i)%ww
    if (allocated(gen(i)%wa)) deallocate(gen(i)%wa)
    allocate(gen(i)%wa(nga))
    gen(i)%wa=cgen(i)%wa

    if (gen(i)%npre>0) then
      if (allocated(gen(i)%pre)) deallocate(gen(i)%pre)
      allocate(gen(i)%pre(cgen(i)%npre))
      gen(i)%pre=cgen(i)%pre
    end if
    if (gen(i)%npost>0) then
      if (allocated(gen(i)%post)) deallocate(gen(i)%post)
      allocate(gen(i)%post(cgen(i)%npost))
      gen(i)%post=cgen(i)%post
    end if
  end do

end subroutine del


subroutine update_view

  use global_widgets
  use gtk, only: gtk_notebook_prev_page, gtk_notebook_next_page, gtk_notebook_get_n_pages, gtk_notebook_get_current_page
 
  integer :: numpg, curpg

  numpg=gtk_notebook_get_n_pages(notebook)
  curpg=gtk_notebook_get_current_page(notebook)

  if(curpg.ge.numpg-1)then
    call gtk_notebook_prev_page(notebook)
    call gtk_notebook_next_page(notebook)
    else
    call gtk_notebook_next_page(notebook)
    call gtk_notebook_prev_page(notebook)
  end if

end subroutine update_view

subroutine undo_redo_ini

use some_widgets

  if((counts==0))then
    cmax=50; cmaxh=25; back=1; nmax=10; togglen=0
    allocate(saveindex(cmax,2),memor2(cmax),memor3(cmax,4))
    allocate(memor5(cmax,4),memor6(cmax,6),memor7(cmax,4),memor4(cmax,3),undone1(cmax+1),memor1(cmax,6))
    allocate(ngback(nmax+1),ngback2(nmax+1))
    saveindex2=0;saveindex=0;ngback=0;ngback2=0
    record(:)=""; recf=0

  end if
  ng_o=ng; w0=-1

end subroutine undo_redo_ini

subroutine undo_redo(quidfaciam,g1,g2,g3,g4,d1,d2) ! we keep track of changes

  integer :: ngb,cmaxh,quidfaciam,g1,g2,g3,g4,m,mm,mmm,mmmm,mmmmm,mmmmmm,m1,m2,stoppp
  real*8 :: d1,d2

  if(quidfaciam.lt.8)then
    counts=counts+1
    undone1(counts+1)=counts
    undone1t=0
  endif
  if(quidfaciam.lt.8)  saveindex2(quidfaciam)=saveindex2(quidfaciam)+1
  if(quidfaciam.lt.8)  m=saveindex2(quidfaciam)

  if (counts>cmax-1)then
    cmax=cmax+50
    if(allocated(mi)) deallocate(mi)
    allocate(mi(cmax,2))
    if(allocated(undone2)) deallocate(undone2)
    allocate(undone2(cmax+1))
    mi=0
    mi(:cmax-50,:)=saveindex(:cmax-50,:)
    undone2=0
    undone2(:cmax-49)=undone1(:cmax-49)
    if(allocated(saveindex)) deallocate(saveindex)
    allocate(saveindex(cmax,2)) 
    if(allocated(undone1))then
      deallocate(undone1)
    end if
    allocate(undone1(cmax+1))
    saveindex=0
    undone1=0
    undone1(:cmax-49)=undone2(:cmax-49)
    saveindex(:cmax-50,:)=mi(:cmax-50,:)
  end if

  if(quidfaciam.lt.8) saveindex(counts,1)=quidfaciam
  if(quidfaciam.lt.8) saveindex(counts,2)=m

   if(quidfaciam==2)then 

    memor2(m)=g1 ! if negative, deletion
    if(size(memor2).eq.m)then
      if(size(m0,1).lt.m)then
        if(allocated(m0)) deallocate(m0)
        allocate(m0(m))
      end if
      m0=0
      m0(:m)=memor2(:m)
      cmaxh=m+25
      if(allocated(memor2)) deallocate(memor2)
      allocate(memor2(cmaxh)) 
      memor2=0
      memor2(:cmaxh-25)=m0(:cmaxh-25)
    end if 

  elseif(quidfaciam==3)then 
    memor3(m,1)=g1; memor3(m,2)=g2; memor3(m,3)=d1; memor3(m,4)=d2
    if(size(memor3,1).eq.m)then
      if(size(mr,1).lt.m)then
        if(allocated(mr)) deallocate(mr)
        allocate(mr(m,6))
      end if
      mr=0d0
      mr(:m,:4)=memor3(:m,:4)
      cmaxh=m+25
      if(allocated(memor3)) deallocate(memor3)
      allocate(memor3(cmaxh,4)) 
      memor3=0d0
      memor3(:cmaxh-25,:4)=mr(:cmaxh-25,:4)
    end if
  
  elseif(quidfaciam==4)then 
    memor4(m,1)=g1; memor4(m,2)=g2; memor4(m,3)=g3
    if(size(memor4,1).eq.m)then
      if(size(mr,1).lt.m)then
        if(allocated(mr)) deallocate(mr)
        allocate(mr(m,6))
      end if
      mr=0d0
      mr(:m,:3)=memor4(:m,:3)
      cmaxh=m+25
      if(allocated(memor4)) deallocate(memor4)
      allocate(memor4(cmaxh,3)) 
      memor4=0d0
      memor4(:cmaxh-25,:3)=mr(:cmaxh-25,:3)
    end if

  elseif(quidfaciam==5)then 
    memor5(m,1)=g1; memor5(m,2)=g2; memor5(m,3)=d1; memor5(m,4)=d2
    if(size(memor5,1).eq.m)then
      if(size(mr,1).lt.m)then
        if(allocated(mr)) deallocate(mr)
        allocate(mr(m,6))
      end if
      mr=0d0
      mr(:m,:4)=memor5(:m,:4)
      cmaxh=m+25
      if(allocated(memor5)) deallocate(memor5)
      allocate(memor5(cmaxh,4)) 
      memor5=0d0
      memor5(:cmaxh-25,:4)=mr(:cmaxh-25,:4)
    end if
  
  elseif(quidfaciam==7)then 
    memor7(m,1)=g1; memor7(m,2)=g2; memor3(m,3)=d1; memor3(m,4)=d2
    if(size(memor7,1).eq.m)then
      if(size(mr,1).lt.m)then
        if(allocated(mr)) deallocate(mr)
        allocate(mr(m,6))
      end if
      mr=0d0
      mr(:m,:4)=memor7(:m,:4)
      cmaxh=m+25
      if(allocated(memor7)) deallocate(memor7)
      allocate(memor7(cmaxh,4)) 
      memor7=0d0
      memor7(:cmaxh-25,:4)=mr(:cmaxh-25,:4)
    end if
  
  elseif(quidfaciam==6)then 
    memor6(m,1)=g1; memor6(m,2)=g2; memor6(m,3)=g3; memor6(m,4)=g4; memor6(m,5)=d1; memor6(m,6)=d2
    if(size(memor6,1).eq.m)then
      if(size(mr,1).lt.m)then
        if(allocated(mr)) deallocate(mr)
        allocate(mr(m,6))
      end if
      mr=0d0
      mr(:m,:6)=memor6(:m,:6)
      cmaxh=m+25
      if(allocated(memor6)) deallocate(memor6)
      allocate(memor6(cmaxh,6)) 
      memor6=0d0
      memor6(:cmaxh-25,:6)=mr(:cmaxh-25,:6)
    end if

  else if(quidfaciam==1)then
    memor1(m,1)=g1; memor1(m,2)=g2; memor1(m,3)=g3; memor1(m,4)=g4; memor1(m,5)=d1; memor1(m,6)=d2
    if(size(memor1,1).eq.m)then
      if(size(mr,1).lt.m)then
        if(allocated(mr)) deallocate(mr)
        allocate(mr(m,6))
      end if
      mr=0d0
      mr(:m,:6)=memor1(:m,:6)
      cmaxh=m+25
      if(allocated(memor1)) deallocate(memor1)
      allocate(memor1(cmaxh,6)) 
      memor1=0d0
      memor1(:cmaxh-25,:6)=mr(:cmaxh-25,:6)
    end if

  else if(quidfaciam==8)then ! retrieve previous state

    if(undone1t==0) back=0
      m=saveindex(undone1(counts+1),2)
      m1=saveindex(undone1(counts+1),1)
    if(undone1(counts+1).gt.1)then
      m2=undone1(counts+1)
      undone1(counts+1)=undone1(m2)
      stoppp=0
    else
      print*, "YOU REACHED THE ORIGINAL STATE. EVERYTHING HAS BEEN UNDONE", m, memor2(m), quidfaciam
    stoppp=1

    end if

  if(m1==2)then

    if(memor2(m).gt.0)then

      call del(ng)

      ng=ng-1
      call update_view

    elseif(memor2(m).lt.0)then ! for the deletion case, we'll need to retrieve externally saved data to restore

      call iniread 
      call readsnap(record(recf))
      recf=recf-1

   endif

  elseif(m1==3)then   ! gene properties
    mm=memor3(m,1)
    if(memor3(m,2)==1)then; gen(mm)%kindof=memor3(m,3)
    elseif(memor3(m,2)==2)then; gen(mm)%mu=memor3(m,3)
    elseif(memor3(m,2)==3)then; gen(mm)%diffu=memor3(m,3)
    elseif(memor3(m,2)==4)then; gen(mm)%idiffu=memor3(m,3)
    endif
  
  elseif(m1==4)then  ! pre/post
    if(memor4(m,3)==1)then ! undoing creation of an arrow
      mm=memor4(m,2)
        do mmm=1,gen(mm)%npre
          if(gen(mm)%pre(mmm)==memor4(m,1))then
            gen(mm)%pre(mmm)=0
            do mmmm=1,ng
              do mmmmm=1,gen(mmmm)%nww
                if((gen(mmmm)%ww(mmmmm,1)==memor4(m,1)).and.(gen(mmmm)%ww(mmmmm,2)==memor4(m,2)))then
                  gen(mmmm)%ww(mmmmm,3)=0.0d0
                  gen(mmmm)%ww(mmmmm,:2)=0 
                  if(mmmmm.lt.gen(mmmm)%nww)then
                    do mmmmmm=mmmmm,gen(mmmm)%nww-1
                      gen(mmmm)%ww(mmmmmm,:)=gen(mmmm)%ww(mmmmmm+1,:)
                    end do
                  end if
                  gen(mmmm)%ww(gen(mmmm)%nww,3)=0d0
                  gen(mmmm)%ww(gen(mmmm)%nww,:2)=0
                  gen(mmmm)%nww=gen(mmmm)%nww-1
                end if
              end do
            end do
            if(mmm.lt.gen(mm)%npre)then
              do mmmm=mmm,gen(mm)%npre-1
                gen(mm)%pre(mmmm)=gen(mm)%pre(mmmm+1)
              end do
            end if
            gen(mm)%pre(gen(mm)%npre)=0d0
            gen(mm)%npre=gen(mm)%npre-1; exit
          end if
        end do
        mm=memor4(m,1)
        do mmm=1,gen(mm)%npost
          if(gen(mm)%post(mmm)==memor4(m,2))then
            gen(mm)%post(mmm)=0
            if(mmm.lt.gen(mm)%npost)then
              do mmmm=mmm,gen(mm)%npost-1
                gen(mm)%post(mmmm)=gen(mm)%post(mmmm+1)
              end do
            end if
            gen(mm)%post(gen(mm)%npost)=0d0
            gen(mm)%npost=gen(mm)%npost-1; exit
          end if
        end do
    else ! undoing deletion of an arrow
      mm=memor4(m,1)
      gen(mm)%npost=gen(mm)%npost+1
      gen(mm)%post(gen(mm)%npost)=memor4(m,2)
      mm=memor4(m,2)
      gen(mm)%npre=gen(mm)%npre+1
      gen(mm)%pre(gen(mm)%npre)=memor4(m,1)
    endif

  elseif(m1==5)then ! w-matrix
    gen(memor5(m,2))%w(memor5(m,1))=memor5(m,3)

  elseif(m1==6)then ! r-matrix
    if(memor6(m,4)==0)then ! only ww strength changed
      do mm=1,ng
        do mmm=1,gen(mm)%nww
          if((gen(mm)%ww(mmm,1)==memor6(m,1)).and.(gen(mm)%ww(mmm,2)==memor6(m,2))) gen(mm)%ww(mmm,3)=memor6(m,5)
        end do
      end do
    else if(memor6(m,4)==1)then ! undo creation of connection
      do mm=1,ng
        do mmm=1,gen(mm)%nww
          if((gen(mm)%ww(mmm,1)==memor6(m,1)).and.(gen(mm)%ww(mmm,2)==memor6(m,2)))then
            gen(mm)%ww(mmm,3)=0.0d0
            gen(mm)%ww(mmm,:2)=0 
            if(mmm.lt.gen(mm)%nww)then
              do mmmm=mmm,gen(mm)%nww-1
                gen(mm)%ww(mmmm,:)=gen(mm)%ww(mmmm+1,:)
              end do
            end if
            gen(mm)%ww(gen(mm)%nww,3)=0d0
            gen(mm)%ww(gen(mm)%nww,:2)=0
            gen(mm)%nww=gen(mm)%nww-1
          end if
        end do
      end do
    endif

  elseif(m1==7)then ! wa-matrix
    gen(memor7(m,1))%wa(memor7(m,2))=memor7(m,3)
  
  elseif(m1==1)then ! kadh-matrix
    if(memor1(m,1)*memor1(m,2).ne.0)then
      gen(memor1(m,1))%wa(1)=memor1(m,3)
      gen(memor1(m,2))%wa(1)=memor1(m,4)
      kadh(memor1(m,3),memor1(m,4))=memor1(m,5)
      kadh(memor1(m,4),memor1(m,3))=memor1(m,6)
    else
      kadh(memor1(m,3),memor1(m,4))=memor1(m,5)
      kadh(memor1(m,4),memor1(m,3))=memor1(m,6)
    end if

  end if

  if(counts-back.gt.0)then
    back=back+1
    undone1t=1
  end if

  call update_view

endif

call update_view

end subroutine undo_redo

subroutine read_data(iin)  ! This initializes the reading of data files

  use global_widgets
  use io
  use genetic
  use general

  character*140 :: fij, dir
  character*140 :: fil
  integer :: diri, iin

  call getarg(1,dir)
print*, "THIS IS THE FILE1", dir
  if(len(trim(dir)).lt.1)then  ! if no file read in (or just one character)

    call system("size=$(pwd); echo ${#size} > path.dat ")  ! read in length of intrinsic PATH

    open(1, file='path.dat')
      read(1,*), pd
    close(1)

    call readsnap0

    if(iin==0) call undo_redo_ini

  else
    
    print*, "I AM READING THIS FILE:", trim(dir), len(trim(dir))
    print*, "If it crashes here, it is not my fault. Check your file and that you put its name correctly."

    call iniread   
    fil="cp "//trim(dir)//" dummyfile.dat"
    call system(trim(fil))
    fil='dummyfile.dat'

    call readsnap(fil)

call system("rm dummyfile.dat ") !R!

    if(iin==0) call undo_redo_ini

  end if

  faster=-1 ! this sets the "draw slow" mode as default
  style1=1 ! white background

    counts=0; back=0; ngc=0
if(allocated(genename)) deallocate(genename)
allocate(genename(ng))
genename=""

if(allocated(npos)) deallocate(npos)
  allocate(npos(ng,6))
npos=0.0d0
if(allocated(kpos)) deallocate(kpos)
  allocate(kpos(ntipusadh,2))
kpos=0.0d0

which2=-1; what2=-1

end subroutine read_data

subroutine read_data_again  ! This initializes the reading of data files

  use global_widgets
  use general
  use genetic
  use io

  character*140 :: fij, dir
  character*140 :: fil
  integer :: diri, iin

  call getarg(1,dir)
print*, "THIS IS THE FILE2", dir
  if(len(trim(dir)).lt.1)then  ! if no file read in (or just one character)

    call readsnap0

    call writesnap
    call system("cp fort.1 file1.dat") 
   ! open(3,file='file1.dat')
   !   write(3,*), "dummyfile.dat"
   ! close(3)

  else
    
    print*, "I AM READING THIS FILE:", trim(dir), len(trim(dir))
    print*, "If it crashes here, it is not my fault. Check your file and that you put its name correctly."

    call iniread

    fil="cp "//trim(dir)//" dummyfile.dat"

    call system(trim(fil))
    fil='dummyfile.dat'

call iniio
call iniarrays

  call readsnap(fil)

  open(3,file='file1.dat')
  dir=trim(dir)
    write(3,*), dir
  print*, dir
  close(3)

  call system("rm dummyfile.dat ") !R!

  end if

  faster=1 ! this sets the "draw slow" mode as default
  counts=0; back=0; ngc=0

if(allocated(genename)) deallocate(genename)
allocate(genename(ng))
genename=""

end subroutine read_data_again
subroutine read_data_again2

    use gtk, only: gtk_notebook_prev_page, gtk_notebook_next_page, gtk_notebook_get_n_pages, gtk_notebook_get_current_page
    use global_widgets
    
    character*140 :: fij,fcp

call system("size=$(pwd); echo ${#size} > path.dat ")  ! read in length of intrinsic PATH

  open(1, file='path.dat')
    read(1,*), pd
  close(1)

    open(3,file='file1.dat')
    read(3,*), fij
    close(3)

    fij=trim(fij)

    call iniread
    call iniio
    call iniarrays
    print*, "THIS IS MY NEW FILE",fij,"TT",len(fij),len(trim(fij))
    fcp="cp "//trim(fij)//" fme.dat"
    call system(trim(fcp))
    fij='fme.dat'

fij='dn.dat'
    call readsnap(fij)

    call update_view

end subroutine read_data_again2

subroutine cell_params_names

    cellparams=""
    cellparams(1)="cell growth"
    cellparams(2)="cell cycle rate"
    cellparams(3)="apoptosis"
    cellparams(4)="secretion rate"
    cellparams(5)="secretability"
    cellparams(6)="repcel secreted"
    cellparams(7)="da/reqcel secreted"
    cellparams(8)="cell polarization"
    cellparams(9)="growth polarization"
    cellparams(10)="division minsize"
    cellparams(11)="division polarizability"
    cellparams(12)="division asymmetry"
    cellparams(13)="E>M Transition"
    cellparams(14)="ECM proteolysis"
    cellparams(15)="division maxsize"
    cellparams(16)="polarization precision"        

end subroutine cell_params_names

function delete_event (widget, event, gdata) result(ret)  bind(c)

  integer(c_int)    :: ret
  type(c_ptr), value :: widget, event, gdata

  run_status = FALSE
  ret = FALSE

end function delete_event

subroutine pending_events ()
   
  do while(IAND(gtk_events_pending(), run_status) /= FALSE)
    boolresult = gtk_main_iteration_do(FALSE) ! False for non-blocking
  end do

end subroutine pending_events

end module basic_handlers

module drawings

  use iso_c_binding, only: c_int
  
  use gtk, only: gtk_container_add, gtk_drawing_area_new, gtk_events_pending, gtk&
  &_main, gtk_main_iteration, gtk_main_iteration_do, gtk_widget_get_window, gtk_w&
  &idget_show, gtk_window_new, gtk_window_set_default, gtk_window_set_default_siz&
  &e, &
  &TRUE, FALSE, c_null_char, GTK_window_TOPLEVEL, gtk_init, g_signal_connect, &
  &CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL, gtk_disable_setlocale
  
  use cairo
  
  use gdk, only: gdk_cairo_create
  use general
  use io
  use basic_handlers
  
  integer(kind=c_int) :: width, height
  
contains

  function draw_network (widget, event, gdata) result(ret)  bind(c)  ! FIRST MATRIX
    use iso_c_binding, only: c_int, c_ptr
    use io
    use global_widgets

    real*8 :: ev = 0.0000001 ! minimal value that we consider.
    integer(c_int)    :: ret
    character*140 :: fii
    character*140 :: fij
    character*3 :: str3
    character*12 :: str12
    character*8 :: str8
    type(c_ptr), value, intent(in) :: widget, event, gdata
    type(c_ptr) :: my_cairo_context, pat
    integer :: cstatus, i3, lx,ly
    integer :: t,t1
    real*8 :: tt1,tt2,tt3,tt4,i1,i2,i4,i5,xx,yy,rr,vv,bl,ta1,ta2,t1x,t1y,t2x,t2y,tr1,tr2


if(rompetout.ne.1)then
if(yetthere==1)then
call gtk_widget_destroy(my_choice_window2)
yetthere=0
endif
if(yetthere1==1)then
call gtk_widget_destroy(my_choice_window)
yetthere1=0
endif
if(yetthere3==1)then
call gtk_widget_destroy(my_choice_window3)
yetthere3=0
endif
endif
rompetout=1
    
    my_cairo_context = gdk_cairo_create (gtk_widget_get_window(widget))

if(togglen==0)then
  if(allocated(kpos)) deallocate(kpos)   !this is just to initialize
  allocate(kpos(ntipusadh,2))
  do t=1,ntipusadh
    kpos(t,1)=350d0+260d0*cos(t*pi/ntipusadh)
    kpos(t,2)=350d0+260d0*sin(t*pi/ntipusadh)
  end do
end if


if(style1==-1)then
   call cairo_set_source_rgb(my_cairo_context, 0d0, 0d0, 0d0)
   call cairo_rectangle(my_cairo_context,0d0,0d0,800d0,800d0)
   call cairo_fill(my_cairo_context)
end if

    ! Text:
    write (str8,"(I8.1)") getot
    call cairo_set_source_rgb(my_cairo_context, 0d0, 0d0, 1d0)
    call cairo_select_font_face(my_cairo_context, "Times"//c_null_char, CAIRO_FONT_SLANT_NORMAL, &
                                 &  CAIRO_FONT_WEIGHT_NORMAL)
    call cairo_set_font_size (my_cairo_context, 30d0)
    call cairo_move_to(my_cairo_context, 10d0, 30d0)
    call cairo_show_text (my_cairo_context, "This is the GRN"//c_null_char)
        call cairo_move_to(my_cairo_context, 300d0, 30d0)

    ! Legend

    lx=540;ly=-620
    call cairo_new_sub_path(my_cairo_context)
if(style1.ge.0) call cairo_set_source_rgb(my_cairo_context, 0d0, 0d0, 0d0)
if(style1.lt.0) call cairo_set_source_rgb(my_cairo_context, 1d0, 1d0, 1d0)
    call cairo_select_font_face(my_cairo_context, "Times"//c_null_char, CAIRO_FONT_SLANT_NORMAL, &
                                 &  CAIRO_FONT_WEIGHT_NORMAL)
    call cairo_set_font_size (my_cairo_context, 18d0)
    call cairo_move_to(my_cairo_context, 10d0+lx, 640d0+ly)
    call cairo_show_text (my_cairo_context, "Interactions"//c_null_char)

    call cairo_new_sub_path(my_cairo_context)
    call cairo_move_to(my_cairo_context, 70d0, 660d0)
            call cairo_set_source_rgb(my_cairo_context, 0.2d0, 0.8d0, 0.2d0)
    call cairo_set_line_width(my_cairo_context, 6.0d0)
    call cairo_move_to(my_cairo_context, 10.0d0+lx, 660.0d0+ly)  
    call cairo_line_to(my_cairo_context, 60.0d0+lx, 660.0d0+ly)    
    call cairo_stroke(my_cairo_context)  
if(style1.ge.0) call cairo_set_source_rgb(my_cairo_context, 0d0, 0d0, 0d0)
if(style1.lt.0) call cairo_set_source_rgb(my_cairo_context, 1d0, 1d0, 1d0)
    call cairo_select_font_face(my_cairo_context, "Times"//c_null_char, CAIRO_FONT_SLANT_NORMAL, &
                                 &  CAIRO_FONT_WEIGHT_NORMAL)
    call cairo_set_font_size (my_cairo_context, 14d0)
    call cairo_move_to(my_cairo_context, 70d0+lx, 663d0+ly)
    call cairo_show_text (my_cairo_context, "Upregulation"//c_null_char)

    call cairo_new_sub_path(my_cairo_context)
    call cairo_move_to(my_cairo_context, 70d0+lx, 680d0+ly)
            call cairo_set_source_rgb(my_cairo_context, 1d0, 0d0, 0d0)
    call cairo_set_line_width(my_cairo_context, 6.0d0)
    call cairo_move_to(my_cairo_context, 10.0d0+lx, 680.0d0+ly)  
    call cairo_line_to(my_cairo_context, 60.0d0+lx, 680.0d0+ly)    
    call cairo_stroke(my_cairo_context)
if(style1.ge.0) call cairo_set_source_rgb(my_cairo_context, 0d0, 0d0, 0d0)
if(style1.lt.0) call cairo_set_source_rgb(my_cairo_context, 1d0, 1d0, 1d0)
    call cairo_select_font_face(my_cairo_context, "Times"//c_null_char, CAIRO_FONT_SLANT_NORMAL, &
                                 &  CAIRO_FONT_WEIGHT_NORMAL)
    call cairo_set_font_size (my_cairo_context, 14d0)
    call cairo_move_to(my_cairo_context, 70d0+lx, 683d0+ly)
    call cairo_show_text (my_cairo_context, "Downregulation"//c_null_char)

    ! Legend:

    call cairo_set_source_rgb(my_cairo_context, 0d0, 0d0, 1d0)

    call cairo_set_font_size (my_cairo_context, 18d0)

if(style1.ge.0) call cairo_set_source_rgb(my_cairo_context, 0d0, 0d0, 0d0)
if(style1.lt.0) call cairo_set_source_rgb(my_cairo_context, 1d0, 1d0, 1d0)
    call cairo_move_to(my_cairo_context, 10d0, 620d0)
    call cairo_show_text (my_cairo_context, "Types of Molecules:"//c_null_char)
             call cairo_set_source_rgb(my_cairo_context, 0.8d0, 0.8d0, 0.8d0) 
    call cairo_move_to(my_cairo_context, 10d0, 650d0)

         xx=15d0; yy=645d0
         pat = cairo_pattern_create_radial (xx-3d0, yy-3d0, 2.3d0, xx-4d0, yy-4d0, 14.0d0)
         call cairo_pattern_add_color_stop_rgba (pat, 0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0)
         call cairo_pattern_add_color_stop_rgba (pat, 1d0, 0.8d0, 0.8d0, 0.8d0, 1.0d0)
         call cairo_set_source (my_cairo_context, pat)
         call cairo_arc(my_cairo_context, xx, yy, 10d0, 0d0, 2*pi)
         call cairo_fill (my_cairo_context)
         call cairo_pattern_destroy (pat)

if(style1.ge.0) call cairo_set_source_rgb(my_cairo_context, 0d0, 0d0, 0d0)
if(style1.lt.0) call cairo_set_source_rgb(my_cairo_context, 1d0, 1d0, 1d0)
    call cairo_move_to(my_cairo_context, 10d0, 650d0)
    call cairo_show_text (my_cairo_context, "1 : Regulatory Molecule "//c_null_char)
             call cairo_set_source_rgb(my_cairo_context, 0.0d0, 0.6d0, 0.0d0)
    call cairo_move_to(my_cairo_context, 10d0, 670d0)

         xx=15d0; yy=665d0
         pat = cairo_pattern_create_radial (xx-3d0, yy-3d0, 2.3d0, xx-4d0, yy-4d0, 14.0d0)
         call cairo_pattern_add_color_stop_rgba (pat, 0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0)
         call cairo_pattern_add_color_stop_rgba (pat, 1d0, 0.0d0, 0.6d0, 0.0d0, 1.0d0)
         call cairo_set_source (my_cairo_context, pat)
         call cairo_arc(my_cairo_context, xx, yy, 10d0, 0d0, 2*pi)
         call cairo_fill (my_cairo_context)
         call cairo_pattern_destroy (pat)

    call cairo_set_source_rgb(my_cairo_context, 0.0d0, 0.6d0, 0.0d0)
    call cairo_move_to(my_cairo_context, 10d0, 670d0)
    call cairo_show_text (my_cairo_context, "2 : Regulatory Molecule, translatable"//c_null_char)
             call cairo_set_source_rgb(my_cairo_context, 0.0d0, 0.0d0, 0.8d0)
    call cairo_move_to(my_cairo_context, 10d0, 690d0)

         xx=15d0; yy=685d0
         pat = cairo_pattern_create_radial (xx-3d0, yy-3d0, 2.3d0, xx-4d0, yy-4d0, 14.0d0)
         call cairo_pattern_add_color_stop_rgba (pat, 0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0)
         call cairo_pattern_add_color_stop_rgba (pat, 1d0, 0.0d0, 0.0d0, 0.8d0, 1.0d0)
         call cairo_set_source (my_cairo_context, pat)
         call cairo_arc(my_cairo_context, xx, yy, 10d0, 0d0, 2*pi)
         call cairo_fill (my_cairo_context)
         call cairo_pattern_destroy (pat)

    call cairo_set_source_rgb(my_cairo_context, 0.0d0, 0.0d0, 0.8d0)
    call cairo_move_to(my_cairo_context, 10d0, 690d0)
    call cairo_show_text (my_cairo_context, "3 : Protein, intracellular"//c_null_char)
             call cairo_set_source_rgb(my_cairo_context, 0.8d0, 0.0d0, 0.0d0)

         xx=15d0; yy=705d0
         pat = cairo_pattern_create_radial (xx-3d0, yy-3d0, 2.3d0, xx-4d0, yy-4d0, 14.0d0)
         call cairo_pattern_add_color_stop_rgba (pat, 0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0)
         call cairo_pattern_add_color_stop_rgba (pat, 1d0, 0.8d0, 0.0d0, 0.0d0, 1.0d0)
         call cairo_set_source (my_cairo_context, pat)
         call cairo_arc(my_cairo_context, xx, yy, 10d0, 0d0, 2*pi)
         call cairo_fill (my_cairo_context)
         call cairo_pattern_destroy (pat)

    call cairo_set_source_rgb(my_cairo_context, 0.8d0, 0.0d0, 0.0d0)
    call cairo_move_to(my_cairo_context, 10d0, 710d0)
    call cairo_show_text (my_cairo_context, "4 : Protein, extracellular"//c_null_char)
             call cairo_set_source_rgb(my_cairo_context, 1.0d0, 0.5d0, 1.0d0)

         xx=445d0; yy=645d0
         pat = cairo_pattern_create_radial (xx-3d0, yy-3d0, 2.3d0, xx-4d0, yy-4d0, 14.0d0)
         call cairo_pattern_add_color_stop_rgba (pat, 0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0)
         call cairo_pattern_add_color_stop_rgba (pat, 1d0, 1.0d0, 0.5d0, 1.0d0, 1.0d0)
         call cairo_set_source (my_cairo_context, pat)
         call cairo_arc(my_cairo_context, xx, yy, 10d0, 0d0, 2*pi)
         call cairo_fill (my_cairo_context)
         call cairo_pattern_destroy (pat)

    call cairo_set_source_rgb(my_cairo_context, 1.0d0, 0.5d0, 1.0d0)
    call cairo_move_to(my_cairo_context, 440d0, 650d0)
    call cairo_show_text (my_cairo_context, "5 : Protein, apically localizing"//c_null_char)
             call cairo_set_source_rgb(my_cairo_context, 0.6d0, 0.6d0, 1.0d0)

         xx=445d0; yy=665d0
         pat = cairo_pattern_create_radial (xx-3d0, yy-3d0, 2.3d0, xx-4d0, yy-4d0, 14.0d0)
         call cairo_pattern_add_color_stop_rgba (pat, 0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0)
         call cairo_pattern_add_color_stop_rgba (pat, 1d0, 0.6d0, 0.6d0, 1.0d0, 1.0d0)
         call cairo_set_source (my_cairo_context, pat)
         call cairo_arc(my_cairo_context, xx, yy, 10d0, 0d0, 2*pi)
         call cairo_fill (my_cairo_context)
         call cairo_pattern_destroy (pat)

    call cairo_set_source_rgb(my_cairo_context, 0.6d0, 0.6d0, 1.0d0)
    call cairo_move_to(my_cairo_context, 440d0, 670d0)
    call cairo_show_text (my_cairo_context, "6 : Protein, basally localizing"//c_null_char)
             call cairo_set_source_rgb(my_cairo_context, 0.2d0, 1.0d0, 0.2d0)

         xx=445d0; yy=685d0
         pat = cairo_pattern_create_radial (xx-3d0, yy-3d0, 2.3d0, xx-4d0, yy-4d0, 14.0d0)
         call cairo_pattern_add_color_stop_rgba (pat, 0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0)
         call cairo_pattern_add_color_stop_rgba (pat, 1d0, 0.2d0, 1.0d0, 0.2d0, 1.0d0)
         call cairo_set_source (my_cairo_context, pat)
         call cairo_arc(my_cairo_context, xx, yy, 10d0, 0d0, 2*pi)
         call cairo_fill (my_cairo_context)
         call cairo_pattern_destroy (pat)

    call cairo_set_source_rgb(my_cairo_context, 0.2d0, 0.9d0, 0.2d0)
    call cairo_move_to(my_cairo_context, 440d0, 690d0)
    call cairo_show_text (my_cairo_context, "7 : Receptor, juxtacrine signalling(7)"//c_null_char)
             call cairo_set_source_rgb(my_cairo_context, 0.8d0, 0.8d0, 0.2d0)

         xx=445d0; yy=705d0
         pat = cairo_pattern_create_radial (xx-3d0, yy-3d0, 2.3d0, xx-4d0, yy-4d0, 14.0d0)
         call cairo_pattern_add_color_stop_rgba (pat, 0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0)
         call cairo_pattern_add_color_stop_rgba (pat, 1d0, 1.0d0, 0.8d0, 0.0d0, 1.0d0)
         call cairo_set_source (my_cairo_context, pat)
         call cairo_arc(my_cairo_context, xx, yy, 10d0, 0d0, 2*pi)
         call cairo_fill (my_cairo_context)
         call cairo_pattern_destroy (pat)

    call cairo_set_source_rgb(my_cairo_context, 1.0d0, 0.8d0, 0.0d0)
    call cairo_move_to(my_cairo_context, 440d0, 710d0)
    call cairo_show_text (my_cairo_context, "8 : Receptor, paracrine signalling(4)"//c_null_char)

if(togglen==0)then
  if(allocated(npos)) deallocate(npos)
  allocate(npos(ng,6))
  do t=1,ng
    npos(t,1)=350d0+260d0*cos(t*2d0*pi/ng)
    npos(t,2)=350d0+260d0*sin(t*2d0*pi/ng)
    npos(t,3)=350d0+260d0*cos(t*pi/ng)
    npos(t,4)=350d0+260d0*sin(t*pi/ng)
    npos(t,5)=50d0+700d0/(ng+1)*t
    npos(t,6)=150d0
  end do
end if

    call cairo_new_sub_path(my_cairo_context)

      do t = 1,ng

            call cairo_move_to(my_cairo_context, npos(t,1)+20, npos(t,2))
                xx=npos(t,1); yy=npos(t,2)

if(wcl==t)then

            pat = cairo_pattern_create_radial (xx, yy, 5.6d0, xx, yy, 45.0d0)

                call cairo_pattern_add_color_stop_rgba (pat,0d0, 1d0, 0.8d0, 0d0, 1.0d0)
                if(style1.ge.0) call cairo_pattern_add_color_stop_rgba (pat,1d0, 0.97d0, 0.97d0, 0.97d0, 1.0d0)
                if(style1.lt.0) call cairo_pattern_add_color_stop_rgba (pat,1d0, 0d0, 0d0, 0d0, 1.0d0)
                call cairo_set_source (my_cairo_context,pat)
                call cairo_arc(my_cairo_context, xx, yy, 45d0, 0d0, 2*pi)
                call cairo_fill (my_cairo_context)
                call cairo_pattern_destroy (pat)
end if
        do t1=1,ng
          if(gen(t)%w(t1).gt.ev)then ! new w positive
            call cairo_set_source_rgb(my_cairo_context, 0.6d0, 0.8d0, 0.6d0)
  call cairo_set_source_rgb(my_cairo_context, 0.2d0, 0.8d0, 0.2d0)

              tt1=0.3d0+2.0d0*sqrt(gen(t)%w(t1))
              if(tt1.gt.10) tt1=10
              if(tt1.lt.0.1) tt1=5
              tt1=5
            call cairo_set_line_width(my_cairo_context, tt1)
            call cairo_move_to(my_cairo_context, npos(t,1), npos(t,2))  
            call cairo_line_to(my_cairo_context, npos(t1,1), npos(t1,2))           
            if(t1==t) then
                call cairo_move_to(my_cairo_context, npos(t,1), npos(t,2)+10d0)
                call cairo_arc(my_cairo_context, npos(t,1)-20d0, npos(t,2)+10d0, 22d0, 0d0, 2*pi)
            end if
            call cairo_stroke(my_cairo_context) 
          else if(gen(t)%w(t1).lt.-ev) then  ! new w negative
            call cairo_set_source_rgb(my_cairo_context, 0.9d0, 0.6d0, 0.6d0)
call cairo_set_source_rgb(my_cairo_context, 1.0d0, 0.2d0, 0.2d0)
              tt1=0.3d0+2.0d0*sqrt(gen(t)%w(t1)*-1)
              if(tt1.gt.10) tt1=10
              if(tt1.lt.0.1) tt1=5
              tt1=5
            call cairo_set_line_width(my_cairo_context, tt1)
            call cairo_move_to(my_cairo_context, npos(t,1), npos(t,2))  
            call cairo_line_to(my_cairo_context, npos(t1,1), npos(t1,2))
            if(t1==t) then
                call cairo_move_to(my_cairo_context, npos(t,1)-0d0, npos(t,2)+10d0)
                call cairo_arc(my_cairo_context, npos(t,1)-20d0, npos(t,2)+10d0, 22d0, 0d0, 2*pi)
            end if
            call cairo_stroke(my_cairo_context) 
          end if
        end do
      end do

      do t=1,ng 
        do t1=1,ng
          if(gen(t)%w(t1).gt.ev)then

            if(t==t1)then
              call cairo_set_source_rgb(my_cairo_context, 0.2d0, 0.8d0, 0.2d0)
              ta1=npos(t,1)-sqrt(200d0)-0d0; ta2=npos(t,2)-sqrt(200d0)+3d0
              call cairo_move_to(my_cairo_context, ta1, ta2)
              call cairo_line_to(my_cairo_context, ta1-10d0, ta2+10d0)
              call cairo_line_to(my_cairo_context, ta1-10d0, ta2-10d0)
              call cairo_close_path(my_cairo_context)
              call cairo_fill(my_cairo_context)
              call cairo_stroke(my_cairo_context)
            endif

call cairo_set_source_rgb(my_cairo_context, 0.2d0, 0.8d0, 0.2d0)
              t1x=npos(t,1);t1y=npos(t,2); t2x=npos(t1,1);t2y=npos(t1,2)
              ta1=t2x-t1x; ta2=t2y-t1y
              tt1=sqrt(ta1**2+ta2**2)
              ta1=ta1*20d0/tt1; ta2=ta2*20d0/tt1 ! this is the x/y that lies behind the sphere
              t1x=t1x+ta1; t1y=t1y+ta2; t2x=t2x-ta1; t2y=t2y-ta2

              tt1=0.3d0+2.0d0*sqrt(abs(gen(t)%w(t1)))
              if(tt1.gt.10) tt1=10
              if(tt1.lt.0.1) tt1=0.1
              tt1=5
              call cairo_set_line_width(my_cairo_context, tt1)
              tt1=0.6d0+4.0d0*sqrt(abs(gen(t)%w(t1)))
              if(tt1.gt.12) tt1=12
              if(tt1.lt.3) tt1=3
              tt1=5
              tt3=(t2x-t1x)
              tt4=(t2y-t1y)
              tt2=sqrt((tt3)**2+(tt4)**2)

              if(tt2.le.70d0)then
                tr1=tt3*0.2; tr2=tt4*0.2
              else
                tr1=30d0/tt2*tt3; tr2=30d0/tt2*tt4
              end if

              call cairo_move_to(my_cairo_context, t1x+tr1, t1y+tr2)  
              call cairo_line_to(my_cairo_context, t1x+tt3*0.5, t1y+tt4*0.5)
              call cairo_stroke(my_cairo_context)
              call cairo_set_source_rgb(my_cairo_context, 0.2d0, 0.8d0, 0.2d0)
              call cairo_set_line_width(my_cairo_context, tt1)
              call cairo_move_to(my_cairo_context, t1x+tr1, t1y+tr2)  
              call cairo_line_to(my_cairo_context, t1x+tt3*0.2, t1y+tt4*0.2)
              call cairo_stroke(my_cairo_context) 
              ! arrowheads

              if(tt2.le.70d0)then
                tr1=tt3*0.2; tr2=tt4*0.2
              else
                tr1=30d0/tt2*tt3; tr2=30d0/tt2*tt4
              end if

              call cairo_move_to(my_cairo_context, t1x+tr1, t1y+tr2)
              ta1=t1x+tr1-(t1y-t2y)
              ta2=t1y+tr2+(t1x-t2x)
              tt1=sqrt((t1x-ta1)**2+(t1y-ta2)**2)
              ta1=t1x+tr1-10d0/tt1*(t1y-t2y)
              ta2=t1y+tr2+10d0/tt1*(t1x-t2x)
              call cairo_move_to(my_cairo_context, ta1, ta2)

              call cairo_line_to(my_cairo_context, t1x, t1y)
              ta1=t1x+tr1-(t1y-t2y)
              ta2=t1y+tr2+(t1x-t2x)
              tt1=sqrt((t1x-ta1)**2+(t1y-ta2)**2)
              ta1=t1x+tr1+10d0/tt1*(t1y-t2y)
              ta2=t1y+tr2-10d0/tt1*(t1x-t2x)
              call cairo_line_to(my_cairo_context, ta1, ta2)
              call cairo_close_path(my_cairo_context)
              call cairo_fill(my_cairo_context)
              call cairo_stroke(my_cairo_context)

          elseif(gen(t)%w(t1).lt.-ev)then

            if(t==t1)then
              call cairo_set_source_rgb(my_cairo_context, 1.0d0, 0.2d0, 0.2d0)
              ta1=npos(t,1)+sqrt(200d0); ta2=npos(t,2)-sqrt(200d0)
              call cairo_move_to(my_cairo_context, ta1-36d0, ta2+14d0)
              call cairo_line_to(my_cairo_context, ta1-26d0, ta2-8d0)
              call cairo_stroke(my_cairo_context)
            endif

call cairo_set_source_rgb(my_cairo_context, 1.0d0, 0.2d0, 0.2d0)
              t1x=npos(t,1);t1y=npos(t,2); t2x=npos(t1,1);t2y=npos(t1,2)
              ta1=t2x-t1x; ta2=t2y-t1y
              tt1=sqrt(ta1**2+ta2**2)
              ta1=ta1*20d0/tt1; ta2=ta2*20d0/tt1 ! this is the x/y that lies behind the sphere
              t1x=t1x+ta1; t1y=t1y+ta2; t2x=t2x-ta1; t2y=t2y-ta2
           
              tt1=0.3d0+2.0d0*sqrt(abs(gen(t)%w(t1)))
              if(tt1.gt.10) tt1=10
              if(tt1.lt.0.1) tt1=0.1
              tt1=5
              call cairo_set_line_width(my_cairo_context, tt1)
              tt1=0.6d0+4.0d0*sqrt(abs(gen(t)%w(t1)))
              if(tt1.gt.12) tt1=12
              if(tt1.lt.3) tt1=3
              tt1=5
              tt3=(t2x-t1x)
              tt4=(t2y-t1y)
              tt2=sqrt((tt3)**2+(tt4)**2)
              call cairo_move_to(my_cairo_context, t1x, t1y)  
              call cairo_line_to(my_cairo_context, t1x+tt3*0.5, t1y+tt4*0.5)
              call cairo_stroke(my_cairo_context)
              call cairo_set_source_rgb(my_cairo_context, 1.0d0, 0.2d0, 0.2d0)
              call cairo_set_line_width(my_cairo_context, tt1)
              call cairo_move_to(my_cairo_context, t1x, t1y)  
              call cairo_line_to(my_cairo_context, t1x+tt3*0.2, t1y+tt4*0.2)
              call cairo_stroke(my_cairo_context)

              ! arrowheads

              tr1=0d0; tr2=0d0
              call cairo_move_to(my_cairo_context, t1x+tr1, t1y+tr2)
              ta1=t1x+tr1-(t1y-t2y)
              ta2=t1y+tr2+(t1x-t2x)
              tt1=sqrt((t1x-ta1)**2+(t1y-ta2)**2)
              ta1=t1x+tr1-10d0/tt1*(t1y-t2y)
              ta2=t1y+tr2+10d0/tt1*(t1x-t2x)
              call cairo_line_to(my_cairo_context, ta1, ta2)
              call cairo_stroke(my_cairo_context)

              call cairo_move_to(my_cairo_context, t1x+tr1, t1y+tr2)
              ta1=t1x+tr1-(t1y-t2y)
              ta2=t1y+tr2+(t1x-t2x)
              tt1=sqrt((t1x-ta1)**2+(t1y-ta2)**2)
              ta1=t1x+tr1+10d0/tt1*(t1y-t2y)
              ta2=t1y+tr2-10d0/tt1*(t1x-t2x)
              call cairo_line_to(my_cairo_context, ta1, ta2)
              call cairo_stroke(my_cairo_context)

          endif

        end do

            call cairo_move_to(my_cairo_context, npos(t,1)+20, npos(t,2))
                xx=npos(t,1); yy=npos(t,2)

       if(gen(t)%kindof==1)then
             rr=0.6d0; vv=0.6d0; bl=0.6d0
             call cairo_set_source_rgb(my_cairo_context, 0.7d0, 0.7d0, 0.7d0)
       elseif(gen(t)%kindof==2)then
               rr=0.0d0; vv=0.5d0; bl=0.0d0
             call cairo_set_source_rgb(my_cairo_context, 0.0d0, 0.4d0, 0.0d0)
       elseif(gen(t)%kindof==3)then
             rr=0.0d0; vv=0.0d0; bl=0.8d0
             call cairo_set_source_rgb(my_cairo_context, 0.0d0, 0.0d0, 0.7d0)
       elseif(gen(t)%kindof==4)then
             rr=0.8d0; vv=0.0d0; bl=0.0d0
             call cairo_set_source_rgb(my_cairo_context, 0.8d0, 0.0d0, 0.0d0)
       elseif(gen(t)%kindof==5)then
             rr=1.0d0; vv=0.5d0; bl=1.0d0
             call cairo_set_source_rgb(my_cairo_context, 1.0d0, 0.5d0, 1.0d0)
       elseif(gen(t)%kindof==6)then
             rr=0.6d0; vv=0.6d0; bl=1.0d0
             call cairo_set_source_rgb(my_cairo_context, 0.4d0, 0.7d0, 1.0d0)
       elseif(gen(t)%kindof==7)then
             rr=0.2d0; vv=1.0d0; bl=0.2d0
             call cairo_set_source_rgb(my_cairo_context, 0.2d0, 1.0d0, 0.2d0)
       elseif(gen(t)%kindof==8)then
             rr=1.0d0; vv=0.8d0; bl=0.0d0
             call cairo_set_source_rgb(my_cairo_context, 1.0d0, 1.0d0, 0.2d0)
       endif

                pat = cairo_pattern_create_radial (xx-6d0, yy-6d0, 5.6d0, xx-8d0, yy-8d0, 28.0d0)
                call cairo_pattern_add_color_stop_rgba (pat, 0d0, 1d0, 1d0, 1d0, 1.0d0)
                call cairo_pattern_add_color_stop_rgba (pat, 1d0, rr, vv, bl, 1.0d0)
                call cairo_set_source (my_cairo_context, pat)
                call cairo_arc(my_cairo_context, xx, yy, 20d0, 0d0, 2*pi)
                call cairo_fill (my_cairo_context)
                call cairo_pattern_destroy (pat)

if(style1.ge.0) call cairo_set_source_rgb(my_cairo_context, 0d0, 0d0, 0d0)
if(style1.lt.0) call cairo_set_source_rgb(my_cairo_context, 1d0, 1d0, 1d0)
        call cairo_set_line_width(my_cairo_context, 3.0d0)
        call cairo_select_font_face(my_cairo_context, "Times"//c_null_char, CAIRO_FONT_SLANT_NORMAL, &
                                 &  CAIRO_FONT_WEIGHT_NORMAL)
        call cairo_set_source_rgb(my_cairo_context, 0d0, 0d0, 0d0)
        call cairo_set_font_size (my_cairo_context, 20d0)
        call cairo_move_to(my_cairo_context, npos(t,1)-20d0, npos(t,2)+5d0)

        if(len(trim(genename(t))).eq.0)then

          write (str3,"(I3.1)") t
          call cairo_show_text (my_cairo_context, str3//c_null_char)
        else
          write (str12,"(A12)") genename(t)
          if(len(trim(str12)).ge.4) call cairo_set_font_size (my_cairo_context, 20d0-len(trim(str12)))
          call cairo_show_text (my_cairo_context, str12//c_null_char)
        end if
        
    end do
    
    call cairo_destroy(my_cairo_context)
    ret = FALSE
  
    call gtk_widget_queue_draw(my_drawing_area)

  end function draw_network


  function draw_network3 (widget, event, gdata) result(ret)  bind(c)
    use iso_c_binding, only: c_int, c_ptr
    use io
    use global_widgets
    integer(c_int)    :: ret
    character*140 :: fii
    character*140 :: fij
    character*3 :: str3
    character*12 :: str12
    character*8 :: str8
    type(c_ptr), value, intent(in) :: widget, event, gdata
    type(c_ptr) :: my_cairo_context, pat
    integer :: cstatus
    integer :: t,t1,t2
    real*8 :: tt1,tt2,tt3,tt4,xx,yy,rr,vv,bl

    my_cairo_context = gdk_cairo_create (gtk_widget_get_window(widget))

if(rompetout.ne.2)then
if(yetthere==1)then
call gtk_widget_destroy(my_choice_window2)
yetthere=0
endif
if(yetthere1==1)then
call gtk_widget_destroy(my_choice_window)
yetthere1=0
endif
if(yetthere3==1)then
call gtk_widget_destroy(my_choice_window3)
yetthere3=0
endif
endif
rompetout=2

if(style1==-1)then
   call cairo_set_source_rgb(my_cairo_context, 0d0, 0d0, 0d0)
   call cairo_rectangle(my_cairo_context,0d0,0d0,800d0,800d0)
   call cairo_fill(my_cairo_context)
end if

    ! Legend:
    write (str8,"(I8.1)") getot
    call cairo_set_source_rgb(my_cairo_context, 0d0, 0d0, 1d0)
    call cairo_select_font_face(my_cairo_context, "Times"//c_null_char, CAIRO_FONT_SLANT_NORMAL, &
                                 &  CAIRO_FONT_WEIGHT_NORMAL)
    call cairo_set_font_size (my_cairo_context, 30d0)
    call cairo_move_to(my_cairo_context, 10d0, 30d0)
    call cairo_show_text (my_cairo_context, "This is the Network of Molecule Modifications"//c_null_char)
        call cairo_move_to(my_cairo_context, 300d0, 30d0)

    call cairo_set_font_size (my_cairo_context, 18d0)

if(style1.ge.0) call cairo_set_source_rgb(my_cairo_context, 0d0, 0d0, 0d0)
if(style1.lt.0) call cairo_set_source_rgb(my_cairo_context, 1d0, 1d0, 1d0)
    call cairo_move_to(my_cairo_context, 10d0, 620d0)
    call cairo_show_text (my_cairo_context, "Types of Molecules:"//c_null_char)
             call cairo_set_source_rgb(my_cairo_context, 0.8d0, 0.8d0, 0.8d0) 
    call cairo_move_to(my_cairo_context, 10d0, 650d0)

         xx=15d0; yy=645d0
         pat = cairo_pattern_create_radial (xx-3d0, yy-3d0, 2.3d0, xx-4d0, yy-4d0, 14.0d0)
         call cairo_pattern_add_color_stop_rgba (pat, 0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0)
         call cairo_pattern_add_color_stop_rgba (pat, 1d0, 0.8d0, 0.8d0, 0.8d0, 1.0d0)
         call cairo_set_source (my_cairo_context, pat)
         call cairo_arc(my_cairo_context, xx, yy, 10d0, 0d0, 2*pi)
         call cairo_fill (my_cairo_context)
         call cairo_pattern_destroy (pat)

if(style1.ge.0) call cairo_set_source_rgb(my_cairo_context, 0d0, 0d0, 0d0)
if(style1.lt.0) call cairo_set_source_rgb(my_cairo_context, 1d0, 1d0, 1d0)
    call cairo_move_to(my_cairo_context, 10d0, 650d0)
    call cairo_show_text (my_cairo_context, "1 : Regulatory Molecule "//c_null_char)
             call cairo_set_source_rgb(my_cairo_context, 0.0d0, 0.6d0, 0.0d0)
    call cairo_move_to(my_cairo_context, 10d0, 670d0)

         xx=15d0; yy=665d0
         pat = cairo_pattern_create_radial (xx-3d0, yy-3d0, 2.3d0, xx-4d0, yy-4d0, 14.0d0)
         call cairo_pattern_add_color_stop_rgba (pat, 0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0)
         call cairo_pattern_add_color_stop_rgba (pat, 1d0, 0.0d0, 0.6d0, 0.0d0, 1.0d0)
         call cairo_set_source (my_cairo_context, pat)
         call cairo_arc(my_cairo_context, xx, yy, 10d0, 0d0, 2*pi)
         call cairo_fill (my_cairo_context)
         call cairo_pattern_destroy (pat)

    call cairo_set_source_rgb(my_cairo_context, 0.0d0, 0.6d0, 0.0d0)
    call cairo_move_to(my_cairo_context, 10d0, 670d0)
    call cairo_show_text (my_cairo_context, "2 : Regulatory Molecule, translatable"//c_null_char)
             call cairo_set_source_rgb(my_cairo_context, 0.0d0, 0.0d0, 0.8d0)
    call cairo_move_to(my_cairo_context, 10d0, 690d0)

         xx=15d0; yy=685d0
         pat = cairo_pattern_create_radial (xx-3d0, yy-3d0, 2.3d0, xx-4d0, yy-4d0, 14.0d0)
         call cairo_pattern_add_color_stop_rgba (pat, 0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0)
         call cairo_pattern_add_color_stop_rgba (pat, 1d0, 0.0d0, 0.0d0, 0.8d0, 1.0d0)
         call cairo_set_source (my_cairo_context, pat)
         call cairo_arc(my_cairo_context, xx, yy, 10d0, 0d0, 2*pi)
         call cairo_fill (my_cairo_context)
         call cairo_pattern_destroy (pat)

    call cairo_set_source_rgb(my_cairo_context, 0.0d0, 0.0d0, 0.8d0)
    call cairo_move_to(my_cairo_context, 10d0, 690d0)
    call cairo_show_text (my_cairo_context, "3 : Protein, intracellular"//c_null_char)
             call cairo_set_source_rgb(my_cairo_context, 0.8d0, 0.0d0, 0.0d0)

         xx=15d0; yy=705d0
         pat = cairo_pattern_create_radial (xx-3d0, yy-3d0, 2.3d0, xx-4d0, yy-4d0, 14.0d0)
         call cairo_pattern_add_color_stop_rgba (pat, 0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0)
         call cairo_pattern_add_color_stop_rgba (pat, 1d0, 0.8d0, 0.0d0, 0.0d0, 1.0d0)
         call cairo_set_source (my_cairo_context, pat)
         call cairo_arc(my_cairo_context, xx, yy, 10d0, 0d0, 2*pi)
         call cairo_fill (my_cairo_context)
         call cairo_pattern_destroy (pat)

    call cairo_set_source_rgb(my_cairo_context, 0.8d0, 0.0d0, 0.0d0)
    call cairo_move_to(my_cairo_context, 10d0, 710d0)
    call cairo_show_text (my_cairo_context, "4 : Protein, extracellular"//c_null_char)
             call cairo_set_source_rgb(my_cairo_context, 1.0d0, 0.5d0, 1.0d0)

         xx=445d0; yy=645d0
         pat = cairo_pattern_create_radial (xx-3d0, yy-3d0, 2.3d0, xx-4d0, yy-4d0, 14.0d0)
         call cairo_pattern_add_color_stop_rgba (pat, 0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0)
         call cairo_pattern_add_color_stop_rgba (pat, 1d0, 1.0d0, 0.5d0, 1.0d0, 1.0d0)
         call cairo_set_source (my_cairo_context, pat)
         call cairo_arc(my_cairo_context, xx, yy, 10d0, 0d0, 2*pi)
         call cairo_fill (my_cairo_context)
         call cairo_pattern_destroy (pat)

    call cairo_set_source_rgb(my_cairo_context, 1.0d0, 0.5d0, 1.0d0)
    call cairo_move_to(my_cairo_context, 440d0, 650d0)
    call cairo_show_text (my_cairo_context, "5 : Protein, apically localizing"//c_null_char)
             call cairo_set_source_rgb(my_cairo_context, 0.6d0, 0.6d0, 1.0d0)

         xx=445d0; yy=665d0
         pat = cairo_pattern_create_radial (xx-3d0, yy-3d0, 2.3d0, xx-4d0, yy-4d0, 14.0d0)
         call cairo_pattern_add_color_stop_rgba (pat, 0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0)
         call cairo_pattern_add_color_stop_rgba (pat, 1d0, 0.6d0, 0.6d0, 1.0d0, 1.0d0)
         call cairo_set_source (my_cairo_context, pat)
         call cairo_arc(my_cairo_context, xx, yy, 10d0, 0d0, 2*pi)
         call cairo_fill (my_cairo_context)
         call cairo_pattern_destroy (pat)

    call cairo_set_source_rgb(my_cairo_context, 0.6d0, 0.6d0, 1.0d0)
    call cairo_move_to(my_cairo_context, 440d0, 670d0)
    call cairo_show_text (my_cairo_context, "6 : Protein, basally localizing"//c_null_char)
             call cairo_set_source_rgb(my_cairo_context, 0.2d0, 1.0d0, 0.2d0)

         xx=445d0; yy=685d0
         pat = cairo_pattern_create_radial (xx-3d0, yy-3d0, 2.3d0, xx-4d0, yy-4d0, 14.0d0)
         call cairo_pattern_add_color_stop_rgba (pat, 0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0)
         call cairo_pattern_add_color_stop_rgba (pat, 1d0, 0.2d0, 1.0d0, 0.2d0, 1.0d0)
         call cairo_set_source (my_cairo_context, pat)
         call cairo_arc(my_cairo_context, xx, yy, 10d0, 0d0, 2*pi)
         call cairo_fill (my_cairo_context)
         call cairo_pattern_destroy (pat)

    call cairo_set_source_rgb(my_cairo_context, 0.2d0, 0.9d0, 0.2d0)
    call cairo_move_to(my_cairo_context, 440d0, 690d0)
    call cairo_show_text (my_cairo_context, "7 : Receptor, juxtacrine signalling(7)"//c_null_char)
             call cairo_set_source_rgb(my_cairo_context, 0.8d0, 0.8d0, 0.2d0)

         xx=445d0; yy=705d0
         pat = cairo_pattern_create_radial (xx-3d0, yy-3d0, 2.3d0, xx-4d0, yy-4d0, 14.0d0)
         call cairo_pattern_add_color_stop_rgba (pat, 0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0)
         call cairo_pattern_add_color_stop_rgba (pat, 1d0, 1.0d0, 0.8d0, 0.0d0, 1.0d0)
         call cairo_set_source (my_cairo_context, pat)
         call cairo_arc(my_cairo_context, xx, yy, 10d0, 0d0, 2*pi)
         call cairo_fill (my_cairo_context)
         call cairo_pattern_destroy (pat)

    call cairo_set_source_rgb(my_cairo_context, 1.0d0, 0.8d0, 0.0d0)
    call cairo_move_to(my_cairo_context, 440d0, 710d0)
    call cairo_show_text (my_cairo_context, "8 : Receptor, paracrine signalling(4)"//c_null_char)

    ! Arrows:
    call cairo_new_sub_path(my_cairo_context)

if(togglen==0)then
  if(allocated(npos)) deallocate(npos)
  allocate(npos(ng,6))
  do t=1,ng
    npos(t,1)=350d0+260d0*cos(t*2d0*pi/ng)
    npos(t,2)=350d0+260d0*sin(t*2d0*pi/ng)
    npos(t,3)=350d0+260d0*cos(t*pi/ng)
    npos(t,4)=350d0+260d0*sin(t*pi/ng)
    npos(t,5)=50d0+700d0/(ng+1)*t
    npos(t,6)=150d0
  end do
end if

    do t = 1, ng    ! now we mark different kindof with different colours
         call cairo_move_to(my_cairo_context, npos(t,1)+20, npos(t,2))
         xx=npos(t,1); yy=npos(t,2)

if(wcl==t)then
            pat = cairo_pattern_create_radial (xx, yy, 5.6d0, xx, yy, 45.0d0)
                call cairo_pattern_add_color_stop_rgba (pat,0d0, 1d0, 0.4d0, 0d0, 1.0d0)
                if(style1.ge.0) call cairo_pattern_add_color_stop_rgba (pat,1d0, 0.97d0, 0.97d0, 0.97d0, 1.0d0)
                if(style1.lt.0) call cairo_pattern_add_color_stop_rgba (pat,1d0, 0d0, 0d0, 0d0, 1.0d0)
                call cairo_set_source (my_cairo_context,pat)
                call cairo_arc(my_cairo_context, xx, yy, 45d0, 0d0, 2*pi)
                call cairo_fill (my_cairo_context)
                call cairo_pattern_destroy (pat)
end if
    end do
    do t=1,ng

       if(gen(t)%kindof==1)then
             rr=0.8d0; vv=0.8d0; bl=0.8d0
             call cairo_set_source_rgb(my_cairo_context, 0.8d0, 0.8d0, 0.8d0) 
       elseif(gen(t)%kindof==2)then
             rr=0.0d0; vv=0.6d0; bl=0.0d0
             call cairo_set_source_rgb(my_cairo_context, 0.0d0, 0.4d0, 0.0d0)
       elseif(gen(t)%kindof==3)then
             rr=0.0d0; vv=0.0d0; bl=0.8d0
             call cairo_set_source_rgb(my_cairo_context, 0.0d0, 0.0d0, 0.7d0)
       elseif(gen(t)%kindof==4)then
             rr=0.8d0; vv=0.0d0; bl=0.0d0
             call cairo_set_source_rgb(my_cairo_context, 0.8d0, 0.0d0, 0.0d0)
       elseif(gen(t)%kindof==5)then
             rr=1.0d0; vv=0.5d0; bl=1.0d0
             call cairo_set_source_rgb(my_cairo_context, 1.0d0, 0.5d0, 1.0d0)
       elseif(gen(t)%kindof==6)then
             rr=0.6d0; vv=0.6d0; bl=1.0d0
             call cairo_set_source_rgb(my_cairo_context, 0.4d0, 0.7d0, 1.0d0)
       elseif(gen(t)%kindof==7)then
             rr=0.2d0; vv=1.0d0; bl=0.2d0
             call cairo_set_source_rgb(my_cairo_context, 0.2d0, 1.0d0, 0.2d0)
       elseif(gen(t)%kindof==8)then
             rr=1.0d0; vv=0.8d0; bl=0.2d0
             call cairo_set_source_rgb(my_cairo_context, 1.0d0, 8.0d0, 0.2d0)
             call cairo_set_source_rgb(my_cairo_context, 1.0d0, 1.0d0, 0.2d0) !,0.2d0z
       endif

! now, we plot the connections

         do t1=1,gen(t)%npost
            t2=gen(t)%post(t1)
if(t.ne.t2)then
            call cairo_set_source_rgb(my_cairo_context, 0.5d0, 0.3d0, 0.0d0)
            tt1=0.3d0+2.0d0*1
            tt3=npos(t,1)-npos(t2,1)
            tt4=npos(t,2)-npos(t2,2)
            tt2=sqrt((tt3)**2+(tt4)**2)
            call cairo_set_line_width(my_cairo_context, tt1)
            call cairo_move_to(my_cairo_context, npos(t,1)+tt3/tt2*20, npos(t,2)+tt4/tt2*20)  
            call cairo_line_to(my_cairo_context, npos(t2,1)+tt3*0.2, npos(t2,2)+tt4*0.2)           
            call cairo_stroke(my_cairo_context) 
            tt1=0.6d0+4.0d0*2
            if(tt1.gt.12) tt1=12
            call cairo_set_line_width(my_cairo_context, tt1)
            tt3=npos(t,1)-npos(t2,1)
            tt4=npos(t,2)-npos(t2,2)
            tt2=sqrt((tt3)**2+(tt4)**2)
            call cairo_move_to(my_cairo_context, npos(t2,1)+tt3/tt2*20, npos(t2,2)+tt4/tt2*20)  
            call cairo_line_to(my_cairo_context, npos(t2,1)+tt3*0.2, npos(t2,2)+tt4*0.2)
            call cairo_stroke(my_cairo_context) 
 end if
         end do

         do t1=1,gen(t)%npre
            t2=gen(t)%pre(t1)
if(t.ne.t2)then
            call cairo_set_source_rgb(my_cairo_context, 0.5d0, 0.3d0, 0.0d0)
            tt1=0.3d0+2.0d0*1
            tt3=npos(t2,1)-npos(t,1)
            tt4=npos(t2,2)-npos(t,2)
            tt2=sqrt((tt3)**2+(tt4)**2)
            call cairo_set_line_width(my_cairo_context, tt1)         
            call cairo_stroke(my_cairo_context) 
            tt1=0.6d0+4.0d0*2
            if(tt1.gt.12) tt1=12
            call cairo_set_line_width(my_cairo_context, tt1)
            tt3=npos(t2,1)-npos(t,1)
            tt4=npos(t2,2)-npos(t,2)
            tt2=sqrt((tt3)**2+(tt4)**2)
            call cairo_stroke(my_cairo_context) 
 end if
         end do

         call cairo_move_to(my_cairo_context, npos(t,1)+20, npos(t,2))
         xx=npos(t,1); yy=npos(t,2)

         pat = cairo_pattern_create_radial (xx-6d0, yy-6d0, 5.6d0, xx-8d0, yy-8d0, 28.0d0)
         call cairo_pattern_add_color_stop_rgba (pat, 0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0)
         call cairo_pattern_add_color_stop_rgba (pat, 1d0, rr, vv, bl, 1.0d0)
         call cairo_set_source (my_cairo_context, pat)
         call cairo_arc(my_cairo_context, xx, yy, 20d0, 0d0, 2*pi)
         call cairo_fill (my_cairo_context)
         call cairo_pattern_destroy (pat)

    
if(style1.ge.0) call cairo_set_source_rgb(my_cairo_context, 0d0, 0d0, 0d0)
if(style1.lt.0) call cairo_set_source_rgb(my_cairo_context, 1d0, 1d0, 1d0)
        call cairo_select_font_face(my_cairo_context, "Times"//c_null_char, CAIRO_FONT_SLANT_NORMAL, &
                                 &  CAIRO_FONT_WEIGHT_NORMAL)
call cairo_set_source_rgb(my_cairo_context, 0d0, 0d0, 0d0)
        call cairo_set_font_size (my_cairo_context, 20d0)
        call cairo_move_to(my_cairo_context, npos(t,1)-15d0, npos(t,2)+5d0)

        if(len(trim(genename(t))).eq.0)then
          write (str3,"(I3.1)") t
          call cairo_show_text (my_cairo_context, str3//c_null_char)
        else
          write (str12,"(A12)") genename(t)
          call cairo_show_text (my_cairo_context, str12//c_null_char)
        end if
    end do
    
    call cairo_destroy(my_cairo_context)
    ret = FALSE
  end function draw_network3




  function draw_network4 (widget, event, gdata) result(ret)  bind(c)
    use iso_c_binding, only: c_int, c_ptr
    use io
    use global_widgets
    use some_widgets

    integer(c_int)    :: ret
    character*140 :: fii
    character*140 :: fij
    character*3 :: str3
    character*12 :: str12
    character*8 :: str8
    type(c_ptr), value, intent(in) :: widget, event, gdata
    type(c_ptr) :: my_cairo_context, pat
    integer :: cstatus
    integer :: t,t1,t2,t3,t4
    real*8 :: tt1,tt2,tt3,tt4,tt5,tt6,xx,yy,rr,vv,bl
    real*8 :: t1x,t1y,t2x,t2y,ta1,ta2,tr1,tr2
    real*8 :: ev=0.000001d0

if(rompetout.ne.3)then
if(yetthere==1)then
call gtk_widget_destroy(my_choice_window2)
yetthere=0
endif
if(yetthere1==1)then
call gtk_widget_destroy(my_choice_window)
yetthere1=0
endif
if(yetthere3==1)then
call gtk_widget_destroy(my_choice_window3)
yetthere3=0
endif
endif
rompetout=3

    my_cairo_context = gdk_cairo_create (gtk_widget_get_window(widget))

if(style1==-1)then
   call cairo_set_source_rgb(my_cairo_context, 0d0, 0d0, 0d0)
   call cairo_rectangle(my_cairo_context,0d0,0d0,800d0,800d0)
   call cairo_fill(my_cairo_context)
end if

    ! Lines:
    call cairo_set_source_rgb(my_cairo_context, 0.5d0, 0.5d0, 0.5d0)
    call cairo_set_line_width(my_cairo_context, 2d0)
    do t = 0, int(height), +20
      call cairo_move_to(my_cairo_context, 0d0, t*1d0)  
      call cairo_line_to(my_cairo_context, t*1d0, height*1d0)
      call cairo_stroke(my_cairo_context) 
    end do
  
    ! Text:
    write (str8,"(I8.1)") getot
    call cairo_set_source_rgb(my_cairo_context, 0d0, 0d0, 1d0)
    call cairo_select_font_face(my_cairo_context, "Times"//c_null_char, CAIRO_FONT_SLANT_NORMAL, &
                                 &  CAIRO_FONT_WEIGHT_NORMAL)
    call cairo_set_font_size (my_cairo_context, 30d0)
    call cairo_move_to(my_cairo_context, 10d0, 30d0)
    call cairo_show_text (my_cairo_context, "This is the Network of catalytic interactions"//c_null_char)
        call cairo_move_to(my_cairo_context, 300d0, 30d0)

    ! Arrows:
    call cairo_new_sub_path(my_cairo_context)

if(togglen==0)then
  if(allocated(npos)) deallocate(npos)
  allocate(npos(ng,6))
  do t=1,ng
    npos(t,1)=350d0+260d0*cos(t*2d0*pi/ng)
    npos(t,2)=350d0+260d0*sin(t*2d0*pi/ng)
    npos(t,3)=350d0+260d0*cos(t*pi/ng)
    npos(t,4)=350d0+260d0*sin(t*pi/ng)
    npos(t,5)=50d0+700d0/(ng+1)*t
    npos(t,6)=150d0
  end do
end if

    do t=1,ng

            call cairo_move_to(my_cairo_context, npos(t,1), npos(t,2))
            xx=npos(t,1); yy=npos(t,2)
if(wcl==t)then
            pat = cairo_pattern_create_radial (xx, yy, 5.6d0, xx, yy, 45.0d0)
if(w0==0) call cairo_pattern_add_color_stop_rgba (pat,0d0, 1d0, 0.8d0, 0d0, 1.0d0)
if(w0==1) call cairo_pattern_add_color_stop_rgba (pat,0d0, 1d0, 0.4d0, 0d0, 1.0d0)
if(w0==2) call cairo_pattern_add_color_stop_rgba (pat,0d0, 0d0, 0.8d0, 0.8d0, 1.0d0)
                if(style1.ge.0) call cairo_pattern_add_color_stop_rgba (pat,1d0, 0.97d0, 0.97d0, 0.97d0, 1.0d0)
                if(style1.lt.0) call cairo_pattern_add_color_stop_rgba (pat,1d0, 0d0, 0d0, 0d0, 1.0d0)
                call cairo_set_source (my_cairo_context,pat)
                call cairo_arc(my_cairo_context, xx, yy, 45d0, 0d0, 2*pi)
                call cairo_fill (my_cairo_context)
                call cairo_pattern_destroy (pat)
end if
end do
  do t=1,ng

       do t1=1,gen(t)%npost
  
            t2=gen(t)%post(t1)
            if(t.ne.t2)then
            if(abs(gen(t)%ww(t2,3)).le.ev) call cairo_set_source_rgb(my_cairo_context, 0.9d0, 0.8d0, 0.7d0)
            tt1=0.3d0+2.0d0*1
              tt1=3
            tt3=npos(t,1)-npos(t2,1)
            tt4=npos(t,2)-npos(t2,2)
            tt2=sqrt((tt3)**2+(tt4)**2)
            call cairo_set_line_width(my_cairo_context, tt1)
            call cairo_move_to(my_cairo_context, npos(t,1)+tt3/tt2*20, npos(t,2)+tt4/tt2*20)  
            call cairo_line_to(my_cairo_context, npos(t2,1)+tt3*0.2, npos(t2,2)+tt4*0.2)           
            call cairo_stroke(my_cairo_context) 
            tt1=0.6d0+4.0d0*2
            if(tt1.gt.12) tt1=12
              tt1=3
            call cairo_set_line_width(my_cairo_context, tt1)
            tt3=npos(t,1)-npos(t2,1)
            tt4=npos(t,2)-npos(t2,2)
            tt2=sqrt((tt3)**2+(tt4)**2)
            call cairo_move_to(my_cairo_context, npos(t2,1)+tt3/tt2*20, npos(t2,2)+tt4/tt2*20)  
            call cairo_line_to(my_cairo_context, npos(t2,1)+tt3*0.2, npos(t2,2)+tt4*0.2)
            call cairo_stroke(my_cairo_context) 

              t2x=npos(t,1); t2y=npos(t,2); t1x=npos(t2,1); t1y=npos(t2,2)
              tt3=t2x-t1x; tt4=t2y-t1y
              tt2=sqrt(tt3**2+tt4**2)
              if(tt2.le.70d0)then
                tr1=tt3*0.7; tr2=tt4*0.7
              else
                tr1=30d0/tt2*tt3; tr2=30d0/tt2*tt4
              end if
              call cairo_move_to(my_cairo_context, t1x+tr1, t1y+tr2)
              ta1=t1x+tr1-(t1y-t2y)
              ta2=t1y+tr2+(t1x-t2x)
              tt1=sqrt((t1x-ta1)**2+(t1y-ta2)**2)
              ta1=t1x+tr1-10d0/tt1*(t1y-t2y)
              ta2=t1y+tr2+10d0/tt1*(t1x-t2x)
              call cairo_move_to(my_cairo_context, ta1, ta2)
              call cairo_line_to(my_cairo_context, t1x, t1y)
              ta1=t1x+tr1-(t1y-t2y)
              ta2=t1y+tr2+(t1x-t2x)
              tt1=sqrt((t1x-ta1)**2+(t1y-ta2)**2)
              ta1=t1x+tr1+10d0/tt1*(t1y-t2y)
              ta2=t1y+tr2-10d0/tt1*(t1x-t2x)
              call cairo_line_to(my_cairo_context, ta1, ta2)
              call cairo_close_path(my_cairo_context)
              call cairo_fill(my_cairo_context)
              call cairo_stroke(my_cairo_context)

            else
              print*, "WARNING! WARNING! WARNING! This gene", t, "is it's own pre!"
              print*, "This doesn't make ANY sense."
              print*, "You'd better correct this nonsense."
              print*, "Seriously!"
            end if
         end do

    end do

    do t = 1, ng 

      do t2 = 1, gen(t)%nww
        t3 = gen(t)%ww(t2,1); t4 = gen(t)%ww(t2,2)    
        
        if(gen(t)%ww(t2,3).gt.ev) then   ! since the arrows of catalytic interactions are pointing towards the end of an arrow of pre>post relationship of forms, we need some calculations.
          tt3=npos(t4,1)-npos(t3,1)
          tt4=npos(t4,2)-npos(t3,2)
          tt2=sqrt((tt3)**2+(tt4)**2)
          tt5=npos(t3,1) + tt3*0.8
          tt6=npos(t3,2) + tt4*0.8
          call cairo_set_line_width(my_cairo_context, 3.0d0)!
          call cairo_move_to(my_cairo_context, tt5, tt6)
          call cairo_line_to(my_cairo_context, npos(t4,1), npos(t4,2))
          call cairo_stroke(my_cairo_context) 
          tt5=npos(t3,1) + tt3*0.7
          tt6=npos(t3,2) + tt4*0.7
          tt2=sqrt((tt5-npos(t,1))**2+(tt6-npos(t,1))**2)
          call cairo_set_line_width(my_cairo_context, 3.0d0)
          call cairo_set_source_rgb(my_cairo_context, 0.2d0, 0.8d0, 0.2d0)
          call cairo_move_to(my_cairo_context, npos(t,1), npos(t,2))
          call cairo_line_to(my_cairo_context, tt5, tt6)
          call cairo_stroke(my_cairo_context) 

          if(gen(t)%ww(t2,1)==t)then
            tt3=npos(t,1)+(tt5-npos(t,1))*0.5
            tt4=npos(t,2)+(tt6-npos(t,2))*0.5
            call cairo_line_to(my_cairo_context, tt3, tt4)
            call cairo_stroke(my_cairo_context) 
            tt2=sqrt((tt3-npos(t,1))**2+(tt4-npos(t,2))**2)
            call cairo_set_line_width(my_cairo_context, 3.0d0)
            tt5=npos(t4,1)-npos(t3,1)
            tt6=npos(t4,2)-npos(t3,2)
            tt5=atan(tt6/tt5)
            call cairo_arc (my_cairo_context, tt3, tt4, tt2, tt5, pi+tt5)
            call cairo_stroke(my_cairo_context) 
            call cairo_set_line_width(my_cairo_context, 9.0d0) 
            if((npos(gen(t)%ww(t2,2),1)>npos(gen(t)%ww(t2,1),1)).or.((npos(gen(t)%ww(t2,2),1)==npos(gen(t)%ww(t2,1),1))&
            &.and.(npos(gen(t)%ww(t2,2),2)>npos(gen(t)%ww(t2,1),2))))then
              call cairo_arc (my_cairo_context, tt3, tt4, tt2, tt5, pi*0.1d0+tt5)

            else

              call cairo_arc (my_cairo_context, tt3, tt4, tt2, 0.9d0*pi+tt5, pi+tt5)

            end if
            call cairo_stroke(my_cairo_context)

              t2x=npos(t,1); t2y=npos(t,2); t1x=tt5; t1y=tt6
            if((npos(gen(t)%ww(t2,2),1)>npos(gen(t)%ww(t2,1),1)).or.&
            &((npos(gen(t)%ww(t2,2),1)==npos(gen(t)%ww(t2,1),1))&
            &.and.(npos(gen(t)%ww(t2,2),2)>npos(gen(t)%ww(t2,1),2))))then
              call cairo_arc (my_cairo_context, tt3, tt4, tt2, tt5, pi*0.1d0+tt5)
            else
              call cairo_arc (my_cairo_context, tt3, tt4, tt2, 0.9d0*pi+tt5, pi+tt5)
            end if

              ta1=t1x+tr1-(t1y-t2y)
              ta2=t1y+tr2+(t1x-t2x)
              tt1=sqrt((t1x-ta1)**2+(t1y-ta2)**2)
              ta1=t1x+tr1-10d0/tt1*(t1y-t2y)
              ta2=t1y+tr2+10d0/tt1*(t1x-t2x)
              call cairo_move_to(my_cairo_context, ta1, ta2)
              call cairo_line_to(my_cairo_context, t1x, t1y)
              ta1=t1x+tr1-(t1y-t2y)
              ta2=t1y+tr2+(t1x-t2x)
              tt1=sqrt((t1x-ta1)**2+(t1y-ta2)**2)
              ta1=t1x+tr1+10d0/tt1*(t1y-t2y)
              ta2=t1y+tr2-10d0/tt1*(t1x-t2x)
              call cairo_line_to(my_cairo_context, ta1, ta2)
              call cairo_close_path(my_cairo_context)
              call cairo_fill(my_cairo_context)
              call cairo_stroke(my_cairo_context)

          elseif(gen(t)%ww(t2,2)==t)then

            tt3=npos(t,1)+(tt5-npos(t,1))*0.5
            tt4=npos(t,2)+(tt6-npos(t,2))*0.5
            call cairo_line_to(my_cairo_context, tt3, tt4)
            call cairo_stroke(my_cairo_context) 
            tt2=sqrt((tt3-npos(t,1))**2+(tt4-npos(t,2))**2)
            call cairo_set_line_width(my_cairo_context, 3.0d0)
            tt5=(npos(t4,1)-(npos(t3,1)))
            tt6=(npos(t4,2)-(npos(t3,2)))
            tt5=atan(tt6/tt5)
            call cairo_arc (my_cairo_context, tt3, tt4, tt2, tt5, pi+tt5)
            call cairo_stroke(my_cairo_context)
            call cairo_set_line_width(my_cairo_context, 9.0d0) 
            if((npos(gen(t)%ww(t2,2),1)>npos(gen(t)%ww(t2,1),1)).or.((npos(gen(t)%ww(t2,2),1)==npos(gen(t)%ww(t2,1),1))&
            &.and.(npos(gen(t)%ww(t2,2),2)>npos(gen(t)%ww(t2,1),2))))then
              call cairo_arc (my_cairo_context, tt3, tt4, tt2, 0.9d0*pi+tt5, pi+tt5)
            else
              call cairo_arc (my_cairo_context, tt3, tt4, tt2, tt5, pi*0.1d0+tt5)
            end if
            call cairo_stroke(my_cairo_context) 
          else
            call cairo_set_line_width(my_cairo_context, 3.0d0)!
            call cairo_move_to(my_cairo_context, npos(t,1)+(tt5-npos(t,1))*0.8,&
           & npos(t,2)+(tt6-npos(t,2))*0.8)
            call cairo_line_to(my_cairo_context, tt5, tt6)
            call cairo_stroke(my_cairo_context) 

              t2x=npos(t,1); t2y=npos(t,2); t1x=tt5; t1y=tt6
              tt3=t2x-t1x; tt4=t2y-t1y
              tt2=sqrt(tt3**2+tt4**2)
              if(tt2.le.70d0)then
                tr1=tt3*0.7; tr2=tt4*0.7
              else
                tr1=30d0/tt2*tt3; tr2=30d0/tt2*tt4
              end if
              call cairo_move_to(my_cairo_context, t1x+tr1, t1y+tr2)
              ta1=t1x+tr1-(t1y-t2y)
              ta2=t1y+tr2+(t1x-t2x)
              tt1=sqrt((t1x-ta1)**2+(t1y-ta2)**2)
              ta1=t1x+tr1-10d0/tt1*(t1y-t2y)
              ta2=t1y+tr2+10d0/tt1*(t1x-t2x)
              call cairo_move_to(my_cairo_context, ta1, ta2)
              call cairo_line_to(my_cairo_context, t1x, t1y)
              ta1=t1x+tr1-(t1y-t2y)
              ta2=t1y+tr2+(t1x-t2x)
              tt1=sqrt((t1x-ta1)**2+(t1y-ta2)**2)
              ta1=t1x+tr1+10d0/tt1*(t1y-t2y)
              ta2=t1y+tr2-10d0/tt1*(t1x-t2x)
              call cairo_line_to(my_cairo_context, ta1, ta2)
              call cairo_close_path(my_cairo_context)
              call cairo_fill(my_cairo_context)
              call cairo_stroke(my_cairo_context)

          end if
            call cairo_set_source_rgb(my_cairo_context, 0.5d0, 0.3d0, 0.0d0)
        elseif(gen(t)%ww(t2,3).lt.-ev) then
          tt3=npos(t4,1)-npos(t3,1)
          tt4=npos(t4,2)-npos(t3,2)
          tt2=sqrt((tt3)**2+(tt4)**2)
          tt5=npos(t3,1) + tt3*0.8
          tt6=npos(t3,2) + tt4*0.8
          call cairo_set_line_width(my_cairo_context, 3.0d0)!
          call cairo_move_to(my_cairo_context, tt5, tt6)
          call cairo_line_to(my_cairo_context, npos(t4,1), npos(t4,2))
          call cairo_stroke(my_cairo_context) 
          tt5=npos(t3,1) + tt3*0.7
          tt6=npos(t3,2) + tt4*0.7
          tt2=sqrt((tt5-npos(t,1))**2+(tt6-npos(t,1))**2)
          call cairo_set_line_width(my_cairo_context, 3.0d0)
          call cairo_set_source_rgb(my_cairo_context, 0.9d0, 0.3d0, 0.3d0)
          call cairo_move_to(my_cairo_context, npos(t,1), npos(t,2))
          call cairo_line_to(my_cairo_context, tt5, tt6)
          call cairo_stroke(my_cairo_context)


       if(gen(t)%ww(t2,1)==t)then
          tt3=npos(t,1)+(tt5-npos(t,1))*0.5
          tt4=npos(t,2)+(tt6-npos(t,2))*0.5
          call cairo_line_to(my_cairo_context, tt3, tt4)
          call cairo_stroke(my_cairo_context) 
          tt2=sqrt((tt3-npos(t,1))**2+(tt4-npos(t,2))**2)
            call cairo_set_line_width(my_cairo_context, 3.0d0)
          tt5=npos(t4,1)-npos(t3,1)
          tt6=npos(t4,2)-npos(t3,2)
            tt5=atan(tt6/tt5)
            call cairo_arc (my_cairo_context, tt3, tt4, tt2, tt5, pi+tt5)
            call cairo_stroke(my_cairo_context)
            call cairo_set_line_width(my_cairo_context, 9.0d0) 
            if((npos(gen(t)%ww(t2,2),1)>npos(gen(t)%ww(t2,1),1)).or.((npos(gen(t)%ww(t2,2),1)==npos(gen(t)%ww(t2,1),1))&
            &.and.(npos(gen(t)%ww(t2,2),2)>npos(gen(t)%ww(t2,1),2))))then 
              call cairo_arc (my_cairo_context, tt3, tt4, tt2, tt5, pi*0.1d0+tt5)
            else
              call cairo_arc (my_cairo_context, tt3, tt4, tt2, 0.9d0*pi+tt5, pi+tt5)
            end if
            call cairo_stroke(my_cairo_context)
        elseif(gen(t)%ww(t2,2)==t)then
          tt3=npos(t,1)+(tt5-npos(t,1))*0.5
          tt4=npos(t,2)+(tt6-npos(t,2))*0.5
          call cairo_line_to(my_cairo_context, tt3, tt4)
          call cairo_stroke(my_cairo_context) 
          tt2=sqrt((tt3-npos(t,1))**2+(tt4-npos(t,2))**2)
            call cairo_set_line_width(my_cairo_context, 3.0d0)
          tt5=npos(t4,1)-npos(t3,1)
          tt6=npos(t4,2)-npos(t3,2)
            tt5=atan(tt6/tt5)
            call cairo_arc (my_cairo_context, tt3, tt4, tt2, tt5, pi+tt5)
            call cairo_stroke(my_cairo_context) 
            call cairo_set_line_width(my_cairo_context, 9.0d0) 
            if((npos(gen(t)%ww(t2,2),1)>npos(gen(t)%ww(t2,1),1)).or.((npos(gen(t)%ww(t2,2),1)==npos(gen(t)%ww(t2,1),1))&
            &.and.(npos(gen(t)%ww(t2,2),2)>npos(gen(t)%ww(t2,1),2))))then 
              call cairo_arc (my_cairo_context, tt3, tt4, tt2, 0.9d0*pi+tt5, pi+tt5)
            else
              call cairo_arc (my_cairo_context, tt3, tt4, tt2, tt5, pi*0.1d0+tt5)
            end if
            call cairo_stroke(my_cairo_context)

        else
          call cairo_set_line_width(my_cairo_context, 3.0d0)!
          call cairo_move_to(my_cairo_context, npos(t,1)+(tt5-npos(t,1))*0.8,&
         & npos(t,2)+(tt6-npos(t,2))*0.8)
          call cairo_line_to(my_cairo_context, tt5, tt6)
          call cairo_stroke(my_cairo_context)


              t2x=npos(t,1); t2y=npos(t,2); t1x=tt5; t1y=tt6
              tt3=t2x-t1x; tt4=t2y-t1y
              tt2=sqrt(tt3**2+tt4**2)
              if(tt2.le.70d0)then
                tr1=tt3*0.7; tr2=tt4*0.7
              else
                tr1=30d0/tt2*tt3; tr2=30d0/tt2*tt4
              end if
              call cairo_move_to(my_cairo_context, t1x+tr1, t1y+tr2)
              ta1=t1x+tr1-(t1y-t2y)
              ta2=t1y+tr2+(t1x-t2x)
              tt1=sqrt((t1x-ta1)**2+(t1y-ta2)**2)
              ta1=t1x+tr1-10d0/tt1*(t1y-t2y)
              ta2=t1y+tr2+10d0/tt1*(t1x-t2x)
              call cairo_move_to(my_cairo_context, ta1, ta2)
              call cairo_line_to(my_cairo_context, t1x, t1y)
              ta1=t1x+tr1-(t1y-t2y)
              ta2=t1y+tr2+(t1x-t2x)
              tt1=sqrt((t1x-ta1)**2+(t1y-ta2)**2)
              ta1=t1x+tr1+10d0/tt1*(t1y-t2y)
              ta2=t1y+tr2-10d0/tt1*(t1x-t2x)
              call cairo_line_to(my_cairo_context, ta1, ta2)
              call cairo_close_path(my_cairo_context)
              call cairo_fill(my_cairo_context)
              call cairo_stroke(my_cairo_context)
 
        end if

        call cairo_set_source_rgb(my_cairo_context, 0.5d0, 0.3d0, 0.0d0)
        end if
        
        if((abs(gen(t)%ww(t2,3)).gt.ev).and.(t.ne.gen(t)%ww(t2,1)).and.(t.ne.gen(t)%ww(t2,2)))then
          call cairo_set_source_rgb(my_cairo_context, 0.5d0, 0.3d0, 0.0d0)
        end if
       
        tt1=0.3d0+2.0d0*sqrt(gen(t)%ww(t2,3))
        if(tt1.gt.10) tt1=10
        call cairo_set_line_width(my_cairo_context, 3.0d0)
        call cairo_move_to(my_cairo_context, npos(t3,1), npos(t3,2))  
        call cairo_line_to(my_cairo_context, npos(t4,1), npos(t4,2))  
        call cairo_stroke(my_cairo_context)  

              t2x=npos(t3,1); t2y=npos(t3,2); t1x=npos(t4,1); t1y=npos(t4,2)
              tt3=t2x-t1x; tt4=t2y-t1y
              tt2=sqrt(tt3**2+tt4**2)
              if(tt2.le.70d0)then
                tr1=tt3*0.7; tr2=tt4*0.7
              else
                tr1=30d0/tt2*tt3; tr2=30d0/tt2*tt4
              end if
              call cairo_move_to(my_cairo_context, t1x+tr1, t1y+tr2)
              ta1=t1x+tr1-(t1y-t2y)
              ta2=t1y+tr2+(t1x-t2x)
              tt1=sqrt((t1x-ta1)**2+(t1y-ta2)**2)
              ta1=t1x+tr1-10d0/tt1*(t1y-t2y)
              ta2=t1y+tr2+10d0/tt1*(t1x-t2x)
              call cairo_move_to(my_cairo_context, ta1, ta2)
              call cairo_line_to(my_cairo_context, t1x, t1y)
              ta1=t1x+tr1-(t1y-t2y)
              ta2=t1y+tr2+(t1x-t2x)
              tt1=sqrt((t1x-ta1)**2+(t1y-ta2)**2)
              ta1=t1x+tr1+10d0/tt1*(t1y-t2y)
              ta2=t1y+tr2-10d0/tt1*(t1x-t2x)
              call cairo_line_to(my_cairo_context, ta1, ta2)
              call cairo_close_path(my_cairo_context)
              call cairo_fill(my_cairo_context)
              call cairo_stroke(my_cairo_context) 

      end do

            call cairo_move_to(my_cairo_context, npos(t,1), npos(t,2))
            xx=npos(t,1); yy=npos(t,2)

       if(gen(t)%kindof==1)then
             rr=0.8d0; vv=0.8d0; bl=0.8d0
             call cairo_set_source_rgb(my_cairo_context, 0.8d0, 0.8d0, 0.8d0) 
       elseif(gen(t)%kindof==2)then
             rr=0.0d0; vv=0.5d0; bl=0.0d0
             call cairo_set_source_rgb(my_cairo_context, 0.0d0, 0.4d0, 0.0d0)
       elseif(gen(t)%kindof==3)then
             rr=0.0d0; vv=0.0d0; bl=0.8d0
             call cairo_set_source_rgb(my_cairo_context, 0.0d0, 0.0d0, 0.7d0)
       elseif(gen(t)%kindof==4)then
             rr=0.8d0; vv=0.0d0; bl=0.0d0
             call cairo_set_source_rgb(my_cairo_context, 0.8d0, 0.0d0, 0.0d0)
       elseif(gen(t)%kindof==5)then
             rr=1.0d0; vv=0.5d0; bl=1.0d0
             call cairo_set_source_rgb(my_cairo_context, 1.0d0, 0.5d0, 1.0d0)
       elseif(gen(t)%kindof==6)then
             rr=0.6d0; vv=0.6d0; bl=1.0d0
             call cairo_set_source_rgb(my_cairo_context, 0.4d0, 0.7d0, 1.0d0)
       elseif(gen(t)%kindof==7)then
             rr=0.2d0; vv=1.0d0; bl=0.2d0
             call cairo_set_source_rgb(my_cairo_context, 0.2d0, 1.0d0, 0.2d0)
       elseif(gen(t)%kindof==8)then
             rr=1.0d0; vv=0.8d0; bl=0.2d0
             call cairo_set_source_rgb(my_cairo_context, 1.0d0, 8.0d0, 0.2d0)
             call cairo_set_source_rgb(my_cairo_context, 1.0d0, 1.0d0, 0.2d0) !,0.2d0z
       endif

            pat = cairo_pattern_create_radial (xx-6d0, yy-6d0, 5.6d0, xx-8d0, yy-8d0, 28.0d0)
                call cairo_pattern_add_color_stop_rgba (pat, 0d0, 1d0, 1d0, 1d0, 1.0d0)
                call cairo_pattern_add_color_stop_rgba (pat, 1d0, rr, vv, bl, 1.0d0)
               
                call cairo_set_source (my_cairo_context, pat)
                call cairo_arc(my_cairo_context, xx, yy, 20d0, 0d0, 2*pi)
                call cairo_fill (my_cairo_context)
                call cairo_pattern_destroy (pat)

if(style1.ge.0) call cairo_set_source_rgb(my_cairo_context, 0d0, 0d0, 0d0)
if(style1.lt.0) call cairo_set_source_rgb(my_cairo_context, 1d0, 1d0, 1d0)
      call cairo_set_line_width(my_cairo_context, 3.0d0)

      call cairo_select_font_face(my_cairo_context, "Times"//c_null_char, CAIRO_FONT_SLANT_NORMAL, &
                                 &  CAIRO_FONT_WEIGHT_NORMAL)
call cairo_set_source_rgb(my_cairo_context, 0d0, 0d0, 0d0)
      call cairo_set_font_size (my_cairo_context, 20d0)
      call cairo_move_to(my_cairo_context, npos(t,1)-15d0, npos(t,2)+5d0)

        if(len(trim(genename(t))).eq.0)then
          write (str3,"(I3.1)") t
          call cairo_show_text (my_cairo_context, str3//c_null_char)
        else
          write (str12,"(A12)") genename(t)
          call cairo_show_text (my_cairo_context, str12//c_null_char)
        end if
  end do
    
    ! Legend

    lx=540;ly=-620
    call cairo_new_sub_path(my_cairo_context)
if(style1.ge.0) call cairo_set_source_rgb(my_cairo_context, 0d0, 0d0, 0d0)
if(style1.lt.0) call cairo_set_source_rgb(my_cairo_context, 1d0, 1d0, 1d0)
    call cairo_select_font_face(my_cairo_context, "Times"//c_null_char, CAIRO_FONT_SLANT_NORMAL, &
                                 &  CAIRO_FONT_WEIGHT_NORMAL)
    call cairo_set_font_size (my_cairo_context, 18d0)
    call cairo_move_to(my_cairo_context, 10d0+lx, 640d0+ly)
    call cairo_show_text (my_cairo_context, "Catalytic interactions"//c_null_char)

    call cairo_new_sub_path(my_cairo_context)
    call cairo_move_to(my_cairo_context, 70d0, 660d0)
            call cairo_set_source_rgb(my_cairo_context, 0.2d0, 0.8d0, 0.2d0)
    call cairo_set_line_width(my_cairo_context, 6.0d0)
    call cairo_move_to(my_cairo_context, 10.0d0+lx, 660.0d0+ly)  
    call cairo_line_to(my_cairo_context, 60.0d0+lx, 660.0d0+ly)    
    call cairo_stroke(my_cairo_context)  
if(style1.ge.0) call cairo_set_source_rgb(my_cairo_context, 0d0, 0d0, 0d0)
if(style1.lt.0) call cairo_set_source_rgb(my_cairo_context, 1d0, 1d0, 1d0)
    call cairo_select_font_face(my_cairo_context, "Times"//c_null_char, CAIRO_FONT_SLANT_NORMAL, &
                                 &  CAIRO_FONT_WEIGHT_NORMAL)
    call cairo_set_font_size (my_cairo_context, 14d0)
    call cairo_move_to(my_cairo_context, 70d0+lx, 663d0+ly)
    call cairo_show_text (my_cairo_context, "Upregulation"//c_null_char)

    call cairo_new_sub_path(my_cairo_context)
    call cairo_move_to(my_cairo_context, 70d0+lx, 680d0+ly)
            call cairo_set_source_rgb(my_cairo_context, 1d0, 0d0, 0d0)
    call cairo_set_line_width(my_cairo_context, 6.0d0)
    call cairo_move_to(my_cairo_context, 10.0d0+lx, 680.0d0+ly)  
    call cairo_line_to(my_cairo_context, 60.0d0+lx, 680.0d0+ly)    
    call cairo_stroke(my_cairo_context)
if(style1.ge.0) call cairo_set_source_rgb(my_cairo_context, 0d0, 0d0, 0d0)
if(style1.lt.0) call cairo_set_source_rgb(my_cairo_context, 1d0, 1d0, 1d0)
    call cairo_select_font_face(my_cairo_context, "Times"//c_null_char, CAIRO_FONT_SLANT_NORMAL, &
                                 &  CAIRO_FONT_WEIGHT_NORMAL)
    call cairo_set_font_size (my_cairo_context, 14d0)
    call cairo_move_to(my_cairo_context, 70d0+lx, 683d0+ly)
    call cairo_show_text (my_cairo_context, "Downregulation"//c_null_char)

    call cairo_set_font_size (my_cairo_context, 18d0)
    call cairo_move_to(my_cairo_context, 70d0+lx, 713d0+ly)

    if(w0==-1) call cairo_show_text (my_cairo_context, "CHOOSE A SUBSTRATE."//c_null_char)
    if(w0==0) call cairo_show_text (my_cairo_context, "CHOOSE A PRODUCT."//c_null_char)
    if(w0==1) call cairo_show_text (my_cairo_context, "CHOOSE A CATALYST."//c_null_char)
    if(w0==2) call cairo_show_text (my_cairo_context, "CHOOSE A SUBSTRATE."//c_null_char)

    ! Legend:

    call cairo_set_source_rgb(my_cairo_context, 0d0, 0d0, 1d0)

    call cairo_set_font_size (my_cairo_context, 18d0)

if(style1.ge.0) call cairo_set_source_rgb(my_cairo_context, 0d0, 0d0, 0d0)
if(style1.lt.0) call cairo_set_source_rgb(my_cairo_context, 1d0, 1d0, 1d0)
    call cairo_move_to(my_cairo_context, 10d0, 620d0)
    call cairo_show_text (my_cairo_context, "Types of Molecules:"//c_null_char)
             call cairo_set_source_rgb(my_cairo_context, 0.8d0, 0.8d0, 0.8d0) 
    call cairo_move_to(my_cairo_context, 10d0, 650d0)

         xx=15d0; yy=645d0
         pat = cairo_pattern_create_radial (xx-3d0, yy-3d0, 2.3d0, xx-4d0, yy-4d0, 14.0d0)
         call cairo_pattern_add_color_stop_rgba (pat, 0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0)
         call cairo_pattern_add_color_stop_rgba (pat, 1d0, 0.8d0, 0.8d0, 0.8d0, 1.0d0)
         call cairo_set_source (my_cairo_context, pat)
         call cairo_arc(my_cairo_context, xx, yy, 10d0, 0d0, 2*pi)
         call cairo_fill (my_cairo_context)
         call cairo_pattern_destroy (pat)

if(style1.ge.0) call cairo_set_source_rgb(my_cairo_context, 0d0, 0d0, 0d0)
if(style1.lt.0) call cairo_set_source_rgb(my_cairo_context, 1d0, 1d0, 1d0)
    call cairo_move_to(my_cairo_context, 10d0, 650d0)
    call cairo_show_text (my_cairo_context, "1 : Regulatory Molecule "//c_null_char)
             call cairo_set_source_rgb(my_cairo_context, 0.0d0, 0.6d0, 0.0d0)
    call cairo_move_to(my_cairo_context, 10d0, 670d0)

         xx=15d0; yy=665d0
         pat = cairo_pattern_create_radial (xx-3d0, yy-3d0, 2.3d0, xx-4d0, yy-4d0, 14.0d0)
         call cairo_pattern_add_color_stop_rgba (pat, 0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0)
         call cairo_pattern_add_color_stop_rgba (pat, 1d0, 0.0d0, 0.6d0, 0.0d0, 1.0d0)
         call cairo_set_source (my_cairo_context, pat)
         call cairo_arc(my_cairo_context, xx, yy, 10d0, 0d0, 2*pi)
         call cairo_fill (my_cairo_context)
         call cairo_pattern_destroy (pat)

    call cairo_set_source_rgb(my_cairo_context, 0.0d0, 0.6d0, 0.0d0)
    call cairo_move_to(my_cairo_context, 10d0, 670d0)
    call cairo_show_text (my_cairo_context, "2 : Regulatory Molecule, translatable"//c_null_char)
             call cairo_set_source_rgb(my_cairo_context, 0.0d0, 0.0d0, 0.8d0)
    call cairo_move_to(my_cairo_context, 10d0, 690d0)

         xx=15d0; yy=685d0
         pat = cairo_pattern_create_radial (xx-3d0, yy-3d0, 2.3d0, xx-4d0, yy-4d0, 14.0d0)
         call cairo_pattern_add_color_stop_rgba (pat, 0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0)
         call cairo_pattern_add_color_stop_rgba (pat, 1d0, 0.0d0, 0.0d0, 0.8d0, 1.0d0)
         call cairo_set_source (my_cairo_context, pat)
         call cairo_arc(my_cairo_context, xx, yy, 10d0, 0d0, 2*pi)
         call cairo_fill (my_cairo_context)
         call cairo_pattern_destroy (pat)

    call cairo_set_source_rgb(my_cairo_context, 0.0d0, 0.0d0, 0.8d0)
    call cairo_move_to(my_cairo_context, 10d0, 690d0)
    call cairo_show_text (my_cairo_context, "3 : Protein, intracellular"//c_null_char)
             call cairo_set_source_rgb(my_cairo_context, 0.8d0, 0.0d0, 0.0d0)

         xx=15d0; yy=705d0
         pat = cairo_pattern_create_radial (xx-3d0, yy-3d0, 2.3d0, xx-4d0, yy-4d0, 14.0d0)
         call cairo_pattern_add_color_stop_rgba (pat, 0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0)
         call cairo_pattern_add_color_stop_rgba (pat, 1d0, 0.8d0, 0.0d0, 0.0d0, 1.0d0)
         call cairo_set_source (my_cairo_context, pat)
         call cairo_arc(my_cairo_context, xx, yy, 10d0, 0d0, 2*pi)
         call cairo_fill (my_cairo_context)
         call cairo_pattern_destroy (pat)

    call cairo_set_source_rgb(my_cairo_context, 0.8d0, 0.0d0, 0.0d0)
    call cairo_move_to(my_cairo_context, 10d0, 710d0)
    call cairo_show_text (my_cairo_context, "4 : Protein, extracellular"//c_null_char)
             call cairo_set_source_rgb(my_cairo_context, 1.0d0, 0.5d0, 1.0d0)

         xx=445d0; yy=645d0
         pat = cairo_pattern_create_radial (xx-3d0, yy-3d0, 2.3d0, xx-4d0, yy-4d0, 14.0d0)
         call cairo_pattern_add_color_stop_rgba (pat, 0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0)
         call cairo_pattern_add_color_stop_rgba (pat, 1d0, 1.0d0, 0.5d0, 1.0d0, 1.0d0)
         call cairo_set_source (my_cairo_context, pat)
         call cairo_arc(my_cairo_context, xx, yy, 10d0, 0d0, 2*pi)
         call cairo_fill (my_cairo_context)
         call cairo_pattern_destroy (pat)
!create connection
    call cairo_set_source_rgb(my_cairo_context, 1.0d0, 0.5d0, 1.0d0)
    call cairo_move_to(my_cairo_context, 440d0, 650d0)
    call cairo_show_text (my_cairo_context, "5 : Protein, apically localizing"//c_null_char)
             call cairo_set_source_rgb(my_cairo_context, 0.6d0, 0.6d0, 1.0d0)

         xx=445d0; yy=665d0
         pat = cairo_pattern_create_radial (xx-3d0, yy-3d0, 2.3d0, xx-4d0, yy-4d0, 14.0d0)
         call cairo_pattern_add_color_stop_rgba (pat, 0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0)
         call cairo_pattern_add_color_stop_rgba (pat, 1d0, 0.6d0, 0.6d0, 1.0d0, 1.0d0)
         call cairo_set_source (my_cairo_context, pat)
         call cairo_arc(my_cairo_context, xx, yy, 10d0, 0d0, 2*pi)
         call cairo_fill (my_cairo_context)
         call cairo_pattern_destroy (pat)

    call cairo_set_source_rgb(my_cairo_context, 0.6d0, 0.6d0, 1.0d0)
    call cairo_move_to(my_cairo_context, 440d0, 670d0)
    call cairo_show_text (my_cairo_context, "6 : Protein, basally localizing"//c_null_char)
             call cairo_set_source_rgb(my_cairo_context, 0.2d0, 1.0d0, 0.2d0)

         xx=445d0; yy=685d0
         pat = cairo_pattern_create_radial (xx-3d0, yy-3d0, 2.3d0, xx-4d0, yy-4d0, 14.0d0)
         call cairo_pattern_add_color_stop_rgba (pat, 0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0)
         call cairo_pattern_add_color_stop_rgba (pat, 1d0, 0.2d0, 1.0d0, 0.2d0, 1.0d0)
         call cairo_set_source (my_cairo_context, pat)
         call cairo_arc(my_cairo_context, xx, yy, 10d0, 0d0, 2*pi)
         call cairo_fill (my_cairo_context)
         call cairo_pattern_destroy (pat)

    call cairo_set_source_rgb(my_cairo_context, 0.2d0, 0.9d0, 0.2d0)
    call cairo_move_to(my_cairo_context, 440d0, 690d0)
    call cairo_show_text (my_cairo_context, "7 : Receptor, juxtacrine signalling(7)"//c_null_char)
             call cairo_set_source_rgb(my_cairo_context, 0.8d0, 0.8d0, 0.2d0)

         xx=445d0; yy=705d0
         pat = cairo_pattern_create_radial (xx-3d0, yy-3d0, 2.3d0, xx-4d0, yy-4d0, 14.0d0)
         call cairo_pattern_add_color_stop_rgba (pat, 0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0)
         call cairo_pattern_add_color_stop_rgba (pat, 1d0, 1.0d0, 0.8d0, 0.0d0, 1.0d0)
         call cairo_set_source (my_cairo_context, pat)
         call cairo_arc(my_cairo_context, xx, yy, 10d0, 0d0, 2*pi)
         call cairo_fill (my_cairo_context)
         call cairo_pattern_destroy (pat)

    call cairo_set_source_rgb(my_cairo_context, 1.0d0, 0.8d0, 0.0d0)
    call cairo_move_to(my_cairo_context, 440d0, 710d0)
    call cairo_show_text (my_cairo_context, "8 : Receptor, paracrine signalling(4)"//c_null_char)
    
    call cairo_destroy(my_cairo_context)
    ret = FALSE
  end function draw_network4



  function draw_network5 (widget, event, gdata) result(ret)  bind(c)
    use iso_c_binding, only: c_int, c_ptr
    use io
    use global_widgets

    integer(c_int)    :: ret
    character*140 :: fii
    character*140 :: fij
    character*3 :: str3
    character*12 :: str12
    character*8 :: str8
    type(c_ptr), value, intent(in) :: widget, event, gdata
    type(c_ptr) :: my_cairo_context, pat
    integer :: cstatus
    integer :: t,t1,tot
    real*8 :: tt1,tt2,tt3,tt4,xx,yy,rr,vv,bl

if(rompetout.ne.4)then
if(yetthere==1)then
call gtk_widget_destroy(my_choice_window2)
yetthere=0
endif
if(yetthere1==1)then
call gtk_widget_destroy(my_choice_window)
yetthere1=0
endif
if(yetthere3==1)then
call gtk_widget_destroy(my_choice_window3)
yetthere3=0
endif; endif
rompetout=4

    my_cairo_context = gdk_cairo_create (gtk_widget_get_window(widget))

    tot=nparam_per_node

if(style1==-1)then
   call cairo_set_source_rgb(my_cairo_context, 0d0, 0d0, 0d0)
   call cairo_rectangle(my_cairo_context,0d0,0d0,800d0,800d0)
   call cairo_fill(my_cairo_context)
end if

    ! Lines:
    call cairo_set_source_rgb(my_cairo_context, 0.9d0, 0.6d0, 0.3d0)
    call cairo_set_line_width(my_cairo_context, 2d0)
    do t = 0, int(height), +20
      call cairo_move_to(my_cairo_context, 0d0, t*1d0)  
      call cairo_line_to(my_cairo_context, t*1d0, height*1d0)
      call cairo_stroke(my_cairo_context) 
    end do
  
    ! Text:
    write (str8,"(I8.1)") getot
    call cairo_set_source_rgb(my_cairo_context, 0d0, 0d0, 1d0)
    call cairo_select_font_face(my_cairo_context, "Times"//c_null_char, CAIRO_FONT_SLANT_NORMAL, &
                                 &  CAIRO_FONT_WEIGHT_NORMAL)
    call cairo_set_font_size (my_cairo_context, 30d0)
    call cairo_move_to(my_cairo_context, 10d0, 30d0)
    call cairo_show_text (my_cairo_context, "GRN and Cell behaviours"//c_null_char)
        call cairo_move_to(my_cairo_context, 350d0, 30d0)

if(togglen==0)then
  if(allocated(npos)) deallocate(npos)
  allocate(npos(ng,6))
  do t=1,ng
    npos(t,1)=350d0+260d0*cos(t*2d0*pi/ng)
    npos(t,2)=350d0+260d0*sin(t*2d0*pi/ng)
    npos(t,3)=350d0+260d0*cos(t*pi/ng)
    npos(t,4)=350d0+260d0*sin(t*pi/ng)
    npos(t,5)=50d0+700d0/(ng+1)*t
    npos(t,6)=150d0
  end do
end if
do t = 1, ng
         xx=npos(t,3); yy=npos(t,4)
if(wcl==t)then
            pat = cairo_pattern_create_radial (xx, yy, 5.6d0, xx, yy, 45.0d0)
                call cairo_pattern_add_color_stop_rgba (pat,0d0, 1d0, 0.8d0, 0d0, 1.0d0)
                if(style1.ge.0) call cairo_pattern_add_color_stop_rgba (pat,1d0, 0.97d0, 0.97d0, 0.97d0, 1.0d0)
                if(style1.lt.0) call cairo_pattern_add_color_stop_rgba (pat,1d0, 0d0, 0d0, 0d0, 1.0d0)
                call cairo_set_source (my_cairo_context,pat)
                call cairo_arc(my_cairo_context, xx, yy, 45d0, 0d0, 2*pi)
                call cairo_fill (my_cairo_context)
                call cairo_pattern_destroy (pat)
end if
end do

    !Arrows:
    call cairo_new_sub_path(my_cairo_context)
    do t = 1, ng
       do t1=1,16
          tot=nparam_per_node+t1
          if(gen(t)%wa(tot).gt.0) then
            call cairo_set_source_rgb(my_cairo_context, 0.2d0, 0.8d0, 0.2d0)
            tt1=0.3d0+2.0d0*sqrt(gen(t)%wa(tot))
            if(tt1.gt.10) tt1=10
            call cairo_set_line_width(my_cairo_context, tt1)
            call cairo_move_to(my_cairo_context, npos(t,3), npos(t,4))   
            call cairo_line_to(my_cairo_context, 80d0+t1*500/16, 250d0)         
            call cairo_stroke(my_cairo_context) 
          else if(gen(t)%wa(tot).lt.0) then
            call cairo_set_source_rgb(my_cairo_context, 1d0, 0d0, 0d0)
              tt1=0.3d0+2.0d0*sqrt(gen(t)%wa(tot)*-1)
              if(tt1.gt.10) tt1=10
            call cairo_set_line_width(my_cairo_context, tt1)
            call cairo_move_to(my_cairo_context, npos(t,3), npos(t,4))     
            call cairo_line_to(my_cairo_context, 80d0+t1*500/16, 250d0) 

            call cairo_stroke(my_cairo_context) 

         end if

      end do

         xx=npos(t,3); yy=npos(t,4)


       if(gen(t)%kindof==1)then
             rr=0.8d0; vv=0.8d0; bl=0.8d0
             call cairo_set_source_rgb(my_cairo_context, 0.8d0, 0.8d0, 0.8d0) 
       elseif(gen(t)%kindof==2)then
             rr=0.0d0; vv=0.5d0; bl=0.0d0
             call cairo_set_source_rgb(my_cairo_context, 0.0d0, 0.4d0, 0.0d0)
       elseif(gen(t)%kindof==3)then
             rr=0.0d0; vv=0.0d0; bl=0.8d0
             call cairo_set_source_rgb(my_cairo_context, 0.0d0, 0.0d0, 0.7d0)
       elseif(gen(t)%kindof==4)then
             rr=0.8d0; vv=0.0d0; bl=0.0d0
             call cairo_set_source_rgb(my_cairo_context, 0.8d0, 0.0d0, 0.0d0)
       elseif(gen(t)%kindof==5)then
             rr=1.0d0; vv=0.5d0; bl=1.0d0
             call cairo_set_source_rgb(my_cairo_context, 1.0d0, 0.5d0, 1.0d0)
       elseif(gen(t)%kindof==6)then
             rr=0.6d0; vv=0.6d0; bl=1.0d0
             call cairo_set_source_rgb(my_cairo_context, 0.4d0, 0.7d0, 1.0d0)
       elseif(gen(t)%kindof==7)then
             rr=0.2d0; vv=1.0d0; bl=0.2d0
             call cairo_set_source_rgb(my_cairo_context, 0.2d0, 1.0d0, 0.2d0)
       elseif(gen(t)%kindof==8)then
             rr=1.0d0; vv=0.8d0; bl=0.2d0
             call cairo_set_source_rgb(my_cairo_context, 1.0d0, 8.0d0, 0.2d0)
             call cairo_set_source_rgb(my_cairo_context, 1.0d0, 1.0d0, 0.2d0) !,0.2d0z
       endif

         pat = cairo_pattern_create_radial (xx-6d0, yy-6d0, 5.6d0, xx-8d0, yy-8d0, 28.0d0)
         call cairo_pattern_add_color_stop_rgba (pat, 0d0, 1d0, 1d0, 1d0, 1d0)
         call cairo_pattern_add_color_stop_rgba (pat, 1d0, rr, vv, bl, 1d0)
         call cairo_set_source (my_cairo_context, pat)
         call cairo_arc(my_cairo_context, xx, yy, 20d0, 0d0, 2*pi)
         call cairo_fill (my_cairo_context)
         call cairo_pattern_destroy (pat)
          call cairo_set_line_width(my_cairo_context, 3.0d0)
          call cairo_stroke(my_cairo_context)

! Now, write the numbers (Genes)

if(style1.ge.0) call cairo_set_source_rgb(my_cairo_context, 0d0, 0d0, 0d0)
if(style1.lt.0) call cairo_set_source_rgb(my_cairo_context, 1d0, 1d0, 1d0)
      call cairo_set_line_width(my_cairo_context, 3.0d0)
      call cairo_select_font_face(my_cairo_context, "Times"//c_null_char, CAIRO_FONT_SLANT_NORMAL, &
                                 &  CAIRO_FONT_WEIGHT_NORMAL)
call cairo_set_source_rgb(my_cairo_context, 0d0, 0d0, 0d0)
      call cairo_set_font_size (my_cairo_context, 20d0)
      call cairo_move_to(my_cairo_context, npos(t,3)-20d0, npos(t,4))

        if(len(trim(genename(t))).eq.0)then
          write (str3,"(I3.1)") t
          call cairo_show_text (my_cairo_context, str3//c_null_char)
        else
          write (str12,"(A12)") genename(t)
          call cairo_show_text (my_cairo_context, str12//c_null_char)
        end if
    end do

    call cairo_move_to(my_cairo_context, 10d0, 660d0)

    ! Legend:

    call cairo_set_source_rgb(my_cairo_context, 0d0, 0d0, 1d0)

    call cairo_set_font_size (my_cairo_context, 18d0)

if(style1.ge.0) call cairo_set_source_rgb(my_cairo_context, 0d0, 0d0, 0d0)
if(style1.lt.0) call cairo_set_source_rgb(my_cairo_context, 1d0, 1d0, 1d0)
    call cairo_move_to(my_cairo_context, 10d0, 620d0)
    call cairo_show_text (my_cairo_context, "Types of Molecules:"//c_null_char)
             call cairo_set_source_rgb(my_cairo_context, 0.8d0, 0.8d0, 0.8d0) 
    call cairo_move_to(my_cairo_context, 10d0, 650d0)

         xx=15d0; yy=645d0
         pat = cairo_pattern_create_radial (xx-3d0, yy-3d0, 2.3d0, xx-4d0, yy-4d0, 14.0d0)
         call cairo_pattern_add_color_stop_rgba (pat, 0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0)
         call cairo_pattern_add_color_stop_rgba (pat, 1d0, 0.8d0, 0.8d0, 0.8d0, 1.0d0)
         call cairo_set_source (my_cairo_context, pat)
         call cairo_arc(my_cairo_context, xx, yy, 10d0, 0d0, 2*pi)
         call cairo_fill (my_cairo_context)
         call cairo_pattern_destroy (pat)

if(style1.ge.0) call cairo_set_source_rgb(my_cairo_context, 0d0, 0d0, 0d0)
if(style1.lt.0) call cairo_set_source_rgb(my_cairo_context, 1d0, 1d0, 1d0)
    call cairo_move_to(my_cairo_context, 10d0, 650d0)
    call cairo_show_text (my_cairo_context, "1 : Regulatory Molecule "//c_null_char)
             call cairo_set_source_rgb(my_cairo_context, 0.0d0, 0.6d0, 0.0d0)
    call cairo_move_to(my_cairo_context, 10d0, 670d0)

         xx=15d0; yy=665d0
         pat = cairo_pattern_create_radial (xx-3d0, yy-3d0, 2.3d0, xx-4d0, yy-4d0, 14.0d0)
         call cairo_pattern_add_color_stop_rgba (pat, 0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0)
         call cairo_pattern_add_color_stop_rgba (pat, 1d0, 0.0d0, 0.6d0, 0.0d0, 1.0d0)
         call cairo_set_source (my_cairo_context, pat)
         call cairo_arc(my_cairo_context, xx, yy, 10d0, 0d0, 2*pi)
         call cairo_fill (my_cairo_context)
         call cairo_pattern_destroy (pat)

    call cairo_set_source_rgb(my_cairo_context, 0.0d0, 0.6d0, 0.0d0)
    call cairo_move_to(my_cairo_context, 10d0, 670d0)
    call cairo_show_text (my_cairo_context, "2 : Regulatory Molecule, translatable"//c_null_char)
             call cairo_set_source_rgb(my_cairo_context, 0.0d0, 0.0d0, 0.8d0)
    call cairo_move_to(my_cairo_context, 10d0, 690d0)

         xx=15d0; yy=685d0
         pat = cairo_pattern_create_radial (xx-3d0, yy-3d0, 2.3d0, xx-4d0, yy-4d0, 14.0d0)
         call cairo_pattern_add_color_stop_rgba (pat, 0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0)
         call cairo_pattern_add_color_stop_rgba (pat, 1d0, 0.0d0, 0.0d0, 0.8d0, 1.0d0)
         call cairo_set_source (my_cairo_context, pat)
         call cairo_arc(my_cairo_context, xx, yy, 10d0, 0d0, 2*pi)
         call cairo_fill (my_cairo_context)
         call cairo_pattern_destroy (pat)

    call cairo_set_source_rgb(my_cairo_context, 0.0d0, 0.0d0, 0.8d0)
    call cairo_move_to(my_cairo_context, 10d0, 690d0)
    call cairo_show_text (my_cairo_context, "3 : Protein, intracellular"//c_null_char)
             call cairo_set_source_rgb(my_cairo_context, 0.8d0, 0.0d0, 0.0d0)

         xx=15d0; yy=705d0
         pat = cairo_pattern_create_radial (xx-3d0, yy-3d0, 2.3d0, xx-4d0, yy-4d0, 14.0d0)
         call cairo_pattern_add_color_stop_rgba (pat, 0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0)
         call cairo_pattern_add_color_stop_rgba (pat, 1d0, 0.8d0, 0.0d0, 0.0d0, 1.0d0)
         call cairo_set_source (my_cairo_context, pat)
         call cairo_arc(my_cairo_context, xx, yy, 10d0, 0d0, 2*pi)
         call cairo_fill (my_cairo_context)
         call cairo_pattern_destroy (pat)

    call cairo_set_source_rgb(my_cairo_context, 0.8d0, 0.0d0, 0.0d0)
    call cairo_move_to(my_cairo_context, 10d0, 710d0)
    call cairo_show_text (my_cairo_context, "4 : Protein, extracellular"//c_null_char)
             call cairo_set_source_rgb(my_cairo_context, 1.0d0, 0.5d0, 1.0d0)

         xx=445d0; yy=645d0
         pat = cairo_pattern_create_radial (xx-3d0, yy-3d0, 2.3d0, xx-4d0, yy-4d0, 14.0d0)
         call cairo_pattern_add_color_stop_rgba (pat, 0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0)
         call cairo_pattern_add_color_stop_rgba (pat, 1d0, 1.0d0, 0.5d0, 1.0d0, 1.0d0)
         call cairo_set_source (my_cairo_context, pat)
         call cairo_arc(my_cairo_context, xx, yy, 10d0, 0d0, 2*pi)
         call cairo_fill (my_cairo_context)
         call cairo_pattern_destroy (pat)

    call cairo_set_source_rgb(my_cairo_context, 1.0d0, 0.5d0, 1.0d0)
    call cairo_move_to(my_cairo_context, 440d0, 650d0)
    call cairo_show_text (my_cairo_context, "5 : Protein, apically localizing"//c_null_char)
             call cairo_set_source_rgb(my_cairo_context, 0.6d0, 0.6d0, 1.0d0)

         xx=445d0; yy=665d0
         pat = cairo_pattern_create_radial (xx-3d0, yy-3d0, 2.3d0, xx-4d0, yy-4d0, 14.0d0)
         call cairo_pattern_add_color_stop_rgba (pat, 0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0)
         call cairo_pattern_add_color_stop_rgba (pat, 1d0, 0.6d0, 0.6d0, 1.0d0, 1.0d0)
         call cairo_set_source (my_cairo_context, pat)
         call cairo_arc(my_cairo_context, xx, yy, 10d0, 0d0, 2*pi)
         call cairo_fill (my_cairo_context)
         call cairo_pattern_destroy (pat)

    call cairo_set_source_rgb(my_cairo_context, 0.6d0, 0.6d0, 1.0d0)
    call cairo_move_to(my_cairo_context, 440d0, 670d0)
    call cairo_show_text (my_cairo_context, "6 : Protein, basally localizing"//c_null_char)
             call cairo_set_source_rgb(my_cairo_context, 0.2d0, 1.0d0, 0.2d0)

         xx=445d0; yy=685d0
         pat = cairo_pattern_create_radial (xx-3d0, yy-3d0, 2.3d0, xx-4d0, yy-4d0, 14.0d0)
         call cairo_pattern_add_color_stop_rgba (pat, 0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0)
         call cairo_pattern_add_color_stop_rgba (pat, 1d0, 0.2d0, 1.0d0, 0.2d0, 1.0d0)
         call cairo_set_source (my_cairo_context, pat)
         call cairo_arc(my_cairo_context, xx, yy, 10d0, 0d0, 2*pi)
         call cairo_fill (my_cairo_context)
         call cairo_pattern_destroy (pat)

    call cairo_set_source_rgb(my_cairo_context, 0.2d0, 0.9d0, 0.2d0)
    call cairo_move_to(my_cairo_context, 440d0, 690d0)
    call cairo_show_text (my_cairo_context, "7 : Receptor, juxtacrine signalling(7)"//c_null_char)
             call cairo_set_source_rgb(my_cairo_context, 0.8d0, 0.8d0, 0.2d0)

         xx=445d0; yy=705d0
         pat = cairo_pattern_create_radial (xx-3d0, yy-3d0, 2.3d0, xx-4d0, yy-4d0, 14.0d0)
         call cairo_pattern_add_color_stop_rgba (pat, 0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0)
         call cairo_pattern_add_color_stop_rgba (pat, 1d0, 1.0d0, 0.8d0, 0.0d0, 1.0d0)
         call cairo_set_source (my_cairo_context, pat)
         call cairo_arc(my_cairo_context, xx, yy, 10d0, 0d0, 2*pi)
         call cairo_fill (my_cairo_context)
         call cairo_pattern_destroy (pat)

    call cairo_set_source_rgb(my_cairo_context, 1.0d0, 0.8d0, 0.0d0)
    call cairo_move_to(my_cairo_context, 440d0, 710d0)
    call cairo_show_text (my_cairo_context, "8 : Receptor, paracrine signalling(4)"//c_null_char)

! Write the numbers (Cell behaviours)
    do t=1,16
       call cairo_set_source_rgb(my_cairo_context, 0.2d0, 0.6d0, 0.9d0)
        call cairo_set_line_width(my_cairo_context, 3.0d0)
        call cairo_select_font_face(my_cairo_context, "Arial"//c_null_char, CAIRO_FONT_SLANT_NORMAL, &
                                 &  CAIRO_FONT_WEIGHT_NORMAL)
        call cairo_set_font_size (my_cairo_context, 20d0)
        call cairo_move_to(my_cairo_context, 80d0+t*500/16, 240d0)

        call cairo_rotate (my_cairo_context, 67.6d0)

        call cairo_show_text (my_cairo_context, trim(cellparams(t))//c_null_char)
        call cairo_rotate (my_cairo_context, -67.6d0)
    end do

    call cairo_move_to(my_cairo_context, 10d0, 80d0)
    call cairo_show_text (my_cairo_context, "Cell behaviours"//c_null_char)

    call cairo_destroy(my_cairo_context)
    ret = FALSE
  end function draw_network5



  function draw_network7 (widget, event, gdata) result(ret)  bind(c)
    use iso_c_binding, only: c_int, c_ptr
    use io
    use global_widgets

    integer(c_int)    :: ret
    character*140 :: fii
    character*140 :: fij
    character*3 :: str3
    character*12 :: str12
    character*8 :: str8
    type(c_ptr), value, intent(in) :: widget, event, gdata
    type(c_ptr) :: my_cairo_context, pat
    integer :: cstatus
    integer :: t,t1,tot,s
    real*8 :: tt1,tt2,tt3,tt4,xx,yy,rr,vv,bl
    integer :: npars(23)

    my_cairo_context = gdk_cairo_create (gtk_widget_get_window(widget))

if(rompetout.ne.5)then
if(yetthere==1)then
call gtk_widget_destroy(my_choice_window2)
yetthere=0
endif
if(yetthere1==1)then
call gtk_widget_destroy(my_choice_window)
yetthere1=0
endif
if(yetthere3==1)then
call gtk_widget_destroy(my_choice_window3)
yetthere3=0
endif; endif
rompetout=5

 npars=(/5,6,7,8,9,10,11,12,13,14,15,16,20,21,22,23,24,25,26,27,28,34,35/)

    tot=23!nparam_per_node

if(style1==-1)then
   call cairo_set_source_rgb(my_cairo_context, 0d0, 0d0, 0d0)
   call cairo_rectangle(my_cairo_context,0d0,0d0,800d0,800d0)
   call cairo_fill(my_cairo_context)
end if

    ! Lines:
    call cairo_set_source_rgb(my_cairo_context, 0.9d0, 0.6d0, 0.3d0)
    call cairo_set_line_width(my_cairo_context, 2d0)
    do t = 0, int(height), +20
      call cairo_move_to(my_cairo_context, 0d0, t*1d0)  
      call cairo_line_to(my_cairo_context, t*1d0, height*1d0)
      call cairo_stroke(my_cairo_context) 
    end do
  
    ! Text:
    write (str8,"(I8.1)") getot
    call cairo_set_source_rgb(my_cairo_context, 0d0, 0d0, 1d0)
    call cairo_select_font_face(my_cairo_context, "Times"//c_null_char, CAIRO_FONT_SLANT_NORMAL, &
                                 &  CAIRO_FONT_WEIGHT_NORMAL)
    call cairo_set_font_size (my_cairo_context, 30d0)
    call cairo_move_to(my_cairo_context, 10d0, 30d0)
    call cairo_show_text (my_cairo_context, "GRN and Node properties"//c_null_char)
        call cairo_move_to(my_cairo_context, 350d0, 30d0)


if(togglen==0)then
  if(allocated(npos)) deallocate(npos)
  allocate(npos(ng,6))
  do t=1,ng
    npos(t,1)=350d0+260d0*cos(t*2d0*pi/ng)
    npos(t,2)=350d0+260d0*sin(t*2d0*pi/ng)
    npos(t,3)=350d0+260d0*cos(t*pi/ng)
    npos(t,4)=350d0+260d0*sin(t*pi/ng)
    npos(t,5)=50d0+700d0/(ng+1)*t
    npos(t,6)=150d0
  end do
end if

do t=1,ng
         xx=npos(t,3); yy=npos(t,4)
if(wcl==t)then
            pat = cairo_pattern_create_radial (xx, yy, 5.6d0, xx, yy, 45.0d0)
                call cairo_pattern_add_color_stop_rgba (pat,0d0, 1d0, 0.8d0, 0d0, 1.0d0)
                if(style1.ge.0) call cairo_pattern_add_color_stop_rgba (pat,1d0, 0.97d0, 0.97d0, 0.97d0, 1.0d0)
                if(style1.lt.0) call cairo_pattern_add_color_stop_rgba (pat,1d0, 0d0, 0d0, 0d0, 1.0d0)
                call cairo_set_source (my_cairo_context,pat)
                call cairo_arc(my_cairo_context, xx, yy, 45d0, 0d0, 2*pi)
                call cairo_fill (my_cairo_context)
                call cairo_pattern_destroy (pat)
end if
end do

    !Arrows:
    call cairo_new_sub_path(my_cairo_context)

    do t = 1, ng
       do s=1,23
       t1=npars(s)
          if(gen(t)%wa(t1).gt.0) then
            call cairo_set_source_rgb(my_cairo_context, 0.2d0, 0.8d0, 0.2d0)
            tt1=0.3d0+2.0d0*sqrt(gen(t)%wa(t1))
            if(tt1.gt.10) tt1=10
            call cairo_set_line_width(my_cairo_context, tt1)
            call cairo_move_to(my_cairo_context, npos(t,3), npos(t,4))  
            call cairo_line_to(my_cairo_context, 80d0+s*500/20, 200d0)        
            call cairo_stroke(my_cairo_context) 
          else if(gen(t)%wa(t1).lt.0) then
            call cairo_set_source_rgb(my_cairo_context, 1d0, 0d0, 0d0)
              tt1=0.3d0+2.0d0*sqrt(gen(t)%wa(t1)*-1)
              if(tt1.gt.10) tt1=10
            call cairo_set_line_width(my_cairo_context, tt1)
            call cairo_move_to(my_cairo_context, npos(t,3), npos(t,4))  

             call cairo_line_to(my_cairo_context, 80d0+s*500/20, 200d0)    

            call cairo_stroke(my_cairo_context) 

         end if

      end do

         xx=npos(t,3); yy=npos(t,4)

       if(gen(t)%kindof==1)then
             rr=0.8d0; vv=0.8d0; bl=0.8d0
             call cairo_set_source_rgb(my_cairo_context, 0.8d0, 0.8d0, 0.8d0) 
       elseif(gen(t)%kindof==2)then
             rr=0.0d0; vv=0.5d0; bl=0.0d0
             call cairo_set_source_rgb(my_cairo_context, 0.0d0, 0.4d0, 0.0d0)
       elseif(gen(t)%kindof==3)then
             rr=0.0d0; vv=0.0d0; bl=0.8d0
             call cairo_set_source_rgb(my_cairo_context, 0.0d0, 0.0d0, 0.7d0)
       elseif(gen(t)%kindof==4)then
             rr=0.8d0; vv=0.0d0; bl=0.0d0
             call cairo_set_source_rgb(my_cairo_context, 0.8d0, 0.0d0, 0.0d0)
       elseif(gen(t)%kindof==5)then
             rr=1.0d0; vv=0.5d0; bl=1.0d0
             call cairo_set_source_rgb(my_cairo_context, 1.0d0, 0.5d0, 1.0d0)
       elseif(gen(t)%kindof==6)then
             rr=0.6d0; vv=0.6d0; bl=1.0d0
             call cairo_set_source_rgb(my_cairo_context, 0.4d0, 0.7d0, 1.0d0)
       elseif(gen(t)%kindof==7)then
             rr=0.2d0; vv=1.0d0; bl=0.2d0
             call cairo_set_source_rgb(my_cairo_context, 0.2d0, 1.0d0, 0.2d0)
       elseif(gen(t)%kindof==8)then
             rr=1.0d0; vv=0.8d0; bl=0.0d0
             call cairo_set_source_rgb(my_cairo_context, 1.0d0, 8.0d0, 0.0d0)
             call cairo_set_source_rgb(my_cairo_context, 1.0d0, 1.0d0, 0.0d0) !,0.2d0z
       endif

         pat = cairo_pattern_create_radial (xx-6d0, yy-6d0, 5.6d0, xx-8d0, yy-8d0, 28.0d0)
         call cairo_pattern_add_color_stop_rgba (pat, 0d0, 1d0, 1d0, 1d0, 1d0)
         call cairo_pattern_add_color_stop_rgba (pat, 1d0, rr, vv, bl, 1d0)
         call cairo_set_source (my_cairo_context, pat)
         call cairo_arc(my_cairo_context, xx, yy, 20d0, 0d0, 2*pi)
         call cairo_fill (my_cairo_context)
         call cairo_pattern_destroy (pat)
          call cairo_set_line_width(my_cairo_context, 3.0d0)
          call cairo_stroke(my_cairo_context)

! Now, write the numbers (Genes)
call cairo_set_source_rgb(my_cairo_context, 0d0, 0d0, 0d0)
      call cairo_set_line_width(my_cairo_context, 3.0d0)
      call cairo_select_font_face(my_cairo_context, "Times"//c_null_char, CAIRO_FONT_SLANT_NORMAL, &
                                 &  CAIRO_FONT_WEIGHT_NORMAL)
      call cairo_set_font_size (my_cairo_context, 20d0)
      call cairo_move_to(my_cairo_context, npos(t,3)-20d0, npos(t,4))

        if(len(trim(genename(t))).eq.0)then
          write (str3,"(I3.1)") t
          call cairo_show_text (my_cairo_context, str3//c_null_char)
        else
          write (str12,"(A12)") genename(t)
          call cairo_show_text (my_cairo_context, str12//c_null_char)
        end if
    end do

    call cairo_move_to(my_cairo_context, 10d0, 660d0)

    ! Legend:

    call cairo_set_source_rgb(my_cairo_context, 0d0, 0d0, 1d0)

    call cairo_set_font_size (my_cairo_context, 18d0)

if(style1.ge.0) call cairo_set_source_rgb(my_cairo_context, 0d0, 0d0, 0d0)
if(style1.lt.0) call cairo_set_source_rgb(my_cairo_context, 1d0, 1d0, 1d0)
    call cairo_move_to(my_cairo_context, 10d0, 620d0)
    call cairo_show_text (my_cairo_context, "Types of Molecules:"//c_null_char)
             call cairo_set_source_rgb(my_cairo_context, 0.8d0, 0.8d0, 0.8d0) 
    call cairo_move_to(my_cairo_context, 10d0, 650d0)

         xx=15d0; yy=645d0
         pat = cairo_pattern_create_radial (xx-3d0, yy-3d0, 2.3d0, xx-4d0, yy-4d0, 14.0d0)
         call cairo_pattern_add_color_stop_rgba (pat, 0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0)
         call cairo_pattern_add_color_stop_rgba (pat, 1d0, 0.8d0, 0.8d0, 0.8d0, 1.0d0)
         call cairo_set_source (my_cairo_context, pat)
         call cairo_arc(my_cairo_context, xx, yy, 10d0, 0d0, 2*pi)
         call cairo_fill (my_cairo_context)
         call cairo_pattern_destroy (pat)

if(style1.ge.0) call cairo_set_source_rgb(my_cairo_context, 0d0, 0d0, 0d0)
if(style1.lt.0) call cairo_set_source_rgb(my_cairo_context, 1d0, 1d0, 1d0)
    call cairo_move_to(my_cairo_context, 10d0, 650d0)
    call cairo_show_text (my_cairo_context, "1 : Regulatory Molecule "//c_null_char)
             call cairo_set_source_rgb(my_cairo_context, 0.0d0, 0.6d0, 0.0d0)
    call cairo_move_to(my_cairo_context, 10d0, 670d0)

         xx=15d0; yy=665d0
         pat = cairo_pattern_create_radial (xx-3d0, yy-3d0, 2.3d0, xx-4d0, yy-4d0, 14.0d0)
         call cairo_pattern_add_color_stop_rgba (pat, 0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0)
         call cairo_pattern_add_color_stop_rgba (pat, 1d0, 0.0d0, 0.6d0, 0.0d0, 1.0d0)
         call cairo_set_source (my_cairo_context, pat)
         call cairo_arc(my_cairo_context, xx, yy, 10d0, 0d0, 2*pi)
         call cairo_fill (my_cairo_context)
         call cairo_pattern_destroy (pat)

    call cairo_set_source_rgb(my_cairo_context, 0.0d0, 0.6d0, 0.0d0)
    call cairo_move_to(my_cairo_context, 10d0, 670d0)
    call cairo_show_text (my_cairo_context, "2 : Regulatory Molecule, translatable"//c_null_char)
             call cairo_set_source_rgb(my_cairo_context, 0.0d0, 0.0d0, 0.8d0)
    call cairo_move_to(my_cairo_context, 10d0, 690d0)

         xx=15d0; yy=685d0
         pat = cairo_pattern_create_radial (xx-3d0, yy-3d0, 2.3d0, xx-4d0, yy-4d0, 14.0d0)
         call cairo_pattern_add_color_stop_rgba (pat, 0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0)
         call cairo_pattern_add_color_stop_rgba (pat, 1d0, 0.0d0, 0.0d0, 0.8d0, 1.0d0)
         call cairo_set_source (my_cairo_context, pat)
         call cairo_arc(my_cairo_context, xx, yy, 10d0, 0d0, 2*pi)
         call cairo_fill (my_cairo_context)
         call cairo_pattern_destroy (pat)

    call cairo_set_source_rgb(my_cairo_context, 0.0d0, 0.0d0, 0.8d0)
    call cairo_move_to(my_cairo_context, 10d0, 690d0)
    call cairo_show_text (my_cairo_context, "3 : Protein, intracellular"//c_null_char)
             call cairo_set_source_rgb(my_cairo_context, 0.8d0, 0.0d0, 0.0d0)

         xx=15d0; yy=705d0
         pat = cairo_pattern_create_radial (xx-3d0, yy-3d0, 2.3d0, xx-4d0, yy-4d0, 14.0d0)
         call cairo_pattern_add_color_stop_rgba (pat, 0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0)
         call cairo_pattern_add_color_stop_rgba (pat, 1d0, 0.8d0, 0.0d0, 0.0d0, 1.0d0)
         call cairo_set_source (my_cairo_context, pat)
         call cairo_arc(my_cairo_context, xx, yy, 10d0, 0d0, 2*pi)
         call cairo_fill (my_cairo_context)
         call cairo_pattern_destroy (pat)

    call cairo_set_source_rgb(my_cairo_context, 0.8d0, 0.0d0, 0.0d0)
    call cairo_move_to(my_cairo_context, 10d0, 710d0)
    call cairo_show_text (my_cairo_context, "4 : Protein, extracellular"//c_null_char)
             call cairo_set_source_rgb(my_cairo_context, 1.0d0, 0.5d0, 1.0d0)

         xx=445d0; yy=645d0
         pat = cairo_pattern_create_radial (xx-3d0, yy-3d0, 2.3d0, xx-4d0, yy-4d0, 14.0d0)
         call cairo_pattern_add_color_stop_rgba (pat, 0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0)
         call cairo_pattern_add_color_stop_rgba (pat, 1d0, 1.0d0, 0.5d0, 1.0d0, 1.0d0)
         call cairo_set_source (my_cairo_context, pat)
         call cairo_arc(my_cairo_context, xx, yy, 10d0, 0d0, 2*pi)
         call cairo_fill (my_cairo_context)
         call cairo_pattern_destroy (pat)

    call cairo_set_source_rgb(my_cairo_context, 1.0d0, 0.5d0, 1.0d0)
    call cairo_move_to(my_cairo_context, 440d0, 650d0)
    call cairo_show_text (my_cairo_context, "5 : Protein, apically localizing"//c_null_char)
             call cairo_set_source_rgb(my_cairo_context, 0.6d0, 0.6d0, 1.0d0)

         xx=445d0; yy=665d0
         pat = cairo_pattern_create_radial (xx-3d0, yy-3d0, 2.3d0, xx-4d0, yy-4d0, 14.0d0)
         call cairo_pattern_add_color_stop_rgba (pat, 0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0)
         call cairo_pattern_add_color_stop_rgba (pat, 1d0, 0.6d0, 0.6d0, 1.0d0, 1.0d0)
         call cairo_set_source (my_cairo_context, pat)
         call cairo_arc(my_cairo_context, xx, yy, 10d0, 0d0, 2*pi)
         call cairo_fill (my_cairo_context)
         call cairo_pattern_destroy (pat)

    call cairo_set_source_rgb(my_cairo_context, 0.6d0, 0.6d0, 1.0d0)
    call cairo_move_to(my_cairo_context, 440d0, 670d0)
    call cairo_show_text (my_cairo_context, "6 : Protein, basally localizing"//c_null_char)
             call cairo_set_source_rgb(my_cairo_context, 0.2d0, 1.0d0, 0.2d0)

         xx=445d0; yy=685d0
         pat = cairo_pattern_create_radial (xx-3d0, yy-3d0, 2.3d0, xx-4d0, yy-4d0, 14.0d0)
         call cairo_pattern_add_color_stop_rgba (pat, 0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0)
         call cairo_pattern_add_color_stop_rgba (pat, 1d0, 0.2d0, 1.0d0, 0.2d0, 1.0d0)
         call cairo_set_source (my_cairo_context, pat)
         call cairo_arc(my_cairo_context, xx, yy, 10d0, 0d0, 2*pi)
         call cairo_fill (my_cairo_context)
         call cairo_pattern_destroy (pat)

    call cairo_set_source_rgb(my_cairo_context, 0.2d0, 0.9d0, 0.2d0)
    call cairo_move_to(my_cairo_context, 440d0, 690d0)
    call cairo_show_text (my_cairo_context, "7 : Receptor, juxtacrine signalling(7)"//c_null_char)
             call cairo_set_source_rgb(my_cairo_context, 0.8d0, 0.8d0, 0.2d0)

         xx=445d0; yy=705d0
         pat = cairo_pattern_create_radial (xx-3d0, yy-3d0, 2.3d0, xx-4d0, yy-4d0, 14.0d0)
         call cairo_pattern_add_color_stop_rgba (pat, 0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0)
         call cairo_pattern_add_color_stop_rgba (pat, 1d0, 1.0d0, 0.8d0, 0.0d0, 1.0d0)
         call cairo_set_source (my_cairo_context, pat)
         call cairo_arc(my_cairo_context, xx, yy, 10d0, 0d0, 2*pi)
         call cairo_fill (my_cairo_context)
         call cairo_pattern_destroy (pat)

    call cairo_set_source_rgb(my_cairo_context, 1.0d0, 0.8d0, 0.0d0)
    call cairo_move_to(my_cairo_context, 440d0, 710d0)
    call cairo_show_text (my_cairo_context, "8 : Receptor, paracrine signalling(4)"//c_null_char)

! Write the numbers (Cell behaviours)

    do s=1,23
       t=npars(s)
       call cairo_set_font_size (my_cairo_context, 20d0)
       call cairo_set_source_rgb(my_cairo_context, 0.6d0, 0.5d0, 0.9d0)
        call cairo_set_line_width(my_cairo_context, 3.0d0)
        call cairo_select_font_face(my_cairo_context, "Arial"//c_null_char, CAIRO_FONT_SLANT_NORMAL, &
                                 &  CAIRO_FONT_WEIGHT_NORMAL)
        call cairo_move_to(my_cairo_context, 85d0+s*500/20, 190d0)

        call cairo_rotate (my_cairo_context, 67.6d0)

        call cairo_show_text (my_cairo_context, trim(nodeparams(t))//c_null_char)
        call cairo_rotate (my_cairo_context, -67.6d0)
    end do

    call cairo_move_to(my_cairo_context, 10d0, 80d0)
    call cairo_set_font_size (my_cairo_context, 20d0)
    call cairo_show_text (my_cairo_context, "Node parameters"//c_null_char)

    call cairo_destroy(my_cairo_context)
    ret = FALSE
  end function draw_network7


  function draw_network8 (widget, event, gdata) result(ret)  bind(c)
    use iso_c_binding, only: c_int, c_ptr
    use io; use genetic; use general
    use global_widgets

    integer(c_int)    :: ret
    character*140 :: fii
    character*140 :: fij
    character*3 :: str3
    character*12 :: str12
    character*8 :: str8
    type(c_ptr), value, intent(in) :: widget, event, gdata
    type(c_ptr) :: my_cairo_context, pat
    integer :: cstatus
    integer :: t,t1,tot,s
    real*8 :: tt1,tt2,tt3,tt4,xx,yy,rr,vv,bl

if(rompetout.ne.6)then
if(yetthere==1)then
call gtk_widget_destroy(my_choice_window2)
yetthere=0
endif
if(yetthere1==1)then
call gtk_widget_destroy(my_choice_window)
yetthere1=0
endif
if(yetthere3==1)then
call gtk_widget_destroy(my_choice_window3)
yetthere3=0
endif
endif
rompetout=6

    my_cairo_context = gdk_cairo_create (gtk_widget_get_window(widget))

if(style1==-1)then
   call cairo_set_source_rgb(my_cairo_context, 0d0, 0d0, 0d0)
   call cairo_rectangle(my_cairo_context,0d0,0d0,800d0,800d0)
   call cairo_fill(my_cairo_context)
end if

    ! Lines:
    call cairo_set_source_rgb(my_cairo_context, 0.9d0, 0.6d0, 0.3d0)
    call cairo_set_line_width(my_cairo_context, 2d0)
    do t = 0, int(height), +20
      call cairo_move_to(my_cairo_context, 0d0, t*1d0)  
      call cairo_line_to(my_cairo_context, t*1d0, height*1d0)
      call cairo_stroke(my_cairo_context) 
    end do
  
    ! Text:
    write (str8,"(I8.1)") getot
    call cairo_set_source_rgb(my_cairo_context, 0d0, 0d0, 1d0)
    call cairo_select_font_face(my_cairo_context, "Times"//c_null_char, CAIRO_FONT_SLANT_NORMAL, &
                                 &  CAIRO_FONT_WEIGHT_NORMAL)
    call cairo_set_font_size (my_cairo_context, 30d0)
    call cairo_move_to(my_cairo_context, 10d0, 30d0)
    call cairo_show_text (my_cairo_context, "Modifying adhesions"//c_null_char)
    call cairo_move_to(my_cairo_context, 445d0, 70d0)
    call cairo_set_source_rgb(my_cairo_context, 0d0, 0d0, 0d0)
    call cairo_show_text (my_cairo_context, "Add a new Bij element"//c_null_char)
        call cairo_move_to(my_cairo_context, 350d0, 30d0)
! kadhs
if(togglen==0)then
  if(allocated(kpos)) deallocate(kpos)   ! more npos
  allocate(kpos(ntipusadh,2))
  do t=1,ntipusadh
    kpos(t,1)=350d0+260d0*cos(t*pi/ntipusadh)
    kpos(t,2)=350d0+260d0*sin(t*pi/ntipusadh)
  end do
end if

if(togglen==0)then
  if(allocated(npos)) deallocate(npos)   ! more npos
  allocate(npos(ng,6))
  do t=1,ng
    npos(t,1)=350d0+260d0*cos(t*2d0*pi/ng)
    npos(t,2)=350d0+260d0*sin(t*2d0*pi/ng)
    npos(t,3)=350d0+260d0*cos(t*pi/ng)
    npos(t,4)=350d0+260d0*sin(t*pi/ng)
    npos(t,5)=50d0+700d0/(ng+1d0)*t
    npos(t,6)=150d0
  end do
end if

do t=1,ng
         xx=npos(t,5); yy=npos(t,6)
if(wcl==t)then
            pat = cairo_pattern_create_radial (xx, yy, 5.6d0, xx, yy, 45.0d0)
                call cairo_pattern_add_color_stop_rgba (pat,0d0, 1d0, 0.8d0, 0d0, 1.0d0)
                if(style1.ge.0) call cairo_pattern_add_color_stop_rgba (pat,1d0, 0.97d0, 0.97d0, 0.97d0, 1.0d0)
                if(style1.lt.0) call cairo_pattern_add_color_stop_rgba (pat,1d0, 0d0, 0d0, 0d0, 1.0d0)
                call cairo_set_source (my_cairo_context,pat)
                call cairo_arc(my_cairo_context, xx, yy, 45d0, 0d0, 2*pi)
                call cairo_fill (my_cairo_context)
                call cairo_pattern_destroy (pat)
end if
end do

    !Arrows:
  call cairo_new_sub_path(my_cairo_context)

    do t = 1,ng

      s = gen(t)%wa(1) ! its kadh  
      if(s.ne.0)then
        call cairo_set_source_rgb(my_cairo_context, 0.5d0, 0.5d0, 0.9d0)
        call cairo_set_line_width(my_cairo_context, 1d0)
        call cairo_move_to(my_cairo_context, npos(t,5), npos(t,6))  
        call cairo_line_to(my_cairo_context, kpos(s,1), kpos(s,2))        
        call cairo_stroke(my_cairo_context)
      end if

        xx=npos(t,5); yy=npos(t,6)

       if(gen(t)%kindof==1)then
             rr=0.8d0; vv=0.8d0; bl=0.8d0
             call cairo_set_source_rgb(my_cairo_context, 0.8d0, 0.8d0, 0.8d0) 
       elseif(gen(t)%kindof==2)then
             rr=0.0d0; vv=0.5d0; bl=0.0d0
             call cairo_set_source_rgb(my_cairo_context, 0.0d0, 0.4d0, 0.0d0)
       elseif(gen(t)%kindof==3)then
             rr=0.0d0; vv=0.0d0; bl=0.8d0
             call cairo_set_source_rgb(my_cairo_context, 0.0d0, 0.0d0, 0.7d0)
       elseif(gen(t)%kindof==4)then
             rr=0.8d0; vv=0.0d0; bl=0.0d0
             call cairo_set_source_rgb(my_cairo_context, 0.8d0, 0.0d0, 0.0d0)
       elseif(gen(t)%kindof==5)then
             rr=1.0d0; vv=0.5d0; bl=1.0d0
             call cairo_set_source_rgb(my_cairo_context, 1.0d0, 0.5d0, 1.0d0)
       elseif(gen(t)%kindof==6)then
             rr=0.6d0; vv=0.6d0; bl=1.0d0
             call cairo_set_source_rgb(my_cairo_context, 0.4d0, 0.7d0, 1.0d0)
       elseif(gen(t)%kindof==7)then
             rr=0.2d0; vv=1.0d0; bl=0.2d0
             call cairo_set_source_rgb(my_cairo_context, 0.2d0, 1.0d0, 0.2d0)
       elseif(gen(t)%kindof==8)then
             rr=1.0d0; vv=0.8d0; bl=0.0d0
             call cairo_set_source_rgb(my_cairo_context, 1.0d0, 8.0d0, 0.2d0)
             call cairo_set_source_rgb(my_cairo_context, 1.0d0, 1.0d0, 0.2d0) !,0.2d0z
       endif

       pat = cairo_pattern_create_radial (xx-6d0, yy-6d0, 5.6d0, xx-8d0, yy-8d0, 28.0d0)
       call cairo_pattern_add_color_stop_rgba (pat, 0d0, 1d0, 1d0, 1d0, 1d0)
       call cairo_pattern_add_color_stop_rgba (pat, 1d0, rr, vv, bl, 1d0)
       call cairo_set_source (my_cairo_context, pat)
       call cairo_arc(my_cairo_context, xx, yy, 20d0, 0d0, 2*pi)
       call cairo_fill (my_cairo_context)
       call cairo_pattern_destroy (pat)
       call cairo_set_line_width(my_cairo_context, 3.0d0)
       call cairo_stroke(my_cairo_context)

! Now, write the numbers (Genes)
if(style1.ge.0) call cairo_set_source_rgb(my_cairo_context, 0d0, 0d0, 0d0)
if(style1.lt.0) call cairo_set_source_rgb(my_cairo_context, 1d0, 1d0, 1d0)
       call cairo_set_line_width(my_cairo_context, 3.0d0)
       call cairo_select_font_face(my_cairo_context, "Times"//c_null_char, CAIRO_FONT_SLANT_NORMAL, &
                                 &  CAIRO_FONT_WEIGHT_NORMAL)
       call cairo_set_font_size (my_cairo_context, 20d0)
       call cairo_move_to(my_cairo_context, npos(t,5)-20d0, npos(t,6))
       call cairo_set_source_rgb(my_cairo_context, 0d0, 0d0, 0d0)
       if(len(trim(genename(t))).eq.0)then
         write (str3,"(I3.1)") t
         call cairo_show_text (my_cairo_context, str3//c_null_char)
       else
         write (str12,"(A12)") genename(t)
         call cairo_show_text (my_cairo_context, str12//c_null_char)
       end if
    end do

! ARROWS between kadh
    if(ntipusadh.ne.0)then
    do s = 1,ntipusadh
      do t = 1,ntipusadh
        if(kadh(s,t).gt.ev) then
          call cairo_set_source_rgb(my_cairo_context, 0.2d0, 0.8d0, 0.2d0)
          if(s.ne.t)then
            tt1=0.3d0+2.0d0*sqrt(kadh(s,t))
            if(tt1.gt.10) tt1=10
            call cairo_set_line_width(my_cairo_context, tt1)
            call cairo_move_to(my_cairo_context, kpos(s,1), kpos(s,2))  
            call cairo_line_to(my_cairo_context, kpos(t,1), kpos(t,2))        
          else
            call cairo_move_to(my_cairo_context, kpos(s,1), kpos(s,2)+10d0)
            call cairo_arc(my_cairo_context, kpos(t,1)-20d0, kpos(t,2)+10d0, 22d0, 0d0, 2*pi)
          end if
          call cairo_stroke(my_cairo_context) 
        else if(kadh(s,t).lt.-ev) then
          call cairo_set_source_rgb(my_cairo_context, 1d0, 0d0, 0d0)
          if(s.ne.t)then
            tt1=0.3d0+2.0d0*sqrt(kadh(s,t)*-1)
            if(tt1.gt.10) tt1=10
            call cairo_set_line_width(my_cairo_context, tt1)
            call cairo_move_to(my_cairo_context, kpos(s,1), kpos(s,2))  
            call cairo_line_to(my_cairo_context, kpos(t,1), kpos(t,2))        
          else
            call cairo_move_to(my_cairo_context, kpos(s,1), kpos(s,2)+10d0)
            call cairo_arc(my_cairo_context, kpos(t,1)-20d0, kpos(t,2)+10d0, 22d0, 0d0, 2*pi)
          end if
          call cairo_stroke(my_cairo_context)
        end if
      end do

       xx=kpos(s,1); yy=kpos(s,2)

       rr=0.7d0; vv=0.7d0; bl=0.2d0
      ! call cairo_set_source_rgb(my_cairo_context, 1.0d0, 8.0d0, 0.2d0)

       pat = cairo_pattern_create_radial (xx-6d0, yy-6d0, 5.6d0, xx-8d0, yy-8d0, 28.0d0)
       call cairo_pattern_add_color_stop_rgba (pat, 0d0, 1d0, 1d0, 1d0, 1d0)
       call cairo_pattern_add_color_stop_rgba (pat, 1d0, rr, vv, bl, 1d0)
       call cairo_set_source (my_cairo_context, pat)
       call cairo_arc(my_cairo_context, xx, yy, 20d0, 0d0, 2*pi)
       call cairo_fill (my_cairo_context)
       call cairo_pattern_destroy (pat)
       call cairo_set_line_width(my_cairo_context, 3.0d0)
       call cairo_stroke(my_cairo_context)

! Now, write the numbers (kadh)
if(style1.ge.0) call cairo_set_source_rgb(my_cairo_context, 0d0, 0d0, 0d0)
if(style1.lt.0) call cairo_set_source_rgb(my_cairo_context, 1d0, 1d0, 1d0)
call cairo_set_source_rgb(my_cairo_context, 0d0, 0d0, 0d0)
       call cairo_set_line_width(my_cairo_context, 3.0d0)
       call cairo_select_font_face(my_cairo_context, "Times"//c_null_char, CAIRO_FONT_SLANT_NORMAL, &
                                 &  CAIRO_FONT_WEIGHT_NORMAL)
       call cairo_set_font_size (my_cairo_context, 20d0)
       call cairo_move_to(my_cairo_context, kpos(s,1)-20d0, kpos(s,2))
         write (str3,"(I3.1)") s
         call cairo_show_text (my_cairo_context, str3//c_null_char)
       call cairo_move_to(my_cairo_context, 10d0, 660d0)
     end do
     end if

    ! Legend:

    call cairo_set_source_rgb(my_cairo_context, 0d0, 0d0, 1d0)

    call cairo_set_font_size (my_cairo_context, 18d0)

if(style1.ge.0) call cairo_set_source_rgb(my_cairo_context, 0d0, 0d0, 0d0)
if(style1.lt.0) call cairo_set_source_rgb(my_cairo_context, 1d0, 1d0, 1d0)
    call cairo_move_to(my_cairo_context, 10d0, 620d0)
    call cairo_show_text (my_cairo_context, "Types of Molecules:"//c_null_char)
             call cairo_set_source_rgb(my_cairo_context, 0.8d0, 0.8d0, 0.8d0) 
    call cairo_move_to(my_cairo_context, 10d0, 650d0)

         xx=15d0; yy=645d0
         pat = cairo_pattern_create_radial (xx-3d0, yy-3d0, 2.3d0, xx-4d0, yy-4d0, 14.0d0)
         call cairo_pattern_add_color_stop_rgba (pat, 0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0)
         call cairo_pattern_add_color_stop_rgba (pat, 1d0, 0.8d0, 0.8d0, 0.8d0, 1.0d0)
         call cairo_set_source (my_cairo_context, pat)
         call cairo_arc(my_cairo_context, xx, yy, 10d0, 0d0, 2*pi)
         call cairo_fill (my_cairo_context)
         call cairo_pattern_destroy (pat)

if(style1.ge.0) call cairo_set_source_rgb(my_cairo_context, 0d0, 0d0, 0d0)
if(style1.lt.0) call cairo_set_source_rgb(my_cairo_context, 1d0, 1d0, 1d0)
    call cairo_move_to(my_cairo_context, 10d0, 650d0)
    call cairo_show_text (my_cairo_context, "1 : Regulatory Molecule "//c_null_char)
             call cairo_set_source_rgb(my_cairo_context, 0.0d0, 0.6d0, 0.0d0)
    call cairo_move_to(my_cairo_context, 10d0, 670d0)

         xx=15d0; yy=665d0
         pat = cairo_pattern_create_radial (xx-3d0, yy-3d0, 2.3d0, xx-4d0, yy-4d0, 14.0d0)
         call cairo_pattern_add_color_stop_rgba (pat, 0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0)
         call cairo_pattern_add_color_stop_rgba (pat, 1d0, 0.0d0, 0.6d0, 0.0d0, 1.0d0)
         call cairo_set_source (my_cairo_context, pat)
         call cairo_arc(my_cairo_context, xx, yy, 10d0, 0d0, 2*pi)
         call cairo_fill (my_cairo_context)
         call cairo_pattern_destroy (pat)

    call cairo_set_source_rgb(my_cairo_context, 0.0d0, 0.6d0, 0.0d0)
    call cairo_move_to(my_cairo_context, 10d0, 670d0)
    call cairo_show_text (my_cairo_context, "2 : Regulatory Molecule, translatable"//c_null_char)
             call cairo_set_source_rgb(my_cairo_context, 0.0d0, 0.0d0, 0.8d0)
    call cairo_move_to(my_cairo_context, 10d0, 690d0)

         xx=15d0; yy=685d0
         pat = cairo_pattern_create_radial (xx-3d0, yy-3d0, 2.3d0, xx-4d0, yy-4d0, 14.0d0)
         call cairo_pattern_add_color_stop_rgba (pat, 0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0)
         call cairo_pattern_add_color_stop_rgba (pat, 1d0, 0.0d0, 0.0d0, 0.8d0, 1.0d0)
         call cairo_set_source (my_cairo_context, pat)
         call cairo_arc(my_cairo_context, xx, yy, 10d0, 0d0, 2*pi)
         call cairo_fill (my_cairo_context)
         call cairo_pattern_destroy (pat)

    call cairo_set_source_rgb(my_cairo_context, 0.0d0, 0.0d0, 0.8d0)
    call cairo_move_to(my_cairo_context, 10d0, 690d0)
    call cairo_show_text (my_cairo_context, "3 : Protein, intracellular"//c_null_char)
             call cairo_set_source_rgb(my_cairo_context, 0.8d0, 0.0d0, 0.0d0)

         xx=15d0; yy=705d0
         pat = cairo_pattern_create_radial (xx-3d0, yy-3d0, 2.3d0, xx-4d0, yy-4d0, 14.0d0)
         call cairo_pattern_add_color_stop_rgba (pat, 0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0)
         call cairo_pattern_add_color_stop_rgba (pat, 1d0, 0.8d0, 0.0d0, 0.0d0, 1.0d0)
         call cairo_set_source (my_cairo_context, pat)
         call cairo_arc(my_cairo_context, xx, yy, 10d0, 0d0, 2*pi)
         call cairo_fill (my_cairo_context)
         call cairo_pattern_destroy (pat)

    call cairo_set_source_rgb(my_cairo_context, 0.8d0, 0.0d0, 0.0d0)
    call cairo_move_to(my_cairo_context, 10d0, 710d0)
    call cairo_show_text (my_cairo_context, "4 : Protein, extracellular"//c_null_char)
             call cairo_set_source_rgb(my_cairo_context, 1.0d0, 0.5d0, 1.0d0)

         xx=445d0; yy=645d0
         pat = cairo_pattern_create_radial (xx-3d0, yy-3d0, 2.3d0, xx-4d0, yy-4d0, 14.0d0)
         call cairo_pattern_add_color_stop_rgba (pat, 0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0)
         call cairo_pattern_add_color_stop_rgba (pat, 1d0, 1.0d0, 0.5d0, 1.0d0, 1.0d0)
         call cairo_set_source (my_cairo_context, pat)
         call cairo_arc(my_cairo_context, xx, yy, 10d0, 0d0, 2*pi)
         call cairo_fill (my_cairo_context)
         call cairo_pattern_destroy (pat)

    call cairo_set_source_rgb(my_cairo_context, 1.0d0, 0.5d0, 1.0d0)
    call cairo_move_to(my_cairo_context, 440d0, 650d0)
    call cairo_show_text (my_cairo_context, "5 : Protein, apically localizing"//c_null_char)
             call cairo_set_source_rgb(my_cairo_context, 0.6d0, 0.6d0, 1.0d0)

         xx=445d0; yy=665d0
         pat = cairo_pattern_create_radial (xx-3d0, yy-3d0, 2.3d0, xx-4d0, yy-4d0, 14.0d0)
         call cairo_pattern_add_color_stop_rgba (pat, 0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0)
         call cairo_pattern_add_color_stop_rgba (pat, 1d0, 0.6d0, 0.6d0, 1.0d0, 1.0d0)
         call cairo_set_source (my_cairo_context, pat)
         call cairo_arc(my_cairo_context, xx, yy, 10d0, 0d0, 2*pi)
         call cairo_fill (my_cairo_context)
         call cairo_pattern_destroy (pat)

    call cairo_set_source_rgb(my_cairo_context, 0.6d0, 0.6d0, 1.0d0)
    call cairo_move_to(my_cairo_context, 440d0, 670d0)
    call cairo_show_text (my_cairo_context, "6 : Protein, basally localizing"//c_null_char)
             call cairo_set_source_rgb(my_cairo_context, 0.2d0, 1.0d0, 0.2d0)

         xx=445d0; yy=685d0
         pat = cairo_pattern_create_radial (xx-3d0, yy-3d0, 2.3d0, xx-4d0, yy-4d0, 14.0d0)
         call cairo_pattern_add_color_stop_rgba (pat, 0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0)
         call cairo_pattern_add_color_stop_rgba (pat, 1d0, 0.2d0, 1.0d0, 0.2d0, 1.0d0)
         call cairo_set_source (my_cairo_context, pat)
         call cairo_arc(my_cairo_context, xx, yy, 10d0, 0d0, 2*pi)
         call cairo_fill (my_cairo_context)
         call cairo_pattern_destroy (pat)

    call cairo_set_source_rgb(my_cairo_context, 0.2d0, 0.9d0, 0.2d0)
    call cairo_move_to(my_cairo_context, 440d0, 690d0)
    call cairo_show_text (my_cairo_context, "7 : Receptor, juxtacrine signalling(7)"//c_null_char)
             call cairo_set_source_rgb(my_cairo_context, 0.8d0, 0.8d0, 0.2d0)

         xx=445d0; yy=705d0
         pat = cairo_pattern_create_radial (xx-3d0, yy-3d0, 2.3d0, xx-4d0, yy-4d0, 14.0d0)
         call cairo_pattern_add_color_stop_rgba (pat, 0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0)
         call cairo_pattern_add_color_stop_rgba (pat, 1d0, 1.0d0, 0.8d0, 0.0d0, 1.0d0)
         call cairo_set_source (my_cairo_context, pat)
         call cairo_arc(my_cairo_context, xx, yy, 10d0, 0d0, 2*pi)
         call cairo_fill (my_cairo_context)
         call cairo_pattern_destroy (pat)

    call cairo_set_source_rgb(my_cairo_context, 1.0d0, 0.8d0, 0.0d0)
    call cairo_move_to(my_cairo_context, 440d0, 710d0)
    call cairo_show_text (my_cairo_context, "8 : Receptor, paracrine signalling(4)"//c_null_char)
             call cairo_set_source_rgb(my_cairo_context, 0.7d0, 0.7d0, 0.2d0)

         xx=440d0; yy=25d0
         pat = cairo_pattern_create_radial (xx-3d0, yy-3d0, 2.3d0, xx-4d0, yy-4d0, 14.0d0)
         call cairo_pattern_add_color_stop_rgba (pat, 0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0)
         call cairo_pattern_add_color_stop_rgba (pat, 1d0, 0.7d0, 0.7d0, 0.2d0, 1.0d0)
         call cairo_set_source (my_cairo_context, pat)
         call cairo_arc(my_cairo_context, xx, yy, 10d0, 0d0, 2*pi)
         call cairo_fill (my_cairo_context)
         call cairo_pattern_destroy (pat)

if(style1.ge.0) call cairo_set_source_rgb(my_cairo_context, 0d0, 0d0, 0d0)
if(style1.lt.0) call cairo_set_source_rgb(my_cairo_context, 1d0, 1d0, 1d0)
    call cairo_move_to(my_cairo_context, 445d0, 30d0)
    call cairo_show_text (my_cairo_context, "Adhesion molecule types"//c_null_char)

    call cairo_destroy(my_cairo_context)
    ret = FALSE

  end function draw_network8

end module drawings

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module on_display_handlers

use iso_c_binding, only: c_ptr, c_funloc
use drawings
use genetic

type(c_ptr) :: spinButtonc1,spinButtonc2,spinButtonc3,spinButtonc4,spinButtonc5,spinButtonc6,spinButtonc7
type(c_ptr) :: spinButtonc8,spinButtonc9
type(c_ptr) :: spinButtoncc1, spinButtoncc2, entryc
 
contains

subroutine GRN_plot(plotid)

use global_widgets

  integer::plotid,widthm,heightm
  type(c_ptr) :: my_window_n
  type(c_ptr) :: my_drawing_area_n

  call gtk_disable_setlocale ()
  call gtk_init ()
  
  ! Properties of the main window :
  widthm = 700
  heightm = 700
  
  if(plotid.gt.10)then
    my_window_n = gtk_window_new (GTK_window_TOPLEVEL)
    call gtk_window_set_default_size(my_window_n, widthm, heightm)
    call gtk_window_set_title(my_window_n, "GRN: genetic interactions / T-Matrix"//c_null_char)
    my_drawing_area_n= gtk_drawing_area_new()
  end if
  call g_signal_connect (my_window, "delete-event"//c_null_char, c_funloc(delete_event))

  if(plotid==1) call g_signal_connect (my_drawing_area, "expose-event"//c_null_char, c_funloc(draw_network))
  if(plotid==2) call g_signal_connect (my_drawing_area, "expose-event"//c_null_char, c_funloc(draw_network3))
  if(plotid==3) call g_signal_connect (my_drawing_area, "expose-event"//c_null_char, c_funloc(draw_network4))
  if(plotid==4) call g_signal_connect (my_drawing_area, "expose-event"//c_null_char, c_funloc(draw_network5))
  if(plotid==5) call g_signal_connect (my_drawing_area, "expose-event"//c_null_char, c_funloc(draw_network7))
  if(plotid==6) call g_signal_connect (my_drawing_area, "expose-event"//c_null_char, c_funloc(draw_network8))
  if(plotid==11) call g_signal_connect (my_drawing_area_n, "expose-event"//c_null_char, c_funloc(draw_network))
  if(plotid==12) call g_signal_connect (my_drawing_area_n, "expose-event"//c_null_char, c_funloc(draw_network3))
  if(plotid==13) call g_signal_connect (my_drawing_area_n, "expose-event"//c_null_char, c_funloc(draw_network4))
  if(plotid==14) call g_signal_connect (my_drawing_area_n, "expose-event"//c_null_char, c_funloc(draw_network5))
  if(plotid==15) call g_signal_connect (my_drawing_area_n, "expose-event"//c_null_char, c_funloc(draw_network7))
  if(plotid==16) call g_signal_connect (my_drawing_area_n, "expose-event"//c_null_char, c_funloc(draw_network8))

  if(plotid.gt.10)then
    call gtk_container_add(my_window_n, my_drawing_area_n)
    call gtk_widget_show (my_drawing_area_n)
    call gtk_widget_show (my_window_n)
  end if

end subroutine GRN_plot

function name_all_connect(widget, event, gdata) result(ret)  bind(c)

  use global_widgets
  use genetic

  type(c_ptr), value :: widget, gdata, event

  integer(c_int)    :: ret

  type(c_ptr) :: boxc,tablen
  type(c_ptr) :: butt_all_names
  integer :: leng,choice,whichi,checky
  character*12 :: entryda
  character*3 :: genenu
  character*2 :: strng1,strng2

  allocate(labelna(ng)) 

  if(yetthere2==1) call gtk_widget_destroy(my_choice_windown) !XX!

  my_choice_windown = gtk_window_new (GTK_WINDOW_TOPLEVEL)
  yetthere2=1
  if(ng.le.25)then
    leng=75+25*ng
  else
    leng=700
  end if

  call g_signal_connect (my_choice_windown, "delete-event"//c_null_char, c_funloc(delete_event))
  call gtk_window_set_default_size(my_choice_windown, leng, 100)
  call gtk_window_set_title(my_choice_windown, "Name the genes"//c_null_char)

    butt_all_names= gtk_button_new_with_mnemonic ("Apply all"//c_null_char)
    tablen = gtk_table_new (6_c_int, 2_c_int, TRUE)

  if(allocated(entryna)) deallocate(entryna)
  allocate(entryna(ng))

  do i=1,ng
    write (genenu,"(I3.1)") i
    labelna(i) = gtk_label_new(genenu//c_null_char)
    entryna(i) = gtk_entry_new()
    entryda = genename(i)
    call gtk_entry_set_text(entryna(i),adjustl(entryda))
 
    write (strng1,"(I2.1)") i-1
    strng1=trim(strng1)//"_c_int"
    strng1=trim(strng1)
    write (strng2,"(I2.1)") i
    strng2=trim(strng2)//"_c_int"
    strng2=trim(strng2)

    call gtk_table_attach_defaults(tablen, labelna(i), 0_c_int, 1_c_int, i-1, i)
    call gtk_table_attach_defaults(tablen, entryna(i), 1_c_int, 2_c_int, i-1, i)
  end do

    write (strng1,"(I2.1)") i-1
    strng1=trim(strng1)//"_c_int"
    strng1=trim(strng1)
    write (strng2,"(I2.1)") i
    strng2=trim(strng2)//"_c_int"
    strng2=trim(strng2)

  call g_signal_connect (butt_all_names, "clicked"//c_null_char, c_funloc(my_names_accept))
  call gtk_table_attach_defaults(tablen, butt_all_names, 0_c_int, 2_c_int, i-1, i)

  call gtk_container_add (my_choice_windown, tablen)

  call gtk_widget_show_all (my_choice_windown)


end function name_all_connect

subroutine choice_prop_box_build(whichi)

  use global_widgets

  type(c_ptr) :: boxc
  type(c_ptr) :: labelc1, labelc2, labelc3, labelc4, labelc5, labelc6, buttonc, buttond, buttoncn, tablec
  integer :: choice,whichi,checky
  character*12 :: entryd
  real*8 :: whichr,ngr,p1,p2,p3,p4

if(rb_arrange_tog.ne.1)then
  whichr=whichi; ngr=ng
  p1=gen(whichr)%kindof
  p2=gen(whichr)%diffu
  p3=gen(whichr)%idiffu
  p4=gen(whichr)%mu

  my_choice_window2 = gtk_window_new (GTK_WINDOW_TOPLEVEL)
  yetthere=1

  call g_signal_connect (my_choice_window2, "delete-event"//c_null_char, c_funloc(delete_event))
  call gtk_window_set_default_size(my_choice_window2, 250, 150)
  call gtk_window_set_title(my_choice_window2, "Gene properties"//c_null_char)
  labelc1 = gtk_label_new("Gene"//c_null_char)
  labelc2 = gtk_label_new("kindof"//c_null_char)
  labelc3 = gtk_label_new("Diffusivity"//c_null_char)
  labelc5 = gtk_label_new("Degradation rate"//c_null_char)
  labelc6 = gtk_label_new("Name this gene"//c_null_char)
  spinButtonc1 = gtk_spin_button_new (gtk_adjustment_new(whichr,whichr,whichr,&
       & 1.0d0,0.5d0,0d0),0.05d0, 0_c_int)
  spinButtonc2 = gtk_spin_button_new (gtk_adjustment_new(p1,1d0,8d0,&
       & 1.0d0,0.5d0,0d0),0.05d0, 0_c_int)
  spinButtonc3 = gtk_spin_button_new (gtk_adjustment_new(p2,-10d6,+10d6,&
       & 0.1d0,0.5d0,0d0),0.05d0, 7_c_int)
  spinButtonc5 = gtk_spin_button_new (gtk_adjustment_new(p4,-10d6,+10d6,&
       & 0.1d0,0.5d0,0d0),0.05d0, 7_c_int)
  entryc = gtk_entry_new()
  entryd = genename(int(whichr))

  call gtk_entry_set_text(entryc,trim(entryd))
  buttonc = gtk_button_new_with_mnemonic ("Accept"//c_null_char)
  buttoncn = gtk_button_new_with_mnemonic ("Apply Name"//c_null_char)
  tablec = gtk_table_new (6_c_int, 2_c_int, TRUE)
  call g_signal_connect (buttonc, "clicked"//c_null_char, c_funloc(my_p_accept))
  call g_signal_connect (buttoncn, "clicked"//c_null_char, c_funloc(my_p_accept))
 
  call gtk_table_attach_defaults(tablec, labelc1, 0_c_int, 2_c_int, 0_c_int, 1_c_int)
  call gtk_table_attach_defaults(tablec, labelc2, 0_c_int, 2_c_int, 2_c_int, 3_c_int)
  call gtk_table_attach_defaults(tablec, labelc3, 0_c_int, 2_c_int, 4_c_int, 5_c_int)
  call gtk_table_attach_defaults(tablec, labelc5, 0_c_int, 2_c_int, 6_c_int, 7_c_int)
  call gtk_table_attach_defaults(tablec, spinButtonc1, 0_c_int, 2_c_int, 1_c_int, 2_c_int)
  call gtk_table_attach_defaults(tablec, spinButtonc2, 0_c_int, 2_c_int, 3_c_int, 4_c_int)
  call gtk_table_attach_defaults(tablec, spinButtonc3, 0_c_int, 2_c_int, 5_c_int, 6_c_int)
  call gtk_table_attach_defaults(tablec, spinButtonc5, 0_c_int, 2_c_int, 7_c_int, 8_c_int)
  call gtk_table_attach_defaults(tablec, labelc6, 0_c_int, 2_c_int, 8_c_int, 9_c_int)
  call gtk_table_attach_defaults(tablec, entryc, 0_c_int, 1_c_int, 9_c_int, 10_c_int)
  call gtk_table_attach_defaults(tablec, buttonc, 0_c_int, 2_c_int, 10_c_int, 11_c_int)
  call gtk_table_attach_defaults(tablec, buttoncn, 1_c_int, 2_c_int, 9_c_int, 10_c_int)

  call gtk_container_add (my_choice_window2, tablec)

  call gtk_widget_show_all (my_choice_window2)

  test0=1

end if

end subroutine choice_prop_box_build

subroutine choice_box_build(whichsaves,which,cuals,buttt)

  use global_widgets
  use some_widgets

  type(c_ptr) :: boxc,checkd
  type(c_ptr) :: labelc1, labelc2, labelc3, labelc4, buttonc, buttond, tablec
  integer :: choice,which,whichsaves,cuals,buttt,mm1,mm2,mi,mm3,mm4
  real*8 :: whichsaver,whichr,ngr,interac,nkad

  if(yetthere1==1) call gtk_widget_destroy(my_choice_window) !XX!
  whichsaver=whichsaves; whichr=which; ngr=ng; nkad=ntipusadh
  choiceboxthere=1
  if((cuals==0))then

    my_choice_window = gtk_window_new (GTK_WINDOW_TOPLEVEL)
    yetthere1=1
    call gtk_window_set_default_size(my_choice_window, 150, 250)
    call gtk_window_set_title(my_choice_window, "Gene interactions"//c_null_char)
    labelc1 = gtk_label_new("Gen 1"//c_null_char)

    spinButtonc4 = gtk_spin_button_new (gtk_adjustment_new(whichsaver,whichsaver,whichsaver,&
       & 1.0d0,0.5d0,0d0),0.05d0, 0_c_int)
    labelc2 = gtk_label_new("Gen 2"//c_null_char)
    spinButtonc6 = gtk_spin_button_new (gtk_adjustment_new(whichr,whichr,whichr,&
       & 1.0d0,0.5d0,0d0),0.05d0, 0_c_int)
    labelc3 = gtk_label_new("Interaction strength"//c_null_char)

    interac=gen(which)%w(whichsaves)

    if(interac==0) interac=0.1d0
    if((faster==-1).and.(buttt.eq.3).and.(interac.eq.0)) interac=-0.1d0
    if((faster==-1).and.(buttt.eq.1).and.(interac.eq.0)) interac=0.1d0
    if((faster==-1).and.(buttt.eq.3).and.(interac.gt.0).and.(modify.ne.1)) interac=-1*interac
    if((faster==-1).and.(buttt.eq.1).and.(interac.lt.0).and.(modify.ne.1)) interac=-1*interac
    if((faster==-1).and.(buttt.eq.2)) interac=0d0
    if(quot.ne.0) interac=quot
    if(deleteif==1) interac=0d0
    if(faster==-1)then 

      m=gtk_spin_button_get_value(spinButtonc4)
      mm=gtk_spin_button_get_value(spinButtonc6)
      call undo_redo(5,m,mm,0,0,gen(m)%w(mm),interac)

      gen(which)%w(whichsaves)=interac
    end if

    spinButtonc7 = gtk_spin_button_new (gtk_adjustment_new(interac,-10d6,+10d6,&
       & 0.1d0,0.5d0,0d0),0.05d0, 7_c_int)
    tablec = gtk_table_new (7_c_int, 2_c_int, TRUE)
    buttonc = gtk_button_new_with_mnemonic ("Accept Changes"//c_null_char)

    call g_signal_connect (buttonc, "clicked"//c_null_char, c_funloc(my_accept))
    call gtk_table_attach_defaults(tablec, labelc1, 0_c_int, 2_c_int, 0_c_int, 1_c_int)
    call gtk_table_attach_defaults(tablec, labelc2, 0_c_int, 2_c_int, 2_c_int, 3_c_int)
    call gtk_table_attach_defaults(tablec, labelc3, 0_c_int, 2_c_int, 4_c_int, 5_c_int)
    call gtk_table_attach_defaults(tablec, spinButtonc4, 0_c_int, 2_c_int, 1_c_int, 2_c_int)
    call gtk_table_attach_defaults(tablec, spinButtonc6, 0_c_int, 2_c_int, 3_c_int, 4_c_int)
    call gtk_table_attach_defaults(tablec, spinButtonc7, 0_c_int, 2_c_int, 5_c_int, 6_c_int)
    call gtk_table_attach_defaults(tablec, buttonc, 0_c_int, 2_c_int, 6_c_int, 7_c_int)
    call gtk_container_add (my_choice_window, tablec)

  else if(cuals==1)then
    my_choice_window = gtk_window_new (GTK_WINDOW_TOPLEVEL)
    yetthere1=1
    call gtk_window_set_default_size(my_choice_window, 150, 150)
    call gtk_window_set_title(my_choice_window, "Relationships of Forms"//c_null_char)
    labelc1 = gtk_label_new("Gen 1"//c_null_char)
    spinButtoncc1 = gtk_spin_button_new (gtk_adjustment_new(whichsaver,0d0,ngr,&
       & 1.0d0,0.5d0,0d0),0.05d0, 0_c_int)
    labelc2 = gtk_label_new("Gen 2"//c_null_char)
    spinButtoncc2 = gtk_spin_button_new (gtk_adjustment_new(whichr,0d0,ngr,&
       & 1.0d0,0.5d0,0d0),0.05d0, 0_c_int)
    buttonc = gtk_button_new_with_mnemonic ("Draw an Arrow"//c_null_char)
    buttond = gtk_button_new_with_mnemonic ("Undraw the Arrow"//c_null_char)
    tablec = gtk_table_new (6_c_int, 2_c_int, TRUE)

    call g_signal_connect (buttonc, "clicked"//c_null_char, c_funloc(my_v_accept))
    call g_signal_connect (buttond, "clicked"//c_null_char, c_funloc(my_v_delete))
    call gtk_table_attach_defaults(tablec, labelc1, 0_c_int, 2_c_int, 0_c_int, 1_c_int)
    call gtk_table_attach_defaults(tablec, labelc2, 0_c_int, 2_c_int, 2_c_int, 3_c_int)
    call gtk_table_attach_defaults(tablec, spinButtoncc1, 0_c_int, 2_c_int, 1_c_int, 2_c_int)
    call gtk_table_attach_defaults(tablec, spinButtoncc2, 0_c_int, 2_c_int, 3_c_int, 4_c_int)
    call gtk_table_attach_defaults(tablec, buttonc, 0_c_int, 2_c_int, 4_c_int, 5_c_int)
    call gtk_table_attach_defaults(tablec, buttond, 0_c_int, 2_c_int, 5_c_int, 6_c_int)

    call gtk_container_add (my_choice_window, tablec)
    if(faster==-1) call my_v_accept_s

  else if(cuals==4)then

    my_choice_window = gtk_window_new (GTK_WINDOW_TOPLEVEL)
    yetthere1=1
    call gtk_window_set_default_size(my_choice_window, 150, 250)
    call gtk_window_set_title(my_choice_window, "Modifying adhesions"//c_null_char)

    tablec = gtk_table_new (7_c_int, 2_c_int, TRUE)
    buttonc = gtk_button_new_with_mnemonic ("Accept Changes"//c_null_char)
    buttond = gtk_button_new_with_mnemonic ("Create a new B_ij element"//c_null_char)

if(ntipusadh==0)then  ! 1234

    labelc1 = gtk_label_new("Do you want to create a B_ij element?"//c_null_char)
    call g_signal_connect (buttond, "clicked"//c_null_char, c_funloc(my_adh_add))
    call gtk_table_attach_defaults(tablec, labelc1, 0_c_int, 2_c_int, 0_c_int, 1_c_int)
    call gtk_table_attach_defaults(tablec, buttond, 0_c_int, 2_c_int, 1_c_int, 2_c_int)

    call gtk_container_add (my_choice_window, tablec)

else !1234

if(buttt==111)then  ! GG ! DEPRECIATED

    if(whichr.le.0) whichr=whichsaver
    if(whichsaver.le.0) whichsaver=whichr
    labelc1 = gtk_label_new("Gen 1"//c_null_char)
    spinButtonc4 = gtk_spin_button_new (gtk_adjustment_new(whichsaver,0d0,ngr,&
       & 1.0d0,0.5d0,0d0),0.05d0, 0_c_int)
    labelc2 = gtk_label_new("Gen 2"//c_null_char)
    spinButtonc6 = gtk_spin_button_new (gtk_adjustment_new(whichr,0d0,ngr,&
       & 1.0d0,0.5d0,0d0),0.05d0, 0_c_int)
    labelc4 = gtk_label_new("Kadh #"//c_null_char)
    spinButtonc8 = gtk_spin_button_new (gtk_adjustment_new(gen(whichsaver)%wa(1),0d0,nkad,&
       & 1.0d0,0.5d0,0d0),0.05d0, 0_c_int)
    spinButtonc9 = gtk_spin_button_new (gtk_adjustment_new(gen(whichr)%wa(1),0d0,nkad,&
       & 1.0d0,0.5d0,0d0),0.05d0, 0_c_int)
    labelc3 = gtk_label_new("Interaction strength"//c_null_char)

    if(faster==-1)then 
      if(interac==0) interac=0.1d0
      if(quot.ne.0) interac=quot
      if(deleteif==1) interac=0d0
      m=gtk_spin_button_get_value(spinButtonc4)
      mm=gtk_spin_button_get_value(spinButtonc6)
      mm1=gen(m)%wa(1); mm2=gen(mm)%wa(1)
      call undo_redo(1,m,mm,mm1,mm2,kadh(mm1,mm2),kadh(mm2,mm1))
      gen(which)%w(whichsaves)=interac
      call update_view
    end if

    spinButtonc7 = gtk_spin_button_new (gtk_adjustment_new(interac,-10d6,+10d6,&
       & 0.1d0,0.5d0,0d0),0.05d0, 7_c_int)
    call gtk_table_attach_defaults(tablec, spinButtonc4, 0_c_int, 1_c_int, 1_c_int, 2_c_int)
    call gtk_table_attach_defaults(tablec, spinButtonc6, 1_c_int, 2_c_int, 1_c_int, 2_c_int)
    call gtk_table_attach_defaults(tablec, labelc1, 0_c_int, 1_c_int, 0_c_int, 1_c_int)
    call gtk_table_attach_defaults(tablec, labelc2, 1_c_int, 2_c_int, 0_c_int, 1_c_int)
    call g_signal_connect (buttonc, "clicked"//c_null_char, c_funloc(my_adh_accept))

else if(buttt==2)then ! A-G

    labelc1 = gtk_label_new("Regulatory Molecule"//c_null_char)
    spinButtonc4 = gtk_spin_button_new (gtk_adjustment_new(whichsaver,whichsaver,whichsaver,&
       & 0.0d0,0.5d0,0d0),0.05d0, 0_c_int)

    labelc4 = gtk_label_new("Bij (Adhesion type)"//c_null_char)
    spinButtonc8 = gtk_spin_button_new (gtk_adjustment_new(whichr,0d0,nkad,&
       & 1.0d0,0.5d0,0d0),0.05d0, 0_c_int)

    labelc3 = gtk_label_new("Interaction strength"//c_null_char)
    interac=kadh(whichr,whichr)
    
    if(faster==-1)then 
      if(interac==0) interac=0.1d0
      if(quot.ne.0) interac=quot
      if(deleteif==1) interac=0d0
      m=gtk_spin_button_get_value(spinButtonc4)
      mm=gtk_spin_button_get_value(spinButtonc4)!spinButtonc6
      mm1=gen(m)%wa(1); mm2=gen(mm)%wa(1)
      call undo_redo(1,m,mm,mm1,mm2,kadh(mm1,mm2),kadh(mm2,mm1))
      gen(whichsaves)%wa(1)=which
      if(deleteif==1) gen(whichsaves)%wa(1)=0
      call update_view
    end if
    spinButtonc7 = gtk_spin_button_new (gtk_adjustment_new(interac,-10d6,+10d6,&
       & 0.1d0,0.5d0,0d0),0.05d0, 7_c_int)
    call gtk_table_attach_defaults(tablec, spinButtonc4, 0_c_int, 1_c_int, 1_c_int, 2_c_int)
    call gtk_table_attach_defaults(tablec, labelc1, 0_c_int, 1_c_int, 0_c_int, 1_c_int)
    call g_signal_connect (buttonc, "clicked"//c_null_char, c_funloc(my_adh_accept))

else if(buttt==3)then ! A-A

    labelc4 = gtk_label_new("Bij (Adhesion type)"//c_null_char)
    spinButtonc8 = gtk_spin_button_new (gtk_adjustment_new(whichsaver,0d0,nkad,&
       & 1.0d0,0.5d0,0d0),0.05d0, 0_c_int)
    spinButtonc9 = gtk_spin_button_new (gtk_adjustment_new(whichr,0d0,nkad,&
       & 1.0d0,0.5d0,0d0),0.05d0, 0_c_int)
    labelc3 = gtk_label_new("Interaction strength"//c_null_char)
    interac=kadh(whichsaver,whichr)
    if(faster==-1)then
      if(interac==0) interac=0.1d0 
      if(quot.ne.0) interac=quot
      if(deleteif==1) interac=0d0
      m=0;mm=0
      do mm1=1,ng
        if(gen(mm1)%wa(1)==int(whichsaver))then
        m=mm1
        exit
        endif
      end do
      do mm2=1,ng
        if(gen(mm2)%wa(1)==int(whichr))then
        mm=mm2
        exit
        endif
      end do
      mm3=whichsaver; mm4=which
      call undo_redo(1,m,mm,mm3,mm4,kadh(whichsaver,whichr),kadh(whichr,whichsaver))

      kadh(whichsaver,whichr)=interac
      kadh(whichr,whichsaver)=interac

      call update_view
    end if
    spinButtonc7 = gtk_spin_button_new (gtk_adjustment_new(interac,-10d6,+10d6,&
       & 0.1d0,0.5d0,0d0),0.05d0, 7_c_int)
    call gtk_table_attach_defaults(tablec, spinButtonc9, 1_c_int, 2_c_int, 3_c_int, 4_c_int)
   call g_signal_connect (buttonc, "clicked"//c_null_char, c_funloc(my_adh_accept2))

  end if

    call g_signal_connect (buttond, "clicked"//c_null_char, c_funloc(my_adh_add))
    call gtk_table_attach_defaults(tablec, labelc3, 0_c_int, 2_c_int, 4_c_int, 5_c_int)
    call gtk_table_attach_defaults(tablec, labelc4, 0_c_int, 1_c_int, 2_c_int, 3_c_int)
    call gtk_table_attach_defaults(tablec, labelc4, 1_c_int, 2_c_int, 2_c_int, 3_c_int)

    call gtk_table_attach_defaults(tablec, spinButtonc8, 0_c_int, 1_c_int, 3_c_int, 4_c_int)
    call gtk_table_attach_defaults(tablec, spinButtonc7, 0_c_int, 2_c_int, 5_c_int, 6_c_int)
    call gtk_table_attach_defaults(tablec, buttonc, 0_c_int, 2_c_int, 6_c_int, 7_c_int)
    call gtk_table_attach_defaults(tablec, buttond, 0_c_int, 2_c_int, 7_c_int, 8_c_int)

    call gtk_container_add (my_choice_window, tablec)
 
  end if !1234

  end if

  if(cuals.ne.1) call gtk_widget_show_all (my_choice_window)

  call update_view

end subroutine choice_box_build

subroutine choice_box_build_ww(whichsave,whichsave2,which,cualw)

use global_widgets

  type(c_ptr) :: boxc,checkd
  type(c_ptr) :: labelc1, labelc2, labelc3, labelc4, buttonc, buttond, tablec
  integer :: choice,which,whichsave,whichsave2,cualw,m,mm,mmm,mmmm
  real*8 :: whichsaver,whichsave2r,whichr,ngr,interac,interac2

  choiceboxthere=1
  if(yetthere3==1) call gtk_widget_destroy(my_choice_window3)

  whichsaver=whichsave; whichsave2r=whichsave2; whichr=which; ngr=ng

  if(cualw==1)then
    my_choice_window3 = gtk_window_new (GTK_WINDOW_TOPLEVEL)
      yetthere3=1
    call gtk_window_set_default_size(my_choice_window3, 150, 250)
    call gtk_window_set_title(my_choice_window3, "Catalytic interactions"//c_null_char)
    labelc1 = gtk_label_new("Pre-form"//c_null_char)
    spinButtonc1 = gtk_spin_button_new (gtk_adjustment_new(whichsaver,whichsaver,whichsaver,&
       & 1.0d0,0.5d0,0d0),0.05d0, 0_c_int)
    labelc2 = gtk_label_new("Post-form"//c_null_char)
    spinButtonc2 = gtk_spin_button_new (gtk_adjustment_new(whichsave2r,whichsave2r,whichsave2r,&
       & 1.0d0,0.5d0,0d0),0.05d0, 0_c_int)
    spinButtonc3 = gtk_spin_button_new (gtk_adjustment_new(whichr,whichr,whichr,&
       & 1.0d0,0.5d0,0d0),0.05d0, 0_c_int)
    labelc3 = gtk_label_new("Catalysed by"//c_null_char)

    if(gen(which)%nww.ge.1)then
      mmmm=0
      do mmm=1,gen(which)%nww
        if(gen(which)%ww(mmm,1)==whichsave)then

          mmmm=1
          interac=gen(which)%ww(mmm,3)
          if(abs(interac).lt.0.0000001) interac=0.1d0
          if(quot.ne.0) interac=quot
          if(deleteif==1) interac=0d0
          interac2=interac
    spinButtonc4 = gtk_spin_button_new (gtk_adjustment_new(interac2,-10d6,+10d6,&
       & 0.1d0,0.5d0,0d0),0.05d0, 7_c_int)
          if(faster==-1)then
          m=gtk_spin_button_get_value(spinButtonc1)
          mm=gtk_spin_button_get_value(spinButtonc2)
            gen(which)%ww(mmm,3)=interac

          endif
        end if
      end do
      elseif(mmmm.ne.0)then

        interac=0.1d0
        if(quot.ne.0) interac=quot
        if(deleteif==1) interac=0d0
        interac2=interac
    spinButtonc4 = gtk_spin_button_new (gtk_adjustment_new(interac2,-10d6,+10d6,&
       & 0.1d0,0.5d0,0d0),0.05d0, 7_c_int)
        m=gtk_spin_button_get_value(spinButtonc1)
          mm=gtk_spin_button_get_value(spinButtonc2)
        if(faster==-1) call my_ww_accept_s
        if(quot.ne.0) interac=quot
        if(deleteif==1) interac=0d0
        if(faster==-1) gen(which)%ww(gen(which)%nww,3)=interac
    else
        interac=0.1d0
        if(quot.ne.0) interac=quot
        if(deleteif==1) interac=0d0
        interac2=interac
            spinButtonc4 = gtk_spin_button_new (gtk_adjustment_new(interac2,-10d6,+10d6,&
       & 0.1d0,0.5d0,0d0),0.05d0, 7_c_int)
        if(faster==-1) gen(which)%ww(1,3)=interac
        if(faster==-1)then
          m=gtk_spin_button_get_value(spinButtonc1)
          mm=gtk_spin_button_get_value(spinButtonc2)

          if(faster==-1) call my_ww_accept_s

          gen(which)%ww(1,3)=interac
        end if
      end if

    if(interac2.lt.0.000001) interac2=0.1d0
    if(quot.ne.0) interac=quot
    if(deleteif==1) interac=0d0
        if(quot.ne.0) interac2=quot
        if(deleteif==1) interac2=0d0
    spinButtonc3 = gtk_spin_button_new (gtk_adjustment_new(whichr,whichr,whichr,&
       & 1.0d0,0.5d0,0d0),0.05d0, 0_c_int)
    labelc4 = gtk_label_new("Catalytic activity"//c_null_char)
    tablec = gtk_table_new (7_c_int, 2_c_int, TRUE)
    spinButtonc4 = gtk_spin_button_new (gtk_adjustment_new(interac2,-10d6,+10d6,&
       & 0.1d0,0.5d0,0d0),0.05d0, 7_c_int)
    buttonc = gtk_button_new_with_mnemonic ("Accept Changes"//c_null_char)

    call g_signal_connect (buttonc, "clicked"//c_null_char, c_funloc(my_ww_accept))
    call gtk_table_attach_defaults(tablec, labelc1, 0_c_int, 2_c_int, 0_c_int, 1_c_int)
    call gtk_table_attach_defaults(tablec, labelc2, 0_c_int, 2_c_int, 2_c_int, 3_c_int)
    call gtk_table_attach_defaults(tablec, labelc3, 0_c_int, 2_c_int, 4_c_int, 5_c_int)
    call gtk_table_attach_defaults(tablec, labelc4, 0_c_int, 2_c_int, 6_c_int, 7_c_int)
    call gtk_table_attach_defaults(tablec, spinButtonc1, 0_c_int, 2_c_int, 1_c_int, 2_c_int)
    call gtk_table_attach_defaults(tablec, spinButtonc2, 0_c_int, 2_c_int, 3_c_int, 4_c_int)
    call gtk_table_attach_defaults(tablec, spinButtonc3, 0_c_int, 2_c_int, 5_c_int, 6_c_int)
    call gtk_table_attach_defaults(tablec, spinButtonc4, 0_c_int, 2_c_int, 7_c_int, 8_c_int)
    call gtk_table_attach_defaults(tablec, buttonc, 0_c_int, 2_c_int, 8_c_int, 9_c_int)
    call gtk_container_add (my_choice_window3, tablec)

  end if

  if(faster==-1) call my_ww_accept_s
  call gtk_widget_show_all (my_choice_window3)

end subroutine choice_box_build_ww

subroutine choice_box_param_build(which,par,cuals,buttt)

  use gtk, only: gtk_notebook_prev_page, gtk_notebook_next_page, gtk_notebook_get_n_pages, gtk_notebook_get_current_page
  use global_widgets

  type(c_ptr) :: boxc,checkd
  type(c_ptr) :: labelc1, labelc2, labelc3, labelc4, buttonc, buttond, tablec
  integer :: choice,which,whichsaves,cuals,par,m,mm,buttt
  real*8 :: parr,whichr,ngr,interac,interaco
  character*24 strs

  parr=par; whichr=which; ngr=ng
  choiceboxthere=1
  if(yetthere==1) call gtk_widget_destroy(my_choice_window2)
  my_choice_window2 = gtk_window_new (GTK_WINDOW_TOPLEVEL)
  yetthere=1
  call gtk_window_set_default_size(my_choice_window2, 150, 250)
  if(cuals==2) call gtk_window_set_title(my_choice_window2, "Node Property Regulation"//c_null_char)
  if(cuals==3) call gtk_window_set_title(my_choice_window2, "Cell behaviour Regulation"//c_null_char)
  labelc1 = gtk_label_new("Gene"//c_null_char)
  spinButtonc1 = gtk_spin_button_new (gtk_adjustment_new(whichr,whichr,whichr,&
       & 1.0d0,0.5d0,0d0),0.05d0, 0_c_int)
  if(cuals==2)then
    labelc2 = gtk_label_new("Node Parameter"//c_null_char)
    labelc4 = gtk_label_new(trim(nodeparams(parr))//c_null_char)
    interac=gen(which)%wa(parr)
    spinButtonc2 = gtk_spin_button_new (gtk_adjustment_new(parr,parr,parr,&
       & 1.0d0,0.5d0,0d0),0.05d0, 0_c_int)
  else
    labelc2 = gtk_label_new("Cell Parameter"//c_null_char)
    labelc4 = gtk_label_new(trim(cellparams(parr))//c_null_char)
    interac=gen(which)%wa(parr+nparam_per_node)
    spinButtonc2 = gtk_spin_button_new (gtk_adjustment_new(parr,parr,parr,&
       & 1.0d0,0.5d0,0d0),0.05d0, 0_c_int)
  endif
  if(abs(interac).lt.0.0000001)then
    interac=0.1d0
    if(buttt==3) interac=-0.1d0
    if(buttt==2) interac=0d0
    if(quot.ne.0) interac=quot
    if(deleteif==1) interac=0d0
  end if
  if(buttt==2) interac=0d0
  if(deleteif==1) interac=0d0
  if((faster==-1).and.(cuals==2))then
    interaco=gen(which)%wa(parr)
    gen(which)%wa(parr)=interac
    m=gtk_spin_button_get_value(spinButtonc1)
    mm=gtk_spin_button_get_value(spinButtonc2)
    call undo_redo(7,m,mm,0,0,interaco,interac)
  end if
  if((faster==-1).and.(cuals==3))then
    interaco=gen(which)%wa(parr+nparam_per_node)
    gen(which)%wa(parr+nparam_per_node)=interac
    m=gtk_spin_button_get_value(spinButtonc1)
    mm=gtk_spin_button_get_value(spinButtonc2)+nparam_per_node
    call undo_redo(7,m,mm,0,0,interaco,interac)
  end if
  labelc3 = gtk_label_new("Interaction strength"//c_null_char)
  spinButtonc3 = gtk_spin_button_new (gtk_adjustment_new(interac,-10d6,+10d6,&
       & 0.1d0,0.5d0,0d0),0.05d0, 7_c_int)
  tablec = gtk_table_new (7_c_int, 2_c_int, TRUE)
  buttonc = gtk_button_new_with_mnemonic ("Accept Changes"//c_null_char)
  if(cuals==2) call g_signal_connect (buttonc, "clicked"//c_null_char, c_funloc(my_wa_accept))
  if(cuals==3) call g_signal_connect (buttonc, "clicked"//c_null_char, c_funloc(my_wa2_accept))

  call gtk_table_attach_defaults(tablec, labelc1, 0_c_int, 2_c_int, 0_c_int, 1_c_int)
  call gtk_table_attach_defaults(tablec, labelc2, 0_c_int, 2_c_int, 2_c_int, 3_c_int)
  call gtk_table_attach_defaults(tablec, labelc3, 0_c_int, 2_c_int, 5_c_int, 6_c_int)
  call gtk_table_attach_defaults(tablec, labelc4, 0_c_int, 2_c_int, 3_c_int, 4_c_int)
  call gtk_table_attach_defaults(tablec, spinButtonc1, 0_c_int, 2_c_int, 1_c_int, 2_c_int)
  call gtk_table_attach_defaults(tablec, spinButtonc2, 0_c_int, 2_c_int, 4_c_int, 5_c_int)
  call gtk_table_attach_defaults(tablec, spinButtonc3, 0_c_int, 2_c_int, 6_c_int, 7_c_int)
  call gtk_table_attach_defaults(tablec, buttonc, 0_c_int, 2_c_int, 7_c_int, 8_c_int)

  call gtk_container_add (my_choice_window2, tablec)

  call gtk_widget_show_all (my_choice_window2)

  call update_view

end subroutine choice_box_param_build

function my_accept(widget, event, gdata) result(ret)  bind(c)

  use gtk, only: gtk_notebook_prev_page, gtk_notebook_next_page, gtk_notebook_get_n_pages, gtk_notebook_get_current_page
  use global_widgets
  use genetic

  type(c_ptr), value :: widget, gdata, event

  integer(c_int)    :: ret
  integer :: m, mm, mmm

  m=gtk_spin_button_get_value(spinButtonc4)
  mm=gtk_spin_button_get_value(spinButtonc6)

  call undo_redo(5,m,mm,0,0,gen(m)%w(mm),gtk_spin_button_get_value(spinButtonc7))

  gen(mm)%w(m)= gtk_spin_button_get_value(spinButtonc7)

  call gtk_widget_destroy(my_choice_window)
  choiceboxthere=0
  call update_npag
  call update_view

end function my_accept

subroutine my_accept_s

  use gtk, only: gtk_notebook_prev_page, gtk_notebook_next_page, gtk_notebook_get_n_pages, gtk_notebook_get_current_page
  use global_widgets
  use genetic

  integer :: m, mm, mmm

if(choiceboxthere==1)then
  m=gtk_spin_button_get_value(spinButtonc4)
  mm=gtk_spin_button_get_value(spinButtonc6)
  call undo_redo(5,m,mm,0,0,gen(m)%w(mm),gtk_spin_button_get_value(spinButtonc7))
  gen(mm)%w(m)= gtk_spin_button_get_value(spinButtonc7)
  choiceboxthere=0
end if
  call update_npag
  call update_view

end subroutine my_accept_s

function my_wa_accept(widget, event, gdata) result(ret)  bind(c)

  use gtk, only: gtk_notebook_prev_page, gtk_notebook_next_page, gtk_notebook_get_n_pages, gtk_notebook_get_current_page
  use global_widgets

  type(c_ptr), value :: widget, gdata, event
  integer(c_int)    :: ret
  integer :: m, mm

  call gtk_widget_destroy(my_choice_window)
  m=gtk_spin_button_get_value(spinButtonc1)
  mm=gtk_spin_button_get_value(spinButtonc2)
  call undo_redo(7,m,mm,0,0,gen(m)%wa(mm),gtk_spin_button_get_value(spinButtonc3))
  gen(m)%wa(mm)= gtk_spin_button_get_value(spinButtonc3)

  call gtk_widget_destroy(my_choice_window)

  call update_view

  call gtk_widget_destroy(my_choice_window2)

end function my_wa_accept

subroutine my_wa_accept_s

  use gtk, only: gtk_notebook_prev_page, gtk_notebook_next_page, gtk_notebook_get_n_pages, gtk_notebook_get_current_page
  use global_widgets

  integer :: m, mm

  call gtk_widget_destroy(my_choice_window)
  m=gtk_spin_button_get_value(spinButtonc1)
  mm=gtk_spin_button_get_value(spinButtonc2)
  call undo_redo(7,m,mm,0,0,gen(m)%wa(mm),gtk_spin_button_get_value(spinButtonc3))
  gen(m)%wa(mm)= gtk_spin_button_get_value(spinButtonc3)

  call gtk_widget_destroy(my_choice_window)

  call update_view

  call gtk_widget_destroy(my_choice_window2)

end subroutine my_wa_accept_s

function my_wa2_accept(widget, event, gdata) result(ret)  bind(c)

  use gtk, only: gtk_notebook_prev_page, gtk_notebook_next_page, gtk_notebook_get_n_pages, gtk_notebook_get_current_page
  use global_widgets

  type(c_ptr), value :: widget, gdata, event
  integer(c_int)    :: ret
  integer :: m, mm

  m=gtk_spin_button_get_value(spinButtonc1)
  mm=gtk_spin_button_get_value(spinButtonc2)+nparam_per_node
  call undo_redo(7,m,mm,0,0,gen(m)%wa(mm),gtk_spin_button_get_value(spinButtonc3))

  gen(m)%wa(mm)= gtk_spin_button_get_value(spinButtonc3)

  call update_view

  call gtk_widget_destroy(my_choice_window2)

end function my_wa2_accept

subroutine my_wa2_accept_s

  use gtk, only: gtk_notebook_prev_page, gtk_notebook_next_page, gtk_notebook_get_n_pages, gtk_notebook_get_current_page
  use global_widgets

  integer :: m, mm

  m=gtk_spin_button_get_value(spinButtonc1)
  mm=gtk_spin_button_get_value(spinButtonc2)+nparam_per_node
  call undo_redo(7,m,mm,0,0,gen(m)%wa(mm),gtk_spin_button_get_value(spinButtonc3))

  gen(m)%wa(mm)= gtk_spin_button_get_value(spinButtonc3)

  call update_view

  call gtk_widget_destroy(my_choice_window2)

end subroutine my_wa2_accept_s

function my_ww_accept(widget, event, gdata) result(ret)  bind(c)

  use global_widgets

  type(c_ptr), value :: widget, gdata, event
  integer(c_int)    :: ret
  integer :: ma, mp, mk, m, mm, mn
  real*8 :: mv
  real*8 :: gesww(ng*ng,3)
    
  ma=gtk_spin_button_get_value(spinButtonc1)
  mp=gtk_spin_button_get_value(spinButtonc2)
  mk=gtk_spin_button_get_value(spinButtonc3)
  mv=gtk_spin_button_get_value(spinButtonc4)
  mn=0
  do m=1,gen(ma)%npost
    if(gen(ma)%post(m)==mp) mn=1
  end do
  if(mn==0)then
    do m=1,gen(ma)%npre
      if(gen(ma)%pre(m)==mp) mn=2
    end do
  end if
  if(mn==2)then; mn=ma; ma=mp; mp=mn
  else if(mn==0)then; print*, "These forms do not give rise to another"
  endif
    
  if(mn.ne.0)then
    mn=3
    do m=1,gen(mk)%nww
      if((gen(mk)%ww(m,1)==ma).and.(gen(mk)%ww(m,2)==mp))then ! easiest case: only mv changed
        call undo_redo(6,ma,mp,mk,0,gen(mk)%ww(m,3),mv) 
        gen(mk)%ww(m,3)=mv; mn=4; exit
      end if
    end do
  end if

    call undo_redo(6,ma,mp,mk,0,0d0,mv)     
    gesww=0d0; geswn=0
    do mm=1,gen(mk)%nww
      gesww(mm,:)=gen(mk)%ww(mm,:)
    end do
    gen(mk)%nww=gen(mk)%nww+1
    if(allocated(gen(mk)%ww)) deallocate(gen(mk)%ww)
    allocate(gen(mk)%ww(gen(mk)%nww,3))
    do mm=1,gen(mk)%nww-1
      gen(mk)%ww(mm,:)=gesww(mm,:)
    end do

    gen(mk)%ww(gen(mk)%nww,1)=ma
    gen(mk)%ww(gen(mk)%nww,2)=mp
    gen(mk)%ww(gen(mk)%nww,3)=mv

  call update_view
  call gtk_widget_destroy(my_choice_window3)
  yetthere3=0

end function my_ww_accept

subroutine my_ww_accept_s

  use gtk, only: gtk_notebook_prev_page, gtk_notebook_next_page, gtk_notebook_get_n_pages, gtk_notebook_get_current_page
  use global_widgets

  integer :: ma, mp, mk, m, mm, mn
  real*8 :: mv
  real*8 :: gesww(ng*ng,3)

  ma=gtk_spin_button_get_value(spinButtonc1)

  mp=gtk_spin_button_get_value(spinButtonc2)

  mk=gtk_spin_button_get_value(spinButtonc3)

  mv=gtk_spin_button_get_value(spinButtonc4)

   
  mn=0
  do m=1,gen(ma)%npost
    if(gen(ma)%post(m)==mp) mn=1
  end do
  if(mn==0)then
    do m=1,gen(ma)%npre
      if(gen(ma)%pre(m)==mp) mn=2
    end do
  end if
  if(mn==2)then; mn=ma; ma=mp; mp=mn
  else if(mn==0)then; print*, "These forms do not give rise to another"
  endif
    
  if(mn.ne.0)then
    mn=3
    do m=1,gen(mk)%nww
      if((gen(mk)%ww(m,1)==ma).and.(gen(mk)%ww(m,2)==mp))then ! easiest case: only mv changed
        call undo_redo(6,ma,mp,mk,0,gen(mk)%ww(m,3),mv) 
        gen(mk)%ww(m,3)=mv; mn=4; exit
      end if
    end do
  end if
  if(mn==3)then ! Draw a new ww. More fun. 
    call undo_redo(6,ma,mp,mk,0,0d0,mv)      
    gesww=0d0; geswn=0
    do mm=1,gen(mk)%nww
      gesww(mm,:)=gen(mk)%ww(mm,:)
    end do
    gen(mk)%nww=gen(mk)%nww+1
    if(allocated(gen(mk)%ww)) deallocate(gen(mk)%ww)
    allocate(gen(mk)%ww(gen(mk)%nww,3))
    do mm=1,gen(mk)%nww-1
      gen(mk)%ww(mm,:)=gesww(mm,:)
    end do
    gen(mk)%ww(gen(mk)%nww,1)=ma
    gen(mk)%ww(gen(mk)%nww,2)=mp
    gen(mk)%ww(gen(mk)%nww,3)=mv
  end if    

  call update_view
        

end subroutine my_ww_accept_s


function my_v_accept(widget, event, gdata) result(ret)  bind(c)

  use gtk, only: gtk_notebook_prev_page, gtk_notebook_next_page, gtk_notebook_get_n_pages, gtk_notebook_get_current_page
  use global_widgets

  type(c_ptr), value :: widget, gdata, event
  integer(c_int)    :: ret
  integer :: mn, m, mm, npo, npr, wpo, wpr, gf
  integer :: ges3(ng,ng,2)

  m=gtk_spin_button_get_value(spinButtoncc1)
  mm=gtk_spin_button_get_value(spinButtoncc2)
  wpo=0; wpr=0; ges3=0
  call undo_redo(4,m,mm,1,0,0d0,0d0) 
  do mn=1,gen(m)%npost
    if(gen(m)%post(mn)==mm) wpo=mn
  end do
  do mn=1,gen(mm)%npre
    if(gen(mm)%pre(mn)==m) wpr=mn
  end do
  if((wpr==0).and.(wpo==0))then
    if(gen(m)%kindof==1)then
      gen(m)%kindof=2
      print*, "ATTENTION! The kindof of gene ", m, "has been changed to 2."
    end if
    if(gen(mm)%kindof<3)then
      gen(mm)%kindof=3
      print*, "ATTENTION! The kindof of gene ", mm, "has been changed to 3."
    end if
      
    gen(m)%npost=gen(m)%npost+1; gen(mm)%npre=gen(mm)%npre+1
    npo=gen(m)%npost; npr=gen(mm)%npre
    if(gen(mm)%npre.ge.1)then
      do gf=1,gen(mm)%npre-1
        ges3(mm,gf,1)=gen(mm)%pre(gf)
      end do
    end if
    if(gen(m)%npost.ge.1)then
      do gf=1,gen(m)%npost-1
        ges3(m,gf,2)=gen(m)%post(gf)
      end do
    end if

    if(allocated(gen(m)%post)) deallocate(gen(m)%post)
    allocate(gen(m)%post(npo))
    if(allocated(gen(mm)%pre)) deallocate(gen(mm)%pre)
    allocate(gen(mm)%pre(npr))

    if(gen(mm)%npre.ge.1)then
      do gf=1,gen(mm)%npre-1
         gen(mm)%pre(gf)=ges3(mm,gf,1)
      end do
    end if
    if(gen(m)%npost.ge.1)then
      do gf=1,gen(m)%npost-1
        gen(m)%post(gf)=ges3(m,gf,2)
      end do
    end if

    gen(m)%post(npo)=mm
    gen(mm)%pre(npr)=m

    print*,  mm, " is now generated by ", m, "."
    print*, gen(m)%npost, gen(mm)%npre

  else
    print*, "The arrow has already been drawn."
  end if

  call update_view

end function my_v_accept

function my_v_accept_n(widget, event, gdata) result(ret)  bind(c)

  use gtk, only: gtk_notebook_prev_page, gtk_notebook_next_page, gtk_notebook_get_n_pages, gtk_notebook_get_current_page
  use global_widgets

  type(c_ptr), value :: widget, gdata, event
  integer(c_int)    :: ret
  integer :: mn, m, mm, npo, npr, wpo, wpr, gf
  integer :: ges3(ng,ng,2)

if(open_tab==0)then

elseif(open_tab==1)then
  m=gtk_spin_button_get_value(spinButtoncc1)
  mm=gtk_spin_button_get_value(spinButtoncc2)
  wpo=0; wpr=0; ges3=0
  call undo_redo(4,m,mm,1,0,0d0,0d0) 
  do mn=1,gen(m)%npost
    if(gen(m)%post(mn)==mm) wpo=mn
  end do
  do mn=1,gen(mm)%npre
    if(gen(mm)%pre(mn)==m) wpr=mn
  end do
  if((wpr==0).and.(wpo==0))then
    if(gen(m)%kindof==1)then
      gen(m)%kindof=2
      print*, "ATTENTION! The kindof of gene ", m, "has been changed to 2."
    end if
    if(gen(mm)%kindof<3)then
      gen(mm)%kindof=3
      print*, "ATTENTION! The kindof of gene ", mm, "has been changed to 3."
    end if
      
    gen(m)%npost=gen(m)%npost+1; gen(mm)%npre=gen(mm)%npre+1
    npo=gen(m)%npost; npr=gen(mm)%npre
    if(gen(mm)%npre.ge.1)then
      do gf=1,gen(mm)%npre-1
        ges3(mm,gf,1)=gen(mm)%pre(gf)
      end do
    end if
    if(gen(m)%npost.ge.1)then
      do gf=1,gen(m)%npost-1
        ges3(m,gf,2)=gen(m)%post(gf)
      end do
    end if

    if(allocated(gen(m)%post)) deallocate(gen(m)%post)
    allocate(gen(m)%post(npo))
    if(allocated(gen(mm)%pre)) deallocate(gen(mm)%pre)
    allocate(gen(mm)%pre(npr))

    if(gen(mm)%npre.ge.1)then
      do gf=1,gen(mm)%npre-1
         gen(mm)%pre(gf)=ges3(mm,gf,1)
      end do
    end if
    if(gen(m)%npost.ge.1)then
      do gf=1,gen(m)%npost-1
        gen(m)%post(gf)=ges3(m,gf,2)
      end do
    end if

    gen(m)%post(npo)=mm
    gen(mm)%pre(npr)=m

    print*, mm, " is now generated by ", m, "."
    print*, gen(m)%npost, gen(mm)%npre

  else
    print*, "The arrow has already been drawn."
  end if
endif

  call update_view

end function my_v_accept_n


subroutine my_v_accept_s

  use gtk, only: gtk_notebook_prev_page, gtk_notebook_next_page, gtk_notebook_get_n_pages, gtk_notebook_get_current_page
  use global_widgets

  integer :: mn, m, mm, npo, npr, wpo, wpr, gf
  integer :: ges3(ng,ng,2)

  m=gtk_spin_button_get_value(spinButtoncc1)
  mm=gtk_spin_button_get_value(spinButtoncc2)
  wpo=0; wpr=0; ges3=0
  call undo_redo(4,m,mm,1,0,0d0,0d0) 
  do mn=1,gen(m)%npost
    if(gen(m)%post(mn)==mm) wpo=mn
  end do
  do mn=1,gen(mm)%npre
    if(gen(mm)%pre(mn)==m) wpr=mn
  end do
  if((wpr==0).and.(wpo==0))then
    if(gen(m)%kindof==1)then
      gen(m)%kindof=2
      print*, "ATTENTION! The kindof of gene ", m, "has been changed to 2."
    end if
    if(gen(mm)%kindof<3)then
      gen(mm)%kindof=3
      print*, "ATTENTION! The kindof of gene ", mm, "has been changed to 3."
    end if
      
    gen(m)%npost=gen(m)%npost+1; gen(mm)%npre=gen(mm)%npre+1
    npo=gen(m)%npost; npr=gen(mm)%npre
    if(gen(mm)%npre.ge.1)then
      do gf=1,gen(mm)%npre-1
         ges3(mm,gf,1)=gen(mm)%pre(gf)
      end do
    end if
    if(gen(m)%npost.ge.1)then
      do gf=1,gen(m)%npost-1
        ges3(m,gf,2)=gen(m)%post(gf)
      end do
    end if

    if(allocated(gen(m)%post)) deallocate(gen(m)%post)
    allocate(gen(m)%post(npo))
    if(allocated(gen(mm)%pre)) deallocate(gen(mm)%pre)
    allocate(gen(mm)%pre(npr))

    if(gen(mm)%npre.ge.1)then
      do gf=1,gen(mm)%npre-1
        gen(mm)%pre(gf)=ges3(mm,gf,1)
      end do
    end if
    if(gen(m)%npost.ge.1)then
      do gf=1,gen(m)%npost-1
        gen(m)%post(gf)=ges3(m,gf,2)
      end do
    end if

    gen(m)%post(npo)=mm
    gen(mm)%pre(npr)=m

    print*, mm, " is now generated by ", m, "."
    print*, gen(m)%npost, gen(mm)%npre

  else
    print*, "The arrow has already been drawn."
  end if

  call update_view

end subroutine my_v_accept_s


function my_p_accept(widget, event, gdata) result(ret)  bind(c)

  use global_widgets
  use gtk, only: gtk_notebook_prev_page, gtk_notebook_next_page, gtk_notebook_get_n_pages, gtk_notebook_get_current_page
  use basic_handlers
  use gtk_sup

  type(c_ptr), value :: widget, gdata, event
  type(c_ptr) :: genenamec
  character*12 :: genenamen
  integer(c_int)    :: ret
  integer :: m, mk
  real*8 :: mmk

  m=gtk_spin_button_get_value(spinButtonc1)
  mk=gtk_spin_button_get_value(spinButtonc2)

  genenamec=(gtk_entry_get_text(entryc))
  call c_f_string(genenamec, genenamen)

  if(len(trim(genenamen)).gt.0) genename(m)=trim(genenamen)

  if((mk==1).and.(gen(m)%npre+gen(m)%npost.ne.0))then
    print*, "Cannot change to kindof 1, as pres or posts are there."
  else if((mk==2).and.(gen(m)%npre.ne.0))then
    print*, "Cannot change to kindof 2, as pres are there."

  else
    if(gen(m)%kindof.ne.mk)then
      mmk=mk
      call undo_redo(3,m,1,0,0,gen(m)%kindof,mmk) 
      gen(m)%kindof = mk
    end if
  end if   

  if(gen(m)%diffu.ne.gtk_spin_button_get_value(spinButtonc3))then
    call undo_redo(3,m,3,0,0,gen(m)%diffu,gtk_spin_button_get_value(spinButtonc3)) 
    gen(m)%diffu = gtk_spin_button_get_value(spinButtonc3)
  end if

  if(gen(m)%mu.ne.gtk_spin_button_get_value(spinButtonc5))then
    call undo_redo(3,m,2,0,0,gen(m)%mu,gtk_spin_button_get_value(spinButtonc5)) 
    gen(m)%mu = gtk_spin_button_get_value(spinButtonc5)
  end if
   call gtk_widget_destroy(my_choice_window2)

  call update_view

end function my_p_accept

function my_names_accept(widget, event, gdata) result(ret)  bind(c)

  use global_widgets
  use gtk, only: gtk_notebook_prev_page, gtk_notebook_next_page, gtk_notebook_get_n_pages, gtk_notebook_get_current_page
  use basic_handlers
  use gtk_sup

  type(c_ptr), value :: widget, gdata, event
  type(c_ptr) :: genenamec
  character*12 :: genenamen
  integer(c_int)    :: ret

do i=1,ng
  genenamec=(gtk_entry_get_text(entryna(i)))
  call c_f_string(genenamec, genenamen)
print*, "HERE?", i, ng
print*, "HERE??", i, genenamen
  if(len(trim(genenamen)).gt.0) genename(i)=trim(genenamen)
end do

   call gtk_widget_destroy(my_choice_windown)

  call update_view

end function my_names_accept


subroutine my_v_delete_s

  use gtk, only: gtk_notebook_prev_page, gtk_notebook_next_page, gtk_notebook_get_n_pages, gtk_notebook_get_current_page
  use global_widgets
  use genetic

  integer :: mn, mnm, m, mm, mc, npo, npr, wpo, wpr, cto
  real*8 :: mr

  m=gtk_spin_button_get_value(spinButtonc1)
  mm=gtk_spin_button_get_value(spinButtonc2)
  mc=gtk_spin_button_get_value(spinButtonc3)
  mr=gtk_spin_button_get_value(spinButtonc4)
  if(deleteif==1)then
    m=aaa0
    mm=bbb0
  endif
  if((gen(m)%npost.ne.0).or.(gen(mm)%npre.ne.0))then  ! actually, this should be an .and.
    wpo=0; wpr=0
   cto=0
do mn=1,ng
  do mnm=1,gen(mn)%nww
    if((gen(mn)%ww(mnm,1)==m).and.(gen(mn)%ww(mnm,2)==mm)) cto=cto+1
  end do
end do
   if(cto.ne.0) call undo_redo(6,m,mm,mc,0,mr,0d0)  
   if(cto.le.1)then
   call undo_redo(4,m,mm,0,0,0d0,0d0)
    do mn=1,gen(m)%npost
      if(gen(m)%post(mn)==mm)then; wpo=mn; exit; end if
    end do
    do mn=1,gen(mm)%npre
      if(gen(mm)%pre(mn)==m)then; wpr=mn; exit; end if
    end do

    if((wpo.ne.0).or.(wpr.ne.0))then
      if(wpo.lt.gen(m)%npost)then
        do mn=wpo,gen(m)%npost-1
          gen(m)%post(mn)=gen(m)%post(mn+1)
        end do
      end if
      if(wpr.lt.gen(mm)%npre)then
        do mn=wpr,gen(mm)%npre-1
          gen(mm)%pre(mn)=gen(mm)%pre(mn+1)
        end do
      end if

      gen(m)%npost=gen(m)%npost-1; gen(mm)%npre=gen(mm)%npre-1
      npo=gen(m)%npost; npr=gen(mm)%npre
      do mn=1,ng
        wpo=0
        do mnm=1,gen(mn)%nww
          if((gen(mn)%ww(mnm,1)==m).and.(gen(mn)%ww(mnm,2)==mm))then
            gen(mn)%ww(mnm,3)=0d0
            wpo=1
          end if
          if(wpo==1)then
            if(mnm.lt.gen(mn)%nww)then
              gen(mn)%ww(mnm,:)=gen(mn)%ww(mnm+1,:)
            else
              gen(mn)%nww=gen(mn)%nww-1
            endif
          endif
        end do
      end do

      print*, mm, " is not any longer produced by ", m, "."
    else
      print*, "No pres or posts."
    end if
  else ! selective catalyst 
    do mn=1,gen(mc)%nww
      if((gen(mc)%ww(mn,1)==m).and.(gen(mc)%ww(mn,2)==mm))then
        gen(mc)%ww(mn,3)=0d0
        gen(mc)%ww(mn,1)=0; gen(mc)%ww(mn,2)=0
        goto 476
      endif
    end do
    476 if(mn.ne.gen(mc)%nww)then
     do mnm=mn,gen(mc)%nww-1
      gen(mc)%ww(mnm,:)=gen(mc)%ww(mnm+1,:)
     end do
    end if
    gen(mc)%nww=gen(mc)%nww-1
    endif    
  endif
call update_npag
  call update_view

end subroutine my_v_delete_s

function my_v_delete(widget, event, gdata) result(ret)  bind(c)

  use gtk, only: gtk_notebook_prev_page, gtk_notebook_next_page, gtk_notebook_get_n_pages, gtk_notebook_get_current_page
  use global_widgets

  type(c_ptr), value :: widget, gdata, event
  integer(c_int)    :: ret
  integer :: mn, mnm, m, mm, npo, npr, wpo, wpr
  real*8 :: mr

  m=gtk_spin_button_get_value(spinButtonc1)
  mm=gtk_spin_button_get_value(spinButtonc2)
  mc=gtk_spin_button_get_value(spinButtonc3)
  mr=gtk_spin_button_get_value(spinButtonc4)

  if(deleteif==1)then
    m=aaa0
    mm=bbb0
  endif 
  if((gen(m)%npost.ne.0).or.(gen(mm)%npre.ne.0))then  ! actually, this should be an .and.
    wpo=0; wpr=0
   cto=0
do mn=1,ng
  do mnm=1,gen(mn)%nww
    if((gen(mn)%ww(mnm,1)==m).and.(gen(mn)%ww(mnm,2)==mm)) cto=cto+1
  end do
end do
   if(cto.ne.0) call undo_redo(6,m,mm,mc,0,mr,0d0)  
   if(cto.le.1)then
    call undo_redo(4,m,mm,0,0,0d0,0d0)
    do mn=1,gen(m)%npost
      if(gen(m)%post(mn)==mm)then; wpo=mn; exit; end if
    end do
    do mn=1,gen(mm)%npre
      if(gen(mm)%pre(mn)==m)then; wpr=mn; exit; end if
    end do

    if((wpo.ne.0).or.(wpr.ne.0))then
      if(wpo.lt.gen(m)%npost)then
        do mn=wpo,gen(m)%npost-1
          gen(m)%post(mn)=gen(m)%post(mn+1)
        end do
      end if
      if(wpr.lt.gen(mm)%npre)then
        do mn=wpr,gen(mm)%npre-1
          gen(mm)%pre(mn)=gen(mm)%pre(mn+1)
        end do
      end if

      gen(m)%npost=gen(m)%npost-1; gen(mm)%npre=gen(mm)%npre-1
      npo=gen(m)%npost; npr=gen(mm)%npre
      do mn=1,ng
        wpo=0
        do mnm=1,gen(mn)%nww
          if((gen(mn)%ww(mnm,1)==m).and.(gen(mn)%ww(mnm,2)==mm))then
            gen(mn)%ww(mnm,3)=0d0
            wpo=1
          end if
          if(wpo==1)then
            if(mnm.lt.gen(mn)%nww)then
              gen(mn)%ww(mnm,:)=gen(mn)%ww(mnm+1,:)
            else
              gen(mn)%nww=gen(mn)%nww-1
            endif
          endif
        end do
      end do

      print*, mm, " is not any longer produced by ", m, "."
    else
      print*, "No pres or posts."
    end if
  else ! selective catalyst
    do mn=1,gen(mc)%nww
      if((gen(mc)%ww(mn,1)==m).and.(gen(mc)%ww(mn,2)==mm))then
        gen(mc)%ww(mn,3)=0d0
        gen(mc)%ww(mn,1)=0; gen(mc)%ww(mn,2)=0
        goto 476
      endif
    end do
    476 if(mn.ne.gen(mc)%nww)then
     do mnm=mn,gen(mc)%nww-1
      gen(mc)%ww(mnm,:)=gen(mc)%ww(mnm+1,:)
     end do
    end if
    gen(mc)%nww=gen(mc)%nww-1
    endif    
  endif
call update_npag
  call update_view
end function my_v_delete

function my_v_connect(widget, event, gdata) result(ret)  bind(c)

  use gtk, only: gtk_notebook_prev_page, gtk_notebook_next_page, gtk_notebook_get_n_pages, gtk_notebook_get_current_page
  use global_widgets

  type(c_ptr), value :: widget, gdata, event
  integer(c_int)    :: ret

if(w4==1)then
  call my_ww_accept_s
if(w1*w2*w3.ne.0)then
  call choice_box_build_ww(w1,w2,w3,w4)
endif
end if

end function my_v_connect

function my_adh_accept(widget, event, gdata) result(ret)  bind(c)

  use global_widgets
  use gtk, only: gtk_notebook_prev_page, gtk_notebook_next_page, gtk_notebook_get_n_pages, gtk_notebook_get_current_page
  use basic_handlers
  use gtk_sup

  type(c_ptr), value :: widget, gdata, event
  type(c_ptr) :: genenamec
  character*12 :: genenamen
  integer(c_int)    :: ret
  integer :: m, mk, ma, mka, mn, m00, mk0
  real*8 :: mmk, mmk0, mmk00

  m=gtk_spin_button_get_value(spinButtonc4)
  mk=gtk_spin_button_get_value(spinButtonc4) ! 6
  ma=gtk_spin_button_get_value(spinButtonc8)
  mka=gtk_spin_button_get_value(spinButtonc8)!9
  mmk=gtk_spin_button_get_value(spinButtonc7)

  m00=gen(m)%wa(1); mk0=gen(mk)%wa(1)
  mmk0=kadh(m00,mk0); mmk00=kadh(mk0,m00)

if(ma*mka.ne.0)then
if(m*mk.ne.0)then
  gen(m)%wa(1)=ma
  gen(mk)%wa(1)=mka
end if

  kadh(ma,mka)=mmk
  kadh(mka,ma)=mmk

if(m*mk.ne.0)then
  print*, "The adhesion between protein", m, "and", mk, "has been changed."
  print*, "Protein", m, "has got adhesivity type", ma, "."
  print*, "Protein", mk, "has got adhesivity type", mka, "."
else
  print*, "The adhesion between adhesivity type", ma,"and ", mka,"has been changed."
end if
  print*, "Their adhesivity is now:", mmk, "."
end if

if(ma==0) gen(m)%wa(1)=0
if(mka==0) gen(mk)%wa(1)=0

  call undo_redo(1,m,mk,m00,mk0,mmk0,mmk00)

  call update_view

end function my_adh_accept

function my_adh_accept2(widget, event, gdata) result(ret)  bind(c)

  use global_widgets
  use gtk, only: gtk_notebook_prev_page, gtk_notebook_next_page, gtk_notebook_get_n_pages, gtk_notebook_get_current_page
  use basic_handlers
  use gtk_sup

  type(c_ptr), value :: widget, gdata, event
  type(c_ptr) :: genenamec
  character*12 :: genenamen
  integer(c_int)    :: ret
  integer :: m, mk, ma, mka, mn, m00, mk0
  real*8 :: mmk, mmk0, mmk00

  ma=gtk_spin_button_get_value(spinButtonc8)
  mka=gtk_spin_button_get_value(spinButtonc9)
  mmk=gtk_spin_button_get_value(spinButtonc7)

  mmk0=kadh(ma,mka)
  mmk00=kadh(mka,ma)

if(ma*mka.ne.0)then
  kadh(ma,mka)=mmk
  kadh(mka,ma)=mmk

  print*, "The adhesion between adhesion type", ma, "and", mka, "has been changed."
  print*, "Their adhesivity is now:", mmk, "."

do m=1,ntipusadh
  do mn=1,ntipusadh

end do; end do

end if

  m=0; mk=0
  do mn=1,ng
    if(gen(mn)%wa(1)==ma)then
      m=mn
      exit
    end if
  end do
  do mn=1,ng
    if(gen(mn)%wa(1)==mka)then
      mk=mn
      exit
    end if
  end do
  call undo_redo(1,m,mk,ma,mka,mmk0,mmk00)

  call update_view

end function my_adh_accept2

function my_adh_add(widget, event, gdata) result(ret)  bind(c)

  use global_widgets
  use gtk, only: gtk_notebook_prev_page, gtk_notebook_next_page, gtk_notebook_get_n_pages, gtk_notebook_get_current_page
  use basic_handlers
  use gtk_sup

  type(c_ptr), value :: widget, gdata, event
  type(c_ptr) :: genenamec
  character*12 :: genenamen
  integer(c_int)    :: ret
  integer :: m, mk, ma, mka
  real*8 :: mmk
  real*8,allocatable :: kadhs(:,:),kpos2(:,:)

if(ntipusadh.gt.0)then
  if(allocated(kadhs)) deallocate(kadhs)
  allocate(kadhs(ntipusadh,ntipusadh))
  kadhs=0d0
  kadhs=kadh
end if

  if(allocated(kadh)) deallocate(kadh)

if(ntipusadh.gt.0)then
  do i=1,ntipusadh
    kpos2(i,:)=kpos(i,:)
    do m=1,ntipusadh
      kadh(i,m)=kadhs(i,m)
    end do
  end do
end if

  ntipusadh=ntipusadh+1

  if(allocated(kpos)) deallocate(kpos)
  allocate(kpos(ntipusadh,2))
  kpos=kpos2

  do i=1,ntipusadh
    kadh(ntipusadh,i)=0d0
    kadh(i,ntipusadh)=0d0
  end do


  call update_view

end function my_adh_add

subroutine my_adh_add_s

  use global_widgets
  use gtk, only: gtk_notebook_prev_page, gtk_notebook_next_page, gtk_notebook_get_n_pages, gtk_notebook_get_current_page
  use basic_handlers
  use gtk_sup

  type(c_ptr) :: genenamec
  character*12 :: genenamen
  integer(c_int)    :: ret
  integer :: m, mk, ma, mka
  real*8 :: mmk
  real*8,allocatable :: kadhs(:,:),kpos2(:,:)

if(ntipusadh.gt.0)then
  if(allocated(kadhs)) deallocate(kadhs)
  allocate(kadhs(ntipusadh,ntipusadh))
  kadhs=0d0
  kadhs=kadh
end if

  if(allocated(kadh)) deallocate(kadh)
  allocate(kadh(ntipusadh+1,ntipusadh+1))
  kadh=0d0

  if(allocated(kpos2)) deallocate(kpos2)
  allocate(kpos2(ntipusadh+1,2))
  kpos2=0d0

if(ntipusadh.gt.0)then
  do i=1,ntipusadh
    kpos2(i,:)=kpos(i,:)
    do m=1,ntipusadh
      kadh(i,m)=kadhs(i,m)
    end do
  end do
end if

  ntipusadh=ntipusadh+1

  if(allocated(kpos)) deallocate(kpos)
  allocate(kpos(ntipusadh,2))
  kpos=kpos2

  do i=1,ntipusadh
    kadh(ntipusadh,i)=0d0
    kadh(i,ntipusadh)=0d0
  end do


  call update_view

end subroutine my_adh_add_s

subroutine move_nodes(cualn,xpos,ypos,cuall)

  use gtk, only: gtk_notebook_prev_page, gtk_notebook_next_page, gtk_notebook_get_n_pages, gtk_notebook_get_current_page
  use global_widgets

  integer:: cualn, cuall
  integer :: xpos,ypos

  type(c_ptr) :: boxc, checkd
  type(c_ptr) :: labelc1, labelc2, labelc3, labelc4, buttonc, buttond, tablec

  if(cuall.lt.2)then
    npos(cualn,1)=xpos
    npos(cualn,2)=ypos
  elseif((cuall.ge.2).and.(cuall.lt.4))then
    npos(cualn,3)=xpos
    npos(cualn,4)=ypos
  elseif(cuall==4)then
    npos(cualn,5)=xpos
    npos(cualn,6)=ypos
  elseif(cuall==5)then
    kpos(cualn,1)=xpos
    kpos(cualn,2)=ypos
  end if

  call update_view
     
end subroutine move_nodes

subroutine info_prop(cualp,pr)

  integer :: cualp,pr

  if(cualp==3)then
    if(pr==1)then
      print*, "This defines the equilibrium distance for nodes in the same face and cell."
    elseif(pr==2)then
      print*, "This defines the maximum distance for interaction between nodes."
    elseif(pr==3)then
      print*, "This defines the elasticity for nodes in the same cell (cohesion)."
    elseif(pr==4)then
      print*, "This defines the unspecific adhesion (adhesion) between nodes in different cells."
    elseif(pr==5)then
      print*, "This defines the strength of the repulsion between nodes in the same cell."
    elseif(pr==6)then
      print*, "This defines the strength of the repulsion between nodes in different cells."
    elseif(pr==7)then
      print*, "The strength of the lateral epithelial surface tension (torsion)."
    elseif(pr==8)then
      print*, "The strength of the apico-basal epithelial surface tension (torsion)."
    elseif(pr==9)then
      print*, "This defines the equilibrium distance for the spring of the ellipse (epithelial sides)."
    elseif(pr==10)then
      print*, "This defines the spring constant for the ellipse (epithelial sides)."
    elseif(pr==11)then
      print*, "This determines the motility of the nodes."
    elseif(pr==12)then
      print*, "This determines the maximal node displacement that may occur by noise."
    elseif(pr==13)then
      print*, "This is the equilibrium radius component due to internal contraction and external deformations."
    elseif(pr==14)then
      print*, "This is the equilibrium radius component due to node growth."
    elseif(pr==15)then
      print*, "This is the equilibrium radius component due to plasticity."
    elseif(pr==16)then
      print*, "This is the equilibrium radius component correcting for volume conservation."
    elseif(pr==17)then
      print*, "This is the differentiation state of the node, 0 is no differentiation, 1 is total differentiation."
    elseif(pr==18)then
      print*, "This defines the elastic constant for locally fixed nodes."
    elseif(pr==19)then
      print*, "This defines the plasticity constant for epithelial plastic deformation."
    elseif(pr==20)then
      print*, "This is the volume conservation constant for epithelial plastic deformation."
    end if
  elseif(cualp==4)then
    if(pr==1)then
      print*, "This is simply cell growth."
    elseif(pr==2)then
      print*, "This lets the cell cycle increase, when==1 the cell can divide if it has the right size."
    elseif(pr==3)then
      print*, "This mediates apoptosis."
    elseif(pr==4)then
      print*, "This is the rate of secretion the gene promotes."
    elseif(pr==5)then
      print*, "A value of 1 means that this gene itself is secretable."
    elseif(pr==6)then
      print*, "This is the intercellular repulsion of the secreted node."
    elseif(pr==7)then
      print*, "The maximal interaction distance (as proportion of reqcel) of the secreted node."
    elseif(pr==8)then
      print*, "This gene is polarizing the cell [affecting polarization vector of the cell]."
    elseif(pr==9)then
      print*, "This gene causes the cells to grow into the direction of polarization with a noise of 1-THIS."
    elseif(pr==10)then
      print*, "This gene affects the minimal number of nodes required for a cell to divide."
    elseif(pr==11)then
      print*, "The dependence of the plane of division from chemical gradients over the Hertwig vector."
    elseif(pr==12)then
      print*, "This gene activates assymetric division."
    elseif(pr==13)then
      print*, "This gene activates epitelial-to-mesenchymal transition (EMT)."
    elseif(pr==14)then
      print*, "This gene promotes ECM proteolysis."
    elseif(pr==15)then
      print*, "This gene changes the max number of nodes that are allowed within a cell without division."
    elseif(pr==16)then
      print*, "This biases noise by the polarization vector."
    end if
  end if

end subroutine info_prop

subroutine GRN_search(plotid,xpos,ypos,butt,which,cual)

  use some_widgets

  integer::plotid,which,butt,cual,mm,mn,nn,whichx,which3,parx,whichi,wos,what,xpos,ypos
  real*8 ::radiu,wo,won
  integer ::npars(23)
  showup=1  !!!RZ
  choiceboxthere=0
  if(cual.le.1)then ! CASE 1
    radiu=sqrt((xpos-350d0)**2+(ypos-350d0)**2)

  if(butt==2)then
    call gtk_widget_destroy(my_choice_window)
  end if
  if(togglen==0)then !ONLY PASS OVER HERE IF YOU HAVE NOT MOVED ANY NODE
    if(abs(radiu-250d0).lt.40)then
      which=-1
      if(ypos.gt.335)then
        wo=acos((xpos-350d0)/radiu)*0.5*ng/pi
        if(abs(wo-nint(wo)).lt.0.2) whichi=int(wo+0.5)
      else
        wo=acos((xpos-350d0)/radiu)*0.5*ng/pi
        if(abs(wo-nint(wo)).lt.0.2) whichi=ng-int(wo+0.5)
      endif
      if(whichi==0) whichi=ng
      if((whichi.ge.0).and.(whichi.le.ng+1))then
      call update_view
      wcl=whichi
        print*, "************************************", xpos, ypos
        print*, "You have chosen this form: ", whichi
        if(gen(whichi)%kindof==1) print*, "This is a gene with only transcription"
        if(gen(whichi)%kindof==2) print*, "This is a gene with transcription and translation"
        if(gen(whichi)%kindof==3) print*, "This is an intracellular gene product"
        if(gen(whichi)%kindof==4) print*, "This is an extracellularly diffusive protein"
        if(gen(whichi)%kindof==5) print*, "This is a basally localizing protein"
        if(gen(whichi)%kindof==6) print*, "This is an apically localizing protein"
        if(gen(whichi)%kindof==7) print*, "This is a membrane receptor for juxtacrine signalling"
        if(gen(whichi)%kindof==8) print*, "This is a membrane receptor for diffusive extracellular ligands"
        print*, "kindof: ", gen(whichi)%kindof
        print*, "Diffusion rate: ", gen(whichi)%diffu
        print*, "Degradation rate: ", gen(whichi)%mu
        if(gen(whichi)%npre.ne.0) print*, "arisen from: ", gen(whichi)%pre(1:gen(whichi)%npre)
        if(gen(whichi)%npost.ne.0) print*, "giving rise to: ", gen(whichi)%post(1:gen(whichi)%npost)
        print*, "affecting the following mechanical parameters/cell behaviours: "
        k=0
        do i=1,nparam
          if(gen(whichi)%wa(i).ne.0)then
            print*, "parameter ", i, " : ", gen(whichi)%wa(i)
            k=k+1
          end if
        end do
        if(k==0) print*, "none"
        print*, "___________________________________", whichi, whichii
        if(((butt.eq.1).or.((butt.eq.3).and.(cual.eq.0))).and.(whichii.gt.0))then
          print*, "Actions of ", whichii, " on ", whichi
          print*, "T-Matrix: ", gen(whichi)%w(whichii)
          if(gen(whichi)%w(whichii).gt.0) print*, "Form ", whichii, " is a transcriptional activator of ", whichi, " ."
          if(gen(whichi)%w(whichii).lt.0) print*, "Form ", whichii, " is a transcriptional repressor of ", whichi, " ."
          print*, "___________________________________"
         if(cual==0) aaa0=whichii; bbb0=whichi; ccc0=cual; ddd0=butt
         showup=-1 !!!RZ
         if((cual==0).and.(faster==-1)) call choice_box_build(whichii,whichi,cual,butt) ! we modify their interactions
          if(whichi.gt.0)then; whichii=0; whichi=0; endif
        endif
        if(butt.eq.1)then
          whichii=whichi
        end if
      end if
    end if
  else ! PASS OVER HERE IF YOU HAVE MOVED NODES
    do mm=1,ng
      if((abs(ypos-npos(mm,2)).lt.20).and.(abs(xpos-npos(mm,1)).lt.20))then; whichi=mm; endif ! RZ ; which=whichi
      if(whichi==-1)then; which=0;whichii=0; endif
    end do
    if((whichi.gt.0).and.(whichi.le.ng+1))then
    call update_view
    wcl=whichi
    which=whichi !!!!!
      print*, "************************************", xpos, ypos, "***"
      print*, "You have chosen this form: ", which
      if(gen(whichi)%kindof==1) print*, "This is a gene with only transcription"
      if(gen(whichi)%kindof==2) print*, "This is a gene with transcription and translation"
      if(gen(whichi)%kindof==3) print*, "This is an intracellular gene product"
      if(gen(whichi)%kindof==4) print*, "This is an extracellularly diffusive protein"
      if(gen(whichi)%kindof==5) print*, "This is a basally localizing protein"
      if(gen(whichi)%kindof==6) print*, "This is an apically localizing protein"
      if(gen(whichi)%kindof==7) print*, "This is a membrane receptor for juxtacrine signalling"
      if(gen(whichi)%kindof==8) print*, "This is a membrane receptor for diffusive extracellular ligands"
      print*, "kindof: ", gen(whichi)%kindof
      print*, "Diffusion rate: ", gen(whichi)%diffu
      print*, "Degradation rate: ", gen(whichi)%mu
      if(gen(whichi)%npre.ne.0) print*, "arisen from: ", gen(whichi)%pre(1:gen(whichi)%npre)
      if(gen(whichi)%npost.ne.0) print*, "giving rise to: ", gen(whichi)%post(1:gen(whichi)%npost)
      print*, "affecting the following mechanical parameters/cell behaviours: "
      k=0
      do i=1,nparam
        if(gen(whichi)%wa(i).ne.0)then
          print*, "parameter ", i, " : ", gen(whichi)%wa(i)
          k=k+1
        end if
      end do
      if(k==0) print*, "none"
      print*, "___________________________________"
    end if
if(whichii.ge.0)then; whichi=which; else; whichi=0; endif


    if(togglen==-1)then
    whichi=wcl !!!!!
    whichii=wcl !!!!
      if(((butt.eq.1).or.((butt.eq.3).and.(cual.eq.0))).and.(whichii.gt.0))then
        do mm=1,ng
          if(whichi==-1) whichi=0
          if((abs(ypos-npos(mm,2)).lt.20).and.(abs(xpos-npos(mm,1)).lt.20)) whichii=mm
        end do
        print*, "LETS MOVE IT:", whichii, " TO:", xpos,ypos
        call move_nodes(whichii,xpos,ypos,cual)       
        print*, "___________________________________"      
      endif
    elseif(togglen==1)then
      if(((butt.eq.1).or.((butt.eq.3).and.(cual.eq.0))).and.(whichii.gt.0).and.(whichi.gt.0))then
      if((whichii.le.ng).and.(whichii.gt.0))then
        print*, "Actions of ", whichii, " on ", whichi
        print*, "T-Matrix: ", gen(whichi)%w(whichii)
        if(gen(whichi)%w(whichii).gt.0) print*, "Form ", whichii, " is a transcriptional activator of ", whichi, " ."
        if(gen(whichi)%w(whichii).lt.0) print*, "Form ", whichii, " is a transcriptional repressor of ", whichi, " ."
       print*, "___________________________________"
         if(cual==0) aaa0=whichii; bbb0=whichi; ccc0=cual; ddd0=butt
         showup=-1
         if((cual==0).and.(faster==-1)) call choice_box_build(whichii,whichi,cual,butt) ! we modify their interactions
         if(whichi.gt.0)then; whichii=0; whichi=0; endif
               if(butt.eq.1)then
        if(whichi.gt.0)then; whichii=0; whichi=0; endif
      end if; endif
      end if
      if(whichi.gt.0) whichii=whichi
    end if
    end if
  end if
  if((cual==1))then
    
if(which3.gt.ng) which3=0 !!!
showup=-1
    if(togglen==0)then !ONLY PASS OVER HERE IF YOU HAVE NOT MOVED ANY NODE
  
      radiu=sqrt((xpos-350d0)**2+(ypos-350d0)**2)

      if(butt==2)then
        call gtk_widget_destroy(my_choice_window)
      end if
      if(abs(radiu-250d0).lt.40)then
        if(ypos.gt.335)then
          wo=acos((xpos-350d0)/radiu)*0.5*ng/pi
          if(abs(wo-nint(wo)).lt.0.2) which=int(wo+0.5)
        else
          wo=acos((xpos-350d0)/radiu)*0.5*ng/pi
          if(abs(wo-nint(wo)).lt.0.2) which=ng-int(wo+0.5)
        endif
        if(which==0) which=ng
        if((which.gt.0).and.(which.le.ng+1))then
call update_view
        wcl=which
        print*, "************************************"
        print*, "You have chosen this form: ", which
        if(gen(which)%kindof==1) print*, "This is a gene with only transcription"
        if(gen(which)%kindof==2) print*, "This is a gene with transcription and translation"
        if(gen(which)%kindof==3) print*, "This is an intracellular gene product"
        if(gen(which)%kindof==4) print*, "This is an extracellularly diffusive protein"
        if(gen(which)%kindof==5) print*, "This is a basally localizing protein"
        if(gen(which)%kindof==6) print*, "This is an apically localizing protein"
        if(gen(which)%kindof==7) print*, "This is a membrane receptor for juxtacrine signalling"
        if(gen(which)%kindof==8) print*, "This is a membrane receptor for diffusive extracellular ligands"
        print*, "kindof: ", gen(which)%kindof
        print*, "Diffusion rate: ", gen(which)%diffu
        print*, "Degradation rate: ", gen(which)%mu
        if(gen(which)%npre.ne.0) print*, "arisen from: ", gen(which)%pre(1:gen(which)%npre)
        if(gen(which)%npost.ne.0) print*, "giving rise to: ", gen(which)%post(1:gen(which)%npost)
        print*, "affecting the following mechanical parameters/cell behaviours: "
        k=0
        do i=1,nparam
          if(gen(which)%wa(i).ne.0)then
            print*, "parameter ", i, " : ", gen(which)%wa(i)
            k=k+1
          end if
        end do
        if(k==0) print*, "none"
        print*, "___________________________________"
        if(butt==1)then
          w0=0
          which3=which2; which2=whichsave; whichsave=which
        if((butt==1).and.(which2.gt.0).and.(which3.eq.0).and.(cual==1))then
          w0=1
          aaa0=which2; bbb0=whichsave; ccc0=cual; ddd0=butt
         if(faster==-1) call choice_box_build(which2,whichsave,cual,butt) ! we modify their interactions
        endif
        if((butt==1).and.(whichsave.gt.0).and.(which3.gt.0).and.(cual==1))then
          w0=2
          print*, "Catalyzation of the reaction of ", which3, " into ", which2, " by ", whichsave, ".", w0
          k=0
          do i=1,gen(whichsave)%nww
            if(gen(whichsave)%ww(i,1)==which3)then
              print*, "The catalytic activity is: ", gen(whichsave)%ww(i,3), "."; k=1; exit
            end if
          enddo
          if(k==0) print*, "There is no catalytic interaction so far."
          print*, "___________________________________"
          w1=which3; w2=which2; w3=whichsave; w4=cual
          !if(faster==-1) call choice_box_build_ww(which3,which2,whichsave,cual) ! X
           call choice_box_build_ww(which3,which2,whichsave,cual) 
          if(which3.ne.0)then
            which3=0; which=-1; whichsave=0; which2=0
          end if
        end if; endif
      end if
    end if
  else ! IF MOVED
    which=-1; showup=-1
  do mm=1,ng
    if(which==-1) which=0
    if((abs(ypos-npos(mm,2)).lt.20).and.(abs(xpos-npos(mm,1)).lt.20)) which=mm
  end do
    if(which==0) which=ng
    if((which.gt.0).and.(which.le.ng+1))then
     call update_view
      wcl=which
      print*, "************************************"
      print*, "You have chosen this form: ", which
      if(gen(which)%kindof==1) print*, "This is a gene with only transcription"
      if(gen(which)%kindof==2) print*, "This is a gene with transcription and translation"
      if(gen(which)%kindof==3) print*, "This is an intracellular gene product"
      if(gen(which)%kindof==4) print*, "This is an extracellularly diffusive protein"
      if(gen(which)%kindof==5) print*, "This is a basally localizing protein"
      if(gen(which)%kindof==6) print*, "This is an apically localizing protein"
      if(gen(which)%kindof==7) print*, "This is a membrane receptor for juxtacrine signalling"
      if(gen(which)%kindof==8) print*, "This is a membrane receptor for diffusive extracellular ligands"
      print*, "kindof: ", gen(which)%kindof
      print*, "Diffusion rate: ", gen(which)%diffu
      print*, "Degradation rate: ", gen(which)%mu
      if(gen(which)%npre.ne.0) print*, "arisen from: ", gen(which)%pre(1:gen(which)%npre)
      if(gen(which)%npost.ne.0) print*, "giving rise to: ", gen(which)%post(1:gen(which)%npost)
      print*, "affecting the following mechanical parameters/cell behaviours: "
      k=0
      do i=1,nparam
        if(gen(which)%wa(i).ne.0)then
          print*, "parameter ", i, " : ", gen(which)%wa(i)
          k=k+1
        end if
      end do
      if(k==0) print*, "none"
      print*, "___________________________________"
      if(togglen==-1)then
        if((butt==1).and.(which.gt.0))then
call update_view
          wcl=which
          print*, "LETS MOVE IT:", which, " TO:", xpos,ypos
          which3=which2; which2=whichsave; whichsave=which
        elseif((butt==1).and.(whichsave.gt.0).and.(which2.gt.0).and.(cual==1))then
          call move_nodes(whichsave,xpos,ypos,cual)        
          print*, "___________________________________"      
        endif
      elseif(togglen==1)then
          if(butt==1)then
            w0=0
            which3=which2; which2=whichsave; whichsave=which
          if((butt==1).and.(which2.gt.0).and.(which3.eq.0).and.(cual==1))then
            w0=1
          aaa0=which2; bbb0=whichsave; ccc0=cual; ddd0=butt
          if(faster==-1) call choice_box_build(which2,whichsave,cual,butt) ! we modify their interactions ! XY
          elseif((butt==1).and.(whichsave.gt.0).and.(which3.gt.0).and.(cual==1))then
            w0=2
            print*, "Catalyzation of the reaction of ", which2, " into ", whichsave, " by ", which, "."
            k=0
            do i=1,gen(which)%nww
              if(gen(which)%ww(i,1)==which2)then
                print*, "The catalytic activity is: ", gen(which)%ww(i,3), "."; k=1; exit
              end if
            enddo
            if(k==0) print*, "There is no catalytic interaction so far."
            print*, "___________________________________"
            w1=which3; w2=which2; w3=whichsave; w4=cual
            !if(faster==-1) call choice_box_build_ww(which3,which2,whichsave,cual) ! X
            call choice_box_build_ww(which3,which2,whichsave,cual) 
          if(which3.ne.0)then
            which3=0; which=-1; whichsave=0; which2=0
          end if
          end if
        end if; endif
      end if
    endif


  elseif((cual.eq.2).or.(cual.eq.3))then
showup=-1
    npars=(/5,6,7,8,9,10,11,12,13,14,15,16,20,21,22,23,24,25,26,27,28,34,35/)

    if(togglen==0)then !ONLY PASS OVER HERE IF YOU HAVE NOT MOVED ANY NODE

      radiu=sqrt((xpos-350d0)**2+(ypos-350d0)**2)

      if(butt==2)then
        call gtk_widget_destroy(my_choice_window)
      end if
      if(abs(radiu-250d0).lt.40)then
        which=-1
        if(ypos.gt.335)then
          wo=acos((xpos-350d0)/radiu)*ng/pi
          if(abs(wo-nint(wo)).lt.0.2) which=int(wo+0.5)
        endif
        if(which==0) which=ng
        if((which.gt.0).and.(which.le.ng+1))then
call update_view
          wcl=which
          print*, "************************************"
          print*, "You have chosen this form: ", which
          if(cual==2)then
            print*, "It regulates the following node properties:"
            nn=0
            do mm=1,23
              if(gen(which)%wa(mm).ne.0)then
                nn=nn+1
                print*, mm, nodeparams(mm), gen(which)%wa(mm)
              end if
            end do
            if(nn==0) print*, "none"
          else
            print*, "It regulates the following cell properties:"
            nn=0
            do mm=1,16
              if(gen(which)%wa(nparam_per_node+mm).ne.0)then
                nn=nn+1
                print*, nparam_per_node+mm, cellparams(nparam_per_node+mm), gen(which)%wa(nparam_per_node+mm)
              end if
            end do
            if(nn==0) print*, "none"
          end if
          print*, "___________________________________"
          if(butt.gt.0)then
            whichsave=which
          end if
        end if
      end if
        if((xpos.le.663.5).and.(xpos.ge.100))then
          if(cual==2)then
            if((ypos.le.200).and.(ypos.ge.150))then
              wo=int((xpos-100)/24.5)+1
              wos=wo
              wo= npars(wo)
              if((butt.gt.0).and.(whichsave.gt.0))then; parx=wo; endif
                print*, "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"
                print*, "This is node property", int(wo), ", i.e. ", nodeparams(int(wo)), "."
                call info_prop(cual,wos)
                print*, "It is regulated by the following forms."
                mn=0
                do mm=1,ng
                  if(gen(mm)%wa(wo).ne.0)then
                    mn=mn+1
                    print*, mm, gen(mm)%wa(wo)
                  end if
                end do
                if(mn==0) print*, "none"
              end if
            else if(cual==3)then
            if((ypos.le.250).and.(ypos.ge.100))then
              wo=int((xpos-100)/30.625)+1
              if((butt.gt.0).and.(whichsave.gt.0))then; parx=wo; endif
                print*, "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"
                print*, "This is cell property", int(wo), ", i.e. ", cellparams(int(wo)), "."
                call info_prop(cual,int(wo))
                print*, "It is regulated by the following forms."
                mn=0
                do mm=1,ng
                  if(gen(mm)%wa(nparam_per_node+wo).ne.0)then
                    mn=mn+1
                    print*, mm, gen(mm)%wa(nparam_per_node+wo)
                  end if
                end do
                if(mn==0) print*, "none"
              end if
            endif
          end if
          if((parx.gt.0).and.(parx.le.35))then
            print*, whichsave, parx, cual, wo
            print*, "_________________________________"
            aaa0=whichsave; bbb0=parx; ccc0=cual; ddd0=butt
            if(faster==-1) call choice_box_param_build(whichsave,parx,cual,butt) !X
          end if
        else  ! NODES HAVE BEEN MOVED

        if(butt==2)then
          call gtk_widget_destroy(my_choice_window)
        end if
        which=-1
        do mm=1,ng
          if(which==-1) which=0
          if((abs(ypos-npos(mm,4)).lt.20).and.(abs(xpos-npos(mm,3)).lt.20))then
            which=mm
          end if
          print*, ypos, xpos, which, npos(mm,3)
        end do
        if(butt.ge.1)then
          if(which.gt.0) whichsave=which
        endif
        if((which.gt.0).and.(which.le.ng+1))then
call update_view
          wcl=which
          print*, "************************************"
          print*, "You have chosen this form: ", which
          if(cual==2)then
            print*, "It regulates the following node properties:"
            nn=0
            do mm=1,23
              if(gen(which)%wa(mm).ne.0)then
                nn=nn+1
                print*, mm, nodeparams(mm), gen(which)%wa(mm)
              end if
            end do
            if(nn==0) print*, "none"
          else if(cual==3)then
            print*, "It regulates the following cell properties:"
            nn=0
            do mm=1,16
              if(gen(which)%wa(nparam_per_node+mm).ne.0)then
                nn=nn+1
                print*, nparam_per_node+mm, cellparams(nparam_per_node+mm), gen(which)%wa(nparam_per_node+mm)
              end if
            end do
            if(nn==0) print*, "none"
          end if
          print*, "___________________________________"
          !call choice_prop_box_build(which)
        end if
    
        if(togglen==-1)then
          if((butt==1).and.(which.gt.0))then
          elseif((butt==1).and.(whichsave.gt.0))then
call update_view
            wcl=whichsave
            print*, "LETS MOVE IT:", whichsave, " TO:", xpos,ypos
            call move_nodes(whichsave,xpos,ypos,cual)      
            print*, "___________________________________"      
          endif
        elseif(togglen==1)then
            if((xpos.le.663.5).and.(xpos.ge.100))then
              if(cual==2)then
                if((ypos.le.200).and.(ypos.ge.150))then
                  wo=int((xpos-100)/24.5)+1
                  wos=wo
                  wo= npars(wo)
                  if((butt.gt.0).and.(whichsave.gt.0))then; parx=wo; endif
                  print*, "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"
                  print*, "This is node property", int(wo), ", i.e. ", nodeparams(int(wo)), "."
                  call info_prop(cual,wos)
                  print*, "It is regulated by the following forms."
                  mn=0
                  do mm=1,ng
                    if(gen(mm)%wa(wo).ne.0)then
                      mn=mn+1
                      print*, mm, gen(mm)%wa(wo)
                    end if
                  end do
                  if(mn==0) print*, "none"
                end if
              elseif(cual==3)then
                if((ypos.le.250).and.(ypos.ge.100))then
                  wo=int((xpos-100)/30.625)+1
                  if((butt.gt.0).and.(whichsave.gt.0))then; parx=wo; endif
                  print*, "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"
                  print*, "This is cell property", int(wo), ", i.e. ", cellparams(int(wo)), "."
                  call info_prop(cual,int(wo))
                  print*, "It is regulated by the following forms."
                  mn=0
                  do mm=1,ng
                    if(gen(mm)%wa(nparam_per_node+wo).ne.0)then
                      mn=mn+1
                   end if
                 end do
                 if(mn==0) print*, "none"
               end if
             end if
           end if
        endif
      end if

      if((parx.gt.0).and.(parx.le.35))then
      print*, whichsave, parx, cual, wo
      print*, "_________________________________"
      aaa0=whichsave; bbb0=parx; ccc0=cual; ddd0=butt
      if(faster==-1) call choice_box_param_build(whichsave,parx,cual,butt) !X
    end if

elseif((cual==4).or.(cual==5))then ! CASE 3
showup=-1
      if((abs(590-xpos).le.140).and.(abs(60-ypos).le.10)) call my_adh_add_s
      if(butt==2)then
        call gtk_widget_destroy(my_choice_window)
      end if
      which=-1
      what=-1
      do mm=1,ng
        radiu=sqrt((xpos-npos(mm,5))**2+(ypos-npos(mm,6))**2)
        if(radiu.lt.20) which=mm
        if((which.gt.0).and.(which.le.ng))then
call update_view
          wcl=which
          print*, "************************************"
          print*, "You have chosen this form: ", which
          if(gen(which)%wa(1).ne.0)then
            print*, "It mediates adhesion type:",gen(which)%wa
          else
            print*, "It doesn't mediate any kind of adhesion."
          end if
          if((butt.gt.0).and.(which.gt.0))then
            which2=whichsave
            whichsave=which
            what2=0
            exit
          end if
        end if
      end do
      do mm=1,ntipusadh
        radiu=sqrt((xpos-kpos(mm,1))**2+(ypos-kpos(mm,2))**2)
        if(radiu.lt.20) what=mm
        if((butt.gt.0).and.(what.gt.0))then
          print*, "This is adhesion molecule", what
          what2=whatsave
          whatsave=what
          which2=0
          exit
        end if
      end do

        if((which2.gt.0).or.((what.eq.-1).and.(what2.eq.0).and.(whatsave.eq.0)))then ! This will be depreciated
        whatsave=0; what=0; what2=0
        wcl=whichsave
        if(togglen.ge.0)then
        else
          print*, "LETS MOVE IT:", whichsave, " TO:", xpos,ypos
          call move_nodes(whichsave,xpos,ypos,4)      
          print*, "___________________________________" 
        end if     

       elseif((what2.gt.0).or.((which.eq.-1).and.(which2.eq.0).and.(whichsave.eq.0)))then
        whichsave=0; which=0; which2=0
        if(togglen.ge.0)then
          aaa0=what2; bbb0=whatsave; ccc0=4; ddd0=3
          if(faster==-1) call choice_box_build(what2,whatsave,4,3) !X
        else
          whichsave=0
          print*, "LETS MOVE IT:", whatsave, " TO:", xpos,ypos
          call move_nodes(whatsave,xpos,ypos,5)      
          print*, "___________________________________" 
        end if  

      elseif((what.eq.-1).and.(which2.gt.0))then
        wcl=which2
        if(togglen.ge.0)then
          aaa0=which2; bbb0=whatsave; ccc0=4; ddd0=2
          if(faster==-1) call choice_box_build(which2,whatsave,4,2) !X
        else
          print*, "LETS MOVE IT:", which2, " TO:", xpos,ypos
          call move_nodes(which2,xpos,ypos,4)   
          whatsave=0; what=0; what2=0
          print*, "___________________________________" 
        end if     
      elseif((what.eq.-1).and.(which2.eq.0))then
        wcl=whichsave
        if(togglen.ge.0)then
          aaa0=whichsave; bbb0=whatsave; ccc0=4; ddd0=2
          if(faster==-1) call choice_box_build(whichsave,whatsave,4,2)!X
        else
          print*, "LETS MOVE IT:", whichsave, " TO:", xpos,ypos
          call move_nodes(whichsave,xpos,ypos,4)   
          whatsave=0; what=0; what2=0
          print*, "___________________________________" 
        end if     
       elseif((which.eq.-1).and.(what2.gt.0))then
        wcl=whichsave
        if(togglen.ge.0)then
          aaa0=whichsave; bbb0=what2; ccc0=4; ddd0=2
          if(faster==-1) call choice_box_build(whichsave,what2,4,2)!X
        else
          whichsave=0; which=0; which2=0
          print*, "LETS MOVE IT:", what2, " TO:", xpos,ypos
          call move_nodes(what2,xpos,ypos,5)      
          print*, "___________________________________" 
        end if     
      elseif((which.eq.-1).and.(what2.eq.0))then
        wcl=whichsave
        if(togglen.ge.0)then
          aaa0=whichsave; bbb0=whatsave; ccc0=4; ddd0=2
          if(faster==-1) call choice_box_build(whichsave,whatsave,4,2) !X
        else
          whichsave=0; which=0; which2=0
          print*, "LETS MOVE IT:", whatsave, " TO:", xpos,ypos
          call move_nodes(whatsave,xpos,ypos,5)      
          print*, "___________________________________" 
        end if        
      end if

  end if
      
end subroutine GRN_search

end module on_display_handlers


module table_handlers

  use gtk_hl
  use gtk, only: gtk_button_new, gtk_check_button_new, gtk_container_add, gtk_ent&
       &ry_get_text, gtk_entry_get_text_length, gtk_entry_new, gtk_entry_set_text, gtk&
       &_main, gtk_main_quit, gtk_widget_destroy, gtk_toggle_button_get_active, gtk_to&
       &ggle_button_set_active, gtk_widget_show, gtk_widget_show_all, gtk_window_new, &
       & gtk_init, gtk_tree_path_new, gtk_tree_view_get_path_at_pos, gtk_vbox_new,&
       &gtk_tree_view_convert_widget_to_bin_window_coords, gtk_tree_selection_select_path, &
       &gtk_menu_popup, gtk_table_new, gtk_disable_setlocale
  use g, only: alloca
  use gdk_events
  use genetic
  use general
  use io
  use basic_handlers
  
  implicit none

  type(c_ptr) :: ihwin,ihwin2,ihscrollcontain,ihlist, base,&
       &  qbut, dbut, lbl
  real*8,allocatable  :: neww(:,:),oldw(:,:),oldwa(:,:) 
  integer :: open_tab,whichtab

contains

subroutine my_destroy_tab(widget, gdata) bind(c)

  type(c_ptr), value :: widget, gdata

  call gtk_widget_destroy(ihwin)
  call gtk_main_quit ()

end subroutine my_destroy_tab

subroutine my_destroy_tab2(widget, gdata) bind(c)

  type(c_ptr), value :: widget, gdata

  call gtk_widget_destroy(ihwin2)
  call gtk_main_quit ()

end subroutine my_destroy_tab2

subroutine my_accept_changes_tab (widget, gdata) bind(c)

  type(c_ptr), value :: widget, gdata
  integer :: m,mm

if(open_tab==0)then
if(whichtab==1)then
if(allocated(neww))then
  do m=1,ng; do mm=1,ng
  if(gen(m)%w(mm).ne.neww(m,mm))then
                       call undo_redo(5,m,mm,0,0,gen(m)%w(mm),neww(m,mm))  ! UNDO REDO SAVE
    gen(m)%w(mm)=neww(m,mm)
  end if
  end do
    call undo_redo(3,m,1,0,0,gen(m)%kindof,neww(m,ng+3)) 
    call undo_redo(3,m,3,0,0,gen(m)%diffu,neww(m,ng+2))
    call undo_redo(3,m,2,0,0,gen(m)%mu,neww(m,ng+1)) 
    gen(m)%diffu=neww(m,ng+1)
    gen(m)%mu=neww(m,ng+2)
    gen(m)%kindof=int(neww(m,ng+3))
  end do
endif
!elseif(open_tab==2)then
else
if(allocated(oldw))then
  do m=1,ng; do mm=1,ng
  if(gen(m)%w(mm).ne.oldw(m,mm))then
                       call undo_redo(5,m,mm,0,0,gen(m)%w(mm),oldw(m,mm))  ! UNDO REDO SAVE
    gen(m)%w(mm)=oldw(m,mm)
  end if
  end do
    call undo_redo(3,m,1,0,0,gen(m)%kindof,oldw(m,ng+3)) 
    call undo_redo(3,m,3,0,0,gen(m)%diffu,oldw(m,ng+2))
    call undo_redo(3,m,2,0,0,gen(m)%mu,oldw(m,ng+1))
    gen(m)%diffu=oldw(m,ng+1)
    gen(m)%mu=oldw(m,ng+2)
    gen(m)%kindof=int(oldw(m,ng+3))
  end do
endif
endif
elseif((open_tab==2))then
if(whichtab==2)then
if(allocated(oldwa))then
  do m=1,ng; do mm=1,23
  if(gen(m)%wa(mm).ne.oldwa(m,mm))then
    call undo_redo(7,m,mm,0,0,gen(m)%wa(mm),oldwa(m,mm))  ! UNDO REDO SAVE
    gen(m)%wa(mm)=oldwa(m,mm)
  end if
  end do
  end do
endif
else
  do m=1,ng; do mm=1,23
  if(gen(m)%wa(mm).ne.oldwa(m,mm))then
    call undo_redo(7,m,mm,0,0,gen(m)%wa(mm),neww(m,mm))  ! UNDO REDO SAVE
  end if
  end do
  end do
endif
elseif((open_tab==3))then
if(whichtab==2)then
if(allocated(oldwa))then
  do m=1,ng; do mm=1,16
  if(gen(m)%wa(mm+23).ne.oldwa(m,mm))then
    call undo_redo(7,m,mm+23,0,0,gen(m)%wa(mm+23),oldwa(m,mm))  ! UNDO REDO SAVE
    gen(m)%wa(mm+23)=oldwa(m,mm)
  end if
  end do
  end do
endif
else
  do m=1,ng; do mm=1,16
  if(gen(m)%wa(mm+23).ne.oldwa(m,mm))then
    call undo_redo(7,m,mm+23,0,0,gen(m)%wa(mm+23),neww(m,mm))  ! UNDO REDO SAVE
  end if
  end do
  end do
endif
endif

  call update_view

if(whichtab==2)then
    call gtk_widget_destroy(ihwin2)
  call gtk_main_quit ()
    call gtk_widget_destroy(ihwin)
  call gtk_main_quit ()
endif

end subroutine my_accept_changes_tab
  
  
end module table_handlers

module Input_handlers

use gtk_hl

use gtk
use gtk_os_dependent

implicit none

type(c_ptr) :: window

logical, private :: file_is_changed = .FALSE.
character(len=240) :: filename=''   ! this is the chosen file

contains

subroutine my_destroy2(widget, gdata) bind(c)

  type(c_ptr), value :: widget, gdata
  integer(kind=c_int) :: ok

  call gtk_widget_destroy(window)
  call gtk_main_quit ()

end subroutine my_destroy2

subroutine do_open(widget, gdata) bind(c)
    
    use basic_handlers
    use general
    use genetic
    use io

    type(c_ptr), value :: widget, gdata

    type(c_ptr) :: c_string
    character(len=200) :: inln
    integer :: ios
    integer :: idxs
    character(len=120) :: filenamen=''

    call gtk_window_set_title(window, "Choose an input file"//c_null_char)
    c_string = gtk_file_chooser_get_filename(widget)
    call convert_c_string(c_string, filenamen)
    call g_free(c_string)

    idxs = index(filenamen, '/', .true.)+1
    call gtk_window_set_title(window, trim(filenamen(idxs:))//c_null_char)
    print*, filenamen
    filename=trim(filenamen)

  end subroutine do_open

subroutine display_file(widget, gdata) bind(c)

  type(c_ptr), value :: widget, gdata
  character(len=240) :: f1

  open(1, file='file1.dat')
    read(1,*), f1
  close(1)

  print*, "This is your file: ", f1

end subroutine display_file

end module Input_handlers


module event_handlers

use gtk, only: gtk_container_add, gtk_drawing_area_new, gtk_events_pending, gtk&
  &_main, gtk_main_iteration, gtk_main_iteration_do, gtk_widget_get_window, gtk_w&
  &idget_queue_draw, gtk_widget_show, gtk_window_new, gtk_window_set_default, gtk&
  &_window_set_default_size, gtk_window_set_title, TRUE, FALSE, c_null_ptr, c_null_char, &
  &GDK_COLORSPACE_RGB, GTK_WINDOW_TOPLEVEL, gtk_init, g_signal_connect, &
  &gtk_table_new, gtk_table_attach_defaults, gtk_container_add, gtk_button_new_with_label,&
  &gtk_widget_show_all, gtk_vbox_new, gtk_hbox_new, gtk_box_pack_start, gtk_spin_button_new,&
  &gtk_adjustment_new, gtk_spin_button_get_value, gtk_label_new, &
  &gtk_expander_new_with_mnemonic, gtk_expander_set_expanded, gtk_main_quit, &
  &gtk_toggle_button_new_with_label, gtk_toggle_button_get_active, gtk_notebook_new,&
  &gtk_notebook_append_page, gtk_text_view_new, gtk_text_view_get_buffer, gtk_text_buffer_set_text,&
  &gtk_scrolled_window_new, C_NEW_LINE, gtk_text_buffer_insert_at_cursor, gtk_statusbar_new,&
  &gtk_statusbar_push, gtk_statusbar_get_context_id, gtk_handle_box_new,&
  &CAIRO_STATUS_SUCCESS, CAIRO_STATUS_NO_MEMORY, CAIRO_STATUS_SURFACE_TYPE_MISMATCH,&
  &CAIRO_STATUS_WRITE_ERROR, gtk_button_new_with_mnemonic, gtk_link_button_new_with_label,&
  &gtk_toggle_button_new_with_mnemonic, gtk_label_new_with_mnemonic, &
  &gtk_window_set_mnemonics_visible, gtk_combo_box_text_new, &
  &gtk_combo_box_text_append_text, gtk_combo_box_text_get_active_text, &
  &gtk_combo_box_text_insert_text, gtk_spin_button_set_value, gtk_spin_button_update, &
! these functions are necessary for the cursor choices
       &gtk_button_new, gtk_check_button_new, gtk_container_add, gtk_ent&
       &ry_get_text, gtk_entry_get_text_length, gtk_entry_new, gtk_entry_set_text, gtk&
       &_main, gtk_main_quit, gtk_widget_destroy, gtk_toggle_button_get_active, gtk_to&
       &ggle_button_set_active, gtk_widget_show, gtk_widget_show_all, gtk_window_new, &
       &gtk_init, gtk_tree_path_new, gtk_tree_view_get_path_at_pos, gtk_notebook_get_current_page, &
       &gtk_tree_view_convert_widget_to_bin_window_coords, gtk_tree_selection_select_path, &
       &gtk_menu_popup, gtk_event_box_new, gtk_bin_get_child, gtk_widget_realize, &
       &gtk_notebook_prev_page, gtk_notebook_next_page, gtk_notebook_get_n_pages, &
       &gtk_check_button_new_with_label, gtk_toggle_button_toggled, gtk_toggle_button_set_active, &
       &gtk_toggle_button_get_active, gtk_label_set_justify, GTK_JUSTIFY_LEFT, gtk_misc_set_alignment, &
       &gtk_disable_setlocale

use cairo, only: cairo_create, cairo_destroy, cairo_paint, cairo_set_source, &
  &cairo_surface_write_to_png, cairo_get_target
  
use gdk, only: gdk_cairo_create, gdk_cairo_set_source_pixbuf
  
use gdk_pixbuf, only: gdk_pixbuf_get_n_channels, gdk_pixbuf_get_pixels, gdk_pix&
  &buf_get_rowstride, gdk_pixbuf_new

use g, only: g_usleep, alloca
use gdk_events
use iso_c_binding, only: c_int, c_ptr, c_char
  

implicit none
integer(c_int) :: run_status = TRUE
integer(c_int) :: boolresult
logical :: boolevent 
integer :: choice
  
contains

function delete_event_1 (widget, event, gdata) result(ret)  bind(c)

  use iso_c_binding, only: c_ptr, c_int
  integer(c_int)    :: ret
  type(c_ptr), value :: widget, event, gdata

  run_status = FALSE
  ret = FALSE
  call gtk_main_quit()

end function delete_event_1

subroutine pending_events ()
    
  do while(IAND(gtk_events_pending(), run_status) /= FALSE)
    boolresult = gtk_main_iteration_do(FALSE)
  end do

end subroutine pending_events

function draw_network_1 (widget, event, gdata) result(ret)  bind(c)

  use iso_c_binding, only: c_int, c_ptr
  use global_widgets
  use on_display_handlers
  implicit none
  integer(c_int)    :: ret
  type(c_ptr), value, intent(in) :: widget, event, gdata
  type(c_ptr) :: my_cairo_context
    
  my_cairo_context = gdk_cairo_create (gtk_widget_get_window(widget))
  call gdk_cairo_set_source_pixbuf(my_cairo_context, my_pixbuf, 0d0, 0d0) 
  call cairo_paint(my_cairo_context)    
  call cairo_destroy(my_cairo_context)
  ret = FALSE

end function draw_network_1

  ! Cursor position
function cursor_choice(widget, event, gdata) result(ret) bind(c)

  use on_display_handlers
  use global_widgets

  type(c_ptr), value :: widget, event, gdata
  integer(c_int)     :: ret
  type(c_ptr) :: drawing_area1
  type(gdkeventbutton), pointer :: cursor
  integer(kind=c_int), pointer :: x,y,button
  type(c_ptr) :: popupmenu, item, selection, bx, by
  type(c_ptr),target :: path
  integer :: choice1,which,whichelse,cpx,cpy,cpb


  cual=gtk_notebook_get_current_page(notebook)

  print*,"You have clicked here:", cursor%x, cursor%y, cual
  ret = FALSE  
  if (c_associated(event)) then

    call c_f_pointer(event,cursor)
    cpx=cursor%x;cpy=cursor%y;cpb=cursor%button
    call GRN_search(choice1,cpx,cpy,cpb,which,cual)

  end if

end function cursor_choice

subroutine choosers(filename3)   ! choice of the input file

  use Input_handlers
  implicit none

  type(c_ptr) :: base, jb, junk

  character(len=30), dimension(3) :: filters
  character(len=30), dimension(3) :: filtnames
  character(len=240) :: filename3
  filters(3) = "text/plain"
  filters(2) = "*.f90"
  filters(1) = "*"
  filtnames(3) = "Text files"
  filtnames(2) = "Fortran code"
  filtnames(1) = "All files"

  ! Initialize GTK
  call gtk_disable_setlocale ()
  call gtk_init()

  ! Create a window and a column boxi

  window = hl_gtk_window_new("Choose an input file from this or a child folder!"//c_null_char, &
       & destroy=c_funloc(my_destroy2))

  base = hl_gtk_box_new()
  call gtk_container_add(window, base)
  jb = hl_gtk_box_new(horizontal=TRUE, homogeneous=TRUE)
  call hl_gtk_box_pack(base, jb)
  junk = hl_gtk_file_chooser_button_new(title="Choose an input file"//c_null_char, &
       & filter=filters, filter_name=filtnames, file_set=c_funloc(do_open))

  call hl_gtk_box_pack(jb, junk)
  junk = hl_gtk_button_new("Print choices"//c_null_char, clicked=c_funloc(display_file))
  call hl_gtk_box_pack(base, junk)
  junk = hl_gtk_button_new("Display choices"//c_null_char, clicked=c_funloc(read_data_again_ini))
  call hl_gtk_box_pack(base, junk)
  junk = hl_gtk_button_new("Accept new choice"//c_null_char, clicked=c_funloc(my_destroy2))
  call hl_gtk_box_pack(base, junk)

  call gtk_widget_show_all(window)

  call gtk_main()

  filename3=filename

end subroutine choosers

subroutine read_data_again_ini(widget, gdata) bind(c)

  use basic_handlers
  type(c_ptr), value :: widget, gdata

  call read_data_again

end subroutine read_data_again_ini

subroutine quit_the_shit(widget, gdata) bind(c)

  use basic_handlers
  use gtk

  type(c_ptr), value :: widget, gdata

  call gtk_main_quit()

end subroutine quit_the_shit

subroutine name_enter(widget, gdata) bind(c)  ! THIS IS OUT OF USE?

  use basic_handlers
  use gtk

  type(c_ptr), value :: widget, gdata
  type(c_ptr) :: junk, jb, algn, eboxi, eboxi2, pbar, wini3
  character(len=240) :: filename2
  integer :: pd2

  call choosers(filename2)
    
  open(2, file='path.dat')
    read(2,*), pd2
  close(2)

  open(3, file='file01.dat')  ! save name of file for all applications. This way, you only need to choose anew when the input file is being changed
    write(3,*) "'",filename2,"'"
  close(3)

  open(1, file='file1.dat')  ! save name of file for all applications. This way, you only need to choose anew when the input file is being changed
    write(1,*) "'",filename2(pd2+2:len(filename2)),"'"
  close(1)
 
  if(len(filename2).le.0)then
    print*, "No file or file from outside the folder chosen. This won't work"
    print*, "We will use the default file instead"
    filename2= '1694390492_192954elli.e___May18_19__0.dat__________15____________70_________10000.dat' ! default file
  end if

end subroutine name_enter

function GRN_replot (widget, gdata ) result(ret)  bind(c)

  use iso_c_binding, only: c_ptr
  use global_widgets
  use on_display_handlers
  use basic_handlers
  use io
  use general
  use genetic

  implicit none
    
  integer(c_int)    :: ret, message_id
  type(c_ptr), value :: widget, gdata
  double complex :: cxc
  integer :: iterations
    
  cxc = gtk_spin_button_get_value (spinButton1) + &
      & (0d0, 1d0)*gtk_spin_button_get_value (spinButton2)
  iterations = INT(gtk_spin_button_get_value (spinButton3))
    
  write(string, '("c=",F8.6,"+i*",F8.6,"   ", I6, " iterations")') cxc, iterations
  call gtk_text_buffer_insert_at_cursor (buffer, &
         & string//C_NEW_LINE//c_null_char, -1_c_int)

  call read_data(1)!_again
  call update_view
  call system("./gne.e")

end function GRN_replot

function GRN_style (widget, gdata ) result(ret)  bind(c)

  use iso_c_binding, only: c_ptr
  use global_widgets
  use on_display_handlers
  use basic_handlers
  use io
  use general
  use genetic

  implicit none
    
  integer(c_int)    :: ret, message_id
  type(c_ptr), value :: widget, gdata

if(style1.le.0)then
  style1=1
  print*, "BACKGROUND IS SET WHITE"
else
  style1=-1
  print*, "BACKGROUND IS SET BLACK"
end if
print*, style1

call update_view

end function GRN_style

subroutine GRN_table_editable_connect (widget, gdata) bind(c)

  use iso_c_binding, only: c_ptr
  use global_widgets

  type(c_ptr), value :: widget, gdata
  call GRN_table_editable

end subroutine GRN_table_editable_connect

subroutine GRN_table_editable_connect2 (widget, gdata) bind(c)

  use iso_c_binding, only: c_ptr
  use global_widgets

  type(c_ptr), value :: widget, gdata
  call GRN_table_editable2

end subroutine GRN_table_editable_connect2

subroutine GRN_table_add_connect (widget, gdata) bind(c)

  use iso_c_binding, only: c_ptr
  use global_widgets

  type(c_ptr), value :: widget, gdata

  call GRN_add

end subroutine GRN_table_add_connect

subroutine move_nodes_connect (widget, gdata) bind(c)

  use iso_c_binding, only: c_ptr
  use global_widgets

  type(c_ptr), value :: widget, gdata

  print*, "LET'S REARRANGE THE NETWORK."
  print*, "TO RETURN TO GRN MODIFICATION, PRESS THE 'RETURN...' BUTTON."

  call move_nodes_control

end subroutine move_nodes_connect

subroutine GRN_table_editable  ! mainly taken from the gtk-fortran project examples; modified

  use table_handlers
  use basic_handlers
  use some_widgets
  use genetic
  use general
  use io
  use gtk, only: gtk_notebook_prev_page, gtk_notebook_next_page, gtk_notebook_get_n_pages, gtk_notebook_get_current_page

  implicit none

  type(c_ptr) :: tabtable, tabbox
  character(len=35) :: line
  character*3 :: strg
  integer :: ij, ltr, jk, p, kl, m, mm
  integer(kind=type_kind),dimension(ng+4) :: ctypes
  integer(kind=type_kind),dimension(24) :: c2types
  integer(kind=type_kind),dimension(17) :: c3types
  integer(kind=type_kind),dimension(17) :: c1types
  character(len=20), dimension(ng+4) :: titles
  character(len=20), dimension(24) :: titles2
  character(len=20), dimension(17) :: titles3
  character(len=20), dimension(17) :: titles1
  integer(kind=c_int), dimension(ng+4) :: sortable, editable
  integer(kind=c_int), dimension(17) :: sortable1, editable1
  integer(kind=c_int), dimension(24) :: sortable2, editable2
  integer(kind=c_int), dimension(17) :: sortable3, editable3
  integer(kind=c_int), allocatable, dimension(:) :: colnos
  character(len=12*ng) :: typlist
  integer(kind=c_int) :: dofi
  real*4 :: dofl
  integer :: coln
  integer :: nopars(23)

  whichtab=1
  nopars=(/5,6,7,8,9,10,11,12,13,14,15,16,20,21,22,23,24,25,26,27,28,34,35/)

  call gtk_disable_setlocale ()
  call gtk_init()

open_tab=gtk_notebook_get_current_page(notebook)


if(open_tab==0)then
  ihwin=hl_gtk_window_new("TF interaction matrix: change if needed!"//c_null_char, destroy=c_funloc(my_destroy_tab))
elseif(open_tab==1)then
  ihwin=hl_gtk_window_new("R catalytic interaction matrix: change if needed!"//c_null_char, destroy=c_funloc(my_destroy_tab))
elseif(open_tab==2)then
  ihwin=hl_gtk_window_new("E biomechanical regulation matrix: change if needed!"//c_null_char, destroy=c_funloc(my_destroy_tab))
elseif(open_tab==3)then
  ihwin=hl_gtk_window_new("C cell behaviour regulation matrix: change if needed!"//c_null_char, destroy=c_funloc(my_destroy_tab))
end if

  base = hl_gtk_box_new()
  call gtk_container_add(ihwin, base)

  editable=1 ! the values can be edited
  sortable=0
  editable1=1
  sortable1=0
  editable2=1
  sortable2=0
  editable3=1
  sortable3=0

if(open_tab==0)then
  ctypes(1)=G_TYPE_INT
  do p=2,ng+4
    ctypes(p)=G_TYPE_FLOAT
  end do

  titles(1)= "GENE"
  do p=2,ng+1
    write (strg,"(I3.1)") p-1
    titles(p) = strg
  end do
  titles(ng+2)="Diffusion rate"
  titles(ng+3)="Degradation rate"
  titles(ng+4)="KINDOF"
  coln=ng+5  

  ihlist = hl_gtk_tree_new(ihscrollcontain, types=ctypes, &
       & changed=c_funloc(list_select),&
       &  multiple=TRUE, height=25*(ng+1), swidth=50*(coln), titles=titles, editable=editable, sortable=sortable )

elseif(open_tab==1)then
  c1types(1)=G_TYPE_INT
  do p=2,ng+1
    c1types(p)=G_TYPE_FLOAT
  end do

  titles1(1)= "GENE"
  do p=2,ng+1
    write (strg,"(I3.1)") p-1
    titles1(p) = strg
  end do

  coln=ng+2  

  ihlist = hl_gtk_tree_new(ihscrollcontain, types=c1types, &
       & changed=c_funloc(list_select),&
       &  multiple=TRUE, height=25*(ng+1), swidth=50*(coln), titles=titles1, editable=editable1, sortable=sortable1 )

elseif(open_tab==2)then

  c2types(1)=G_TYPE_INT
  do p=2,24
    c2types(p)=G_TYPE_FLOAT
  end do

  titles2(1)= "Form"
  do p=2,24
    titles2(p) = nodeparams(nopars(p-1))!strg
  end do
  coln=25

  ihlist = hl_gtk_tree_new(ihscrollcontain, types=c2types, &
       & changed=c_funloc(list_select),&
       &  multiple=TRUE, height=25*(ng+1), swidth=50*(coln), titles=titles2, editable=editable2, sortable=sortable2 )

elseif(open_tab==3)then

  c3types(1)=G_TYPE_INT
  do p=2,17
    c3types(p)=G_TYPE_FLOAT
  end do

  titles3(1)= "Form"
  do p=2,17
    titles3(p) = cellparams(p-1)
  end do
  coln=18

  ihlist = hl_gtk_tree_new(ihscrollcontain, types=c3types, &
       & changed=c_funloc(list_select),&
       &  multiple=TRUE, height=25*(ng+1), swidth=50*(coln), titles=titles3, editable=editable3, sortable=sortable3 )

endif

if(open_tab==0)then
  do ij=1,ng
    call hl_gtk_tree_ins(ihlist, row = (/ -1 /))
      write(line,"('List entry number ',I0)") ij
    ltr=len_trim(line)+1
    line(ltr:ltr)=c_null_char
    call hl_gtk_tree_set_cell(ihlist, absrow=ij-1, col=0, ivalue=ij)
    do kl=1,ng+1
      dofl=gen(ij)%w(kl)
      call hl_gtk_tree_set_cell(ihlist, absrow=ij-1, col=kl, fvalue=dofl)
    end do
    dofl=gen(ij)%diffu
    kl=ng+1
    call hl_gtk_tree_set_cell(ihlist, absrow=ij-1, col=kl, fvalue=dofl)
    dofl=gen(ij)%mu
    kl=ng+2
    call hl_gtk_tree_set_cell(ihlist, absrow=ij-1, col=kl, fvalue=dofl)
    dofl=int(gen(ij)%kindof)
    kl=ng+3
    call hl_gtk_tree_set_cell(ihlist, absrow=ij-1, col=kl, fvalue=dofl)
  end do
elseif(open_tab==1)then
  do ij=1,ng
    call hl_gtk_tree_ins(ihlist, row = (/ -1 /))
      write(line,"('List entry number ',I0)") ij
    ltr=len_trim(line)+1
    line(ltr:ltr)=c_null_char
    call hl_gtk_tree_set_cell(ihlist, absrow=ij-1, col=0, ivalue=ij)
    do kl=2,ng+1
      dofl=gen(ij)%w(kl)
      call hl_gtk_tree_set_cell(ihlist, absrow=ij-1, col=kl, fvalue=dofl)
    end do
  end do
elseif(open_tab==2)then

  do ij=1,ng
    call hl_gtk_tree_ins(ihlist, row = (/ -1 /))
      write(line,"('List entry number ',I0)") ij
    ltr=len_trim(line)+1
    line(ltr:ltr)=c_null_char
    call hl_gtk_tree_set_cell(ihlist, absrow=ij-1, col=0, ivalue=ij)
    do kl=1,23
      dofl=gen(ij)%wa(nopars(kl))
      call hl_gtk_tree_set_cell(ihlist, absrow=ij-1, col=kl, fvalue=dofl)
    end do
  end do
elseif(open_tab==3)then

  do ij=1,ng
    call hl_gtk_tree_ins(ihlist, row = (/ -1 /))
      write(line,"('List entry number ',I0)") ij
    ltr=len_trim(line)+1
    line(ltr:ltr)=c_null_char
    call hl_gtk_tree_set_cell(ihlist, absrow=ij-1, col=0, ivalue=ij)
    do kl=1,16
      dofl=gen(ij)%wa(nparam_per_node+kl)
      call hl_gtk_tree_set_cell(ihlist, absrow=ij-1, col=kl, fvalue=dofl)
    end do
  end do
endif

  call hl_gtk_box_pack(base, ihscrollcontain)

  lbl = gtk_label_new("Edit any field to change the interaction of COLUMN on ROW"//c_null_char)
  call hl_gtk_box_pack(base, lbl)

  tabtable = gtk_table_new (1_c_int, 3_c_int, TRUE)
  tabbox = gtk_vbox_new (FALSE, 10_c_int);
  call gtk_box_pack_start (tabbox, tabtable, FALSE, FALSE, 0_c_int)

  qbut = hl_gtk_button_new("Accept Changes"//c_null_char, clicked=c_funloc(my_accept_changes_tab))
  call gtk_table_attach_defaults(tabtable, qbut, 0_c_int, 1_c_int, 0_c_int, 1_c_int)
  call gtk_container_add (tabtable, qbut)

  qbut = gtk_button_new_with_mnemonic ("_Show original_"//c_null_char)
  call g_signal_connect (qbut, "clicked"//c_null_char, c_funloc(GRN_table_editable_connect2))
  call gtk_table_attach_defaults(tabtable, qbut, 1_c_int, 2_c_int, 0_c_int, 1_c_int)
  call gtk_container_add (tabtable, qbut)

  qbut = hl_gtk_button_new("Quit"//c_null_char, clicked=c_funloc(my_destroy_tab))
  call gtk_table_attach_defaults(tabtable, qbut, 2_c_int, 3_c_int, 0_c_int, 1_c_int)
  call gtk_container_add (tabtable, qbut)
  call gtk_container_add (base, tabbox)
  call gtk_widget_show_all(ihwin)

  call gtk_main()

end subroutine GRN_table_editable

subroutine GRN_table_editable2  ! mainly taken from the gtk-fortran project examples; modified

  use table_handlers
  use basic_handlers
  use some_widgets
  use genetic
  use general
  use io
  use gtk, only: gtk_notebook_prev_page, gtk_notebook_next_page, gtk_notebook_get_n_pages, gtk_notebook_get_current_page

  implicit none

  type(c_ptr) :: tabtable, tabbox
  character(len=35) :: line
  character*3 :: strg
  integer :: ij, ltr, jk, p, kl, m, mm
  integer(kind=type_kind),dimension(ng+4) :: ctypes
  integer(kind=type_kind),dimension(24) :: c2types
  integer(kind=type_kind),dimension(17) :: c3types
  integer(kind=type_kind),dimension(17) :: c1types
  character(len=20), dimension(ng+4) :: titles
  character(len=20), dimension(24) :: titles2
  character(len=20), dimension(17) :: titles3
  character(len=20), dimension(17) :: titles1
  integer(kind=c_int), dimension(ng+4) :: sortable, editable
  integer(kind=c_int), dimension(17) :: sortable1, editable1
  integer(kind=c_int), dimension(24) :: sortable2, editable2
  integer(kind=c_int), dimension(17) :: sortable3, editable3
  integer(kind=c_int), allocatable, dimension(:) :: colnos
  character(len=12*ng) :: typlist
  integer(kind=c_int) :: dofi
  real*4 :: dofl
  integer :: coln
  integer :: nopars(23)

  whichtab=2
  nopars=(/5,6,7,8,9,10,11,12,13,14,15,16,20,21,22,23,24,25,26,27,28,34,35/)

  call gtk_disable_setlocale ()
  call gtk_init()

open_tab=gtk_notebook_get_current_page(notebook)


if(open_tab==0)then
  ihwin2=hl_gtk_window_new("TF interaction matrix: change if needed!"//c_null_char, destroy=c_funloc(my_destroy_tab2))
elseif(open_tab==1)then
  ihwin2=hl_gtk_window_new("R catalytic interaction matrix: change if needed!"//c_null_char, destroy=c_funloc(my_destroy_tab2))
elseif(open_tab==2)then
  ihwin2=hl_gtk_window_new("E biomechanical regulation matrix: change if needed!"//c_null_char, destroy=c_funloc(my_destroy_tab2))
elseif(open_tab==3)then
  ihwin2=hl_gtk_window_new("C cell behaviour regulation matrix: change if needed!"//c_null_char, destroy=c_funloc(my_destroy_tab2))
end if

  base = hl_gtk_box_new()
  call gtk_container_add(ihwin2, base)

  editable=1 ! the values can be edited
  sortable=0
  editable1=1
  sortable1=0
  editable2=1
  sortable2=0
  editable3=1
  sortable3=0

if(open_tab==0)then
  ctypes(1)=G_TYPE_INT
  do p=2,ng+4
    ctypes(p)=G_TYPE_FLOAT
  end do

  titles(1)= "GENE"
  do p=2,ng+1
    write (strg,"(I3.1)") p-1
    titles(p) = strg
  end do
  titles(ng+2)="Diffusion rate"
  titles(ng+3)="Degradation rate"
  titles(ng+4)="KINDOF"
  coln=ng+5  

  ihlist = hl_gtk_tree_new(ihscrollcontain, types=ctypes, &
       & changed=c_funloc(list_select),&
       &  multiple=TRUE, height=25*(ng+1), swidth=50*(coln), titles=titles, editable=editable, sortable=sortable )

elseif(open_tab==1)then
  c1types(1)=G_TYPE_INT
  do p=2,ng+1
    c1types(p)=G_TYPE_FLOAT
  end do

  titles1(1)= "GENE"
  do p=2,ng+1
    write (strg,"(I3.1)") p-1
    titles1(p) = strg
  end do

  coln=ng+2  

  ihlist = hl_gtk_tree_new(ihscrollcontain, types=c1types, &
       & changed=c_funloc(list_select),&
       &  multiple=TRUE, height=25*(ng+1), swidth=50*(coln), titles=titles1, editable=editable1, sortable=sortable1 )

elseif(open_tab==2)then

  c2types(1)=G_TYPE_INT
  do p=2,24
    c2types(p)=G_TYPE_FLOAT
  end do

  titles2(1)= "Form"
  do p=2,24
    titles2(p) = nodeparams(nopars(p-1))!strg
  end do
  coln=25

  ihlist = hl_gtk_tree_new(ihscrollcontain, types=c2types, &
       & changed=c_funloc(list_select),&
       &  multiple=TRUE, height=25*(ng+1), swidth=50*(coln), titles=titles2, editable=editable2, sortable=sortable2 )

elseif(open_tab==3)then

  c3types(1)=G_TYPE_INT
  do p=2,17
    c3types(p)=G_TYPE_FLOAT
  end do

  titles3(1)= "Form"
  do p=2,17
    titles3(p) = cellparams(p-1)
  end do
  coln=18

  ihlist = hl_gtk_tree_new(ihscrollcontain, types=c3types, &
       & changed=c_funloc(list_select),&
       &  multiple=TRUE, height=25*(ng+1), swidth=50*(coln), titles=titles3, editable=editable3, sortable=sortable3 )

endif

if(open_tab==0)then
  do ij=1,ng
    call hl_gtk_tree_ins(ihlist, row = (/ -1 /))
      write(line,"('List entry number ',I0)") ij
    ltr=len_trim(line)+1
    line(ltr:ltr)=c_null_char
    call hl_gtk_tree_set_cell(ihlist, absrow=ij-1, col=0, ivalue=ij)
    do kl=1,ng+1
      dofl=oldw(ij,kl)!gen(ij)%w(kl)
      call hl_gtk_tree_set_cell(ihlist, absrow=ij-1, col=kl, fvalue=dofl)
    end do
    dofl=oldw(ij,kl+1)!gen(ij)%diffu
    kl=ng+1
    call hl_gtk_tree_set_cell(ihlist, absrow=ij-1, col=kl, fvalue=dofl)
    dofl=oldw(ij,kl+2)!gen(ij)%mu
    kl=ng+2
    call hl_gtk_tree_set_cell(ihlist, absrow=ij-1, col=kl, fvalue=dofl)
    dofl=oldw(ij,kl+3)!int(gen(ij)%kindof)
    kl=ng+3
    call hl_gtk_tree_set_cell(ihlist, absrow=ij-1, col=kl, fvalue=dofl)
  end do
elseif(open_tab==1)then
  do ij=1,ng
    call hl_gtk_tree_ins(ihlist, row = (/ -1 /))
      write(line,"('List entry number ',I0)") ij
    ltr=len_trim(line)+1
    line(ltr:ltr)=c_null_char
    call hl_gtk_tree_set_cell(ihlist, absrow=ij-1, col=0, ivalue=ij)
    do kl=1,ng+1
      dofl=oldw(ij,kl)!dofl=gen(ij)%w(kl)
      call hl_gtk_tree_set_cell(ihlist, absrow=ij-1, col=kl, fvalue=dofl)
    end do
  end do
elseif(open_tab==2)then

  do ij=1,ng
    call hl_gtk_tree_ins(ihlist, row = (/ -1 /))
      write(line,"('List entry number ',I0)") ij
    ltr=len_trim(line)+1
    line(ltr:ltr)=c_null_char
    call hl_gtk_tree_set_cell(ihlist, absrow=ij-1, col=0, ivalue=ij)
    do kl=1,23
      dofl=oldwa(ij,kl)!dofl=gen(ij)%wa(nopars(kl))
      call hl_gtk_tree_set_cell(ihlist, absrow=ij-1, col=kl, fvalue=dofl)
    end do
  end do
elseif(open_tab==3)then

  do ij=1,ng
    call hl_gtk_tree_ins(ihlist, row = (/ -1 /))
      write(line,"('List entry number ',I0)") ij
    ltr=len_trim(line)+1
    line(ltr:ltr)=c_null_char
    call hl_gtk_tree_set_cell(ihlist, absrow=ij-1, col=0, ivalue=ij)
    do kl=1,16
      dofl=oldwa(ij,kl)!dofl=gen(ij)%wa(nparam_per_node+kl)
      call hl_gtk_tree_set_cell(ihlist, absrow=ij-1, col=kl, fvalue=dofl)
    end do
  end do
endif

  call hl_gtk_box_pack(base, ihscrollcontain)

  lbl = gtk_label_new("Edit any field to change the interaction of COLUMN on ROW"//c_null_char)
  call hl_gtk_box_pack(base, lbl)

  tabtable = gtk_table_new (1_c_int, 2_c_int, TRUE)
  tabbox = gtk_vbox_new (FALSE, 10_c_int);
  call gtk_box_pack_start (tabbox, tabtable, FALSE, FALSE, 0_c_int)

  qbut = hl_gtk_button_new("Accept Changes"//c_null_char, clicked=c_funloc(my_accept_changes_tab))
  call gtk_table_attach_defaults(tabtable, qbut, 0_c_int, 1_c_int, 0_c_int, 1_c_int)
  call gtk_container_add (tabtable, qbut)

  qbut = hl_gtk_button_new("Quit"//c_null_char, clicked=c_funloc(my_destroy_tab2))
  call gtk_table_attach_defaults(tabtable, qbut, 1_c_int, 2_c_int, 0_c_int, 1_c_int)
  call gtk_container_add (tabtable, qbut)
  call gtk_container_add (base, tabbox)
  call gtk_widget_show_all(ihwin2)

  call gtk_main()

end subroutine GRN_table_editable2


function my_accept_add (widget, gdata ) result(ret)  bind(c)

  use gtk, only: gtk_notebook_prev_page, gtk_notebook_next_page, gtk_notebook_get_n_pages, gtk_notebook_get_current_page
  use iso_c_binding, only: c_ptr, c_long
  use global_widgets
  use basic_handlers
  use genetic

  implicit none

  integer(c_int)    :: ret
  type(c_ptr), value :: widget, gdata
  integer(kind=c_int) :: ag, m, mm, mmm
  real*8 :: ges1(ng+1,ng+4), geswa1(ng+1,28+16), gesw(ng,ng*ng,3), gespos(ng+1,6), geskpos(ntipusadh,2)
  real*8 :: gex2(nda,ng+1),agex2(nda,ng+1)
  integer :: ges2(ng+1,4), geswa2(ng+1,7),npag2(nga+1),whonpag2(nga+1,ng+1)
  integer :: ges3(ng+1,ng+1,2)
  character*12 :: genenamesave(ng+1)
  type(genes),allocatable :: ges(:)
  integer :: mn
  integer,parameter :: seed = 1
  integer, dimension(8) :: dt

  ges1=0.0d0
  ges2=0
  ges2(ng+1,4)=1
  ges3=0
  geswa1=0.0d0
  geswa2=0
  gesw=0.0d0
  gespos=0.0d0; geskpos=0.0d0
  gex2=0.0d0; agex2=0.0d0
  npag2=0; whonpag2=0
  genenamesave=""

  do m=1,ng
    do mm=1,ng
      ges1(m,mm)=gen(m)%w(mm)
    end do
    do mm=1,nga
      npag2(mm)=npag(mm); whonpag2(mm,m)=whonpag(mm,m)
    end do

    ges1(m,ng+2)=gen(m)%diffu
    ges1(m,ng+3)=gen(m)%idiffu
    ges1(m,ng+4)=gen(m)%mu
    ges2(m,1)=gen(m)%nww
    ges2(m,2)=gen(m)%npre
    ges2(m,3)=gen(m)%npost
    ges2(m,4)=gen(m)%kindof
    genenamesave(m)=genename(m)
    gex2(:,m)=gex(:,m); agex2(:,m)=agex2(:,m)
    if(gen(m)%npre.ne.0)then
      do mm=1,gen(m)%npre
        ges3(m,mm,1)=gen(m)%pre(mm)
      end do
    end if
    if(gen(m)%npost.ne.0)then
      do mm=1,gen(m)%npost
        ges3(m,mm,2)=gen(m)%post(mm)
      end do
    end if

    do mm=1,28
      geswa1(m,mm)=gen(m)%wa(mm)
    end do
    do mm=29,35
      geswa2(m,mm-28)=gen(m)%wa(mm)
    end do
    do mm=1,16
      geswa1(m,mm+28)=gen(m)%wa(mm+nparam_per_node)
    end do
    do mm=1,gen(m)%nww
      gesw(m,mm,:)=gen(m)%ww(mm,:)
    end do

    gespos(m,:)=npos(m,:)

  end do

  ag = int(gtk_spin_button_get_value(spinButtona))
  ng=ng+1
  param(19)=ng

  call initiate_gene

  do m=1,ng
    if (allocated(gen(m)%w)) deallocate(gen(m)%w)
    allocate(gen(m)%w(ng))
    if (allocated(gen(m)%pre)) deallocate(gen(m)%pre)
    allocate(gen(m)%pre(ges2(m,2)))
    if (allocated(gen(m)%post)) deallocate(gen(m)%post)
    allocate(gen(m)%post(ges2(m,3)))
    if (allocated(gen(m)%ww)) deallocate(gen(m)%ww)
    allocate(gen(m)%ww(ng*ng,3))
    if (allocated(npos)) deallocate(npos)
    allocate(npos(ng,6))
    if (allocated(genename)) deallocate(genename)
    allocate(genename(ng))
    genename(ng)=""
print*, "UNDONE", ng
    if(allocated(gex)) deallocate(gex)
    allocate(gex(nda,ng))
    if(allocated(agex)) deallocate(agex)
    allocate(agex(nda,ng))
    if (allocated(npag)) deallocate(npag)
    allocate(npag(nga))
    npag=0
    if (allocated(whonpag)) deallocate(whonpag)  
    allocate(whonpag(nga,ng))
    whonpag=0
  end do
  gex(:,ng)=0.0d0; agex(:,ng)=0.0d0
  do m=1,ng-1
    do mm=1,ng-1
      gen(m)%w(mm)=ges1(m,mm)
    end do
    do mm=1,nga
      npag(mm)=npag2(mm); whonpag(mm,m)=whonpag2(mm,m)
    end do
    gen(m)%diffu=ges1(m,ng+1)
    gen(m)%idiffu=ges1(m,ng+2)
    gen(m)%mu=ges1(m,ng+3)
    gen(m)%nww=ges2(m,1)
    gen(m)%npre=ges2(m,2)
    gen(m)%npost=ges2(m,3)
    gen(m)%kindof=ges2(m,4)
    npos(m,:)=gespos(m,:)
    gex(:,m)=gex2(:,m)
    agex(:,m)=agex2(:,m)
    genename(m)=genenamesave(m)

    if(gen(m)%npre.ne.0)then
      do mm=1,gen(m)%npre
        gen(m)%pre(mm)=ges3(m,mm,1)
      end do
    end if
    if(gen(m)%npost.ne.0)then
      do mm=1,gen(m)%npost
        mn= ges3(m,mm,2)
        gen(m)%post(mm)=mn   
      end do
    end if

    do mm=1,28
      gen(m)%wa(mm)=geswa1(m,mm)
    end do
    do mm=1,16
      gen(m)%wa(mm+nparam_per_node)=geswa1(m,mm+28)
    end do
    do mm=29,35
      gen(m)%wa(mm)=geswa2(m,mm-28)
    end do
    do mm=1,gen(m)%nww
      gen(m)%ww(mm,:)=gesw(m,mm,:)
    end do
  end do

  if((ag.ne.0).and.(ag.lt.ng)) then
    print*, "Gene ", ag, "has been duplicated."

    do mm=1,ng-1
      gen(ng)%w(mm)=ges1(ag,mm)
    end do
    do mm=1,ng-1
      gen(mm)%w(ng)=ges1(mm,ag)
    end do
    gen(ng)%diffu=ges1(ag,ng+1)
    gen(ng)%idiffu=ges1(ag,ng+2)
    gen(ng)%mu=ges1(ag,ng+3)
    gen(ng)%nww=ges2(ag,1)
    gen(ng)%npre=ges2(ag,2)
    gen(ng)%npost=ges2(ag,3)
    gen(ng)%kindof=ges2(ag,4)
    if(gen(ag)%w(ag).ne.0)then
      gen(ag)%w(ng)=gen(ag)%w(ag)
      gen(ng)%w(ng)=gen(ag)%w(ag)
      gen(ng)%w(ag)=gen(ag)%w(ag)
    endif
    if(gen(ag)%npre.ne.0)then
      do mm=1,gen(ag)%npre
        gen(ng)%pre(mm)=ges3(ag,mm,1)
        do m=1,ng-1
          do mmm=1,gen(m)%nww  ! if ag is a post catalyzed by m
            if((gen(m)%ww(mmm,1)==gen(ng)%pre(mm)).and.(gen(m)%ww(mmm,2)==ag))then
              gen(m)%nww=gen(m)%nww+1
              gen(m)%ww(gen(m)%nww,1)=gen(ng)%pre(mm)
              gen(m)%ww(gen(m)%nww,2)=ng
              gen(m)%ww(gen(m)%nww,3)=gen(m)%ww(mmm,3)
              if(gen(m)%ww(mm,3)==ag)then ! if ag is its own catalyzator
                gen(ng)%nww=gen(ng)%nww+1
                gen(ng)%ww(gen(ng)%nww,1)=gen(m)%ww(mm,1)
                gen(ng)%ww(gen(ng)%nww,2)=gen(m)%ww(mm,2)
                gen(ng)%ww(gen(ng)%nww,3)=ng
              end if
            end if
          end do
        end do
      end do
    end if
    if(gen(ag)%npost.ne.0)then
      do mm=1,gen(ag)%npost
        gen(ng)%post(mm)=ges3(ag,mm,2)
        do m=1,ng-1
          do mmm=1,gen(m)%nww  ! if ag is a pre catalyzed by m
            if((gen(m)%ww(mmm,2)==gen(ng)%post(mm)).and.(gen(m)%ww(mmm,1)==ag))then
              gen(m)%nww=gen(m)%nww+1
              gen(m)%ww(gen(m)%nww,2)=gen(ng)%post(mm)
              gen(m)%ww(gen(m)%nww,1)=ng
              gen(m)%ww(gen(m)%nww,3)=gen(m)%ww(mmm,3)
              if(gen(m)%ww(mm,3)==ag)then ! if ag is its own catalyzator
                gen(ng)%nww=gen(ng)%nww+1
                gen(ng)%ww(gen(ng)%nww,1)=gen(m)%ww(mm,1)
                gen(ng)%ww(gen(ng)%nww,2)=gen(m)%ww(mm,2)
                gen(ng)%ww(gen(ng)%nww,3)=ng
              end if
            end if
          end do
        end do
      end do
    end if
    do mm=1,28
      gen(ng)%wa(mm)=geswa1(ag,mm)
    end do
    do mm=1,16
      gen(ng)%wa(mm+nparam_per_node)=geswa1(ag,mm+28)
    end do
    do mm=29,35
      gen(ng)%wa(mm)=geswa2(ag,mm-28)
    end do
    do mm=1,gen(m)%nww
      gen(ng)%ww(mm,:)=gesw(ag,mm,:)
    end do
    gen(ng)%wa(1)=gen(ag)%wa(1)
      
    print*, "This is your new gene:", ng
    print*, "Have fun with it! :) "
    print*, ng, gen(ng)%npre, gen(ng)%npost

    call undo_redo(2,ag,0,0,0,0d0,0d0)

  else
      
    print*, "No gene chosen"

  end if

  call iniread
  call update_npag

  call update_view

end function my_accept_add

function my_accept_del (widget, gdata ) result(ret)  bind(c)

  use iso_c_binding, only: c_ptr, c_long
  use global_widgets
  use basic_handlers
  use genetic
  use general
  implicit none
    
  integer(c_int)    :: ret
  type(c_ptr), value :: widget, gdata
  integer(kind=c_int) :: ag, m, mm,mmm,mmmm,mmmmm,mmmmmm
  real*8 :: ges1(ng,ng+3)
  integer :: ges2(ng,4)
  integer :: ges3(ng,ng,2)
  integer :: p,mn1, mn2, mt
  character*152 :: name_rec
  character(len=8) :: fmtr
  character*1 :: ep
  fmtr='(I1.1)'

  ag = int(gtk_spin_button_get_value(spinButtona2))
  m=ag

recf=recf+1
write (ep,fmtr) recf !!! RZ
record(recf)="Recordfile_"//trim(ep)
print*, "Record coniguration as: ", record(recf)
if(recf==9)then
  do i=2,8
    record(i)=record(i+1)
  end do
  recf=8
endif
call writesnap
name_rec=trim("cp fort.1 "//adjustl(record(recf))) 
call system(name_rec)

call undo_redo(2,recf*(-1),0,0,0,0d0,0d0)

call del(ag)

call update_view
print*, "DEFINE", p, ng, ag !R!
do p=ag,ng-1
print*, "DEFINE", p, ng, ag !R!
genename(p)=genename(p+1)
enddo

ng=ng-1

call update_view

end function my_accept_del

function my_circular (widget, gdata ) result(ret)  bind(c)

  use iso_c_binding, only: c_ptr, c_long
  use global_widgets
  use basic_handlers
  use genetic
  use general
  implicit none
    
  integer(c_int)    :: ret
  type(c_ptr), value :: widget, gdata

  togglen=0

  call update_view

end function my_circular

function my_arrange_manually (widget, gdata ) result(ret)  bind(c)

  use iso_c_binding, only: c_ptr, c_long
  use global_widgets
  use basic_handlers
  use genetic
  use general
  implicit none
    
  integer(c_int)    :: ret
  type(c_ptr), value :: widget, gdata

  togglen=-1

end function my_arrange_manually

function my_arrange_return (widget, gdata ) result(ret)  bind(c)

  use iso_c_binding, only: c_ptr, c_long
  use global_widgets
  use basic_handlers
  use genetic
  use general
  implicit none
    
  integer(c_int)    :: ret
  type(c_ptr), value :: widget, gdata

  togglen=1
  rb_arrange_tog=0

end function my_arrange_return

function my_arrange_saver (widget, gdata ) result(ret)  bind(c)

  use iso_c_binding, only: c_ptr, c_long
  use global_widgets
  use basic_handlers
  use genetic
  use general
  implicit none
    
  integer(c_int)    :: ret
  type(c_ptr), value :: widget, gdata
  integer:: m

  if(allocated(savepos)) deallocate(savepos)
  allocate(savepos(ng+ntipusadh,6))
  savepos=0d0
  saveposn=ng

  do m=1,saveposn
    savepos(m,:)=npos(m,:)
  end do
  if(ntipusadh.ge.0)then
    do m=1,ntipusadh
      savepos(m+saveposn,1:2)=kpos(m,1:2)
    end do
  end if

end function my_arrange_saver

function my_arrange_retriever (widget, gdata ) result(ret)  bind(c)

  use iso_c_binding, only: c_ptr, c_long
  use global_widgets
  use basic_handlers
  use genetic
  use general
  implicit none
   
  integer(c_int)    :: ret
  type(c_ptr), value :: widget, gdata
  integer:: m

  do m=1,saveposn
    npos(m,:)=savepos(m,:)
  end do
  if(ntipusadh.ge.0)then
    do m=1,ntipusadh
      kpos(m,1:2)=savepos(m+saveposn,1:2)
    end do
  end if

  call update_view

end function my_arrange_retriever

subroutine GRN_add

  use table_handlers
  use genetic
  use global_widgets
  use iso_c_binding, only: c_ptr
  use gtk, only: gtk_notebook_prev_page, gtk_notebook_next_page, gtk_notebook_get_n_pages, gtk_notebook_get_current_page
   
  implicit none

  type(c_ptr) :: labela1, labela2, labela3, buttona, buttona2, tablea
  real*8 :: ngr

  ngr=ng

  if(yetthere5==1) call gtk_widget_destroy(add_window)
  add_window = gtk_window_new (GTK_WINDOW_TOPLEVEL)
  yetthere5=1
  call gtk_window_set_default_size(add_window, 50, 100)
  call gtk_window_set_title(add_window, "Add new genes"//c_null_char)

  labela1 = gtk_label_new("Duplicate an existing gene"//c_null_char)
  labela2 = gtk_label_new("or create a default one with '0'."//c_null_char)
  labela3 = gtk_label_new("Delete a gene."//c_null_char)
  spinButtona = gtk_spin_button_new (gtk_adjustment_new(0d0,0d0,ngr,&
       & 1.0d0,0.5d0,0d0),0.05d0, 0_c_int)
  spinButtona2 = gtk_spin_button_new (gtk_adjustment_new(0d0,1d0,ngr,&
       & 1.0d0,0.5d0,0d0),0.05d0, 0_c_int)
  tablea = gtk_table_new (5_c_int, 2_c_int, TRUE)
  buttona = gtk_button_new_with_mnemonic ("Add it NOW"//c_null_char)
  buttona2 = gtk_button_new_with_mnemonic ("Delete it NOW"//c_null_char)

  call g_signal_connect (buttona, "clicked"//c_null_char, c_funloc(my_accept_add))
  call g_signal_connect (buttona2, "clicked"//c_null_char, c_funloc(my_accept_del))
  call gtk_table_attach_defaults(tablea, labela1, 0_c_int, 2_c_int, 0_c_int, 1_c_int)
  call gtk_table_attach_defaults(tablea, labela2, 0_c_int, 2_c_int, 1_c_int, 2_c_int)
  call gtk_table_attach_defaults(tablea, labela3, 0_c_int, 2_c_int, 4_c_int, 5_c_int)
  call gtk_table_attach_defaults(tablea, spinButtona, 0_c_int, 2_c_int, 2_c_int, 3_c_int)
  call gtk_table_attach_defaults(tablea, spinButtona2, 0_c_int, 2_c_int, 5_c_int, 6_c_int)
  call gtk_table_attach_defaults(tablea, buttona, 0_c_int, 2_c_int, 3_c_int, 4_c_int)
  call gtk_table_attach_defaults(tablea, buttona2, 0_c_int, 2_c_int, 6_c_int, 7_c_int)
  call gtk_container_add (add_window, tablea)

  call gtk_widget_show_all (add_window)

end subroutine GRN_add

subroutine move_nodes_control

  use table_handlers
  use genetic
  use global_widgets
  use iso_c_binding, only: c_ptr

  implicit none

  type(c_ptr) :: buttonn, buttonn1, buttonn2, buttonn3, buttonn4, tablen
  real*8 :: ngr

  ngr=ng
  if(yetthere6==1) call gtk_widget_destroy(add_window)
  add_window = gtk_window_new (GTK_WINDOW_TOPLEVEL)
  yetthere6=1
  call gtk_window_set_default_size(add_window, 50, 100)
  call gtk_window_set_title(add_window, "Arrange your GRN nodes"//c_null_char)

  tablen = gtk_table_new (5_c_int, 2_c_int, TRUE)
  buttonn = gtk_button_new_with_mnemonic ("Arrange circular (default)"//c_null_char)
  buttonn1 = gtk_button_new_with_mnemonic ("Arrange manually"//c_null_char)
  buttonn2 = gtk_button_new_with_mnemonic ("Return to GRN modification"//c_null_char)
  buttonn3 = gtk_button_new_with_mnemonic ("Save arrangement"//c_null_char)
  buttonn4 = gtk_button_new_with_mnemonic ("Retrieve arrangement"//c_null_char)

  call g_signal_connect (buttonn, "clicked"//c_null_char, c_funloc(my_circular))
  call g_signal_connect (buttonn1, "clicked"//c_null_char, c_funloc(my_arrange_manually))
  call g_signal_connect (buttonn2, "clicked"//c_null_char, c_funloc(my_arrange_return))
  call g_signal_connect (buttonn3, "clicked"//c_null_char, c_funloc(my_arrange_saver))
  call g_signal_connect (buttonn4, "clicked"//c_null_char, c_funloc(my_arrange_retriever))
  call gtk_table_attach_defaults(tablen, buttonn, 0_c_int, 2_c_int, 0_c_int, 2_c_int)
  call gtk_table_attach_defaults(tablen, buttonn1, 0_c_int, 2_c_int, 2_c_int, 4_c_int)
  call gtk_table_attach_defaults(tablen, buttonn2, 0_c_int, 2_c_int, 8_c_int, 10_c_int)
  call gtk_table_attach_defaults(tablen, buttonn3, 0_c_int, 2_c_int, 4_c_int, 6_c_int)
  call gtk_table_attach_defaults(tablen, buttonn4, 0_c_int, 2_c_int, 6_c_int, 8_c_int)
  call gtk_container_add (add_window, tablen)
  call gtk_widget_show_all (add_window)

end subroutine move_nodes_control

subroutine update_run 

! This works like an accept-button for global parameters
  use iso_c_binding, only: c_ptr
  use global_widgets
  use basic_handlers
  use on_display_handlers
  use table_handlers

  implicit none

  do i=1,4
    if(ffu(i)==1) call gtk_toggle_button_set_active (ffut(i), TRUE)
    ffu(i)= gtk_toggle_button_get_active(ffut(i))
  end do
  do i=7,8
    if(ffu(i)==1) call gtk_toggle_button_set_active (ffut(i), TRUE)
    ffu(i)= gtk_toggle_button_get_active(ffut(i))
  end do

  do i=10,13
    if(ffu(i)==1) call gtk_toggle_button_set_active (ffut(i), TRUE)
    ffu(i)= gtk_toggle_button_get_active(ffut(i))
  end do
  do i=16,18
    if(ffu(i)==1) call gtk_toggle_button_set_active (ffut(i), TRUE)
    ffu(i)= gtk_toggle_button_get_active(ffut(i))
  end do
  do i=20,nfu
    if(ffu(i)==1) call gtk_toggle_button_set_active (ffut(i), TRUE)
    ffu(i)= gtk_toggle_button_get_active(ffut(i))
  end do

  ffu(9)=gtk_spin_button_get_value (spinbutt9)

  ffu(19)=gtk_spin_button_get_value (spinbutt8)

  do i=1,nparam
    param(i)=gtk_spin_button_get_value (spinbutt(i))

  end do

end subroutine update_run

function save_it_plot (widget, gdata ) result(ret)  bind(c)

  use iso_c_binding, only: c_ptr
  use global_widgets
  use basic_handlers
  use on_display_handlers
  use table_handlers
  use io
  use genetic
  use general

  implicit none
    
  type(c_ptr), value :: widget, gdata
  integer(c_int)    :: ret
  integer :: yetthere4=0

  type(c_ptr) :: boxc
  type(c_ptr) :: tablep, buttonpi, labelpi

  if(yetthere4==1) call gtk_widget_destroy(my_choice_windowp) !XX!

  my_choice_windowp = gtk_window_new (GTK_WINDOW_TOPLEVEL)
  yetthere4=1

  call g_signal_connect (my_choice_windowp, "delete-event"//c_null_char, c_funloc(delete_event))
  call gtk_window_set_default_size(my_choice_windowp, 500, 50)
  call gtk_window_set_title(my_choice_windowp, "Save recent plot as:"//c_null_char)
  labelpi = gtk_label_new("Define the name of your network plot here"//c_null_char)

  entrypi = gtk_entry_new()
  name_png='give_me_a_fancy_name_please.png'
  call gtk_entry_set_text(entrypi,name_png)
  buttonpi = gtk_button_new_with_mnemonic ("Save it"//c_null_char)

  tablep = gtk_table_new (3_c_int, 2_c_int, TRUE)
  call g_signal_connect (buttonpi, "clicked"//c_null_char, c_funloc(save_plot))
 
  call gtk_table_attach_defaults(tablep, labelpi, 0_c_int, 2_c_int, 0_c_int, 1_c_int)
  call gtk_table_attach_defaults(tablep, entrypi, 0_c_int, 2_c_int, 1_c_int, 2_c_int)
  call gtk_table_attach_defaults(tablep, buttonpi, 0_c_int, 2_c_int, 2_c_int, 3_c_int)

  call gtk_container_add (my_choice_windowp, tablep)

  call gtk_widget_show_all (my_choice_windowp)

end function save_it_plot

function save_it (widget, gdata ) result(ret)  bind(c)

  use iso_c_binding, only: c_ptr
  use global_widgets
  use basic_handlers
  use on_display_handlers
  use table_handlers
  use io
  use genetic
  use general

  implicit none
    
  type(c_ptr), value :: widget, gdata
  integer(c_int)    :: ret

  type(c_ptr) :: boxc
  type(c_ptr) :: tables, buttonsi, labelsi

  if(namefieldopen==1) call gtk_widget_destroy(my_choice_windows) !XX!

  my_choice_windows = gtk_window_new (GTK_WINDOW_TOPLEVEL)
  namefieldopen=1

  call g_signal_connect (my_choice_windows, "delete-event"//c_null_char, c_funloc(delete_event))
  call gtk_window_set_default_size(my_choice_windows, 500, 50)
  call gtk_window_set_title(my_choice_windows, "Save input file as:"//c_null_char)
  labelsi = gtk_label_new("Define the name of your new input file here"//c_null_char)

  entrysi = gtk_entry_new()
  name_dat='give_me_a_name_please.dat'
  call gtk_entry_set_text(entrysi,name_dat)
  buttonsi = gtk_button_new_with_mnemonic ("Save it"//c_null_char)

  tables = gtk_table_new (3_c_int, 2_c_int, TRUE)
  call g_signal_connect (buttonsi, "clicked"//c_null_char, c_funloc(my_save_accept))
 
  call gtk_table_attach_defaults(tables, labelsi, 0_c_int, 2_c_int, 0_c_int, 1_c_int)
  call gtk_table_attach_defaults(tables, entrysi, 0_c_int, 2_c_int, 1_c_int, 2_c_int)
  call gtk_table_attach_defaults(tables, buttonsi, 0_c_int, 2_c_int, 2_c_int, 3_c_int)

  call gtk_container_add (my_choice_windows, tables)

  call gtk_widget_show_all (my_choice_windows)

end function save_it


function my_save_accept (widget, gdata ) result(ret)  bind(c)

  use iso_c_binding, only: c_ptr
  use global_widgets
  use basic_handlers
  use on_display_handlers
  use table_handlers
  use io
  use genetic
  use general

  implicit none
    
  type(c_ptr), value :: widget, gdata
  integer(c_int)    :: ret
  character*12 :: name_def
  type(c_ptr) :: name_dati
  character*1 :: pagea
  integer, dimension(8) :: dt
  integer :: t1

name_dat='GRN.dat'

  name_dati=gtk_entry_get_text(entrysi)
  call c_f_string(name_dati, name_dat)

   call gtk_widget_destroy(my_choice_windows)
   namefieldopen=0

  print*, "WRITE THE NAME OF THE FILE AS WHICH IT SHALL BE SAVED. THEN PRESS ENTER."
  print*, name_dat
  if(len(trim(name_dat)).gt.119) print*, "More than 120 characters are too much. This is going to FAIL."
  name_dat=trim('"'//trim(name_dat)//'"')
  if((len(trim(name_dat)).ge.119).or.(len(trim(name_dat)).le.3))then

    call date_and_time(values=dt)
    t1=dt(6)+100*dt(5)+10000*dt(4)+1000000*dt(3)
    write (name_def,"(I12.1)") t1

    name_def=trim(name_def)
    name_dat=trim(adjustr(name_def)//".dat") ! this is the default name

  end if

  call iniread 

  call update_run

  call put_param_to_matrix(param)

  param(19)=ng; param(10)=nd; param(9)=ntipusadh

  call writesnap ! check if writesnap or iniwritesnap
  name_dat=trim("cp fort.1 "//adjustl(name_dat)) 

  call system(name_dat)


end function my_save_accept


function run_it (widget, gdata ) result(ret)  bind(c)

  use iso_c_binding, only: c_ptr
  use global_widgets
  use basic_handlers
  use on_display_handlers
  use table_handlers
  use genetic
  use general
  use io


  implicit none
    
  type(c_ptr), value :: widget, gdata
  integer(c_int)    :: ret
  character*12 :: name_def
  character*120 :: name_datt
  character*1 :: pagea
  integer, dimension(8) :: dt
  integer :: t1, t2, t3

  name_dat='GRN.dat'

  call date_and_time(values=dt)
  t1=dt(6)+100*dt(5)+10000*dt(4)+1000000*dt(3)
  write (name_def,"(I12.1)") t1
  name_def=trim(name_def)
  name_dat=trim(adjustr(name_def)//".dat") ! this is the default name

nomfinal='out.dat'
  t3=ng
  call writesnap ! do we need this if we only run the model?
  t2=ng
do while(t2.ne.t3)

  call writesnap
  t2=ng
end do

  name_datt=trim(name_dat)

  name_datt=trim(name_datt)
  name_dat=trim("cp fort.1 "//adjustl(name_dat)) 

  call system(name_dat)
  name_datt=trim('./EMaker '//trim("fort.1"))
  call system(name_datt)

end function run_it

subroutine list_select(lista, gdata) bind(c)

  use genetic
  use table_handlers

  type(c_ptr), value :: lista, gdata

  integer(kind=c_int) :: nsel
  integer(kind=c_int), dimension(:,:), allocatable :: selections
  integer(kind=c_int), dimension(:), allocatable :: dep
  integer(kind=c_int) :: n, n3
  integer(kind=c_int64_t) :: n4
  real(kind=c_float) :: nlog
  character(len=10) :: nodd
  integer :: il,jl,kl
  integer(kind=c_float) :: name
  real(kind=c_float) :: nn
  integer :: nqpars(23)
      
  nqpars=(/5,6,7,8,9,10,11,12,13,14,15,16,20,21,22,23,24,25,26,27,28,34,35/)

  nsel = hl_gtk_tree_get_selections(C_NULL_PTR, selections, selection=lista, &
         & depths=dep)
  if (nsel == 0) then
    print *, "No selection"
    return
  end if

if(open_tab==0)then

  if (nsel == 1) then

    if(allocated(neww)) deallocate(neww)
    allocate(neww(ng,ng+3)) 
    if(allocated(oldw)) deallocate(oldw)
    allocate(oldw(ng,ng+3)) 
    do il=1,ng
      do jl=1,ng
        neww(il,jl)=gen(il)%w(jl)   
      end do
      neww(il,ng+1)=gen(il)%diffu   
      neww(il,ng+2)=gen(il)%mu 
      neww(il,ng+3)=gen(il)%kindof 
    end do
    oldw=neww
    do il=1,ng+3
      do kl=0,ng-1
      selections(:dep(1),1)=kl
      call hl_gtk_tree_get_cell(ihlist, selections(:dep(1),1), 0, ivalue=name)
      call hl_gtk_tree_get_cell(ihlist, selections(:dep(1),1), il, fvalue=nn)
      neww(name,il)=nn
      end do
    end do

  else

  end if

  do il=1,ng
    do jl=1,ng
      if(whichtab==1) gen(il)%w(jl) = neww(il,jl)  
      if(whichtab==2) gen(il)%w(jl) = oldw(il,jl)  
    end do
      if(whichtab==1) gen(il)%diffu=neww(il,ng+1)   
      if(whichtab==1) gen(il)%mu=neww(il,ng+2) 
      if(whichtab==1) gen(il)%kindof=int(neww(il,ng+3))
      if(whichtab==2) gen(il)%diffu=oldw(il,ng+1)   
      if(whichtab==2) gen(il)%mu=oldw(il,ng+2) 
      if(whichtab==2) gen(il)%kindof=int(oldw(il,ng+3))
  end do

elseif(open_tab==2)then

  if (nsel == 1) then

    if(allocated(neww)) deallocate(neww)
    allocate(neww(ng,23)) 
    if(allocated(oldwa)) deallocate(oldwa)
    allocate(oldwa(ng,23)) 

    do il=1,ng
      do jl=1,23
        neww(il,jl)=gen(il)%wa(nqpars(jl))   
      end do
    end do
    oldwa=neww
    do il=1,ng
      do kl=0,ng-1
      selections(:dep(1),1)=kl
      call hl_gtk_tree_get_cell(ihlist, selections(:dep(1),1), 0, ivalue=name)
      call hl_gtk_tree_get_cell(ihlist, selections(:dep(1),1), il, fvalue=nn)
      neww(name,il)=nn
      end do
    end do

  else

  end if

  do il=1,ng
    do jl=1,23
      if(whichtab==1) gen(il)%wa(nqpars(jl)) = neww(il,jl)  
      if(whichtab==2) gen(il)%wa(nqpars(jl)) = oldwa(il,jl)  
    end do
  end do

elseif(open_tab==3)then

  if (nsel == 1) then

    if(allocated(neww)) deallocate(neww)
    allocate(neww(ng,16)) 
    if(allocated(oldwa)) deallocate(oldwa)
    allocate(oldwa(ng,16)) 

    do il=1,ng
      do jl=1,16
        neww(il,jl)=gen(il)%wa(nparam_per_node+jl)   
      end do
    end do
    oldwa=neww
    do il=1,ng
      do kl=0,ng-1
      selections(:dep(1),1)=kl
      call hl_gtk_tree_get_cell(ihlist, selections(:dep(1),1), 0, ivalue=name)
      call hl_gtk_tree_get_cell(ihlist, selections(:dep(1),1), il, fvalue=nn)
      neww(name,il)=nn
      end do
    end do

  else

  end if

  do il=1,ng
    do jl=1,16
      if(whichtab==1) gen(il)%wa(nparam_per_node+jl) = neww(il,jl)  
      if(whichtab==2) gen(il)%wa(nparam_per_node+jl) = oldwa(il,jl) 
    end do
  end do

end if

  deallocate(selections)
end subroutine list_select

function info_connect (widget, gdata ) result(ret)  bind(c)

  use iso_c_binding, only: c_ptr
  use global_widgets
  use basic_handlers
  use io
  use general
  use genetic

  implicit none
        
  integer(c_int)    :: ret
  type(c_ptr), value :: widget, gdata
  type(c_ptr) :: my_cairo_context
  integer:: who,get,gett
  character*140 uui

uui='dn.dat'

print*, "How can I help you?"
print*, "Press 1 for Information about M params"
print*, "Press 2 for Information about L params"
print*, "Press 3 for Information about the current network view"
print*, "Press 4 for Information about the buttons"
print*, "Press 5 for FAQs"
print*, "Press any other NUMBER to leave this information menu"

read*, get

if(get==1)then
  print*, "TELL WHAT THIS MATRIX IS ABOUT"
  print*, "Some parameters apply to the entire model simulations."
  print*, "These are called MODEL PARAMETERS and are explained in the model description."
  print*, "Here is where you can change most of them. Note that some of them stay invariable."
  print*, "Note also that some of them are continuous, while others are discrete."
elseif(get==2)then
  print*, "TELL WHAT THIS MATRIX IS ABOUT"
  print*, "These are the LOGICAL MODEL PARAMETERS."
  print*, "They represent non-standard options or conditions under which the model should be run."
  print*, "Usually, you can toggle between ON/OFF, in two cases, there are more than two options."
  print*, "Read the model description to learn more about them."
elseif(get==3)then
  who=gtk_notebook_get_current_page(notebook)
    ! print here what may be important to know
  if(who==0)then
    print*, "THIS IS THE NETWORK OF TRANSCRIPTIONAL INTERACTIONS."
    print*, "Regulatory molecules, such as genes and proteins, are symbolized by spherical nodes."
    print*, "By default, all regulatory molecules in the network used are arranged in a circle."
    print*, "Different colors represent some of their particular properties, like localization within cells."
    print*, "A legend is always shown underneath the network."
    print*, "They can be changed or adjusted by pressing the 'Modify gene properties' button."
    print*, " "
    print*, "In this view, positive transcriptional interactions are shown as green arrows,"
    print*, "negative ones as red arrows. Autoregulation arrows are of circular shape."
    print*, "In order to draw a new connection, you have two options:"
    print*, "Either you activate the 'Click on Network' mode in the lower part on the button panel"
    print*, "which allows you to automatically draw connections by clicking on two genes sequentially."
    print*, "With the 'Create + connection' and 'Create - connetion', "
    print*, "you can switch between positive and negative interaction arrows."
    print*, "The default magnitude is always 0.1. With 'Delete connection', you delete the chosen arrow."
    print*, "You can modify interaction strengthes via the 'Modify connection' button."
    print*, "Or you activate the 'Click on Buttons' option, which doesn't draw the arrows directly."
    print*, "Instead, you need to choose two genes, then clock on one of the uppermost four buttons,"
    print*, "(whose functions I just explained)"
    print*, "to open a dialog window that asks you to adjust the interaction strength and"
    print*, "finally accept your decision pressing the ACCEPT button."
    print*, " "
  elseif(who==1)then
    print*, "THIS IS THE NETWORK OF CATALYTIC INTERACTIONS."
    print*, "Regulatory molecules, such as genes and proteins, are symbolized by spherical nodes."
    print*, "By default, all regulatory molecules in the network used are arranged in a circle."
    print*, "Different colors represent some of their particular properties, like localization within cells."
    print*, "A legend is always shown underneath the network."
    print*, "They can be changed or adjusted by pressing the 'Modify gene properties' button."
    print*, " "
    print*, "Some regulatory molecules can be transformed into others in a catalytic reaction."
    print*, "Depending on what kind of networks you are modelling, this can be genes translated into proteins,"
    print*, "proteins into active proteins, phosphorylations, membrane-bound into diffusive protein"
    print*, "configurations, bound and unbound proteins etc."
    print*, "For all these reactions, SUBSTRAT, PRODUCT and CATALYST need to be defined."
    print*, "This network view visualizes and allows modification of such reactions."
    print*, "In this view, brown arrows point from SUBSTRAT to the PRODUCT of a catalytic reaction,"
    print*, "that become darker if you define a CATALYST."
    print*, "Positive catalysations are represented by green arrows, pointing from the CATALYST to the"
    print*, "reaction arrow near the PRODUCT, negative ones as red arrows."
    print*, "Notice that the CATALYST may be the SUBSTRAT or PRODUCT of the reaction."
    print*, " "
    print*, "In order to define a new reaction, you have two options:"
    print*, "Either you activate the 'Click on Network' mode in the lower part on the button panel"
    print*, "which allows you to automatically draw such interactions by clicking on three genes sequentially."
    print*, "The first one, underlain by a yellow halo, represents the SUBSTRATE."
    print*, "The second one, underlain by a red halo, represents the PRODUCT."
    print*, "The first one, underlain by a blue halo, represents the CATALYST."
    print*, "With the 'Create + connection' and 'Create - connetion', "
    print*, "you can switch between positive and negative interaction arrows."
    print*, "The default magnitude is always 0.1. With 'Delete connection', you delete the chosen arrow."
    print*, "You can modify interaction strengthes via the 'Modify connection' button."
    print*, "Or you activate the 'Click on Buttons' option, which doesn't draw the arrows directly."
    print*, "Instead, you need to choose two genes, then clock on one of the uppermost four buttons,"
    print*, "(whose functions I just explained)"
    print*, "to open a dialog window that asks you to adjust the interaction strength and"
    print*, "finally accept your decision pressing the ACCEPT button."
    print*, "Notice that this network gets messed up easily and you may want to rearrange your nodes"
    print*, "to increase readability. You can do that in the 'ARRANGE MODE'."
    print*, "Activate it by clicking on the toggle in the lower right corner."
    print*, " "
  elseif(who==2)then
    print*, "THIS IS THE NETWORK OF NODE PROPERTY REGULATIONS."
    print*, "Regulatory molecules, such as genes and proteins, are symbolized by spherical nodes."
    print*, "By default, all regulatory molecules in the network used are arranged in a circle."
    print*, "Different colors represent some of their particular properties, like localization within cells."
    print*, "A legend is always shown underneath the network."
    print*, "They can be changed or adjusted by pressing the 'Modify gene properties' button."
    print*, "NODE PROPERTIES are mostly mechanical properties that a node can have."
    print*, "E.g. stiffness, adhesivity, size, motility, differentiation state etc."
    print*, "They are regulated by gene expression, as being mediated by (mostly structure and surface) proteins"
    print*, "and define the mechanical properties of a cell."
    print*, "Notice that they may be different between nodes of the same cell."
    print*, "Here, they are written in a short form in the upper half of the display."
    print*, "Check the model description for further information."
    print*, " "
    print*, "In this view, positive regulations are shown as green lines,"
    print*, "negative ones as red lines."
    print*, "In order to draw a new connection, you have two options:"
    print*, "Either you activate the 'Click on network' mode in the lower part on the button panel"
    print*, "which allows you to automatically draw connections by clicking on a gene and a node property."
    print*, "With the 'Create + connection' and 'Create - connetion', "
    print*, "you can switch between positive and negative interactions."
    print*, "The default magnitude is always 0.1. With 'Delete connection', you delete the chosen arrow."
    print*, "You can modify interaction strengthes via the 'Modify connection' button."
    print*, "Or you activate the 'Click on buttons' option, which doesn't draw the arrows directly."
    print*, "Instead, you need to choose two elements, then clock on one of the uppermost four buttons,"
    print*, "(whose functions I just explained)"
    print*, "to open a dialog window that asks you to adjust the interaction strength and"
    print*, "finally accept your decision pressing the ACCEPT button."
    print*, " "
  elseif(who==3)then
    print*, "THIS IS THE NETWORK OF CELL BEHAVIOUR REGULATIONS."
    print*, "Regulatory molecules, such as genes and proteins, are symbolized by spherical nodes."
    print*, "By default, all regulatory molecules in the network used are arranged in a circle."
    print*, "Different colors represent some of their particular properties, like localization within cells."
    print*, "A legend is always shown underneath the network."
    print*, "They can be changed or adjusted by pressing the 'Modify gene properties' button."
    print*, "CELL BEHAVIOURS simply determine what cells are capable of."
    print*, "E.g. growth, division, apoptosis, secretion, polarization etc."
    print*, "These behaviours are regulated by gene expression."
    print*, "Notice that they are the same in all nodes of a given cell."
    print*, "Here, they are written in a short form in the upper half of the display."
    print*, "Check the model description for further information."
    print*, " "
    print*, "In this view, positive regulations are shown as green lines,"
    print*, "negative ones as red lines."
    print*, "In order to draw a new connection, you have two options:"
    print*, "Either you activate the 'Click on Network' mode in the lower part on the button panel"
    print*, "which allows you to automatically draw connections by clicking on a gene and a cell behaviour."
    print*, "With the 'Create + connection' and 'Create - connetion', "
    print*, "you can switch between positive and negative interactions."
    print*, "The default magnitude is always 0.1. With 'Delete connection', you delete the chosen arrow."
    print*, "You can modify interaction strengthes via the 'Modify connection' button."
    print*, "Or you activate the 'Click on buttons' option, which doesn't draw the arrows directly."
    print*, "Instead, you need to choose two elements, then clock on one of the uppermost four buttons,"
    print*, "(whose functions I just explained)"
    print*, "to open a dialog window that asks you to adjust the interaction strength and"
    print*, "finally accept your decision pressing the ACCEPT button."
    print*, " "
  elseif(who==4)then
    print*, "THIS IS THE NETWORK OF SPECIFIC ADHESIONS."
    print*, "Regulatory molecules, such as genes and proteins, are symbolized by spherical nodes."
    print*, "By default, all regulatory molecules in the network used are arranged in a circle."
    print*, "Different colors represent some of their particular properties, like localization within cells."
    print*, "A legend is always shown underneath the network."
    print*, "They can be changed or adjusted by pressing the 'Modify gene properties' button."
    print*, "Genes can encode proteins which change the adhesion of a cell."
    print*, "In other words, the adhesion forces between two contacting cells"
    print*, "depend on the particular set of surface proteins expressed on them."
    print*, "We model these interactions by defining adhesion types that are mediated by a gene/protein"
    print*, "give them a number (called Bij element) and represent them here by a golden sphere."
    print*, "Adhesive interactions are defined between these adhesion types and can be"
    print*, "either positive (Adhesion) or negative (Repulsion)."
    print*, " "
    print*, "In this view, positive heterotypic interactions between Bij elements are shown as green lines,"
    print*, "negative ones as red lines. Homotypic interactions are drawn as circles."
    print*, "Blue lines display which genes they are mediated by."
    print*, "In order to draw a new line, you have two options:"
    print*, "Either you activate the 'Click on Network' mode in the lower part on the button panel"
    print*, "which allows you to automatically draw connections by clicking on two Bij elements"
    print*, "or a Bij element and a gene."
    print*, "With the 'Create + connection' and 'Create - connetion', "
    print*, "you can switch between positive and negative interactions."
    print*, "The default magnitude is always 0.1. With 'Delete connection', you delete the chosen arrow."
    print*, "You can modify interaction strengthes via the 'Modify connection' button."
    print*, "Or you activate the 'Click on buttons' option, which doesn't draw the arrows directly."
    print*, "Instead, you need to choose two elements, then clock on one of the uppermost four buttons,"
    print*, "(whose functions I just explained)"
    print*, "to open a dialog window that asks you to adjust the interaction strength and"
    print*, "finally accept your decision pressing the ACCEPT button."
    print*, " "
  end if
elseif(get==4)then
  print*, "COPY THE MANUAL TEXTS HERE"
elseif(get==5)then
  print*, " ******* F A Qs ******* "
  print*, " Choose one of the following questions by entering its number."
  print*, " Press 0 in order to quit. "
  print*, " "
  print*, " 1 - How can I change the input file? "
  print*, " 2 - How can I change the initial morphology?"
  print*, " 3 - How can I try out the changes I have made?"
  print*, " 4 - How can I save my network?"
  print*, " 5 - The circular arrangement is not practical for me. How can I change that?"
  print*, " 6 - I dont understand the names used here. Where can I read more?"
  print*, " 7 - I am confused. What is the difference between the matrices?"
  print*, " 8 - What are the different types of molecules?"
  print*, " 9 - I found a problem that annoys me or have a suggestion."
  read*, gett
  if(gett==1)then
    print*, " 1 - How can I change the input file? "
    print*, ""
    print*, "There are two ways how you can do that."
    print*, ""
    print*, "1.When you run the program by typing ./NetworkMaker in the command line,"
    print*, "  write the name of the input file after: ./NetworkMaker FILE."
    print*, "  Be aware that the PATH to your FILE must be included."
    print*, "2.After typing ./NetworkMaker without defining an input file,"
    print*, "  a small window opens and asks you if you wish to choose your own file"
    print*, "  or want to START the network editor with a default file."
    print*, "  Press the 'Choose a new input file' button to access the file system"
    print*, "  of your computer and make a choice by double-clicking on an input file."
    print*, ""
  elseif(gett==2)then
    print*, " 2 - How can I change the initial morphology?"
    print*, ""
    print*, " The initial morphology cannot be changed by means of NetworkMaker."
    print*, " However, this is done by the Embryo Editor of EmbryoMaker."
    print*, " Please, read the respective manual in order to learn more."
    print*, ""
  elseif(gett==3)then
    print*, " 3 - How can I try out the changes I have made?"
    print*, ""
    print*, " NetworkMaker includes a funtion to run EmbryoMaker with the current network."
    print*, " Simply press the RUN SIMULATION button at the right side."
    print*, " If you want to keep the input conditions, you can write a new input file"
    print*, " by means of the SAVE INPUT FILE button."
    print*, ""
  elseif(gett==4)then
    print*, " 4 - How can I save my network?"
    print*, ""
    print*, " In order to save the network, press the SAVE INPUT FILE button on the right side."
    print*, " You will need to enter a name for the file right here."
    print*, ""
    print*, " You can also save the current view by pressing SAVE PLOT(PNG)."
    print*, " Again, you will be prompted to enter a name here."
    print*, ""
  elseif(gett==5)then
    print*, " 5 - The circular arrangement is not practical for me. How can I change that?"
    print*, ""
    print*, " NetworkMaker comes with the convenient option of manually re-arranging the view."
    print*, " You can toggle the ARRANGE MODE at the bottom of the list of buttons."
    print*, " Now, you will be able to choose genes by clicking on them,"
    print*, " and displace them by clicking on their new position."
    print*, " The connections between the genes will remain unchanged."
    print*, " Since the arrangements of genes differ between the matrices displayed,"
    print*, " you will need to rearrange them individually."
    print*, ""
    print*, " In addition, a sub-menu will be prompted."
    print*, " Press ARRANGE CIRCULAR in order to restall the original circular arrangement."
    print*, " Press SAVE ARRANGEMENT to save the current positions in a matrix."
    print*, " This allows you to return to the circular arrangement"
    print*, " without losing your arrangement."
    print*, " In order to retrieve it later after returning to the drawing mode,"
    print*, " just press the RETRIEVE ARRANGEMENT button."
    print*, " "
  elseif(gett==6)then
    print*, " 6 - I dont understand the names used here. Where can I read more?"
    print*, ""
    print*, " GNOMO has been created to encompass various aspects of development."
    print*, " This is why the names and abbreviations used here may appear messy."
    print*, " Please, read the general manual FILE to learn more."
    print*, ""
  elseif(gett==7)then
    print*, " 7 - I am confused. What is the difference between the matrices?"
    print*, ""
    print*, " The matrices represent different classes of interactions."
    print*, ""
    print*, " T-Matrix: Transcriptional interactions."
    print*, " Arrows between nodes represent transcriptional interactions between genes."
    print*, " This means that the rate of their production is modified by other genes."
    print*, " "
    print*, " R-Matrix: Catalytic interactions."
    print*, " Molecules such as proteins can act as catalysts in various reactions."
    print*, " This applies mostly to reactions in which a molecule A is converted "
    print*, " into another one B. Both A and B are represented independently in the network."
    print*, " Unlike in transcriptional activations or repressions,"
    print*, " catalytic reactions are stochiometrically conservative, which means"
    print*, " that the amount of B cannot exceed the amount of A."
    print*, " Use the left mouse buttons to choose A and B and the right mouse button"
    print*, " for the catalyst."
    print*, " Note that you may want to include the back-reactions as well."
    print*, " Examples are Phophorylations, protein cleavings and reactor-ligand bindings."
    print*, ""
    print*, " E-Matrix: Interactions with Biomechanics."
    print*, " The regulatory molecules are arrayed on the lower half of the display field."
    print*, " On the upper half, you can see the biomechanical parameters used by GNOMO."
    print*, " Read the general user's manual to learn more about each of them."
    print*, " e.g. interaction radii, unspecific adhesivity, elasticity, motility,..."
    print*, " Those can be regulated by molecules."
    print*, " In order to draw or modify such a regulation, you can use this interface."
    print*, ""
    print*, " C-Matrix: Regulation of Cell Behaviours."
    print*, " Genes and proteins can regulate cell behaviours such as"
    print*, " Cell cycle rate, Apoptosis, ECM secretion, Epithelium-Mesenchyme transitions etc."
    print*, " Like in the E-Matrix, you can modify regulatory relationships between genes and"
    print*, " cell behaviours using this interface."
    print*, ""
    print*, " B-Matrix: Specific adhesions."
    print*, " Surface protein change the adhesive properties of cells in a molecule-specific"
    print*, " manner by interacting with other proteins at the surface."
    print*, " By means of this interface, gene products can be linked with different"
    print*, " adhesion types, as represented by golden spheres."
    print*, " this allows to model, e.g., differential adhesion and repulsion."
    print*, ""
    print*, " Note that some properties such as diffusivity and degradation rate can be modified"
    print*, " in a separate window that opens after clicking on the respective molecule."
    print*, ""
  elseif(gett==8)then
    print*, " 8 - What are the different types of molecules?"
    print*, ""
    print*, " Molecule species differ in their spatial properties (localizations), which is relevant"
    print*, " if it comes to integrating into a spatial multi-scale model."
    print*, " Here, we use different colours to represent these differences. "
    print*, " Grey: Nuclear molecules (Genes) that do not change."
    print*, " Dark green: Nuclear molecules (Genes) that can be translated into molecules"
    print*, " with a different localization (Proteins)"
    print*, " Dark blue: intracellularly diffusing molecules"
    print*, " Red: extracellularly diffusive molecules."
    print*, " light blue and pink: apically and basally localizing molecules (in epithelium)"
    print*, " yellow: active surface receptors; bind extracellulary diffusive molecules."
    print*, " bright green: surface-bound receptors for cell-cell signalling."
    print*, ""
  elseif(gett==9)then
    print*, " Many problems may be related to the use of different versions of the operating system."
    print*, " You can write in case you find some bugs or want to voice suggestions"
    print*, " to: roland.zimm@helsinki.fi "
    print*, " If it is about a bug, please mention which operating system you are working with"
    print*, " and copy the error message, if any."
    print*, " Thank you :)"
    print*, ""
  endif
endif

end function info_connect

function save_plot (widget, gdata ) result(ret)  bind(c)

  use iso_c_binding, only: c_ptr
  use global_widgets
  use basic_handlers
  use on_display_handlers
  use table_handlers

  implicit none
    
  integer(c_int)    :: ret
  type(c_ptr), value :: widget, gdata
  type(c_ptr) :: my_cairo_context, name_pngi
  integer(c_int) :: cstatus, message_id, t1
  character*12 :: name_def
  character*1 :: pagea
  integer, dimension(8) :: dt

  name_png='GRN.png'

  name_pngi=gtk_entry_get_text(entrypi)
  call c_f_string(name_pngi, name_png)

   call gtk_widget_destroy(my_choice_windowp)

  print*, "WRITE THE NAME OF THE FILE AS WHICH IT SHALL BE SAVED"
  print*, name_png
  if(len(trim(name_png)).ge.119) print*, "More than 120 characters are too much. This is going to FAIL."
  name_png=adjustr(trim(name_png))
  if((len(trim(name_png)).ge.119).or.(len(trim(name_png)).eq.0))then
    call date_and_time(values=dt)
    t1=dt(6)+100*dt(5)+10000*dt(4)+1000000*dt(3)
    write (name_def,"(I12.1)") t1
    t1=gtk_notebook_get_current_page(notebook)
    write (pagea,"(I1.1)") t1
    name_png=trim("GRN_"//trim(pagea)//"_"//trim(name_def)//".png") ! this is the default name
  end if

  call gtk_widget_realize(my_drawing_area)
  my_cairo_context = gdk_cairo_create (gtk_widget_get_window(my_drawing_area))
  call gdk_cairo_set_source_pixbuf(my_cairo_context, my_pixbuf, 0d0, 0d0) 

  if (.not. computing) then
    !cstatus = cairo_surface_write_to_png(cairo_get_target(my_cairo_context), trim(name_png)//c_null_char)
    if (cstatus == CAIRO_STATUS_SUCCESS) then
      string = "Successfully saved: GRN.png"//c_null_char
    else if (cstatus == CAIRO_STATUS_NO_MEMORY) then
      string = "Failed: memory allocation"//c_null_char
    else if (cstatus == CAIRO_STATUS_SURFACE_TYPE_MISMATCH) then
      string = "Failed: no pixel content"//c_null_char
    else if (cstatus == CAIRO_STATUS_WRITE_ERROR) then
      string = "Failed: I/O error"//c_null_char
    else
      string = "Failed"
    end if

  end if
    ! FIXME: how to save only the drawing_area
  call cairo_destroy(my_cairo_context)

    ret = FALSE

end function save_plot

function plot_extra (widget, gdata ) result(ret)  bind(c)

  use iso_c_binding, only: c_ptr
  use global_widgets
  use basic_handlers
  use on_display_handlers
  use table_handlers

  implicit none
    
  type(c_ptr), value :: widget, gdata
  integer(c_int)    :: ret
  integer :: choic

  choic=gtk_notebook_get_current_page(notebook)

  if(choic==0)then
    choice = 1
  elseif(choic==1)then
    choice = 3
  elseif(choic==2)then
    choice = 5
  elseif(choic==3)then
    choice = 4
  elseif(choic==4)then
    choice = 6
  else
    choice=choic
  end if

  call GRN_plot(10+choice)

end function plot_extra

function undo_redo_connect (widget, gdata ) result(ret)  bind(c)

  use iso_c_binding, only: c_ptr
  use global_widgets
  use basic_handlers
  use on_display_handlers
  use table_handlers

  implicit none
    
  type(c_ptr), value :: widget, gdata
  integer(c_int)    :: ret
  integer :: choic

  if(counts.gt.0) call undo_redo(8,0,0,0,0,0d0,0d0)

end function undo_redo_connect

function draw_fast_connect (widget, gdata ) result(ret)  bind(c)

  use iso_c_binding, only: c_ptr
  use global_widgets
  use basic_handlers
  use on_display_handlers
  use table_handlers

  implicit none

  type(c_ptr), value :: widget, gdata
  integer(c_int)    :: ret
    
  faster=faster*-1

end function draw_fast_connect

function choice_box_build_connect(widget, event, gdata) result(ret)  bind(c)

  use gtk, only: gtk_notebook_prev_page, gtk_notebook_next_page, gtk_notebook_get_n_pages, gtk_notebook_get_current_page
  use iso_c_binding, only: c_ptr
  use genetic
  use some_widgets ! already in global_widgets
  use global_widgets
  use on_display_handlers, only: choice_box_build, choice_box_param_build, choice_box_build_ww

  type(c_ptr), value :: widget, gdata, event
  integer(c_int)    :: ret
  integer :: aaa1,bbb1,ccc1,ddd1

  modify=1

  quot=0.0d0
  bbb1=bbb0; aaa1=aaa0; ddd1=ddd0
  ccc1 = gtk_notebook_get_current_page(notebook)
if(aaa1*bbb1.ne.0)then
  if(ccc1==2)then
    if(gen(aaa1)%wa(bbb1).ne.0)then
      call choice_box_param_build(aaa1,bbb1,ccc1,ddd1)
    endif
  elseif(ccc1==3)then
    if(gen(aaa1)%wa(nparam_per_node+bbb1).ne.0)then
      call choice_box_param_build(aaa1,bbb1,ccc1,ddd1)
    endif
  elseif(ccc1==1)then
   do i=1,gen(w3)%nww    

      if((gen(w3)%ww(i,1).eq.w1).and.(gen(w3)%ww(i,2).eq.w2).and.(gen(w3)%ww(i,3).ne.0))then

        call choice_box_build_ww(w1,w2,w3,w4)
      endif
    end do
  elseif(ccc1==0)then
    if(gen(bbb1)%w(aaa1).ne.0)then
      call choice_box_build(aaa1,bbb1,ccc1,ddd1)
    endif
  elseif(ccc1==4)then
    call choice_box_build(aaa1,bbb1,ccc1,ddd1)
  endif
endif

end function choice_box_build_connect

function add_plus_connect(widget, event, gdata) result(ret)  bind(c)

  use gtk, only: gtk_notebook_prev_page, gtk_notebook_next_page, gtk_notebook_get_n_pages, gtk_notebook_get_current_page
  use iso_c_binding, only: c_ptr
  use genetic
  use on_display_handlers, only: choice_box_build, choice_box_param_build, choice_box_build_ww, my_accept_s,&
  &my_ww_accept_s, my_wa_accept_s, my_wa2_accept_s, my_adh_add_s
  use global_widgets
  use basic_handlers

  type(c_ptr), value :: widget, gdata, event
  integer(c_int)    :: ret
  integer :: aaa1,bbb1,ccc1,ddd1,abcd

  bbb1=bbb0; aaa1=aaa0; ddd1=ddd0
  abcd=0
  quot=0.1d0
  modify=0
  ccc1 = gtk_notebook_get_current_page(notebook)

if(showup==-1)then
if(aaa1*bbb1*ddd1.ne.0)then

  if(ccc1==2)then
    if(gen(aaa1)%wa(bbb1).eq.0)then
      call choice_box_param_build(aaa1,bbb1,ccc1,ddd1)
    endif
  elseif(ccc1==3)then
    if(gen(aaa1)%wa(nparam_per_node+bbb1).eq.0)then
      call choice_box_param_build(aaa1,bbb1,ccc1,ddd1)
    endif
  elseif(ccc1==1)then
    call choice_box_build(w1,w2,1,w4)
    call choice_box_build_ww(w1,w2,w3,w4)
  elseif(ccc1==0)then

    if(gen(bbb1)%w(aaa1).eq.0)then
      call choice_box_build(aaa1,bbb1,ccc1,ddd1)
    endif
  elseif(ccc1==4)then
    call choice_box_build(aaa1,bbb1,ccc1,ddd1)
  endif

endif
endif
if(choiceboxthere==1)then
if(faster==1)then
if(cual==0) call my_accept_s
if(cual==1) call my_ww_accept_s
if(cual==2) call my_wa_accept_s
if(cual==3) call my_wa2_accept_s
endif
endif

end function add_plus_connect

function add_minus_connect(widget, event, gdata) result(ret)  bind(c)

  use gtk, only: gtk_notebook_prev_page, gtk_notebook_next_page, gtk_notebook_get_n_pages, gtk_notebook_get_current_page
  use iso_c_binding, only: c_ptr
  use genetic
  use on_display_handlers, only: choice_box_build, choice_box_param_build, choice_box_build_ww, my_accept_s,&
&my_ww_accept_s, my_wa_accept_s, my_wa2_accept_s, my_adh_add_s
  use global_widgets
  use basic_handlers

  type(c_ptr), value :: widget, gdata, event
  integer(c_int)    :: ret
  integer :: aaa1,bbb1,ccc1,ddd1,abcd

  bbb1=bbb0; aaa1=aaa0; ddd1=ddd0
  abcd=0
  quot=-0.1d0

  ccc1 = gtk_notebook_get_current_page(notebook)

if(showup==-1)then
if(aaa1*bbb1.ne.0)then
  if(ccc1==2)then
    if(gen(aaa1)%wa(bbb1).eq.0)then
      call choice_box_param_build(aaa1,bbb1,ccc1,ddd1)
    endif
  elseif(ccc1==3)then
    if(gen(aaa1)%wa(nparam_per_node+bbb1).eq.0)then
      call choice_box_param_build(aaa1,bbb1,ccc1,ddd1)
    endif
  elseif(ccc1==1)then

    call choice_box_build(w1,w2,1,w4)
    call choice_box_build_ww(w1,w2,w3,w4)
  elseif(ccc1==0)then
    if(gen(bbb1)%w(aaa1).eq.0)then
      call choice_box_build(aaa1,bbb1,ccc1,ddd1)
    endif
  elseif(ccc1==4)then
    call choice_box_build(aaa1,bbb1,ccc1,ddd1)
  endif
end if
end if
if(choiceboxthere==1)then
if(faster==1)then
!call gtk_widget_destroy(my_choice_window)
!choiceboxthere=0
if(cual==0) call my_accept_s
if(cual==1) call my_ww_accept_s
if(cual==2) call my_wa_accept_s
if(cual==3) call my_wa2_accept_s
end if
endif

end function add_minus_connect

function delete_connect(widget, event, gdata) result(ret)  bind(c)

  use gtk, only: gtk_notebook_prev_page, gtk_notebook_next_page, gtk_notebook_get_n_pages, gtk_notebook_get_current_page
  use iso_c_binding, only: c_ptr
  use genetic
  use some_widgets ! already in global_widgets
  use on_display_handlers, only: choice_box_build, choice_box_param_build, choice_box_build_ww, my_v_delete_s,& 
& my_accept_s, my_ww_accept_s, my_wa_accept_s, my_wa2_accept_s, my_adh_add_s
  use global_widgets

  type(c_ptr), value :: widget, gdata, event
  integer(c_int)    :: ret
  integer :: aaa1,bbb1,ccc1,ddd1,abcd

  bbb1=bbb0; aaa1=aaa0; ddd1=ddd0
  abcd=0
  deleteif=1
  ccc1 = gtk_notebook_get_current_page(notebook)
 
if(aaa1*bbb1*ddd1.ne.0)then
  if(ccc1==2)then
    call choice_box_param_build(aaa1,bbb1,ccc1,ddd1)
  elseif(ccc1==3)then
    call choice_box_param_build(aaa1,bbb1,ccc1,ddd1)
  elseif(ccc1==1)then
    call my_v_delete_s
  elseif(ccc1==0)then
    call choice_box_build(aaa1,bbb1,ccc1,ddd1)
  elseif(ccc1==4)then
    call choice_box_build(aaa1,bbb1,ccc1,ddd1)
  endif
  deleteif=0
endif

if(choiceboxthere==1)then
if(cual==0) call my_accept_s
if(cual==1) call my_ww_accept_s
if(cual==2) call my_wa_accept_s
if(cual==3) call my_wa2_accept_s
if(cual==1) call my_v_delete_s
endif
if(cual==1) call my_ww_accept_s
if(cual==1) call my_v_delete_s

end function delete_connect

function gene_properties_connect(widget, event, gdata) result(ret)  bind(c)

  use gtk, only: gtk_notebook_prev_page, gtk_notebook_next_page, gtk_notebook_get_n_pages, gtk_notebook_get_current_page
  use iso_c_binding, only: c_ptr
  use genetic
  use some_widgets ! already in global_widgets
  use on_display_handlers, only: choice_box_build, choice_prop_box_build
  use global_widgets
  use basic_handlers

  type(c_ptr), value :: widget, gdata, event
  integer(c_int)    :: ret
  integer :: cualq

  cualq=wcl
 
  call choice_prop_box_build(cualq)

end function gene_properties_connect

end module event_handlers

module rb_handlers !>>Miquel8-10-14

  use gtk_hl
  use gtk_hl_button
  use gtk, only: gtk_button_new, gtk_container_add, gtk_main, gtk_main_quit, gtk_&
       &widget_destroy, gtk_radio_button_new, gtk_toggle_button_get_active, gtk_widget&
       &_show, gtk_widget_show_all, gtk_window_new, gtk_init, gtk_disable_setlocale
  use basic_handlers
  use event_handlers
  use global_widgets

  implicit none

  type(c_ptr) :: box, window, qbut, group
  type(c_ptr),dimension(6) :: rbut

contains

  subroutine rb_toggle(widget, gdata) bind(c)
    type(c_ptr), value :: widget, gdata

    integer(kind=c_int), pointer :: fdata
    integer(kind=c_int) :: sdata

    if (gtk_toggle_button_get_active(widget) == FALSE) then
       return
    end if

    if (c_associated(gdata)) then
       call c_f_pointer(gdata, fdata)
    end if
    sdata = hl_gtk_radio_group_get_select(group)

    if(sdata==0)then
      faster=-1 ; print*,"FAST MODE ON"

      if(rb_arrange_tog==1)then
        togglen=1
        rb_arrange_tog=0
      end if
    end if
    if(sdata==1)then
      faster=1 ; print*,"SLOW MODE ON"

      if(rb_arrange_tog==1)then
        togglen=1
        rb_arrange_tog=0
      end if
    end if

    !if(sdata==2)then; rb_arrange_tog=1; togglen=-1

    !if(yetthere1==1)then
    !  call gtk_widget_destroy(my_choice_window); yetthere1=0
    !endif

    !if(yetthere==1)then
    !  call gtk_widget_destroy(my_choice_window2); yetthere=0
    !endif

    !call move_nodes_control ; print*,"ARRANGE MODE ON";end if

  end subroutine rb_toggle

end module rb_handlers




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!  H E R E    S T A R T S    T H E    A C T U A L    P R O G R A M  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! programmed by Roland Zimm 2014 based on publicly accessible sourcecode from the gtk-fortran project
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



program grn_editor

  use event_handlers
  use global_widgets
  use on_display_handlers
  use io
  use basic_handlers
  use rb_handlers !>>Miquel8-10-14

  implicit none
  
  integer(c_int) :: message_id, firstTab, secondTab, thirdTab, fourthTab, fifthTab
  integer :: width_e, height_e
  real*8 :: newval
  real*4 :: xal, yal
  real*8 :: parval
  integer(kind=c_int), dimension(3), target :: isel=[ (i-1,i=1,3) ] !>>Miquel8-10-14
  character(len=20) :: rb_label !>>Miquel8-10-14
  character*2 :: str02
!character*12,allocatable :: genenamex(:)
  character*140 :: inputtt
  integer :: cnt, currentw

  type(c_ptr) :: labeltc, buttonlc, buttonln, buttonld, buttonlm, buttonlp, labell, labelr

  call read_data(0)
  call read_data(1)

  call cell_params_names
  call gtk_disable_setlocale ()
  call gtk_init ()

  width_e = 1250
  height_e = 750

  if(allocated(ffut)) deallocate(ffut)
  allocate(ffut(nfu))

  ! Define all widgets
  my_window = gtk_window_new (GTK_WINDOW_TOPLEVEL)
  button2 = gtk_button_new_with_mnemonic ("_Save as PNG"//c_null_char)
  label1 = gtk_label_new("real(c)"//c_null_char)
  spinButton1 = gtk_spin_button_new (gtk_adjustment_new(-0.835d0,-2d0,+2d0,&
       & 0.05d0,0.5d0,0d0),0.05d0, 7_c_int)
  label2 = gtk_label_new("imag(c) "//c_null_char)
  spinButton2 = gtk_spin_button_new (gtk_adjustment_new(-0.2321d0, -2d0, &
       & +2d0,0.05d0,0.5d0,0d0),0.05d0, 7_c_int)
  label3 = gtk_label_new("iterations"//c_null_char)
  spinButton3 = gtk_spin_button_new (gtk_adjustment_new(1000d0,1d0, &
       & +100000d0,10d0,100d0,0d0),10d0, 0_c_int)
  button02 = gtk_button_new_with_mnemonic ("_EXIT_"//c_null_char)
  button03 = gtk_button_new_with_mnemonic ("_Run simulation_"//c_null_char)
  button04 = gtk_button_new_with_mnemonic ("_UNDO_"//c_null_char)
  button05 = gtk_button_new_with_mnemonic ("_Name genes_"//c_null_char)
  button06 = gtk_button_new_with_mnemonic ("_Information_"//c_null_char)
  button07 = gtk_button_new_with_mnemonic ("_Modify as Matrix_"//c_null_char)
  button08 = gtk_button_new_with_mnemonic ("_Change background"//c_null_char)
  button09 = gtk_button_new_with_mnemonic ("_Plot in new window_"//c_null_char)
  button10 = gtk_button_new_with_mnemonic ("_Save plot (PNG)_"//c_null_char)
  button11 = gtk_button_new_with_mnemonic ("_Add/Delete genes_"//c_null_char)
  button12 = gtk_button_new_with_mnemonic ("_Arrange genes_"//c_null_char)
  button13 = gtk_button_new_with_mnemonic ("_Save input file_"//c_null_char)
  buttonlc = gtk_button_new_with_mnemonic ("_Create + connection_"//c_null_char)
  buttonln = gtk_button_new_with_mnemonic ("_Create - connection_"//c_null_char)
  buttonld = gtk_button_new_with_mnemonic ("_Delete connection_"//c_null_char)
  buttonlm = gtk_button_new_with_mnemonic ("_Modify connection_"//c_null_char)
  buttonlp = gtk_button_new_with_mnemonic ("_Modify gene properties"//c_null_char)
  labell = gtk_label_new(""//c_null_char) 
  labelr = gtk_label_new("-------Node connections-------"//c_null_char) 
  label4 = gtk_label_new("Predefined values:"//c_null_char)           
  table = gtk_table_new (8_c_int, 2_c_int, TRUE)
  table2 = gtk_table_new (20_c_int, 1_c_int, TRUE)
  label00 = gtk_label_new("Number of genes:"//c_null_char)
  expander = gtk_expander_new_with_mnemonic ("_M params"//c_null_char)
  expander2 = gtk_expander_new_with_mnemonic ("_L params"//c_null_char)
  box2 = gtk_hbox_new (FALSE, 10_c_int)

rb_arrange_tog=0 !>>Miquel9-10-14

 !!!!!>>>Miquel8-10-14 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


     rb_label="Click on Network"
     rbut(1) = hl_gtk_radio_button_new(group, trim(rb_label)//c_null_char, &
          & toggled=c_funloc(rb_toggle), data=c_loc(isel(1)))

rb_label="Click on Buttons"
     rbut(2) = hl_gtk_radio_button_new(group, trim(rb_label)//c_null_char, &
          & toggled=c_funloc(rb_toggle), data=c_loc(isel(2)))

     !rb_label="Rearrange nodes"
     !rbut(3) = hl_gtk_radio_button_new(group, trim(rb_label)//c_null_char, &
     !     & toggled=c_funloc(rb_toggle), data=c_loc(isel(3)))

  call hl_gtk_radio_group_set_select(group, 3_c_int)

 !!!!!>>>Miquel8-10-14 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    if (allocated(param)) deallocate(param)
    allocate(param(nparam))

    if (allocated(paramo)) deallocate(paramo)
    allocate(paramo(nparam))

    if (allocated(pparam)) deallocate(pparam)
    allocate(pparam(mamax,nparam))

    if (allocated(names_param)) deallocate(names_param)
    allocate(names_param(nparam))

    if (allocated(varglobal_out)) deallocate(varglobal_out)
    allocate(varglobal_out(nvarglobal_out))

    if (allocated(pvarglobal_out)) deallocate(pvarglobal_out)
    allocate(pvarglobal_out(mamax,nvarglobal_out))

    if (allocated(names_varglobal_out)) deallocate(names_varglobal_out)
    allocate(names_varglobal_out(nvarglobal_out))

    if (allocated(names_fu)) deallocate(names_fu)
    allocate(names_fu(nfu))

    if (allocated(nodeparams)) deallocate(nodeparams)
    allocate(nodeparams(nparam_per_node))
    if (allocated(ffu)) deallocate(ffu)
    allocate(ffu(nfu))
    if (allocated(param)) deallocate(param)
    allocate(param(nparam))
    if (allocated(varglobal_out)) deallocate(varglobal_out)
    allocate(varglobal_out(nvarglobal_out))

    call s_names_fu
    call s_names_param
    call s_names_varglobal_out
    call s_names_nodeparams

param=1

  do i=1,nparam
    labelt(i)=gtk_label_new(adjustl(trim(names_param(i)))//c_null_char)

    parval=param(i)

    if(abs(parval*100).lt.10) parval=0.1d0
    if((i==1).or.(i==10).or.(i==9).or.(i==13).or.(i==18).or.(i==19).or.(i==20))then
      spinbutt(i)=gtk_spin_button_new (gtk_adjustment_new(param(i),param(i),param(i),&
       & 0.0d0,0.5d0,0d0),0.05d0, 0_c_int)
    else if(i==27)then
      spinbutt(i)=gtk_spin_button_new (gtk_adjustment_new(param(i),0d0,+1d0,&
       & 1d0,0.5d0,0d0),0.05d0, 0_c_int)
    else if(i==12)then
      spinbutt(i)=gtk_spin_button_new (gtk_adjustment_new(param(i),-100*parval,+100*parval,&
       & 0.01*parval,0.5d0,0d0),0.05d0, 0_c_int)
    else if(i==26)then  ! not described in the methods
      spinbutt(i)=gtk_spin_button_new (gtk_adjustment_new(param(i),param(i),param(i),&
       & 0.0d0,0.0d0,0d0),0.00d0, 3_c_int)
    else if(i==3)then  ! not described in the methods
      spinbutt(i)=gtk_spin_button_new (gtk_adjustment_new(param(i),param(i),param(i),&
       & 0.0d0,0.0d0,0d0),0.00d0, 0_c_int)
    else if(i==4)then  ! not described in the methods
      spinbutt(i)=gtk_spin_button_new (gtk_adjustment_new(param(i),param(i),param(i),&
       & 0.0d0,0.0d0,0d0),0.00d0, 3_c_int)
    else if(i==7)then  ! not described in the methods
      spinbutt(i)=gtk_spin_button_new (gtk_adjustment_new(param(i),param(i),param(i),&
       & 0.0d0,0.0d0,0d0),0.00d0, 3_c_int)
    else if(i==14)then  ! not described in the methods
      spinbutt(i)=gtk_spin_button_new (gtk_adjustment_new(param(i),param(i),param(i),&
       & 0.0d0,0.0d0,0d0),0.00d0, 0_c_int)
    else
      if(parval.ge.0.001)then
        spinbutt(i)=gtk_spin_button_new (gtk_adjustment_new(param(i),-100*parval,+100*parval,&
         & 0.01*parval,0.5d0,0d0),0.05d0, 3_c_int)
      elseif(parval.le.-0.001)then
        spinbutt(i)=gtk_spin_button_new (gtk_adjustment_new(param(i),+100*parval,-100*parval,&
         & 0.01*parval,0.5d0,0d0),0.05d0, 3_c_int)
      else
        spinbutt(i)=gtk_spin_button_new (gtk_adjustment_new(param(i),-0.1d0,0.1d0,&
         & 0.0001d0,0.5d0,0d0),0.05d0, 3_c_int)
      end if
    end if
  end do

! connect the buttons with their button_press_events
  call gtk_window_set_default_size(my_window, width_e, height_e)
  call gtk_window_set_title(my_window, "NetworkMaker"//c_null_char)
  call g_signal_connect (button02, "clicked"//c_null_char, c_funloc(quit_the_shit))!name_enter
  call g_signal_connect (button03, "clicked"//c_null_char, c_funloc(run_it))
  call g_signal_connect (button04, "clicked"//c_null_char, c_funloc(undo_redo_connect))
  call g_signal_connect (button05, "clicked"//c_null_char, c_funloc(name_all_connect))
  call g_signal_connect (button06, "clicked"//c_null_char, c_funloc(info_connect))
  call g_signal_connect (button07, "clicked"//c_null_char, c_funloc(GRN_table_editable_connect))
  call g_signal_connect (button08, "clicked"//c_null_char, c_funloc(GRN_style))
  call g_signal_connect (button09, "clicked"//c_null_char, c_funloc(plot_extra))
  call g_signal_connect (button10, "clicked"//c_null_char, c_funloc(save_it_plot))
  call g_signal_connect (button11, "clicked"//c_null_char, c_funloc(GRN_table_add_connect))
  call g_signal_connect (buttonlc, "clicked"//c_null_char, c_funloc(add_plus_connect))
  call g_signal_connect (buttonln, "clicked"//c_null_char, c_funloc(add_minus_connect))
  call g_signal_connect (buttonld, "clicked"//c_null_char, c_funloc(delete_connect))
  call g_signal_connect (buttonlm, "clicked"//c_null_char, c_funloc(choice_box_build_connect))
  call g_signal_connect (buttonlp, "clicked"//c_null_char, c_funloc(gene_properties_connect))
  call g_signal_connect (button12, "clicked"//c_null_char, c_funloc(move_nodes_connect))
  call g_signal_connect (button13, "clicked"//c_null_char, c_funloc(save_it))

! display the array of flags
    do i=1,nparam
    write (str,"(I2.1)") i+1
    str=trim(str)//"_c_int"
    str=trim(str)
    write (str2,"(I2.1)") i+2
    str2=trim(str2)//"_c_int"
    str2=trim(str2)
    labeltc=gtk_label_new((trim(names_param(i)))//c_null_char)

    call gtk_table_attach_defaults(table, spinbutt(i), 0_c_int, 1_c_int, i+1, i+2)

    xal=0d0; yal=0d0
    call gtk_misc_set_alignment(labeltc, xal, yal)

    call gtk_table_attach_defaults(table, labeltc, 1_c_int, 6_c_int, i+1, i+2)
  end do

  cnt=0
  do i=1,4
    cnt=cnt+1
    write (str,"(I2.1)") i+1
    str=trim(str)//"_c_int"
    str=trim(str)
    write (str2,"(I2.1)") i+2
    str2=trim(str2)//"_c_int"
    str2=trim(str2)
    write (str4,"(I2.1)") cnt
    str=trim(str4)//"_c_int"
    str=trim(str4)
    ffut(i)=gtk_check_button_new_with_label("L"//trim(str4)//adjustl(trim(": "//names_fu(i)))//c_null_char)

    call gtk_table_attach_defaults(table2, ffut(i), 0_c_int, 1_c_int, i+1, i+2)

  end do
  do i=7,8
    cnt=cnt+1
    write (str,"(I2.1)") i-1
    str=trim(str)//"_c_int"
    str=trim(str)
    write (str2,"(I2.1)") i
    str2=trim(str2)//"_c_int"
    str2=trim(str2)
    write (str4,"(I2.1)") cnt
    str=trim(str4)//"_c_int"
    str=trim(str4)
    ffut(i)=gtk_check_button_new_with_label("L"//trim(str4)//adjustl(": "//trim(names_fu(i)))//c_null_char)
    call gtk_table_attach_defaults(table2, ffut(i), 0_c_int, 1_c_int, i-1, i)
  end do

  cnt=cnt+1
  do i=10,13

    cnt=cnt+1
    write (str,"(I2.1)") i
    str=trim(str)//"_c_int"
    str=trim(str)
    write (str2,"(I2.1)") i+1
    str2=trim(str2)//"_c_int"
    str2=trim(str2)
    write (str4,"(I2.1)") cnt
    str=trim(str4)//"_c_int"
    str=trim(str4)
    ffut(i)=gtk_check_button_new_with_label("L"//trim(str4)//adjustl(": "//trim(names_fu(i)))//c_null_char)
    call gtk_table_attach_defaults(table2, ffut(i), 0_c_int, 1_c_int, i, i+1)
  end do
  do i=16,18
    cnt=cnt+1
    write (str,"(I2.1)") i-2
    str=trim(str)//"_c_int"
    str=trim(str)
    write (str2,"(I2.1)") i-1
    str2=trim(str2)//"_c_int"
    str2=trim(str2)
    write (str4,"(I2.1)") cnt
    str=trim(str4)//"_c_int"
    str=trim(str4)
    ffut(i)=gtk_check_button_new_with_label("L"//trim(str4)//adjustl(": "//trim(names_fu(i)))//c_null_char)
    call gtk_table_attach_defaults(table2, ffut(i), 0_c_int, 1_c_int, i-2, i-1)

  end do
  cnt=cnt+1
  do i=20,nfu
    cnt=cnt+1
    write (str,"(I2.1)") i-1
    str=trim(str)//"_c_int"
    str=trim(str)
    write (str2,"(I2.1)") i
    str2=trim(str2)//"_c_int"
    str2=trim(str2)
    write (str4,"(I2.1)") cnt
    str=trim(str4)//"_c_int"
    str=trim(str4)

    ffut(i)=gtk_check_button_new_with_label("L"//trim(str4)//adjustl(": "//trim(names_fu(i)))//c_null_char)

    call gtk_table_attach_defaults(table2, ffut(i), 0_c_int, 1_c_int, i-1, i)
  end do

  do i=1,4

    if(ffu(i)==1) call gtk_toggle_button_set_active (ffut(i), TRUE)
    ffuval= gtk_toggle_button_get_active(ffut(i))
  end do

  do i=7,8

    if(ffu(i)==1) call gtk_toggle_button_set_active (ffut(i), TRUE)
    ffuval= gtk_toggle_button_get_active(ffut(i))
  end do
  do i=10,13

    if(ffu(i)==1) call gtk_toggle_button_set_active (ffut(i), TRUE)
    ffuval= gtk_toggle_button_get_active(ffut(i))
  end do
  do i=16,18

    if(ffu(i)==1) call gtk_toggle_button_set_active (ffut(i), TRUE)
    ffuval= gtk_toggle_button_get_active(ffut(i))
  end do
  do i=20,nfu

    if(ffu(i)==1) call gtk_toggle_button_set_active (ffut(i), TRUE)
    ffuval= gtk_toggle_button_get_active(ffut(i))
  end do

! now for the non-binary flags (9)
  labeltc=gtk_label_new("L7: "//adjustl(trim(names_fu(9)))//c_null_char)
  spinbutt9=gtk_spin_button_new (gtk_adjustment_new(param(9),0d0,2d0,&
       & 1.0d0,0.5d0,0d0),0.05d0, 0_c_int)
  call gtk_table_attach_defaults(table2, labeltc, 0_c_int, 1_c_int, 8_c_int, 9_c_int)
  call gtk_table_attach_defaults(table2, spinbutt9, 0_c_int, 1_c_int, 9_c_int, 10_c_int)
  labeltc=gtk_label_new("L15: "//adjustl(trim(names_fu(19)))//c_null_char)
  spinbutt8=gtk_spin_button_new (gtk_adjustment_new(param(19),0d0,2d0,&
       & 1.0d0,0.5d0,0d0),0.05d0, 0_c_int)
  call gtk_table_attach_defaults(table2, labeltc, 0_c_int, 1_c_int, 17_c_int, 18_c_int)
  call gtk_table_attach_defaults(table2, spinbutt8, 0_c_int, 1_c_int, 18_c_int, 19_c_int)

! pack the expanders and initialize
  call gtk_container_add (expander, table)
  call gtk_expander_set_expanded(expander, FALSE)
  call gtk_container_add (expander2, table2)
  call gtk_expander_set_expanded(expander2, FALSE)
  call gtk_box_pack_start (box2, expander, FALSE, FALSE, 0_c_int)
  call gtk_box_pack_start (box2, expander2, FALSE, FALSE, 1_c_int)
! now draw the notebook windows individually
! w-Matrix
  my_drawing_area = gtk_drawing_area_new()

  call GRN_plot(1)

  eventabbox1 = gtk_event_box_new()
  call gtk_container_add (eventabbox1, my_drawing_area)
  
  notebook = gtk_notebook_new ()
  call g_signal_connect(eventabbox1, "button-press-event"//c_null_char, c_funloc(cursor_choice))

  notebookLabel1 = gtk_label_new_with_mnemonic("_T-Matrix"//c_null_char)
  firstTab = gtk_notebook_append_page (notebook, eventabbox1, notebookLabel1)
  ! ww-Matrix
  my_drawing_area = gtk_drawing_area_new()
  call GRN_plot(3)

  eventabbox1 = gtk_event_box_new()

  call gtk_container_add (eventabbox1, my_drawing_area)

  call g_signal_connect(eventabbox1, "button-press-event"//c_null_char, c_funloc(cursor_choice))
  
  notebookLabel2 = gtk_label_new_with_mnemonic("_R-Matrix"//c_null_char)
  secondTab = gtk_notebook_append_page (notebook, eventabbox1, notebookLabel2)

  ! E-Matrix
  my_drawing_area = gtk_drawing_area_new()
  call GRN_plot(5)

  eventabbox1 = gtk_event_box_new()

  call gtk_container_add (eventabbox1, my_drawing_area)

  call g_signal_connect(eventabbox1, "button-press-event"//c_null_char, c_funloc(cursor_choice))
  
  notebookLabel4 = gtk_label_new_with_mnemonic("_E-Matrix"//c_null_char)
  fourthTab = gtk_notebook_append_page (notebook, eventabbox1, notebookLabel4)
  
  ! C-Matrix
  my_drawing_area = gtk_drawing_area_new()
  call GRN_plot(4)

  eventabbox1 = gtk_event_box_new()

  call gtk_container_add (eventabbox1, my_drawing_area)

  call g_signal_connect(eventabbox1, "button-press-event"//c_null_char, c_funloc(cursor_choice))
  
  notebookLabel5 = gtk_label_new_with_mnemonic("_C-Matrix"//c_null_char)
  fifthTab = gtk_notebook_append_page (notebook, eventabbox1, notebookLabel5)

  ! kadh-Matrix
  
  my_drawing_area = gtk_drawing_area_new()
  call GRN_plot(6)

  eventabbox1 = gtk_event_box_new()

  call gtk_container_add (eventabbox1, my_drawing_area)
  call g_signal_connect(eventabbox1, "button-press-event"//c_null_char, c_funloc(cursor_choice))
  
  notebookLabel6 = gtk_label_new_with_mnemonic("_B-Matrix"//c_null_char)
  fifthTab = gtk_notebook_append_page (notebook, eventabbox1, notebookLabel6)  
  call gtk_box_pack_start (box2, notebook, TRUE, TRUE, 0_c_int)

! place the right sided buttons

  table2 = gtk_table_new (10_c_int, 2_c_int, TRUE)
  call gtk_table_attach_defaults(table2, button03, 0_c_int, 2_c_int, 6_c_int, 7_c_int)  !>>Miquel8-10-14
  call gtk_table_attach_defaults(table2, button13, 0_c_int, 2_c_int, 7_c_int, 8_c_int) !
  call gtk_table_attach_defaults(table2, button10, 0_c_int, 2_c_int, 8_c_int, 9_c_int) !
  call gtk_table_attach_defaults(table2, button09, 0_c_int, 2_c_int, 9_c_int, 10_c_int) !
  call gtk_table_attach_defaults(table2, button11, 0_c_int, 2_c_int, 10_c_int, 11_c_int) !
  call gtk_table_attach_defaults(table2, button12, 0_c_int, 2_c_int, 11_c_int, 12_c_int) !
  call gtk_table_attach_defaults(table2, button05, 0_c_int, 2_c_int, 12_c_int, 13_c_int)!
  call gtk_table_attach_defaults(table2, button07, 0_c_int, 2_c_int, 13_c_int, 14_c_int) !
  call gtk_table_attach_defaults(table2, button08, 0_c_int, 2_c_int, 14_c_int, 15_c_int) !
  call gtk_table_attach_defaults(table2, button06, 0_c_int, 2_c_int, 15_c_int, 16_c_int)!
  call gtk_table_attach_defaults(table2, labell, 0_c_int, 2_c_int, 4_c_int, 5_c_int)!
  call gtk_table_attach_defaults(table2, buttonlc, 0_c_int, 2_c_int, 0_c_int, 1_c_int) !
  call gtk_table_attach_defaults(table2, buttonln, 0_c_int, 2_c_int, 1_c_int, 2_c_int)!
  call gtk_table_attach_defaults(table2, buttonld, 0_c_int, 2_c_int, 2_c_int, 3_c_int) !
  call gtk_table_attach_defaults(table2, buttonlm, 0_c_int, 2_c_int, 3_c_int, 4_c_int)!
  call gtk_table_attach_defaults(table2, buttonlp, 0_c_int, 2_c_int, 4_c_int, 5_c_int)!
  call gtk_table_attach_defaults(table2, button04, 0_c_int, 2_c_int, 16_c_int, 17_c_int) !
  call gtk_table_attach_defaults(table2, button02, 0_c_int, 2_c_int, 21_c_int, 22_c_int)!
  call gtk_table_attach_defaults(table2, rbut(2), 0_c_int, 2_c_int, 18_c_int, 19_c_int)!
  call gtk_table_attach_defaults(table2, rbut(1), 0_c_int, 2_c_int, 19_c_int, 20_c_int)!
  call gtk_table_attach_defaults(table2, rbut(3), 0_c_int, 2_c_int, 20_c_int, 21_c_int)!
  call gtk_table_attach_defaults(table2, labelr, 0_c_int, 2_c_int, 17_c_int, 18_c_int)

  call gtk_box_pack_start (box2, table2, FALSE, FALSE, 0_c_int)
  call gtk_container_add (my_window, box2)
  call gtk_window_set_mnemonics_visible (my_window, TRUE)

! show all
  call gtk_widget_show_all (my_window)
  
! This is maintaining the visualization for storage into a file
  pixwidth_e  = 500
  pixheight_e = 500
  my_pixbuf = gdk_pixbuf_new(GDK_COLORSPACE_RGB, FALSE, 8_c_int, pixwidth_e, pixheight_e)    
  call c_f_pointer(gdk_pixbuf_get_pixels(my_pixbuf), pixel, (/0/))
  nch = gdk_pixbuf_get_n_channels(my_pixbuf)
  rowstride = gdk_pixbuf_get_rowstride(my_pixbuf)
  pixel = char(0)

! Everything. Call main.
  call gtk_main ()

end program grn_editor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!gfortran -w general.mod.f90 io.mod.f90 genetic.mod.f90 grn_editor_whole.f90 -o grnee.e `pkg-config --cflags --libs /home/rolazimm/Desktop/gtk-fortran/gtk-fortran-master/build/src/gtk-2-fortran.pc gtk+-2.0`

