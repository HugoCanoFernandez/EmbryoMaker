
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


module basics

use gtk
use gdk
use gtk_hl
use iso_c_binding, only: c_ptr, c_char, c_int

implicit none

type(c_ptr) :: window

end module basics


module functions

use basics

contains

subroutine my_destroy2(widget, gdata) bind(c)

  type(c_ptr), value :: widget, gdata
  integer(kind=c_int) :: ok

  call gtk_widget_destroy(window)
  call gtk_main_quit ()

end subroutine my_destroy2

function choose_input (widget, gdata ) result(ret)  bind(c)

  use iso_c_binding, only: c_ptr, c_int
  use basics

  !use Input_handlers
  implicit none

  type(c_ptr) :: base, jb, junk
  type(c_ptr), value :: widget, gdata
  integer(c_int)    :: ret

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
  !junk = hl_gtk_button_new("CHOOSE (FROM THE SAME FOLDER)!"//c_null_char, clicked=c_funloc(open_file))
  call hl_gtk_box_pack(jb, junk)
 ! junk = hl_gtk_button_new("Print choices"//c_null_char, clicked=c_funloc(display_file))
 ! call hl_gtk_box_pack(base, junk)
 ! junk = hl_gtk_button_new("Display choices"//c_null_char, clicked=c_funloc(read_data_again_ini))
 ! call hl_gtk_box_pack(base, junk)
  junk = hl_gtk_button_new("-----------------Quit-----------------"//c_null_char, clicked=c_funloc(my_destroy2))
  call hl_gtk_box_pack(base, junk)

  call gtk_widget_show_all(window)

  call gtk_main()

  !filename3=filename

end function choose_input

subroutine do_open(widget, gdata) bind(c)
    
  use gtk_hl
  use gtk_os_dependent
  use basics
  use iso_c_binding, only: c_ptr, c_int

    type(c_ptr), value :: widget, gdata

    type(c_ptr) :: c_string
    character(len=200) :: inln
    integer :: ios
    integer :: idxs
    character(len=120) :: filename, filenamen=''
    character(len=140) :: input

    call gtk_window_set_title(window, "Choose an input file"//c_null_char)
    c_string = gtk_file_chooser_get_filename(widget)
    call convert_c_string(c_string, filenamen)
    call g_free(c_string)

    idxs = index(filenamen, '/', .true.)+1
    call gtk_window_set_title(window, trim(filenamen(idxs:))//c_null_char)
    print*, filenamen
    filename=trim(filenamen)

    print*, "FILENAME", filename
    input="cp "//trim(filename)//" newchoice.dat" ! saves the file as a local "newchoice.dat"
    call system(input)
    input="./NetworkMaker newchoice.dat" ! This launches the interface directly
    call system(input)

end subroutine do_open

subroutine open_file(widget, gdata) bind(c)

  use g
  use gtk_hl
  use gtk_hl_container
  use gtk_hl_entry
  use gtk_hl_button
  use gtk_hl_chooser
  use gdk_events
  use gdk
  use gtk
  use gtk_os_dependent
  use gtk_sup
  type(c_ptr), value :: widget, gdata
  type(c_ptr) :: window

  integer(kind=c_int) :: isel
  character(len=240), dimension(:), allocatable :: chfile
  character(len=240) :: filename
  character(len=30), dimension(4) :: filters
  character(len=30), dimension(4) :: filtnames
print*, "EINS"
  filters(1) = "*"
  filters(2) = "*.txt,*.lis"
  filters(3) = "*.f90"
  filters(4) = "*.dat,*.log"
  filtnames(1) = "All files"
  filtnames(2) = "Text files"
  filtnames(3) = "Fortran code"
  filtnames(4) = "Output files"
print*, "EINSKOMMAFÜNF"
  isel = hl_gtk_file_chooser_show(chfile, create=FALSE,&
         & title="Select input file"//c_null_char, filter=filters, &
         & filter_name=filtnames, wsize=(/ 600_c_int, 400_c_int /), &
         & edit_filters=TRUE, &
         & parent=window, all=TRUE)
print*, "ZWEI"
  if (isel == FALSE) return   ! No selection made
print*, "DREI"
  filename = chfile(1)
  deallocate(chfile)
print*, "VIER"

  !!!call gtk_widget_destroy(window)

end subroutine open_file

function Start_program (widget, gdata ) result(ret)  bind(c)

  use iso_c_binding, only: c_ptr, c_int

  implicit none
    
  type(c_ptr), value :: widget, gdata
  integer(c_int)    :: ret
  character*140 :: dir
  character*150 :: ddir

  call getarg(1,dir)

  if(len(dir).lt.1)then  
    ddir="./NetworkMaker "//trim(dir)
  else
    ddir="./NetworkMaker "
  endif
  call system(ddir)

end function Start_program

function Start_program_0 (widget, gdata ) result(ret)  bind(c)

  use iso_c_binding, only: c_ptr, c_int

  implicit none
    
  type(c_ptr), value :: widget, gdata
  integer(c_int)    :: ret
  character*140 :: dir
  character*150 :: ddir

  call getarg(1,dir)
  ddir="./NetworkMaker "
  call system(ddir)

end function Start_program_0


function End_program (widget, gdata ) result(ret)  bind(c)

  use iso_c_binding, only: c_ptr, c_int

  implicit none
    
  type(c_ptr), value :: widget, gdata
  integer(c_int)    :: ret

  call gtk_main_quit()

end function End_program

end module functions


program grn_interface

  use iso_c_binding, only: c_ptr, c_char, c_int
  use gtk
  use gdk, only: gdk_cairo_create, gdk_cairo_set_source_pixbuf
  
  use gdk_pixbuf, only: gdk_pixbuf_get_n_channels, gdk_pixbuf_get_pixels, gdk_pix&
  &buf_get_rowstride, gdk_pixbuf_new

  use functions
  use basics

  !use event_handlers
  !use global_widgets
  !use on_display_handlers
  !use io
  !use basic_handlers
  !use rb_handlers !>>Miquel8-10-14

  implicit none
  
  type(c_ptr) :: table, box, my_window, button1, button2, button3, button4, my_pixbuf, base, junk, jb
  integer :: width_e, height_e
  real*8 :: newval
  real*4 :: xal, yal
  real*8 :: parval
  character(kind=c_char), dimension(:), pointer :: pixel
  integer(kind=c_int) :: nch, rowstride, pixwidth_e, pixheight_e
  logical :: computing = .false.
  character(LEN=80) :: string
  character*140 :: di
  character*150 :: ddi

  integer(c_int)    :: ret

  character(len=30), dimension(3) :: filters
  character(len=30), dimension(3) :: filtnames
  character(len=240) :: filename3
  filters(3) = "text/plain"
  filters(2) = "*.f90"
  filters(1) = "*"
  filtnames(3) = "Text files"
  filtnames(2) = "Fortran code"
  filtnames(1) = "All files"

  !type(c_ptr) :: labeltc

  call getarg(1,di)
  ddi="./NetworkMaker "//trim(di)
  print*, "STARTTTTTTTTTTT", di, len(trim(di))

  if(len(trim(di))>1) call system(ddi)
    if(len(trim(di))==1)then
      ddi="./NetworkMaker void"
      call system(ddi)
    end if

  call gtk_init ()
  ! Properties of the main_1 window :
  width_e = 300
  height_e = 150
  ! Define all widgets
  my_window = gtk_window_new (GTK_WINDOW_TOPLEVEL)
  button1 = gtk_label_new (" Choose a new input file "//c_null_char)
  button2 = gtk_button_new_with_mnemonic ("*** START ***"//c_null_char)
  button3 = gtk_button_new_with_mnemonic ("Quit"//c_null_char)
  button4 = gtk_button_new_with_mnemonic ("Create a new GRN"//c_null_char)
  table = gtk_table_new (2_c_int, 2_c_int, TRUE)
  box = gtk_hbox_new (FALSE, 2_c_int)

  call gtk_window_set_default_size(my_window, width_e, height_e)
  call gtk_window_set_title(my_window, "Launch NetworkMaker"//c_null_char)
  call g_signal_connect (button2, "clicked"//c_null_char, c_funloc(Start_program))
  call g_signal_connect (button4, "clicked"//c_null_char, c_funloc(Start_program_0))
  call g_signal_connect (button3, "clicked"//c_null_char, c_funloc(End_program))

  call gtk_table_attach_defaults(table, button1, 0_c_int, 2_c_int, 0_c_int, 1_c_int)  !>>Miquel8-10-14
  call gtk_table_attach_defaults(table, button2, 0_c_int, 4_c_int, 2_c_int, 4_c_int) !
  call gtk_table_attach_defaults(table, button3, 0_c_int, 4_c_int, 4_c_int, 5_c_int) 
  call gtk_table_attach_defaults(table, button4, 0_c_int, 4_c_int, 1_c_int, 2_c_int)

  junk = hl_gtk_file_chooser_button_new(title="Choose an input file"//c_null_char, &
       & filter=filters, filter_name=filtnames, file_set=c_funloc(do_open))

  call gtk_table_attach_defaults(table, junk, 2_c_int, 4_c_int, 0_c_int, 1_c_int)
! show all
  call gtk_box_pack_start (box, table, FALSE, FALSE, 0_c_int)
  call gtk_container_add (my_window, box)
  call gtk_window_set_mnemonics_visible (my_window, TRUE)

  call gtk_widget_show_all (my_window)

  call gtk_main ()

end program grn_interface

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!gfortran -w general.mod.f90 io.mod.f90 genetic.mod.f90 grn_editor_whole.f90 -o grnee.e `pkg-config --cflags --libs /home/rolazimm/Desktop/gtk-fortran/gtk-fortran-master/build/src/gtk-2-fortran.pc gtk+-2.0`


