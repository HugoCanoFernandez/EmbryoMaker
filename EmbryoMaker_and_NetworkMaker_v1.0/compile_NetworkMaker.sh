#!/bin/bash

cp src/net/* bin
cp src/core/general.mod.f90 bin
cp src/core/genetic.mod.f90 bin
cp src/core/io.mod.f90 bin
cp src/core/neighboring.mod.f90 bin
cp src/core/geompack3.f90 bin

cd bin

YOURPATH=$(pwd)
  sed -i 's:MYPATH:'$YOURPATH':g' gtk-2-fortran.pc 
echo $(pwd)

gfortran -w gdk.f90 glib-auto.f90 unixonly-auto.f90 gtk.f90 gtk-sup.f90 gtk-hl-container.f90 gtk-hl-entry.f90 gtk-hl-misc.f90 gtk-hl-accelerator.f90 gtk-hl-button.f90 gdk-pixbuf-auto.f90 gtk-hl-tree.f90 gtk-hl-chooser.f90 gtk-hl.f90 gdkevents-auto2.f90 gnomo-net.f90 -o NetworkMaker-launcher `pkg-config --cflags --libs gtk+-2.0` -L$(pwd) `pkg-config --cflags --libs gtk-2-fortran.pc`

gfortran -w general.mod.f90 geompack3.f90 neighboring.mod.f90 genetic.mod.f90 io.mod.f90 gdk.f90 glib-auto.f90 unixonly-auto.f90 gtk.f90 gtk-sup.f90 gtk-hl-container.f90 gtk-hl-entry.f90 gtk-hl-misc.f90 gtk-hl-accelerator.f90 gtk-hl-button.f90 gdk-pixbuf-auto.f90 gtk-hl-tree.f90 gtk-hl-chooser.f90 gtk-hl.f90 gdkevents-auto2.f90 cairo-auto.f90 grn_editor_whole1.f90 -o NetworkMaker `pkg-config --cflags --libs gtk+-2.0` -L$(pwd) `pkg-config --cflags --libs gtk-2-fortran.pc`

sed -i 's:'$YOURPATH':MYPATH:g' gtk-2-fortran.pc
rm *mod
rm *.f90
rm gtk* lib*

cd ..

echo 'executables installed in bin/'
echo
echo 'Run NetworkMaker-launcher to select an EmbryoMaker' 
echo 'file to edit or start a new network from scratch'
echo
echo 'If you already have a file to edit you can just run '
echo 'NetworkMaker and enter the path to the file as first'
echo 'argument in the command line'

