#!/bin/sh
##########################################
#
# Build Fortran modules and compile main EA2D function as a Python module.
#
# Pedro Mediano, Mar 2020
#
##########################################

####
# Added f2py flags and gfortran flags for shared object, 
# otherwise it would not compile for me.
####
export LDFLAGS=-Wl,-rpath=.
export NPY_DISTUTILS_APPEND_FLAGS=1
gfortran -shared -fPIC variables.f90 -o variables.o
gfortran -shared -fPIC mt95f4.f90 -o mt95f4.o
python -m numpy.f2py -c -m EA2D EA2D.f90 mt95f4.o variables.o

