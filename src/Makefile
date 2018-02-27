# ============================================================================
# Name        : Makefile
# Author      : Kensuke Konishi
# Version     : 0.0.1
# Copyright   : It is Complicated.
# Description : Makefile for Hello World in Fortran
# ============================================================================

.PHONY: all clean

# Change this line if you are using a different Fortran compiler
FC = gfortran
FC = ifort

TARGET = ../bin/copt2d
SRCS = copt2d_cpml.f90 parameter_2dwave.f90 main2d_cpml_v2.f90
OBJS = parameter.o cpml.o main.o
MODS = parameter_2dwave.mod copt2d_cpml.mod
SWITCH = -O3 -parallel -par-report1
#DEBUG = -debug -CB -warn all -g
.SUFFIXES: .f90

all: $(TARGET)

$(TARGET): $(OBJS)
	$(FC) $(DEBUG) $(SWITCH) -o $@ $(OBJS)

%.o: %.f90
	$(FC) $(DEBUG) $(SWITCH) $<
	
parameter.o: parameter_2dwave.f90
	$(FC) $(DEBUG) -o parameter.o -c $(SWITCH) parameter_2dwave.f90

cpml.o: copt2d_cpml.f90
	$(FC) $(DEBUG) -o cpml.o -c $(SWITCH) copt2d_cpml.f90

main.o: $(MODS) main2d_cpml.f90 declaration.f90
	$(FC) $(DEBUG) -o main.o -c $(SWITCH) main2d_cpml.f90
	
parameter_2dwave.mod: parameter.o parameter_2dwave.f90
	
copt2d_cpml.mod: cpml.o copt2d_cpml.f90

clean:
	rm -f $(OBJS) $(MODS) $(TARGET) 
 