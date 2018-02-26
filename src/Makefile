# ============================================================================
# Name        : Makefile
# Author      : Kensuke Konishi
# Version     : 0.0.1
# Copyright   : It is Complicated.
# Description : Makefile for Hello World in Fortran
# ============================================================================

.PHONY: all clean

# Change this line if you are using a different Fortran compiler
FORTRAN_COMPILER = gfortran
#FORTRAN_COMPILER = ifort

TARGET = ../bin/copt2d
SRCS = copt2d_cpml.f90 parameter_2dwave.f90 main2d_cpml_v2.f90
OBJS = parameter.o cpml.o main.o
MODS = parameter_2dwave.mod copt2d_cpml.mod
SWITCH = -O3
#DEBUG = -debug -check bounds -warn interface -g
.SUFFIXES: .f90

#bin/%.o: src/%.f90
#	$(FORTRAN_COMPILER) -o $@ -c $(SWITCH) $< 

all: $(TARGET)

$(TARGET): $(OBJS)
	$(FORTRAN_COMPILER) $(DEBUG) $(SWITCH) -o $@ $(OBJS)


%.o: %.f90
	$(FORTRAN_COMPILER) $(SWITCH) $<
	
parameter.o: parameter_2dwave.f90
	$(FORTRAN_COMPILER) -o parameter.o -c $(SWITCH) parameter_2dwave.f90

cpml.o: copt2d_cpml.f90
	$(FORTRAN_COMPILER) -o cpml.o -c $(SWITCH) copt2d_cpml.f90

main.o: $(MODS) main2d_cpml.f90 declaration.f90
	$(FORTRAN_COMPILER) -o main.o -c $(SWITCH) main2d_cpml.f90
	
parameter_2dwave.mod: parameter.o parameter_2dwave.f90
	
copt2d_cpml.mod: cpml.o copt2d_cpml.f90



clean:
	rm -f $(OBJS) $(MODS) $(TARGET) 
 
