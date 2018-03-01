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
MPIFC = mpif90

MPITARGET = ../bin/mpi_copt2d
TARGET = ../bin/copt2d
SRCS = copt2d_cpml.f90 parameter_2dwave.f90 main2d_cpml_v2.f90
OBJS = parameter.o cpml.o 
MODS = parameter_2dwave.mod copt2d_cpml.mod
SWITCH = -O3
#PARALLEL = -parallel -par-report1
#DEBUG = -debug -CB -warn all -g
FLAGS= $(SWITCH) $(PARALLEL) $(DEBUG)
.SUFFIXES: .f90

all: $(TARGET) $(MPITARGET)

$(TARGET): $(OBJS) main.o
	$(FC) $(FLAGS) -o $@ $(OBJS) main.o

$(MPITARGET): $(OBJS) mpi_main.o
	$(MPIFC) $(FLAGS) -o $@ $(OBJS) mpi_main.o

%.o: %.f90
	$(FC) $(FLAGS) $<
	
parameter.o: parameter_2dwave.f90
	$(FC) $(FLAGS) -fast -o parameter.o -c parameter_2dwave.f90

cpml.o: copt2d_cpml.f90
	$(FC) $(FLAGS) -fast -o cpml.o -c copt2d_cpml.f90

main.o: $(MODS) main2d_cpml.f90 declaration.f90
	$(FC) $(FLAGS) -fast -o main.o -c main2d_cpml.f90

mpi_main.o: main2d_cpml_mpi.f90
	$(MPIFC) $(FLAGS) -o mpi_main.o -c main2d_cpml_mpi.f90
		
parameter_2dwave.mod: parameter.o parameter_2dwave.f90
	
copt2d_cpml.mod: cpml.o copt2d_cpml.f90

clean:
	rm -f $(OBJS) $(MODS) $(TARGET) main.o mpi_main.o
 