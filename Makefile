# ============================================================================
# Name        : Makefile
# Author      : Kei Hasegawa
# Version     :
# Copyright   : It is Complicated.
# Description : Makefile for Hello World in Fortran
# ============================================================================
 
 # ============================================================================
# Name        : Makefile
# Author      : Kei Hasegawa, Kensuke Konishi
# Version     : 0.0.1
# Copyright   : It is Complicated.
# Description : Makefile for Hello World in Fortran
# ============================================================================

.PHONY: all clean

# Change this line if you are using a different Fortran compiler
FORTRAN_COMPILER = gfortran
FORTRAN_COMPILER = ifort

TARGET = bin/copt2d
SRCS = src/copt2d_cpml.f90 src/parameter_copt_run.f90 src/main2d_cpml_v2.f90
OBJS = src/parameter.o src/cpml.o src/main.o 
MODS = src/parameter.mod src/copt2d_cpml.mod 
SWITCH = -O3
#DEBUG = -debug -check bounds -warn interface -g
.SUFFIXES: .f90

#bin/%.o: src/%.f90
#	$(FORTRAN_COMPILER) -o $@ -c $(SWITCH) $< 

all: $(TARGET)
	@:

$(TARGET): $(OBJS)
	$(FORTRAN_COMPILER) $(DEBUG) $(SWITCH) -o $@ $(OBJS)

src/parameter.o: src/parameter_copt_run.f90
	@ cd src
	$(FORTRAN_COMPILER) -o src/parameter.o -c $(SWITCH) src/parameter_copt_run.f90

src/cpml.o: src/copt2d_cpml.f90
	@ cd src
	$(FORTRAN_COMPILER) -o src/cpml.o -c $(SWITCH) src/copt2d_cpml.f90

src/main.o: src/parameter.o src/cpml.o $(MODS) src/main2d_cpml.f90
	$(FORTRAN_COMPILER) -o bin/main.o -c $(SWITCH) src/main2d_cpml.f90
	
share/model1.mod: bin/model.o src/model1.f90
	@:
	
share/copt2d_cpml.mod: bin/cpml.o src/copt2d_cpml.f90
	@:


clean:
	rm -f $(OBJS) $(MODS) $(TARGET) *mod
 