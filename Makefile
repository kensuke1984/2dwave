# ============================================================================
# Name        : Makefile
# Author      : Kensuke Konishi
# Version     : 0.0.1
# Copyright   : It is Complicated.
# Description : Makefile for 2dwave
# ============================================================================
 
.PHONY: all clean

#DEBUG = -debug -check bounds -warn interface -g
.SUFFIXES: .f90

all:bin src 
	
bin:
	@mkdir bin

src::
	@(cd src; make)
	
clean:
	@(cd src; make clean)
 
