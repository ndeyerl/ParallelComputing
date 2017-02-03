#implementation: 
#make //this builds the library libgaussq.a containing the objects
#make quadjac
#
# makefile for the pierreQP package to calculate Gauss quadrature points
#

# compile flags
#CFLAGS = -g
CFLAGS = -O2 

CC = gcc
AR = ar
RANLIB = ranlib

OBJS =  GaussCore.o Hermite.o Jacobi.o Laguerre.o

libgaussq.a: $(OBJS)
	$(AR) rv libgaussq.a $(OBJS)
	$(RANLIB) libgaussq.a

quadjac:  quadjac.o panelIA.o 
	gcc $(CFLAGS) -o quadjac quadjac.o panelIA.o libgaussq.a  -llapack -lm




.c.o:
	gcc -c $(CFLAGS) $<

