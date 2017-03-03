###################################################################
#  Makefile for homework 2
#
#  Daniel R. Reynolds
#  SMU Mathematics
#  Math 4370/6370
#  9 February 2017
###################################################################

# compiler & flags
CC = gcc
CFLAGS = -O2

# makefile targets
all: vec1d_test.exe vec2d_test.exe vec2d_b_test.exe

vec1d_test.exe : vec1d_test.c vec1d.o GramSchmidt1d.o get_time.o
	$(CC) $(CFLAGS) $^ -lm -o $@

vec2d_test.exe : vec2d_test.c vec2d.o GramSchmidt2d.o get_time.o
	$(CC) $(CFLAGS) $^ -lm -o $@

vec2d_b_test.exe : vec2d_b_test.c vec2d_b.o GramSchmidt2d_b.o get_time.o
	$(CC) $(CFLAGS) $^ -lm -o $@

vec1d.o : vec1d.c vec1d.h
	$(CC) $(CFLAGS) -c $< -o $@

vec2d.o : vec2d.c vec2d.h
	$(CC) $(CFLAGS) -c $< -o $@

vec2d_b.o : vec2d_b.c vec2d_b.h
	$(CC) $(CFLAGS) -c $< -o $@

GramSchmidt1d.o : GramSchmidt1d.c vec1d.h
	$(CC) $(CFLAGS) -c $< -o $@

GramSchmidt2d.o : GramSchmidt2d.c vec2d.h
	$(CC) $(CFLAGS) -c $< -o $@

GramSchmidt2d_b.o : GramSchmidt2d_b.c vec2d_b.h
	$(CC) $(CFLAGS) -c $< -o $@

get_time.o : get_time.c get_time.h
	$(CC) $(CFLAGS) -c $< -o $@

clean :
	\rm -f *.o *.out a_data

realclean : clean
	\rm -f *.exe *~


####### End of Makefile #######
