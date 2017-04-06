###################################################################
#  Makefile for Lab on debugging MPI programs
#
#  Daniel R. Reynolds
#  SMU Mathematics
#  Math 4370/6370
#  18 March 2015
###################################################################

MPICC = mpicc
CFLAGS = -O2


quadrature.exe : quadrature.c
	$(MPICC) $(CFLAGS) $^ -lm -o $@


clean :
	\rm -f *.o *.out

realclean : clean
	\rm -f *.exe *~


####### End of Makefile #######
