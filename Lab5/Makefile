###################################################################
#  Makefile for Lab 5
#
#  Daniel R. Reynolds
#  SMU Mathematics
#  Math 4370/6370
#  10 February 2015
###################################################################

CXX = icpc
CFLAGS = -O2
OMP = -openmp

driver.exe : driver.cpp vectors.cpp Makefile
	$(CXX) $(OMP) $(CFLAGS) driver.cpp vectors.cpp -o $@

clean :
	\rm -f *.o *.out

realclean : clean
	\rm -f *.exe *~


####### End of Makefile #######
