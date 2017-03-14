###################################################################
#  Makefile for Lab 7
#
#  Daniel R. Reynolds
#  SMU Mathematics
#  Math 4370/6370
#  17 February 2015
###################################################################

CXX = g++
CXXFLAGS = -O2
OMPFLAG = -fopenmp


quadrature.exe : quadrature.cpp
	$(CXX) $(CXXFLAGS) $(OMPFLAG) $^ -lm -o $@


clean :
	\rm -f *.o *.out

realclean : clean
	\rm -f *.exe *~


####### End of Makefile #######
