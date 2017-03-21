###################################################################
#  Makefile for Lab 7: global minimization
#
#  Daniel R. Reynolds
#  SMU Mathematics
#  Math 6370
#  13 January 2013
###################################################################

# compilers/flags for zeno
CXX = icc
CXXFLAGS = -O2
OMPFLAGS = -openmp


# makefile targets
glob_min.exe : glob_min.cpp
	$(CXX) $(CXXFLAGS) $^ -lm -o $@

glob_min2.exe : glob_min2.cpp
	$(CXX) $(CXXFLAGS) $^ -lm -o $@

glob_min_omp.exe : glob_min_omp.cpp
	$(CXX) $(OMPFLAGS) $(CXXFLAGS) $^ -lm -o $@

all : glob_min.exe glob_min2.exe glob_min_omp.exe 

clean :
	\rm -f *.o *.out

realclean : clean
	\rm -f *.exe *~


####### End of Makefile #######
