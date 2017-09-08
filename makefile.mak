OBJS=teradataMain.o utilities.o
CXX=g++ -std=c++0x
LINKS=-I/apps/rhel6/boost/1.60.0_intel-16.0.1.150/include
LINKER=-L/apps/rhel6/boost/1.60.0_intel-16.0.1.150/lib -Xlinker -rpath -Xlinker /apps/rhel6/boost/1.60.0_intel-16.0.1.150/lib -fopenmp -L -lboost_random -lboost_math_c99 -lgsl -lgslcblas -lm
CXXFLAGS=-Wall -c $(LINKS) $(LINKER)
LFLAGS=-Wall -fopenmp -lgsl -lgslcblas -lm
p1: $(OBJS)
        $(CXX) $(LFLAGS) $(OBJS) -o p1
teradataMain.o: teradataMain.cpp utilities.cpp utilities.hpp normal.hpp
        $(CXX) teradataMain.cpp $(CXXFLAGS)
utilties.o: utilities.cpp utilities.hpp normal.hpp
        $(CXX) utilities.cpp $(CXXFLAGS)
clean:
        \rm *.o
