BOOSTDIR=/boost/1.66.0_impi-2017.1.132_intel-17.0.1.132
GSLDIR=/gsl/2.4_gcc-4.8.5


CXX = g++
CXXFLAGS = -std=c++11 -Wall


INCLUDES=\
	-I$(BOOSTDIR)/include\
	-I$(GSLDIR)/include
LINKLIBS=\
	-L$(BOOSTDIR)/lib -Xlinker -rpath -Xlinker $(BOOSTDIR)/lib\
	-L$(GSLDIR)/lib
LIBS=\
	-lboost_math_c99\
	-lboost_random\
	-lgsl\
	-lgslcblas\
	-lm\
	-fopenmp

SRCS=\
	teradataMain.cpp\
	utilities.cpp
OBJ=$(SRCS:.cpp=.o)

print-%: ; @echo $* = $($*)


GeneticDataSimulator: $(OBJ)
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(LINKLIBS) $(LIBS) $(OBJ) -o GeneticDataSimulator 

teradataMain.o: teradataMain.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(LINKLIBS) $(LIBS) -c teradataMain.cpp

utilities.o: utilities.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(LINKLIBS) $(LIBS) -c utilities.cpp

clean:	
	rm -f *.o *~ output

