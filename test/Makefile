includeDir = ../include

libDir = ../../lib
libIncDir = ../../include


indexHeaders = $(includeDir)/fast_lcp.h $(includeDir)/lf_gbwt.h $(includeDir)/compText.h
utilHeaders = $(includeDir)/testing.h $(includeDir)/ioHelp.h
queryHeaders = $(includeDir)/querySupport.h $(includeDir)/setMaximalMatchQuery.h $(includeDir)/longMatchQuery.h
headers = $(indexHeaders) $(utilHeaders) $(queryHeaders)

LDFLAGS = -L$(libDir)/ 
LDLIBS = -lgbwt -lsdsl 
CXXFLAGS = -std=c++14 -fopenmp -O3 -I$(includeDir) -I$(libIncDir)

.PHONY: all clean rebuild

all: doTest generateIndices verifyStructuresAndQueries

clean: 
	rm doTest generateIndices verifyStructuresAndQueries

rebuild: clean all

doTest: $(headers)
	$(CXX) $(CXXFLAGS) -o doTest doTest.cpp $(LDFLAGS) $(LDLIBS) 

generateIndices: $(headers)
	$(CXX) $(CXXFLAGS) -o generateIndices generateIndices.cpp $(LDFLAGS) $(LDLIBS) 

verifyStructuresAndQueries: $(headers)
	$(CXX) $(CXXFLAGS) -o verifyStructuresAndQueries verifyStructuresAndQueries.cpp $(LDFLAGS) $(LDLIBS) 
