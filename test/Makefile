includeDir = ../include

libDir = ../../lib
libIncDir = ../../include
gbwtDir = ~/GBWTQuery/gbwt/bin


indexHeaders = $(includeDir)/fast_lcp.h $(includeDir)/lf_gbwt.h $(includeDir)/compText.h
utilHeaders = $(includeDir)/testing.h $(includeDir)/ioHelp.h
queryHeaders = $(includeDir)/querySupport.h $(includeDir)/setMaximalMatchQuery.h $(includeDir)/longMatchQuery.h
headers = $(indexHeaders) $(utilHeaders) $(queryHeaders)

LDFLAGS = -L$(libDir)/ 
LDLIBS = -lgbwt -lsdsl 
CXXFLAGS = -std=c++14 -fopenmp -O3 -I$(includeDir) -I$(libIncDir)

dataDir = ../toyData/indexes
toTest = bi_fig1 bi_fig1_repeated128
testSizes = 0 10 100
toRemove = 2
Llength = 4

.PHONY: all clean rebuild test

all: doTest generateIndices verifyStructuresAndQueries

clean: 
	rm doTest generateIndices verifyStructuresAndQueries

rebuild: clean all

test: testverify testGenerateandTest

testverify: verifyStructuresAndQueries
	for base in $(toTest);\
	do\
		for size in $(testSizes);\
		do\
			echo ./verifyStructuresAndQueries $(dataDir)/$$base/$$base.gbwt $(dataDir)/$$base/$$base.ri $$size;\
			./verifyStructuresAndQueries $(dataDir)/$$base/$$base.gbwt $(dataDir)/$$base/$$base.ri $$size;\
		done;\
	done

testGenerateandTest: generateIndices doTest
	for base in $(toTest);\
	do\
		echo ./generateIndices $(dataDir)/$$base/$$base $(toRemove);\
		./generateIndices $(dataDir)/$$base/$$base $(toRemove);\
		echo $(gbwtDir)/build_gbwt -l -R -s 1024 -o $(dataDir)/$$base/$$base\Trimmed_withDA $(dataDir)/$$base/$$base\Trimmed;\
		$(gbwtDir)/build_gbwt -l -R -s 1024 -o $(dataDir)/$$base/$$base\Trimmed_withDA $(dataDir)/$$base/$$base\Trimmed;\
		echo ./doTest $(dataDir)/$$base/$$base\Trimmed $(Llength);\
		./doTest $(dataDir)/$$base/$$base\Trimmed $(Llength);\
	done

doTest: $(headers) doTest.cpp
	$(CXX) $(CXXFLAGS) -o doTest doTest.cpp $(LDFLAGS) $(LDLIBS) 

generateIndices: $(headers) generateIndices.cpp
	$(CXX) $(CXXFLAGS) -o generateIndices generateIndices.cpp $(LDFLAGS) $(LDLIBS) 

verifyStructuresAndQueries: $(headers) verifyStructuresAndQueries.cpp
	$(CXX) $(CXXFLAGS) -o verifyStructuresAndQueries verifyStructuresAndQueries.cpp $(LDFLAGS) $(LDLIBS) 
