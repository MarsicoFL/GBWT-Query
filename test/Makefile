includeDir = ../include

libDir = ../../lib
libIncDir = ../../include
gbwtDir = ../../gbwt/bin
SDSL_DIR= ../../sdsl-lite

indexHeaders = $(includeDir)/fast_lcp.h $(includeDir)/lf_gbwt.h $(includeDir)/compText.h
utilHeaders = $(includeDir)/testing.h $(includeDir)/ioHelp.h
queryHeaders = $(includeDir)/querySupport.h $(includeDir)/setMaximalMatchQuery.h $(includeDir)/longMatchQuery.h
headers = $(indexHeaders) $(utilHeaders) $(queryHeaders)

WARNINGFLAGS = #-Wall -Weffc++ -Wextra -Wconversion -Wsign-conversion
PARALLELFLAGS = -fopenmp -pthread

#####################################################################################################
# Apple Clang support modified from GBWT Makefile, at https://github.com/jltsiren/gbwt
include $(SDSL_DIR)/Make.helper

libDir = $(LIB_DIR)
libIncDir = $(INC_DIR)

# Apple Clang does not support OpenMP directly, so we need special handling.
ifeq ($(shell uname -s), Darwin)
    # The compiler complains about -fopenmp instead of missing input.
    ifeq ($(strip $(shell $(MY_CXX) -fopenmp /dev/null -o/dev/null 2>&1 | grep fopenmp | wc -l)), 1)
        $(info The compiler is Apple Clang that needs libomp for OpenMP support.)

        # The compiler only needs to do the preprocessing.
        PARALLELFLAGS = -Xpreprocessor -fopenmp -pthread

        # Find libomp installed by Homebrew or MacPorts.
        ifeq ($(shell if [ -e $(HOMEBREW_PREFIX)/include/omp.h ]; then echo 1; else echo 0; fi), 1)
            $(info Found libomp installed by Homebrew and linked to $(HOMEBREW_PREFIX).)
            PARALLELFLAGS += -I$(HOMEBREW_PREFIX)/include
            LDFLAGS += -L$(HOMEBREW_PREFIX)/lib
        else ifeq ($(shell if [ -d $(HOMEBREW_PREFIX)/opt/libomp/include ]; then echo 1; else echo 0; fi), 1)
            $(info Found a keg-only libomp installed by Homebrew at $(HOMEBREW_PREFIX)/opt/libomp.)
            PARALLELFLAGS += -I$(HOMEBREW_PREFIX)/opt/libomp/include
            LDFLAGS += -L$(HOMEBREW_PREFIX)/opt/libomp/lib
        else ifeq ($(shell if [ -d /opt/local/lib/libomp ]; then echo 1; else echo 0; fi), 1)
            $(info Found libomp installed by MacPorts at /opt/local.)
            PARALLELFLAGS += -I/opt/local/include/libomp
            LDFLAGS += -L/opt/local/lib/libomp
        else
            $(error Could not find libomp. Please install it using Homebrew or MacPorts.)
        endif

        # We also need to link it.
        LDLIBS += -lomp
    endif
endif
#####################################################################################################

LDFLAGS += -L$(libDir)/ 
LDLIBS += -lgbwt -lsdsl 

CXXFLAGS = -DNDEBUG -std=c++14 $(PARALLELFLAGS) -O3 -I$(includeDir) -I$(libIncDir) $(WARNINGFLAGS)


dataDir = ../toyData/indexes
toTest = bi_fig1 bi_fig1_repeated128
testSizes = 0 10 100
toRemove = 2
Llength = 4

1000G_base = ../toyData/indexes/1000G/chr21
1000G_haplotypes_remove = 100
1000G_L_lengths = 100 500 1000 5000 7500 10000 20000 200000

.PHONY: all clean rebuild test

all: doTest generateIndices verifyStructuresAndQueries

clean: 
	rm doTest generateIndices verifyStructuresAndQueries

rebuild: clean all

test: testverify testGenerateandTest



1000GGenerate: generateIndices
	echo ./generateIndices $(1000G_base) $(1000G_haplotypes_remove)
	./generateIndices $(1000G_base) $(1000G_haplotypes_remove)
	echo $(gbwtDir)/build_gbwt -l -R -s 1024 -o $(1000G_base)Trimmed_withDA $(1000G_base)Trimmed
	$(gbwtDir)/build_gbwt -l -R -s 1024 -o $(1000G_base)Trimmed_withDA $(1000G_base)Trimmed

1000GTest: doTest
	for matchLLength in $(1000G_L_lengths);\
	do\
		echo ./doTest $(1000G_base)Trimmed $$matchLLength;\
		./doTest $(1000G_base)Trimmed $$matchLLength;\
	done




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
