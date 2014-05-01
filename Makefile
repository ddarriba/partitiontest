SRC_DIR = src/
OBJ_DIR = src/
PLL_DIR = /usr/local/include/pll

CLIBRARIES = -lm -lpthread -lpll-sse3-pthreads # -fopenmp
CPP = g++
CXXFLAGS = -O3 -I$(PLL_DIR) -D_PLL -DHAVE_CONFIG_H -D_GNU_SOURCE -g -std=c++0x \
   -Waddress -Warray-bounds -Wc++11-compat -Wchar-subscripts -Wenum-compare -Wcomment -Wformat -Wmain -Wmaybe-uninitialized -Wmissing-braces -Wnonnull -Wparentheses -Wreturn-type -Wsequence-point -Wsign-compare -Wstrict-aliasing -Wstrict-overflow=1 -Wswitch -Wtrigraphs -Wuninitialized -Wunused-function -Wunused-label -Wunused-value -Wunused-variable -Wvolatile-register-var -DPTHREADS  #-fopenmp 


PARTEST_CXXFILES = \
	exe/ModelOptimize.cpp \
	exe/ModelSelector.cpp \
	exe/PartitionSelector.cpp \
	indata/PartitionMap.cpp \
	indata/PartitionElement.cpp \
	indata/PartitioningScheme.cpp \
	model/Model.cpp \
	model/NucleicModel.cpp \
	model/ProteicModel.cpp \
	model/SelectionModel.cpp \
	parser/ArgumentParser.cpp \
	parser/ConfigParser.cpp \
	parser/INIReader.cpp \
	search/SearchAlgorithm.cpp \
	search/ExhaustiveSearchAlgorithm.cpp \
	search/HierarchicalClusteringSearchAlgorithm.cpp \
	search/GreedySearchAlgorithm.cpp \
	util/Utilities.cpp \
	util/GlobalDefs.cpp \
	util/PrintMeta.cpp \
	PartitionTest.cpp
	
GLOBAL_DEPS = pll.h globalVariables.h

PARTEST_OBJ = $(patsubst %.cpp, $(OBJ_DIR)%.o, $(PARTEST_CXXFILES))

all: partest 
	$(CPP) $(CXXFLAGS) \
	$(PARTEST_OBJ) -o partest-pll $(CLIBRARIES)


clean:
	-rm $(PARTEST_OBJ) 
 
partest: $(PARTEST_OBJ)
	@echo "==> Finished compiling PARTEST files"

$(OBJ_DIR)%.o: $(SRC_DIR)%.cpp
	$(CPP) $(CXXFLAGS) -c $+ -I$(SRC_DIR) -o $@


.PHONY: partest
