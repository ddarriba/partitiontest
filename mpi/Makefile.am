CXX = @MPICXX@

if BUILD_MPI 

bin_PROGRAMS = partest-mpi
partest_mpi_SOURCES = \
	../src/exe/ModelOptimize.cpp \
	../src/exe/ModelSelector.cpp \
	../src/exe/PartitionSelector.cpp \
	../src/indata/PartitionMap.cpp \
	../src/indata/PartitionElement.cpp \
	../src/indata/PartitioningScheme.cpp \
	../src/indata/TreeManager.cpp \
	../src/indata/PllTreeManager.cpp \
	../src/model/Model.cpp \
	../src/model/NucleicModel.cpp \
	../src/model/ProteicModel.cpp \
	../src/model/SelectionModel.cpp \
	../src/parser/ArgumentParser.cpp \
	../src/parser/ConfigParser.cpp \
	../src/parser/INIReader.cpp \
	../src/search/SearchAlgorithm.cpp \
	../src/search/ExhaustiveSearchAlgorithm.cpp \
	../src/search/HierarchicalClusteringSearchAlgorithm.cpp \
	../src/search/GreedySearchAlgorithm.cpp \
	../src/search/RandomSearchAlgorithm.cpp \
	../src/util/FileUtilities.cpp \
	../src/util/GlobalDefs.cpp \
	../src/util/PrintMeta.cpp \
	../src/util/Utilities.cpp \
	../src/PartitionTest.cpp
partest_mpi_CPPFLAGS = -I../src -DHAVE_MPI -Wall
partest_mpi_LDFLAGS = -I../src -DHAVE_MPI -Wall
if HAVE_PTHREADS
  partest_mpi_CPPFLAGS += -DHAVE_PTHREADS
  partest_mpi_LDFLAGS += -DHAVE_PTHREADS
endif
if HAVE_AVX
  partest_mpi_CPPFLAGS += -DHAVE_AVX
  partest_mpi_LDFLAGS += -DHAVE_AVX
endif
if HAVE_SSE
  partest_mpi_CPPFLAGS += -DHAVE_SSE
  partest_mpi_LDFLAGS += -DHAVE_SSE
endif

endif
