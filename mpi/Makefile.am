CXX = @MPICXX@

if BUILD_MPI 

bin_PROGRAMS = partest-mpi
partest_mpi_SOURCES = \
	$(top_srcdir)/src/exe/ModelOptimize.cpp \
	$(top_srcdir)/src/exe/ModelSelector.cpp \
	$(top_srcdir)/src/exe/PartitionSelector.cpp \
	$(top_srcdir)/src/indata/PartitionMap.cpp \
	$(top_srcdir)/src/indata/PartitionElement.cpp \
	$(top_srcdir)/src/indata/PartitioningScheme.cpp \
	$(top_srcdir)/src/model/Model.cpp \
	$(top_srcdir)/src/model/NucleicModel.cpp \
	$(top_srcdir)/src/model/ProteicModel.cpp \
	$(top_srcdir)/src/model/SelectionModel.cpp \
	$(top_srcdir)/src/parser/ArgumentParser.cpp \
	$(top_srcdir)/src/parser/ConfigParser.cpp \
	$(top_srcdir)/src/parser/INIReader.cpp \
	$(top_srcdir)/src/search/SearchAlgorithm.cpp \
	$(top_srcdir)/src/search/ExhaustiveSearchAlgorithm.cpp \
	$(top_srcdir)/src/search/HierarchicalClusteringSearchAlgorithm.cpp \
	$(top_srcdir)/src/search/GreedySearchAlgorithm.cpp \
	$(top_srcdir)/src/search/RandomSearchAlgorithm.cpp \
	$(top_srcdir)/src/util/Utilities.cpp \
	$(top_srcdir)/src/util/GlobalDefs.cpp \
	$(top_srcdir)/src/util/PrintMeta.cpp \
	$(top_srcdir)/src/PartitionTest.cpp
partest_mpi_CPPFLAGS = -I/usr/local/include/pll -I$(top_srcdir)/src -DHAVE_MPI
partest_mpi_LDFLAGS = -I/usr/local/include/pll -I$(top_srcdir)/src -DHAVE_MPI
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
