bin_PROGRAMS = partest partest-parser

partest_SOURCES = \
	exe/ModelOptimize.cpp \
	exe/ModelSelector.cpp \
	exe/PartitionSelector.cpp \
	indata/PartitionMap.cpp \
	indata/PartitionElement.cpp \
	indata/PartitioningScheme.cpp \
	indata/TreeManager.cpp \
	indata/PllTreeManager.cpp \
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
	search/RandomSearchAlgorithm.cpp \
	util/FileUtilities.cpp \
	util/GlobalDefs.cpp \
	util/PrintMeta.cpp \
	util/Utilities.cpp \
	PartitionTest.cpp



partest_CPPFLAGS = -g -O3 -Wall -DPTHREADS
partest_LDFLAGS = -g -O3 -Wall -DPTHREADS

  partest_CPPFLAGS += -DHAVE_PTHREADS
  partest_LDFLAGS += -DHAVE_PTHREADS

if HAVE_AVX
  partest_CPPFLAGS += -DHAVE_AVX
  partest_LDFLAGS += -DHAVE_AVX
endif

if HAVE_SSE
  partest_CPPFLAGS += -DHAVE_SSE
  partest_LDFLAGS += -DHAVE_SSE
endif

partest_parser_SOURCES = \
	parser/INIReader.cpp \
	partestParserUtils/PartestParserUtils.cpp \
	PartitionTestParser.cpp

partest_parser_CPPFLAGS = -Wall
partest_parser_LDFLAGS = 

pkgincludedir=${includedir}/partest
pkginclude_HEADERS= \
	exe/ModelOptimize.h \
	exe/ModelSelector.h \
	exe/PartitionSelector.h \
	indata/PartitionMap.h \
	indata/TreeManager.h \
	indata/PllTreeManager.h \
	indata/PartitionElement.h \
	indata/PartitioningScheme.h \
	model/Model.h \
	model/NucleicModel.h \
	model/ProteicModel.h \
	model/SelectionModel.h \
	parser/ArgumentParser.h \
	parser/ConfigParser.h \
	parser/INIReader.h \
	search/SearchAlgorithm.h \
	search/ExhaustiveSearchAlgorithm.h \
	search/HierarchicalClusteringSearchAlgorithm.h \
	search/GreedySearchAlgorithm.h \
	search/RandomSearchAlgorithm.h \
	util/Utilities.h \
	util/GlobalDefs.h \
	util/PrintMeta.h \
	util/FileUtilities.h \
	partestParserUtils/PartestParserUtils.h \
	PartitionTestParser.h \
	PartitionTest.h
