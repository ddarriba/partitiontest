/*  PartitionTest, fast selection of the best fit partitioning scheme for
 *  multi-gene data sets.
 *  Copyright May 2013 by Diego Darriba
 *
 *  This program is free software; you may redistribute it and/or modify its
 *  under the terms of the GNU General Public License as published by the Free
 *  Software Foundation; either version 3 of the License, or (at your option)
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 *  or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 *  for more details.
 *
 *  For any other inquiries send an Email to Diego Darriba
 *  ddarriba@udc.es
 */

/**
 * @file GlobalDefs.h
 * @author Diego Darriba
 * @brief Global definitions for Partition Test
 */

#ifndef GLOBALDEFS_H_
#define GLOBALDEFS_H_

#define PROGRAM_VERSION "1.0"
#define PROGRAM_DATE "15 Jul 2014"

#include <pll/pll.h>

#include <vector>
#include <climits>
#include <string>

#ifdef HAVE_MPI
  #include <mpi.h>
  #define I_AM_ROOT myRank==0
#else
  #define I_AM_ROOT 1
#endif


#define MAX_FILE_LENGTH 250

/* ML optimization parameters */
#define AUTO_EPSILON           0.0f
#define ML_PARAM_CUTOFF        PLL_TRUE
#define ML_PARAM_STEPWIDTH     5
#define ML_PARAM_MAXREARRANGE  11
#define ML_PARAM_BESTTRAV      5
#define ML_PARAM_INITIALSET    PLL_FALSE

namespace partest {

typedef unsigned long int bitMask;
#define MAX_PARTITIONS LONG_MAX
typedef std::vector<size_t> t_partitionElementId;
typedef std::vector<t_partitionElementId> t_partitioningScheme;
typedef std::vector<t_partitioningScheme> t_schemesVector;

typedef struct {
	size_t start;
	size_t end;
	size_t id;
} PEsection;

#define DOUBLE_INF 1e140

#define CKP_DIR "ckpfiles"

#ifdef _WIN32
#define char_separator '\\'
const std::string os_separator("\\");
#else
const std::string os_separator("/");
#define char_separator '/'
#endif

#define NUM_DNA_RATES 6
#define NUM_AA_RATES 190
/** Number of states for nucleotide substitution models */
#define NUM_NUC_FREQS 4
/** Number of states for amino-acid replacement models */
#define NUM_PROT_FREQS 20

#define PROT_MATRIX_SIZE 18 // excluding auto & GTR
enum ProtMatrix {
	PROT_MATRIX_DAYHOFF = PLL_DAYHOFF,
	PROT_MATRIX_DCMUT = PLL_DCMUT,
	PROT_MATRIX_JTT = PLL_JTT,
	PROT_MATRIX_MTREV = PLL_MTREV,
	PROT_MATRIX_WAG = PLL_WAG,
	PROT_MATRIX_RTREV = PLL_RTREV,
	PROT_MATRIX_CPREV = PLL_CPREV,
	PROT_MATRIX_VT = PLL_VT,
	PROT_MATRIX_BLOSUM62 = PLL_BLOSUM62,
	PROT_MATRIX_MTMAM = PLL_MTMAM,
	PROT_MATRIX_LG = PLL_LG,
	PROT_MATRIX_MTART = PLL_MTART,
	PROT_MATRIX_MTZOA = PLL_MTZOA,
	PROT_MATRIX_PMB = PLL_PMB,
	PROT_MATRIX_HIVB = PLL_HIVB,
	PROT_MATRIX_HIVW = PLL_HIVW,
	PROT_MATRIX_JTTDCMUT = PLL_JTTDCMUT,
	PROT_MATRIX_FLU = PLL_FLU,
	PROT_MATRIX_AUTO = PLL_AUTO,
	PROT_MATRIX_GTR = PLL_GTR
};

#define NUC_MATRIX_SIZE 22
enum NucMatrix {
	NUC_MATRIX_JC = 0,
	NUC_MATRIX_F81 = 1,
	NUC_MATRIX_K80 = 2,
	NUC_MATRIX_HKY = 3,
	NUC_MATRIX_TrNef = 4,
	NUC_MATRIX_TrN = 5,
	NUC_MATRIX_TPM1 = 6,
	NUC_MATRIX_TPM1uf = 7,
	NUC_MATRIX_TPM2 = 8,
	NUC_MATRIX_TPM2uf = 9,
	NUC_MATRIX_TPM3 = 10,
	NUC_MATRIX_TPM3uf = 11,
	NUC_MATRIX_TIM1ef = 12,
	NUC_MATRIX_TIM1 = 13,
	NUC_MATRIX_TIM2ef = 14,
	NUC_MATRIX_TIM2 = 15,
	NUC_MATRIX_TIM3ef = 16,
	NUC_MATRIX_TIM3 = 17,
	NUC_MATRIX_TVMef = 18,
	NUC_MATRIX_TVM = 19,
	NUC_MATRIX_SYM = 20,
	NUC_MATRIX_GTR = 21
};

#define EX_OK           0       /* successful termination */
#define EX__BASE        64      /* base value for error messages */
#define EX_USAGE        64      /* command line usage error */
#define EX_DATAERR      65      /* data format error */
#define EX_NOINPUT      66      /* cannot open input */
#define EX_NOUSER       67      /* addressee unknown */
#define EX_NOHOST       68      /* host name unknown */
#define EX_UNAVAILABLE  69      /* service unavailable */
#define EX_SOFTWARE     70      /* internal software error */
#define EX_OSERR        71      /* system error (e.g., can't fork) */
#define EX_OSFILE       72      /* critical OS file missing */
#define EX_CANTCREAT    73      /* can't create (user) output file */
#define EX_IOERR        74      /* input/output error */
#define EX_TEMPFAIL     75      /* temp failure; user is invited to retry */
#define EX_PROTOCOL     76      /* remote error in protocol */
#define EX_NOPERM       77      /* permission denied */
#define EX_CONFIG       78      /* configuration error */

enum RateVar {
	RateVarM = 1, RateVarF = 2, RateVarI = 4, RateVarG = 8
};

enum StartTopo {
	StartTopoMP, StartTopoML, StartTopoFIXED, StartTopoUSER
};

enum SearchAlgo {
	SearchExhaustive,
	SearchRandom,
	SearchGreedy,
	SearchGreedyExtended,
	SearchHCluster,
	SearchAuto
};

enum DataType {
	DT_NUCLEIC = 1, DT_PROTEIC = 2
};

enum InformationCriterion {
	AIC, AICC, BIC, DT
};

enum SampleSize {
	SS_ALIGNMENT, SS_SHANNON, SS_CUSTOM
};

enum OptimizeMode {
	OPT_SEARCH, OPT_GTR, OPT_CUSTOM
};

#define CHECKPOINT_LOADED      0
#define CHECKPOINT_SAVED       1
#define CHECKPOINT_UNAVAILABLE 2
#define CHECKPOINT_UNEXISTENT  3

#define DEFAULT_DATA_TYPE DT_NUCLEIC
#define DEFAULT_STARTING_TOPOLOGY StartTopoMP
#define DEFAULT_SEARCH_ALGO SearchHCluster
#define DEFAULT_IC_TYPE BIC
#define DEFAULT_SAMPLE_SIZE SS_ALIGNMENT
#define DEFAULT_OPTIMIZE OPT_SEARCH
#define DEFAULT_DO_F false
#define DEFAULT_DO_I false
#define DEFAULT_DO_G false

/* checkpointing */
extern bool ckpAvailable;
extern std::string ckpPath;
extern std::string ckpStartingTree;
extern std:: string ckpFinalTree;

/* configuration */
/** Number of threads used for optimization */
extern int number_of_threads;
/** Data type (nucleic or amino acid) */
extern DataType data_type;
/** Starting topology for optimizations */
extern StartTopo starting_topology;
/** Algorithm for searching in the partition space */
extern SearchAlgo search_algo;
/** Epsilon for optimization algorithm */
extern double epsilon;
/** Number of samples in configurable algorithms */
extern int max_samples;
/** Whether to stop or not in local maxima */
extern bool non_stop;
/** Whether to thoroughly optimize the final scheme*/
extern bool compute_final_tree;
/** Algorithm for candidate model selection */
extern OptimizeMode optimize_mode;
/** Information criterion for selecting models/partitions */
extern InformationCriterion ic_type;
/** Model rates to optimize */
extern bitMask do_rate;
/** Models to evaluate (proteins only) */
extern bitMask protModels;
/** Determine if branch lengths are optimized for each partition */
extern bool reoptimize_branch_lengths;
/** Determine whether to estimate per-gene branch lengths */
extern bool pergene_branch_lengths;


/* input/output */
extern std::string * input_file;
extern std::string * config_file;
extern std::string * user_tree;
extern std::string * output_dir;
extern std::string * models_logfile;
extern std::string * schemes_logfile;
extern std::string * results_logfile;
extern std::string * log_logfile;
extern bool outputAvailable;
extern bool force_overriding;

#ifdef HAVE_MPI
	extern int myRank;
	extern int numProcs;
#endif

/* data description */
/** Number of taxa in the alignment */
extern size_t num_taxa;
/** Sequence length of the alignment */
extern size_t seq_len;
/** Number of patterns in the alignment */
extern size_t num_patterns;
/** Number of models to evaluate */
extern size_t number_of_models;
/** Number of genes in the input data */
extern size_t number_of_genes;
/**
 * Number of schemes to evaluate if a fixed set
 * of candidate schemes is specified
 */
extern size_t number_of_schemes;
extern std::string ** singleGeneNames;
extern char * starting_tree;
extern char ** pergene_starting_tree;
extern std::vector<t_partitioningScheme> * schemes;

/* data structures */
extern pllQueue * pllPartsQueue;
extern partitionList * pllPartitions;
extern pllAlignmentData * phylip;
extern pllInstance * tree;

extern time_t start_time;

void exit_partest(int status);
std::string timestamp();

}

#endif /* GLOBALDEFS_H_ */
