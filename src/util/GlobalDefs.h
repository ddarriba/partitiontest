/*
 * GlobalDefs.h
 *
 *  Created on: Jan 8, 2013
 *      Author: diego
 */

#ifndef GLOBALDEFS_H_
#define GLOBALDEFS_H_

#include <vector>
#include <climits>
#include <string>
#ifdef _PLL
extern "C" {
#include "pll.h"
}
#endif

/* this two definitions below must match */
#define MAX_PARTITIONS LONG_MAX
typedef std::vector<unsigned int> t_partitionElementId;

#define DOUBLE_INF 1e140

#ifdef _WIN32
#define char_separator '\\'
const std::string os_separator("\\");
#else
const std::string os_separator("/");
#define char_separator '/'
#endif

#define NUM_DNA_STATES  4
#define NUM_DNA_RATES   6
#define NUM_AA_STATES  20
#define NUM_AA_RATES  190
#define NUM_RATE_CATEGORIES 4

typedef std::vector<t_partitionElementId> t_partitioningScheme;
typedef std::vector<t_partitioningScheme> t_schemesVector;

#ifdef _PLL
#ifndef AXML_H
#define AXML_H
#include "pll.h"
#endif
#endif

namespace partest {

typedef unsigned int bitMask;

enum RateVar {
  RateVarM = 1,
  RateVarF = 2,
  RateVarI = 4,
  RateVarG = 8
};

#ifndef _PLL
enum StartTopo {
	StartTopoBIONJ,
	StartTopoML,
	StartTopoFIXED,
	StartTopoUSER
};
#else
enum StartTopo {
	StartTopoMP,
	StartTopoFIXED,
	StartTopoUSER
};
#endif

enum SearchAlgo {
	SearchExhaustive,
	SearchRandom,
	SearchGreedy,
	SearchGreedyExtended,
	SearchHCluster
};

#ifdef _PLL
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
#else
#define PROT_MATRIX_SIZE 14
enum ProtMatrix {
  PROT_MATRIX_LG = 0,
  PROT_MATRIX_WAG = 1,
  PROT_MATRIX_JTT = 2,
  PROT_MATRIX_MTREV = 3,
  PROT_MATRIX_DAYHOFF = 4,
  PROT_MATRIX_DCMUT = 5,
  PROT_MATRIX_RTREV = 6,
  PROT_MATRIX_CPREV = 7,
  PROT_MATRIX_VT = 8,
  PROT_MATRIX_BLOSUM62 = 9,
  PROT_MATRIX_MTMAM = 10,
  PROT_MATRIX_MTART = 11,
  PROT_MATRIX_HIVW = 12,
  PROT_MATRIX_HIVB = 13
};
#endif

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

#define DEFAULT_DATA_TYPE DT_NUCLEIC
#ifndef _PLL
#define DEFAULT_STARTING_TOPOLOGY StartTopoML
#else
#define DEFAULT_STARTING_TOPOLOGY StartTopoMP
#endif
#define DEFAULT_SEARCH_ALGO SearchGreedy
#define DEFAULT_IC_TYPE BIC
#define DEFAULT_SAMPLE_SIZE SS_ALIGNMENT
#define DEFAULT_OPTIMIZE OPT_SEARCH
#define DEFAULT_DO_F false
#define DEFAULT_DO_I false
#define DEFAULT_DO_G false

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

#define EX__MAX 78      /* maximum listed value */

enum DataType {
  DT_NUCLEIC=1, DT_PROTEIC=2
};

enum InformationCriterion {
	AIC, AICC, BIC, DT
};

enum SampleSize {
	SS_ALIGNMENT, SS_SHANNON, SS_CUSTOM
};

enum OptimizeMode {
	OPT_SEARCH, OPT_GTR
};

} /* namespace partest */

#endif /* GLOBALDEFS_H_ */
