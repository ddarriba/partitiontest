#pragma once
#ifndef ARGUMENT_PARSER_HPP
#define ARGUMENT_PARSER_HPP

#include "../options/ParTestOptions.h"
#include <iomanip>

#define ARG_DT_NUCLEIC "nt"
#define ARG_DT_PROTEIC "aa"
#define ARG_TOPO_BIONJ "bionj"
#define ARG_TOPO_ML    "ml"
#define ARG_TOPO_FIXED "fixed"
#define ARG_TOPO_USER  "user"
#define ARG_IC_AIC     "aic"
#define ARG_IC_BIC     "bic"
#define ARG_IC_AICC    "aicc"
#define ARG_IC_DT      "dt"
#define ARG_AF_PHYLIP_SEQ "phylip-seq"
#define ARG_AF_PHYLIP_INT "phylip-int"
#define ARG_AF_NEXUS      "nexus"
#define ARG_AF_CLUSTAL    "clustal"
#define ARG_AF_FASTA      "fasta"
#define ARG_AF_DCSE       "dcse"
#define ARG_AF_GENBANK    "genbank"
#define ARG_SS_ALIGN   "alignment"
#define ARG_SS_SHANNON "shannon"
#define ARG_SS_CUSTOM  "custom"
#define ARG_SEARCH_EXHAUSTIVE "exhaustive"
#define ARG_SEARCH_RANDOM "random"
#define ARG_SEARCH_GREEDY "greedy"
#define ARG_SEARCH_HIERARCHICAL "hcluster"

enum ArgIndex {
  ARG_NULL,
  ARG_INPUT_FILE,
  ARG_INPUT_FORMAT,
  ARG_USER_TREE,
  ARG_DATA_TYPE,
  ARG_IC_TYPE,
  ARG_FREQUENCIES,
  ARG_INV,
  ARG_GAMMA,
  ARG_TOPOLOGY,
  ARG_CONFIG_FILE,
  ARG_NUM_PROCS,
  ARG_SAMPLE_SIZE,
  ARG_SEARCH_ALGORITHM,
  ARG_HELP,
  ARG_CONFIG_HELP,
  ARG_CONFIG_TEMPLATE,
  ARG_END
};

struct option {
    ArgIndex index;
    char char_code;
    const char *long_code;
    bool required_value;
};

namespace partest {

class ArgumentParser {

  public:
    ArgumentParser();
    ~ArgumentParser();
    ArgIndex get_opt(int argc, char *argv[], char *argument, char *value);
    void fill_options(int argc, char *argv[], ParTestOptions * options);
  private:
    void init();
    int index;
    int subindex;
    option* arguments;
};

} /* namespace partest */

#endif
