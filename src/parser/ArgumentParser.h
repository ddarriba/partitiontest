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
 * @file ArgumentParser.h
 */

#pragma once
#ifndef ARGUMENT_PARSER_HPP
#define ARGUMENT_PARSER_HPP

#include "PartitionTest.h"

#include <iomanip>

#define ARG_DT_NUCLEIC "nt"		/** Argument for nucleotide datatype */
#define ARG_DT_PROTEIC "aa"		/** Argument for amino-acid datatype */
#define ARG_TOPO_MP    "mp"		/** Argument for MP starting topology */
#define ARG_TOPO_FIXED "fixed"	/** Argument for fixed JC/JTT starting topology */
#define ARG_TOPO_USER  "user"	/** Argument for custom starting topology */
#define ARG_IC_AIC     "aic"	/** Argument for using AIC for model selection */
#define ARG_IC_BIC     "bic"	/** Argument for using BIC for model selection */
#define ARG_IC_AICC    "aicc"	/** Argument for using AICc for model selection */
#define ARG_IC_DT      "dt"		/** Argument for using DT for model selection */
#define ARG_AF_PHYLIP_SEQ "phylip-seq"	/** Argument for sequential PHYLIP data format */
#define ARG_AF_PHYLIP_INT "phylip-int"	/** Argument for interleaved PHYLIP data format */
#define ARG_AF_NEXUS      "nexus"		/** Argument for NEXUS data format */
#define ARG_AF_CLUSTAL    "clustal"		/** Argument for CLUSTAL data format */
#define ARG_AF_FASTA      "fasta"		/** Argument for FASTA data format */
#define ARG_AF_DCSE       "dcse"		/** Argument for DCSE data format */
#define ARG_AF_GENBANK    "genbank"		/** Argument for GenBank data format */
#define ARG_SS_ALIGN   "alignment"	/** Argument for NxL sample size */
#define ARG_SS_SHANNON "shannon"	/** Argument for Shannon sample size */
#define ARG_SS_CUSTOM  "custom"		/** Argument for custom sample size */
#define ARG_SEARCH_EXHAUSTIVE "exhaustive"	/** Argument for exhaustive search */
#define ARG_SEARCH_RANDOM "random"			/** Argument for random search */
#define ARG_SEARCH_GREEDY "greedy"			/** Argument for greedy search */
#define ARG_SEARCH_GREEDY_EXT "greedyext"	/** Argument for extended greedy search */
#define ARG_SEARCH_HIERARCHICAL "hcluster"	/** Argument for hierarchical clustering search */
#define ARG_SEARCH_AUTO "auto"	/** Argument for hierarchical clustering search */
#define ARG_OPTIMIZE_BESTMODEL "findmodel"	/** Argument for finding the best model for each partition */
#define ARG_OPTIMIZE_GTR       "gtr"       /** Argument for optimizing GTR among partitions */

enum ArgIndex {
	ARG_NULL, ARG_INPUT_FILE, /** Argument for input data file */
	ARG_INPUT_FORMAT, /** Argument for input data format */
	ARG_USER_TREE, /** Argument for input user tree file */
	ARG_DATA_TYPE, /** Argument for data type (aa/nt) */
	ARG_IC_TYPE, /** Argument for selection criterion */
	ARG_FREQUENCIES, /** Argument for including +F models */
	ARG_INV, /** Argument for including +I models */
	ARG_GAMMA, /** Argument for including +G models */
	ARG_TOPOLOGY, /** Argument for starting topology type */
	ARG_FINAL_TREE, /** Argument for computing final tree */
	ARG_CONFIG_FILE, /** Argument for configuration file name */
	ARG_NUM_PROCS, /** Argument for number of processors */
	ARG_SAMPLE_SIZE, /** Argument for sample size type */
	ARG_SEARCH_ALGORITHM, /** Argument for search algorithm */
	ARG_HCLUSTER_REPS, /** Number of hcluster replicates */
	ARG_OPTIMIZE, /** Argument for search algorithm */
	ARG_OUTPUT, /** Argument for setting the output directory */
	ARG_HELP, /** Argument for show help */
	ARG_CONFIG_HELP, /** Argument for show help about configuration */
	ARG_CONFIG_TEMPLATE, /** Argument for show a configuration template */
	ARG_NON_STOP, /** Search until the end */
	ARG_VERSION, /** Show version */
	ARG_END
};

struct option {
	ArgIndex index; /** Argument index */
	char char_code; /** Short option (single char) */
	const char *long_code; /** Long option */
	bool required_value; /** Whether a value is required or not */
};

namespace partest {

class ArgumentParser {

public:
	ArgumentParser(PartitionTest * ptest);
	~ArgumentParser();

	/**
	 * @brief Gets the next argument.
	 *
	 * @param argc Arguments array.
	 * @param argv Number of arguments.
	 * @param[out] argument Read argument (starting with "-" or "--").
	 * @param[out] value Read value (if any).
	 *
	 * @return The index of the parsed argument.
	 */
	ArgIndex get_opt(int argc, char *argv[], char *argument, char *value);

	void parse(int argc, char *argv[]);

private:
	void init();
	int index; /** Current argument index within the arguments array. */
	int subindex; /** Current artument subindex within an arguments chain. */

	option* arguments; /** Read arguments. */
	PartitionTest * ptest;
};

} /* namespace partest */

#endif
