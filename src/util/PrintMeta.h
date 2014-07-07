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

/*
 * @file PrintMeta.hpp
 * @author Diego Darriba
 * @brief Set of utilities for printing meta-data and results
 */

#ifndef PRINTMETA_HPP_
#define PRINTMETA_HPP_

#include "PartitionTest.h"
#include "indata/PartitioningScheme.h"

#include <iostream>

#define OPT_DESCR_LENGTH 45
#define H_RULE_LENGTH 80

namespace partest {

class PrintMeta {
  public:
    static void print_header(std::ostream& output);
    static void print_options(std::ostream& output);
    static void print_usage(std::ostream& output);
    static void print_results(ostream & ofs, PartitioningScheme * bestScheme);
    static void print_results_xml(ostream & ofs, PartitioningScheme * bestScheme);
};

} /* namespace partest */

#endif /* PRINTMETA_HPP_ */
