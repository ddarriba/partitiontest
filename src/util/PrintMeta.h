/*
 * PrintMeta.hpp
 *
 *  Created on: 01/06/2012
 *      Author: diego
 */

#ifndef PRINTMETA_HPP_
#define PRINTMETA_HPP_

#define OPT_DESCR_LENGTH 45
#define H_RULE_LENGTH 80

#include <iostream>
#include "../options/ParTestOptions.h"

namespace partest {

class PrintMeta {
  public:
    static void print_header(std::ostream& output);
    static void print_options(std::ostream& output, ParTestOptions & options);
    static void print_usage(std::ostream& output);
};

} /* namespace partest */

#endif /* PRINTMETA_HPP_ */
