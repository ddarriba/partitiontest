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
 * @file PartitionTestParser.h
 * @author Diego Darriba
 * @brief Main class of the PartitionTest parser for converting input configuration file formats
 */

#ifndef PARTITIONTESTPARSER_H_
#define PARTITIONTESTPARSER_H_

#include "partestParserUtils/PartestParserUtils.h"

#include <string>
#include <vector>
#include <iostream>

namespace partest_parser
{

  enum ConfigFormat
  {
    cfPartitionFinder, /** PartitionFinder format */
    cfRAxML /** RAxML format */
  };

  class PartitionTestParser
  {
  public:
    PartitionTestParser (int argc, char *argv[]);
    virtual ~PartitionTestParser ();

    int parseConfigFile ();
  private:
    void printHelp (std::ostream & out);
    int getFormat (char * str, ConfigFormat * result);

    std::vector<std::string> * partitionStrings; /** Vector of partitions */
    char inputFile[256];
    char outputFile[256];
    char binName[50];
    ConfigFormat format;
  };

} /* namespace partest */

#endif /* PARTITIONTESTPARSER_H_ */
