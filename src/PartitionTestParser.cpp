/*
 * PartitionTestParser.cpp
 *
 *  Created on: Jul 1, 2014
 *      Author: diego
 */

#include "PartitionTestParser.h"

#include <iomanip>
#include <ostream>
#include <fstream>
#include <string.h>
#include <cstdio>
#include <stdio.h>

namespace partest_parser {

using namespace std;

PartitionTestParser::PartitionTestParser(int argc, char *argv[]) {
	partitions = NULL;

	strcpy(binName, argv[0]);
	if (argc < 3 || getFormat(argv[1], &format)) {
		printHelp(cerr);
		exit(1);
	}

	strcpy(inputFile, argv[2]);
	if (!PartestParserUtils::existsFile(inputFile)) {
		cerr << "[ERROR] File " << inputFile << " does not exist." << endl;
		exit(1);
	}

	if (argc > 3) {
		strcpy(outputFile, argv[3]);
	} else {
		strcpy(outputFile, inputFile);
		strcat(outputFile, ".partest.cfg");
	}
	if (PartestParserUtils::existsFile(outputFile)) {
		cerr << "[ERROR] File " << outputFile << " already exists." << endl;
		exit(1);
	}

	cout << endl;
	cout << setw(10) << "" << setw(41) << setfill('-') << "" << setfill(' ')
			<< endl;
	cout << setw(10) << "" << " PartitionTest Configuration File Parser"
			<< endl;
	cout << setw(10) << "" << setw(41) << setfill('-') << "" << setfill(' ')
			<< endl;
	cout << setw(15) << "" << "Input:  " << inputFile << endl;
	cout << setw(15) << "" << "Output: " << outputFile << endl;
	cout << setw(15) << "" << "Format: ";
	switch (format) {
	case cfPartitionFinder:
		cout << "PartitionFinder" << endl;
		break;
	case cfRAxML:
		cout << "RAxML" << endl;
		break;
	}
	cout << endl;
}

PartitionTestParser::~PartitionTestParser() {
	delete partitions;
}

int PartitionTestParser::getFormat(char * str, ConfigFormat * result) {
	int i = 0;
	for (char c = str[i]; c != '\0'; c = str[i++]) {
		if (c < '0' || c > '9') {
			return 1;
		}
	}
	(*result) = (ConfigFormat) atoi(str);
	return 0;
}

void PartitionTestParser::printHelp(ostream & out) {
	out << endl << "Usage: " << binName << " FORMAT INPUT_FILE [OUTPUT_FILE]"
			<< endl << endl;
	out << "    FORMAT         Configuration input file format. (Required)"
			<< endl;
	out << "                   " << cfPartitionFinder << " PartitionFinder"
			<< endl;
	out << "                   " << cfRAxML << " RAxML" << endl;
	out << endl;
	out << "    INPUT_FILE     Input configuration file. (Required)" << endl;
	out << endl;
	out << "    OUTPUT_FILE    Output PartitionTest configuration file."
			<< endl;
	out
			<< "                   If not specified, out will be {INPUT_FILE}.partest.cfg"
			<< endl;
	out << endl;
}

int PartitionTestParser::parseConfigFile() {
	PartestParserUtils parserUtils(inputFile, outputFile);
	switch (format) {
	case cfPartitionFinder:
		parserUtils.parsePartitionFinderFile(&partitionStrings);
		break;
	case cfRAxML: {
		parserUtils.parseRaxmlFile(&partitions);
		break;
	}
	}

	ofstream ofs(outputFile, ios::out);
	ofs << "; PartitionTest configuration file." << endl;
	ofs << "; ---------------------------------." << endl << endl;
	ofs
			<< "; This file has been automatically generated by PartitionTestParser."
			<< endl;
	ofs
			<< "; Inspect and uncomment the execution parameters according to your preferences."
			<< endl;
	ofs << endl;
	ofs << "[input]" << endl;
	ofs << "; msa = input_file.phy" << endl;
	ofs << "; datatype = nt" << endl;
	ofs << "; tree = fixed" << endl << endl;
	ofs << "[search]" << endl;
	ofs << "; algorithm = auto" << endl << endl;
	ofs << "[models]" << endl;
	ofs << "; include = all" << endl << endl;
	ofs << "[partitions]" << endl;
	if (format == cfPartitionFinder) {
		for (size_t i = 0; i < partitionStrings->size(); i++) {
			cout << partitionStrings->at(i) << endl;
			ofs << partitionStrings->at(i) << endl;
		}
	} else {
		for (size_t i = 0; i < partitions->size(); i++) {
			cout << partitions->at(i).name << " : "
					<< partitions->at(i).start[0] << " to "
					<< partitions->at(i).end[0] << "/"
					<< partitions->at(i).stride[0] << endl;
			ofs << partitions->at(i).name << " = ";
			for (int j = 0; j < partitions->at(i).numberOfSections; j++) {
				if (j)
					ofs << ",";
				ofs << partitions->at(i).start[j] << "-"
						<< partitions->at(i).end[j];
				if (partitions->at(i).stride[j] > 1) {
					ofs << "\\" << partitions->at(i).stride[j];
				}
			}
			ofs << endl;
		}
	}
	ofs << endl;
	ofs << "[output]" << endl;
	ofs << "; path = " << outputFile << ".results" << endl << endl;
	ofs.close();

	return 0;
}

} /* namespace partest_parser */

int main(int argc, char * argv[]) {

	partest_parser::PartitionTestParser parser(argc, argv);
	parser.parseConfigFile();

	return 0;
}
