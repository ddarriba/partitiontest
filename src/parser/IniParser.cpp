/*
 * IniParser.cpp
 *
 *  Created on: Apr 28, 2014
 *      Author: diego
 */

#include "IniParser.h"
#include "util/Utilities.h"

#include <iostream>
#include <fstream>
#include <wctype.h>

#include <string.h>
#include <algorithm>
#include <functional>
#include <cctype>
#include <locale>

namespace partest {

const string & ltrim(std::string &s) {
	s.erase(s.begin(),
			std::find_if(s.begin(), s.end(),
					std::not1(std::ptr_fun<int, int>(std::isspace))));
	return s;
}

bool IniParser::existsSection(const char * section) {
	for (const char * s : sections) {
		if (!strcmp(s,section))
			return true;
	}
	return false;
}

int IniParser::indexOfSection(const char * section) {
	int cur = 0;
	for (const char * s : sections) {
		if (!strcmp(s,section))
			return cur;
		cur++;
	}
	return -1;
}

int IniParser::indexOfKey(const char * section, const char * key) {
	int sectionIndex = indexOfSection(section);
	if (sectionIndex < 0) return sectionIndex;
	int cur = 0;
	for (int i=0; i<keys.at(sectionIndex)->size(); i++) {
		if (!strcmp(keys.at(sectionIndex)->at(i), key)) {
			return cur;
		}
		cur++;
	}
	return -1;
}

vector<const char *> * IniParser::GetAllKeys(const char * section) {
	return keys.at(indexOfSection(section));
}

const char * & IniParser::GetValue(const char * section, const char * key, const char * defaultValue) {
	int sectionIndex = indexOfSection(section);
	int keyIndex = indexOfKey(section, key);
	if (keyIndex < 0) return defaultValue;
	return values.at(sectionIndex)->at(keyIndex);
}

IniParser::IniParser(const char * file) {

	ifstream myfile;
	myfile.open(file);

	int currentSection = -1;
	for (string line; getline(myfile, line);) {
		//line = Utilities::trim(line);
		if (line.length() == 0 || ltrim(line).find_first_of(";") == 0)
			continue;

		if (line.find_first_of("[") != string::npos) {
			int start = line.find_first_of("[");
			int end = line.find_first_of("]");
			string trimedline = line.erase(end, line.length() - end).erase(0,
					start + 1);
			char * sectCSTR = (char *) malloc(trimedline.length() + 1);
			strcpy(sectCSTR, trimedline.c_str());
			for (int i=0; i<trimedline.length(); i++)
				sectCSTR[i] = tolower(sectCSTR[i]);
			sections.push_back(sectCSTR);
			currentSection++;
			keys.push_back(new vector<const char *>());
			values.push_back(new vector<const char *>());
		} else {
			int sepPos = line.find_first_of("=");
			if (sepPos == string::npos) {
				cerr << "ERROR IN CONFIG FILE" << endl;
				exit_partest(EX_IOERR);
			}
			char * keyCSTR = (char *) malloc(sepPos + 1);
			memcpy(keyCSTR, line.c_str(), sepPos);
			char * valCSTR = (char *) malloc(line.length() - sepPos + 1);
			memcpy(valCSTR, line.c_str() + sepPos + 1, line.length() - sepPos);
			for (int i=0; i<sepPos; i++)
				keyCSTR[i] = tolower(keyCSTR[i]);
			keys.at(currentSection)->push_back(keyCSTR);
			for (int i=0; i<line.length()-sepPos; i++)
				valCSTR[i] = tolower(valCSTR[i]);
			values.at(currentSection)->push_back(valCSTR);
		}
	}
}

IniParser::~IniParser() {
	// TODO Auto-generated destructor stub
}

} /* namespace partest */
