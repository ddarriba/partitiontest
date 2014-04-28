/*
 * IniParser.h
 *
 *  Created on: Apr 28, 2014
 *      Author: diego
 */

#ifndef INIPARSER_H_
#define INIPARSER_H_

#include <vector>

namespace partest {

using namespace std;

class IniParser {
public:
	IniParser(const char * file);
	virtual ~IniParser();
	vector<const char *> getKeys;
	bool existsSection(const char * section);
	int indexOfSection(const char * section);
	int indexOfKey(const char * section, const char * key);
	vector<const char *> * GetAllKeys(const char * section);
	const char * & GetValue(const char * section, const char * key, const char * defaultValue);
private:
	vector<const char *> sections;
	vector<vector<const char *> *> keys;
	vector<vector<const char *> *> values;
};

} /* namespace partest */

#endif /* INIPARSER_H_ */
