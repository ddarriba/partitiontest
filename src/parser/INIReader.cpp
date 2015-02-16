/* Read an INI file into easy-to-access name/value pairs.
 * Copyright (c) 2009, Brush Technology
 * All rights reserved.
 */

#include <algorithm>
#include <iostream>
#include <cctype>
#include <cstdlib>
#include <cstdio>
#include <ctype.h>
#include <cstring>
#include <sstream>
#include <iomanip>

#include "INIReader.h"

using namespace std;

INIReader::INIReader(string filename) {
	curgene = 0;
	_genes = 0;
	_error = ini_parse(filename.c_str(), ValueHandler, this);
}

INIReader::~INIReader() {
//	if (_genes) {
//		delete _genes;
//	}
}

int INIReader::ParseError() {
	return _error;
}

string INIReader::Get(string section, string name, string default_value) {
	string key = MakeKey(section, name);
	return _values.count(key) ? _values[key] : default_value;
}

long INIReader::GetInteger(string section, string name, long default_value) {
	string valstr = Get(section, name, "");
	const char* value = valstr.c_str();
	char* end;
	// This parses "1234" (decimal) and also "0x4D2" (hex)
	long n = strtol(value, &end, 0);
	return end > value ? n : default_value;
}

double INIReader::GetReal(string section, string name, double default_value) {
	string valstr = Get(section, name, "");
	const char* value = valstr.c_str();
	char* end;
	double n = strtod(value, &end);
	return end > value ? n : default_value;
}

bool INIReader::GetBoolean(string section, string name, bool default_value) {
	string valstr = Get(section, name, "");
	// Convert to lower case to make string comparisons case-insensitive
	transform(valstr.begin(), valstr.end(), valstr.begin(), ::tolower);
	if (valstr == "true" || valstr == "yes" || valstr == "on" || valstr == "1")
		return true;
	else if (valstr == "false" || valstr == "no" || valstr == "off"
			|| valstr == "0")
		return false;
	else
		return default_value;
}

vector<string> * INIReader::getGenes(string section) {

	if (!_genes) {
		map<string, string>::iterator iter;
		int ngenes = 0;
		for (iter = _values.begin(); iter != _values.end(); iter++) {
			if (!(iter->first.substr(0, section.length() + 1).compare(
					section + "."))) {
				ngenes++;
			}
		}
		_genes = new vector<string>((size_t) ngenes * 2);
		size_t _curgene = 0;
		for (iter = _values.begin(); iter != _values.end(); iter++) {
			if (!(iter->first.substr(0, section.length() + 1).compare(
					section + "."))) {
				_genes->at(_curgene) = iter->first.substr((section.length() + 6),
						iter->first.length() - (section.length() + 6));
				_genes->at(_curgene + 1) = iter->second;
				_curgene += 2;
			}
		}
	}
	return _genes;
}

vector<string> * INIReader::getSchemes(string section) {

	vector<string> * schemeLines = new vector<string>();

	map<string, string>::iterator iter;
	for (iter = _values.begin(); iter != _values.end(); iter++) {
		if (!(iter->first.substr(0, (section.length() + 1)).compare(section + "."))) {
			schemeLines->push_back(iter->second);
		}
	}

	return schemeLines;
}

string INIReader::MakeKey(string section, string name) {
	string key = section + "." + name;
	// Convert to lower case to make section/name lookups case-insensitive
	transform(key.begin(), key.end(), key.begin(), ::tolower);
	return key;
}

int INIReader::ValueHandler(void* user, const char* section, const char* name,
		const char* value) {
	INIReader* reader = (INIReader*) user;
	string key = MakeKey(section, name);
	if (!(key.substr(0, 11).compare("partitions."))) {
		/* PartitionTest gene definition */
		/* A numerical value is added as prefix to the gene name for keeping the original order */
		stringstream ss;
		ss << "partitions.";
		ss << setw(5) << setfill('0') << reader->curgene++ << setfill(' ');
		ss << key.substr(11, key.length() - 11);
		reader->_values[ss.str()] += value;
	} else if (!(key.substr(0, 12).compare("data_blocks."))) {
		/* PartitionFinder gene definition */
		/* A numerical value is added as prefix to the gene name for keeping the original order */
		stringstream ss;
		ss << "data_blocks.";
		ss << setw(5) << setfill('0') << reader->curgene++ << setfill(' ');
		ss << key.substr(12, key.length() - 12);
		reader->_values[ss.str()] += value;
	} else {
		reader->_values[key] += value;
	}

	return 1;
}

/* Strip whitespace chars off end of given string, in place. Return s. */
static char* rstrip(char* s) {
	char* p = s + strlen(s);
	while (p > s && isspace((unsigned char) (*--p)))
		*p = '\0';
	return s;
}

/* Return pointer to first non-whitespace char in given string. */
static char* lskip(const char* s) {
	while (*s && isspace((unsigned char) (*s)))
		s++;
	return (char*) s;
}

/* Return pointer to first char c or ';' comment in given string, or pointer to
 null at end of string if neither found. ';' must be prefixed by a whitespace
 character to register as a comment. */
static char* find_char_or_comment(const char* s, char c) {
	int was_whitespace = 0;
	while (*s && *s != c && !(was_whitespace && *s == ';')) {
		was_whitespace = isspace((unsigned char) (*s));
		s++;
	}
	return (char*) s;
}

/* Version of strncpy that ensures dest (size bytes) is null-terminated. */
static char* strncpy0(char* dest, const char* src, size_t size) {
	strncpy(dest, src, size);
	dest[size - 1] = '\0';
	return dest;
}

/* See documentation in header file. */
int INIReader::ini_parse_file(FILE* file,
		int (*handler)(void*, const char*, const char*, const char*),
		void* user) {
	/* Uses a fair bit of stack (use heap instead if you need to) */
#if INI_USE_STACK
	char line[INI_MAX_LINE];
#else
	char* line;
#endif
	char section[MAX_SECTION] = "";
	char prev_name[MAX_NAME] = "";

	char* start;
	char* end;
	char* name;
	char* value;
	int lineno = 0;
	int error = 0;

#if !INI_USE_STACK
	line = (char*)malloc(INI_MAX_LINE);
	if (!line) {
		return -2;
	}
#endif

	/* Scan through file line by line */
	while (fgets(line, INI_MAX_LINE, file) != NULL) {
		lineno++;
		start = line;
#if INI_ALLOW_BOM
		if (lineno == 1 && (unsigned char) start[0] == 0xEF
				&& (unsigned char) start[1] == 0xBB
				&& (unsigned char) start[2] == 0xBF) {
			start += 3;
		}
#endif
		start = lskip(rstrip(start));

		if (*start == ';' || *start == '#') {
			/* Per Python ConfigParser, allow '#' comments at start of line */
		}
#if INI_ALLOW_MULTILINE
		else if (*prev_name && *start && start > line) {
			/* Non-black line with leading whitespace, treat as continuation
			 of previous name's value (as per Python ConfigParser). */
			if (!handler(user, section, prev_name, start) && !error)
				error = lineno;
		}
#endif
		else if (*start == '[') {
			/* A "[section]" line */
			end = find_char_or_comment(start + 1, ']');
			if (*end == ']') {
				*end = '\0';
				strncpy0(section, start + 1, sizeof(section));
				*prev_name = '\0';
			} else if (!error) {
				/* No ']' found on section line */
				error = lineno;
			}
		} else if (*start && *start != ';') {
			/* Not a comment, must be a name[=:]value pair */
			end = find_char_or_comment(start, '=');
			if (!*end) {
#if INI_ALLOW_MULTILINE
				if (!handler(user, section, prev_name, start) && !error)
					error = lineno;
#else
				cerr << "[ERROR] Line " << lineno
						<< "in the configuration file is too large for the current buffer" << endl;
				exit_partest(EX_IOERR);
#endif
			} else {
				if (*end != '=') {
					end = find_char_or_comment(start, ':');
				}
				if (*end == '=' || *end == ':') {
					*end = '\0';
					name = rstrip(start);
					value = lskip(end + 1);
					end = find_char_or_comment(value, '\0');
					if (*end == ';')
						*end = '\0';
					rstrip(value);

					/* Valid name[=:]value pair found, call handler */
					strncpy0(prev_name, name, sizeof(prev_name));
					if (!handler(user, section, name, value) && !error)
						error = lineno;
				} else if (!error) {
					/* No '=' or ':' found on name[=:]value line */
					error = lineno;
				}
			}
		}

#if INI_STOP_ON_FIRST_ERROR
		if (error)
		break;
#endif
	}

#if !INI_USE_STACK
	free(line);
#endif
	return error;
}

/* See documentation in header file. */
int INIReader::ini_parse(const char* filename,
		int (*handler)(void*, const char*, const char*, const char*),
		void* user) {
	FILE* file;
	int error;

	file = fopen(filename, "r");
	if (!file)
		return -1;
	error = ini_parse_file(file, handler, user);
	fclose(file);
	return error;
}
