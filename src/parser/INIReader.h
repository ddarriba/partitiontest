// Read an INI file into easy-to-access name/value pairs.
// Copyright (c) 2009, Brush Technology
// All rights reserved.
// inih and INIReader are released under the New BSD license (see LICENSE.txt).
// Go to the project home page for more info:
//
// http://code.google.com/p/inih/

#ifndef __INIREADER_H__
#define __INIREADER_H__

#include <map>
#include <vector>
#include <string>

#define MAX_SECTION 50
#define MAX_NAME 50

/* Nonzero to allow multi-line value parsing, in the style of Python's
 ConfigParser. If allowed, ini_parse() will call the handler with the same
 name for each subsequent line parsed. */
#ifndef INI_ALLOW_MULTILINE
#define INI_ALLOW_MULTILINE 1
#endif

/* Nonzero to allow a UTF-8 BOM sequence (0xEF 0xBB 0xBF) at the start of
 the file. See http://code.google.com/p/inih/issues/detail?id=21 */
#ifndef INI_ALLOW_BOM
#define INI_ALLOW_BOM 1
#endif

/* Nonzero to use stack, zero to use heap (malloc/free). */
#ifndef INI_USE_STACK
#define INI_USE_STACK 1
#endif

/* Stop parsing on first error (default is to keep parsing). */
#ifndef INI_STOP_ON_FIRST_ERROR
#define INI_STOP_ON_FIRST_ERROR 0
#endif

/* Maximum line length for any line in INI file. */
#ifndef INI_MAX_LINE
#define INI_MAX_LINE 200
#endif

using namespace std;

// Read an INI file into easy-to-access name/value pairs. (Note that I've gone
// for simplicity here rather than speed, but it should be pretty decent.)
class INIReader {
public:
	// Construct INIReader and parse given filename. See ini.h for more info
	// about the parsing.
	INIReader(string filename);

	// Return the result of ini_parse(), i.e., 0 on success, line number of
	// first error on parse error, or -1 on file open error.
	int ParseError();

	// Get a string value from INI file, returning default_value if not found.
	string Get(string section, string name, string default_value);

	// Get an integer (long) value from INI file, returning default_value if
	// not found or not a valid integer (decimal "1234", "-1234", or hex "0x4d2").
	long GetInteger(string section, string name, long default_value);

	// Get a real (floating point double) value from INI file, returning
	// default_value if not found or not a valid floating point value
	// according to strtod().
	double GetReal(string section, string name, double default_value);

	// Get a boolean value from INI file, returning default_value if not found or if
	// not a valid true/false value. Valid true values are "true", "yes", "on", "1",
	// and valid false values are "false", "no", "off", "0" (not case sensitive).
	bool GetBoolean(string section, string name, bool default_value);

	map<string, string> * getGenes(string section);
	vector<string> * getSchemes(string section);

private:
	int _error;
	map<string, string> _values;
	static string MakeKey(string section, string name);
	static int ValueHandler(void* user, const char* section, const char* name,
			const char* value);

	/* Parse given INI-style file. May have [section]s, name=value pairs
	 (whitespace stripped), and comments starting with ';' (semicolon). Section
	 is "" if name=value pair parsed before any section heading. name:value
	 pairs are also supported as a concession to Python's ConfigParser.

	 For each name=value pair parsed, call handler function with given user
	 pointer as well as section, name, and value (data only valid for duration
	 of handler call). Handler should return nonzero on success, zero on error.

	 Returns 0 on success, line number of first error on parse error (doesn't
	 stop on first error), -1 on file open error, or -2 on memory allocation
	 error (only when INI_USE_STACK is zero).
	 */
	int ini_parse(const char* filename,
			int (*handler)(void* user, const char* section, const char* name,
					const char* value), void* user);

	/* Same as ini_parse(), but takes a FILE* instead of filename. This doesn't
	 close the file when it's finished -- the caller must do that. */
	int ini_parse_file(FILE* file,
			int (*handler)(void* user, const char* section, const char* name,
					const char* value), void* user);
};

#endif  // __INIREADER_H__
