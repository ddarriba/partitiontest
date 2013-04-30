#include "Utilities.h"
#include "parser/ArgumentParser.h"
#include <iomanip>
#include <string.h>
#include <stdlib.h>

namespace partest {

CfgMap Utilities::config;

long Utilities::factorial(unsigned int x) {
	if (x <= 1)
		return 1;
	else
		return x * factorial(x - 1);
}

long Utilities::combinatorial(unsigned int a, unsigned int b) {
	if (a <= 1)
		return 1;

	return factorial(a) / (factorial(b) * factorial(a - b));
}

int Utilities::bell(int n) {
	/* HARDCODED for n <= 15 */
	switch (n) {
	case 0:
		return 1;
	case 1:
		return 1;
	case 2:
		return 2;
	case 3:
		return 5;
	case 4:
		return 15;
	case 5:
		return 52;
	case 6:
		return 203;
	case 7:
		return 877;
	case 8:
		return 4140;
	case 9:
		return 21147;
	case 10:
		return 115975;
	case 11:
		return 678570;
	case 12:
		return 4213597;
	case 13:
		return 27644437;
	case 14:
		return 190899322;
	case 15:
		return 1382958545;
	default:
		return 1382958545;
	}
}

double Utilities::stringToDouble(const string& s) {
	istringstream i(s);
	double x;
	if (!(i >> x))
		return 0;
	return x;
}

int Utilities::stringToInt(const string& s) {
	istringstream i(s);
	int x;
	if (!(i >> x))
		return 0;
	return x;
}

string Utilities::timeToString(time_t time) {
	time_t seconds = time;
	int hours = seconds / 3600;
	seconds -= hours * 3600;
	int minutes = seconds / 60;
	seconds -= minutes * 60;
	stringstream ss;
	ss << setfill('0') << setw(2) << hours << "h:" << setw(2) << minutes << "m:"
			<< setw(2) << seconds << "s";
	return (ss.str());
}

int Utilities::countWords(string strString) {
	int nSpaces = 0;
	unsigned int i = 0;

	// Skip over spaces at the beginning of the word
	while (isspace(strString.at(i)))
		i++;

	for (; i < strString.length(); i++) {
		if (isspace(strString.at(i))) {
			nSpaces++;

			// Skip over duplicate spaces & if a NULL character is found, we're at the end of the string
			while (isspace(strString.at(i + 1)))
				i++;
			if (strString.at(i) == '\0') {
				nSpaces--;
			}
		}
	}

	// The number of words = the number of spaces + 1
	return nSpaces + 1;
}

// generic solution
int Utilities::numDigits(int number) {
	int digits = 0;
	if (number < 0)
		digits = 1;  // remove this line if '-' counts as a digit
	while (number) {
		number /= 10;
		digits++;
	}
	return digits;
}

// variable-precision SWAR algorithm
int Utilities::setbitsCount(unsigned int value) {
	value = value - ((value >> 1) & 0x55555555);
	value = (value & 0x33333333) + ((value >> 2) & 0x33333333);
	return (((value + (value >> 4)) & 0x0F0F0F0F) * 0x01010101) >> 24;
}

bool Utilities::existProperty(const string& property) {
	return !config[property].empty();
}

string Utilities::getValue(const string& property) {
	return config[property];
}

char * Utilities::getTempFilename() {
	bool valid_name = false;
	char * filename = 0;
	while (!valid_name) {
		stringstream chk_filename;
		chk_filename << "/tmp/partest" << rand() << ".tmp";
		FILE *check = fopen(chk_filename.str().c_str(), "r");
		if (!check) {
			filename = strdup(chk_filename.str().c_str());
			valid_name = true;
		} else {
			fclose(check);
		}
	}
	return filename;
}

char * Utilities::getStreamTempFilename() {
	bool valid_name = false;
	char * filename = 0;
	while (!valid_name) {
		stringstream chk_filename;
		chk_filename << "/tmp/partest_stream_" << rand() << ".tmp";
		FILE *check = fopen(chk_filename.str().c_str(), "r");
		if (!check) {
			filename = strdup(chk_filename.str().c_str());
			valid_name = true;
		} else {
			fclose(check);
		}
	}
	return filename;
}

/*!
 * Simple function that copies a file from initialFilePath to outputFilePath
 */
int Utilities::copyFile(string initialFilePath, string outputFilePath) {

	ifstream initialFile(initialFilePath.c_str(), ios::in | ios::binary);
	ofstream outputFile(outputFilePath.c_str(), ios::out | ios::binary);
	//defines the size of the buffer
	initialFile.seekg(0, ios::end);
	long fileSize = initialFile.tellg();
	//Requests the buffer of the predefined size

	//As long as both the input and output files are open...
	if (initialFile.is_open() && outputFile.is_open()) {
		short * buffer = new short[fileSize / 2];
		//Determine the file's size
		//Then starts from the beginning
		initialFile.seekg(0, ios::beg);
		//Then read enough of the file to fill the buffer
		initialFile.read((char*) buffer, fileSize);
		//And then write out all that was read
		outputFile.write((char*) buffer, fileSize);
		delete[] buffer;
	}
	//If there were any problems with the copying process, let the user know
	else if (!outputFile.is_open()) {
		cout << "I couldn't open " << outputFilePath << " for copying!\n";
		return 1;
	} else if (!initialFile.is_open()) {
		cout << "I couldn't open " << initialFilePath << " for copying!\n";
		return 1;
	}

	initialFile.close();
	outputFile.close();

	return 0;
}

void Utilities::exit_partest(int exitValue) {
	switch (exitValue) {
	case EX_USAGE:
		//ArgumentParser::display_usage();
		break;
	}

	exit(exitValue);
}

FILE * Utilities::myfopen(const char *path, const char *mode, bool lazy) {
	FILE *fp = fopen(path, mode);

	if (strcmp(mode, "r") == 0 || strcmp(mode, "rb") == 0) {
		if (fp)
			return fp;
		else {
			if (lazy) {
				return (FILE *) NULL;
			} else {
				printf(
						"\n Error: the file %s you want to open for reading does not exist, exiting ...\n\n",
						path);
				exit(-1);

			}
		}
	} else {
		if (fp)
			return fp;
		else {
			if (lazy) {
				return (FILE *) NULL;
			} else {
				printf(
						"\n Error: the file %s you want to open for writing or appending can not be opened [mode: %s], exiting ...\n\n",
						path, mode);
				exit(-1);
			}
		}
	}
}

int Utilities::myGetline(char **lineptr, int *n, FILE *stream) {
	char *line, *p;
	int size, copy, len;
	int chunkSize = 256 * sizeof(char);

	if (*lineptr == NULL || *n < 2) {
		line = (char *) realloc(*lineptr, chunkSize);
		if (line == NULL)
			return -1;
		*lineptr = line;
		*n = chunkSize;
	}

	line = *lineptr;
	size = *n;

	copy = size;
	p = line;

	while (1) {
		while (--copy > 0) {
			register int c = getc(stream);
			if (c == EOF)
				goto lose;
			else {
				*p++ = c;
				if (c == '\n' || c == '\r')
					goto win;
			}
		}

		/* Need to enlarge the line buffer.  */
		len = p - line;
		size *= 2;
		line = (char *) realloc(line, size);
		if (line == NULL)
			goto lose;
		*lineptr = line;
		*n = size;
		p = line + len;
		copy = size - len;
	}

	lose: if (p == *lineptr)
		return -1;
	/* Return a partial line since we got an error in the middle.  */
	win: *p = '\0';
	return p - *lineptr;
}

} /* namespace partest */
