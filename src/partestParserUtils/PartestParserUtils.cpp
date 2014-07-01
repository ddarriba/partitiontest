/*
 * PartestParserUtils.cpp
 *
 *  Created on: Jul 1, 2014
 *      Author: diego
 */

#include "PartestParserUtils.h"

#include "parser/INIReader.h"

#include <algorithm>
#include <iostream>
#include <string.h>
#include <cstdio>

namespace partest_parser {

using namespace std;

PartestParserUtils::PartestParserUtils(char * inputFile, char * outputFile) :
		inputFile(inputFile), outputFile(outputFile) {
}

PartestParserUtils::~PartestParserUtils() {
}

int PartestParserUtils::parsePartitionFinderFile(
		vector<string> ** partitionInfos) {

	INIReader ini(inputFile);

	/** PARTITIONS **/
	std::map<std::string, string> * keys = ini.getGenes("data_blocks");
	(*partitionInfos) = new vector<string>();

	int partitionId = 0;
	std::map<std::string, std::string>::iterator iter;
	for (iter = keys->begin(); iter != keys->end(); iter++) {
		string partitionLine(iter->first);
		replace( iter->second.begin(), iter->second.end(), '\\', '/');
		replace( iter->second.begin(), iter->second.end(), ';', '\0');
		partitionLine.append(" = ");
		partitionLine.append(iter->second.begin(), iter->second.end());
		(*partitionInfos)->push_back(partitionLine);
		partitionId++;
	}
	delete keys;

	return 0;
}

int PartestParserUtils::parseRaxmlFile(
		vector<partitionInfo> ** partitionInfos) {
	FILE *f;
	int numberOfModels = 0;
	int nbytes = 0;
	char *ch;
	char *cc = (char *) NULL;
	char **p_names;
	int n, i, l;
	int lower, upper, modulo;
	char buf[256];
	int pairsCount;

	f = fopen(inputFile, "rb");

	while (myGetline(&cc, &nbytes, f) > -1) {
		if (!lineContainsOnlyWhiteChars(cc)) {
			numberOfModels++;
		}
		if (cc)
			free(cc);
		cc = (char *) NULL;
	}

	rewind(f);

	p_names = (char **) malloc(sizeof(char *) * numberOfModels);
	(*partitionInfos) = new vector<partitionInfo>(numberOfModels);

	i = 0;
	while (myGetline(&cc, &nbytes, f) > -1) {
		if (!lineContainsOnlyWhiteChars(cc)) {
			n = strlen(cc);
			p_names[i] = (char *) malloc(sizeof(char) * (n + 1));
			strcpy(&(p_names[i][0]), cc);
			i++;
		}
		if (cc)
			free(cc);
		cc = (char *) NULL;
	}

	for (i = 0; i < numberOfModels; i++) {

		ch = p_names[i];
		pairsCount = 0;
		skipWhites(&ch);

		if (*ch == '=') {
			printf("Identifier missing prior to '=' in %s\n", p_names[i]);
			return (-1);
		}

		char * partitionName;
		analyzeIdentifier(&ch, i, &partitionName);
		(*partitionInfos)->at(i).name = string(partitionName);
		ch++;

		numberPairs: pairsCount++;

		(*partitionInfos)->at(i).numberOfSections = pairsCount;

		skipWhites(&ch);

		if (!isNum(*ch)) {
			printf("%c Number expected in %s\n", *ch, p_names[i]);
			return (-1);
		}

		l = 0;
		while (isNum(*ch)) {
			/*printf("%c", *ch);*/
			buf[l] = *ch;
			ch++;
			l++;
		}
		buf[l] = '\0';
		lower = atoi(buf);
		(*partitionInfos)->at(i).start[pairsCount - 1] = lower;

		skipWhites(&ch);

		/* NEW */

		if ((*ch != '-') && (*ch != ',')) {
			if (*ch == '\0' || *ch == '\n' || *ch == '\r') {
				upper = lower;
				goto SINGLE_NUMBER;
			} else {
				printf("'-' or ',' expected in %s\n", p_names[i]);
				return (-1);
			}
		}

		if (*ch == ',') {
			upper = lower;
			goto SINGLE_NUMBER;
		}

		/* END NEW */

		ch++;

		skipWhites(&ch);

		if (!isNum(*ch)) {
			printf("%c Number expected in %s\n", *ch, p_names[i]);
			return (-1);
		}

		l = 0;
		while (isNum(*ch)) {
			buf[l] = *ch;
			ch++;
			l++;
		}
		buf[l] = '\0';
		upper = atoi(buf);
		SINGLE_NUMBER: (*partitionInfos)->at(i).end[pairsCount - 1] = upper;

		if (upper < lower) {
			printf(
					"Upper bound %d smaller than lower bound %d for this partition: %s\n",
					upper, lower, p_names[i]);
			return (-1);
		}

		skipWhites(&ch);

		if (*ch == '\0' || *ch == '\n' || *ch == '\r') /* PC-LINEBREAK*/
		{
			continue;
		}

		if (*ch == ',') {
			ch++;
			goto numberPairs;
		}

		if (*ch == '\\') {
			ch++;
			skipWhites(&ch);

			if (!isNum(*ch)) {
				printf("%c Number expected in %s\n", *ch, p_names[i]);
				return (-1);
			}

			l = 0;
			while (isNum(*ch)) {
				buf[l] = *ch;
				ch++;
				l++;
			}
			buf[l] = '\0';
			modulo = atoi(buf);
			(*partitionInfos)->at(i).stride[pairsCount - 1] = modulo;

			skipWhites(&ch);
			if (*ch == '\0' || *ch == '\n' || *ch == '\r') {
				continue;
			}
			if (*ch == ',') {
				ch++;
				goto numberPairs;
			}
		}

		if (*ch == '/') {
			printf(
					"\nRAxML detected the character \"/\" in your partition file.\n");
			printf(
					"Did you mean to write something similar to this: \"DNA, p1=1-100\\3\" ?\n");
			printf(
					"It's actually a backslash, not a slash, the program will exit now with an error!\n\n");
		} else {
			printf(
					"\nRAxML detected the character \"%c\" in your partition file,\n",
					*ch);
			printf("while it does not belong there!\n");
			printf(
					"\nAre you sure that your partition file complies with the RAxML partition file format?\n");
			printf(
					"\nActually reading the manual, does indeed do help a lot\n\n");
			printf("The program will exit now with an error!\n\n");
		}

		printf(
				"The problematic line in your partition file is this one here:\n\n");

		printf("%s\n\n", p_names[i]);

		return (-1);
	}

	printf("    ...done RAxML parsing\n");

	fclose(f);

	return (0);
}

int PartestParserUtils::lineContainsOnlyWhiteChars(char * str) {
	if (!strlen(str))
		return true;
	int i = 0;
	for (char c = str[i]; c != '\0' && c != '\n'; c = str[i++]) {
		if (c != ' ') {
			return false;
		}
	}
	return true;
}

void PartestParserUtils::skipWhites(char **ch) {
	while (**ch == ' ' || **ch == '\t')
		*ch = *ch + 1;
}

int PartestParserUtils::myGetline(char **lineptr, int *n, FILE *stream) {
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

/**
 * This function has modifications from original RAxML's. It is used to
 * bypass the partition model and get the partition name.
 */
void PartestParserUtils::analyzeIdentifier(char **ch, int modelNumber,
		char **modelName) {
	char ident[2048] = "", model[2048] = "";

	int i = 0, n, j, containsComma = 0;

	while (**ch != '=') {
		if (**ch != ' ' && **ch != '\t') {
			ident[i] = **ch;
			i++;
		}
		*ch = *ch + 1;
	}

	n = i;
	i = 0;

	for (i = 0; i < n; i++)
		if (ident[i] == ',')
			containsComma = 1;

	if (!containsComma) {
		printf(
				"Error, model file must have format: DNA or AA model, then a comma, and then the partition name\n");
		exit(-1);
	} else {
		int openBracket = 0, closeBracket = 0, openPos = 0, closePos = 0;

		i = 0;

		while (ident[i] != ',') {
			if (ident[i] == '[') {
				openPos = i;
				openBracket++;
			}
			if (ident[i] == ']') {
				closePos = i;
				closeBracket++;
			}
			model[i] = ident[i];
			i++;
		}

		if (closeBracket > 0 || openBracket > 0) {
			if (!((closeBracket == 1) && (openBracket == 1)
					&& (openPos < closePos))) {
				printf(
						"\nError: Apparently you want to specify a user-defined protein substitution model that shall be read from file\n");
				printf(
						"It must be enclosed in opening and closing bracktes like this: [fileName]\n\n");
				printf("you specified: %s\n\n", model);
				exit(-1);
			}
		}

		i = 0;
		while (ident[i++] != ',')
			;

		(*modelName) = (char*) malloc((n - i + 1) * sizeof(char));

		j = 0;
		while (i < n)
			(*modelName)[j++] = ident[i++];

		(*modelName)[j] = '\0';
	}
}

} /* namespace partest */
