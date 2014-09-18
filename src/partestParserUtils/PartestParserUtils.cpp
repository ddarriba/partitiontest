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
 * @file PartestParserUtils.cpp
 * @author Diego Darriba
 * @brief
 */

#include "PartestParserUtils.h"

#include "parser/INIReader.h"

#include <algorithm>
#include <iostream>
#include <string.h>
#include <sstream>
#include <cstdio>

namespace partest_parser {

using namespace std;

PartestParserUtils::PartestParserUtils(char * inputFile, char * outputFile) :
		inputFile(inputFile), outputFile(outputFile) {
}

PartestParserUtils::~PartestParserUtils() {
}

int PartestParserUtils::parsePartitionFinderFile(vector<string> ** partitions,
		char ** alignment, char ** models, pfSearch * searchAlgo,
		selection * icSelection) {

	INIReader ini(inputFile);

	const char * value;
	value = ini.Get("", "alignment", "").c_str();
	if (strcmp(value, "")) {
		(*alignment) = (char *) malloc(strlen(value) + 1);
		strcpy((*alignment), value);
	}
	value = ini.Get("", "models", "").c_str();
	if (strcmp(value, "")) {
		(*models) = (char *) malloc(strlen(value) + 1);
		strcpy((*models), value);
	}
	value = ini.Get("", "model_selection", "").c_str();
	if (!strcmp(value, "BIC")) {
		*icSelection = icBIC;
	} else if (!strcmp(value, "AIC")) {
		*icSelection = icAIC;
	} else if (!strcmp(value, "AICc")) {
		*icSelection = icAICc;
	} else {
		*icSelection = icUNDEFINED;
	}

	value = ini.Get("schemes", "search", "").c_str();
	if (!strcmp(value, "greedy")) {
		*searchAlgo = searchGreedy;
	} else if (!strcmp(value, "hcluster")) {
		*searchAlgo = searchHCluster;
	} else if (!strcmp(value, "rcluster")) {
		*searchAlgo = searchRCluster;
	} else {
		*searchAlgo = searchUNDEFINED;
	}

	/** PARTITIONS **/
	vector<string> * keys = ini.getGenes("data_blocks");
	int num_genes = keys->size()*2;
	(*partitions) = new vector<string>(num_genes);

	for (int partitionId = 0; partitionId < num_genes; partitionId++) {
		string partitionLine(keys->at(partitionId * 2));
		replace(keys->at(partitionId * 2 + 1).begin(), keys->at(partitionId * 2 + 1).end(), '\\', '/');
		replace(keys->at(partitionId * 2 + 1).begin(), keys->at(partitionId * 2 + 1).end(), ';', '\0');
		partitionLine.append(" = ");
		partitionLine.append(keys->at(partitionId * 2 + 1).begin(), keys->at(partitionId * 2 + 1).end());
		(*partitions)->push_back(partitionLine);
	}
	delete keys;

	return 0;
}

int PartestParserUtils::parseRaxmlFile(vector<string> ** partitions) {
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

	(*partitions) = new vector<string>(numberOfModels);

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
		partitionInfo pInfo;
		pInfo.name = string(partitionName);
		ch++;

		numberPairs: pairsCount++;

		pInfo.numberOfSections = pairsCount;

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
		pInfo.start[pairsCount - 1] = lower;

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
		SINGLE_NUMBER: pInfo.end[pairsCount - 1] = upper;

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
			pInfo.stride[pairsCount - 1] = modulo;

			stringstream ss;
			ss << pInfo.name << " = ";
			for (int j = 0; j < pInfo.numberOfSections; j++) {
				if (j)
					ss << ",";
				ss << pInfo.start[j] << "-" << pInfo.end[j];
				if (pInfo.stride[j] > 1) {
					ss << "\\" << pInfo.stride[j];
				}
			}
			ss << endl;
			(*partitions)->at(i) = ss.str();

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
