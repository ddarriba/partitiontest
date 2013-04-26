/*
 * PLLAlignment.h
 *
 *  Created on: Jan 21, 2013
 *      Author: diego
 */

#ifndef PLLALIGNMENT_H_
#define PLLALIGNMENT_H_

#include "Alignment.h"
#ifndef AXML_H
#define AXML_H
#include "../../axml.h"
#endif
extern "C" {
	#include "parser/phylip/phylip.h"
	#include "utils.h"
	struct pllPhylip * pllPhylipParse (const char *);
}

namespace partest {

class PLLAlignment: public Alignment {
public:
	PLLAlignment(PLLAlignment * alignment, int firstPosition,
			int lastPosition);
	PLLAlignment(PLLAlignment * alignment, int * firstPosition,
			int * lastPosition, int numberOfSections);
	PLLAlignment(std::string alignmentFile, DataType dataType);
	virtual ~PLLAlignment();
	Alignment * splitAlignment(int firstPosition, int lastPosition);
	Alignment * splitAlignment(int * firstPosition, int * lastPosition,
			int numberOfSections);
	tree * getTree() { return tr; }
private:
	tree * tr;
};

} /* namespace partest */
#endif /* PLLAlignment_H_ */
