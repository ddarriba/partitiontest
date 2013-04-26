/*
 * PLLModelOptimize.h
 *
 *  Created on: Jan 8, 2013
 *      Author: diego
 */

#ifndef PLLMODELOPTIMIZE_H_
#define PLLMODELOPTIMIZE_H_

#include "ModelOptimize.h"
#include "../model/ModelSet.h"

/* external C */
#ifndef AXML_H
#define AXML_H
#include "axml.h"
#endif
//#include "parser/phylip.h"
//extern "C" {
//void read_phylip_msa (tree * tr, const char * filename, int format, int type);
//}

namespace partest {

class PLLModelOptimize: public ModelOptimize {
public:
	PLLModelOptimize(ParTestOptions * options);
	int optimizeModel(Model * model, PartitionElement * partitionElement, int index,
			int groupCount);
	virtual ~PLLModelOptimize();
private:
	tree * tr;
	double **empiricalFrequencies;
};

} /* namespace partest */
#endif /* PLLMODELOPTIMIZE_H_ */
