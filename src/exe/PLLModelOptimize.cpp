/*
 * PLLModelOptimize.cpp
 *
 *  Created on: Jan 8, 2013
 *      Author: diego
 */

#include "PLLModelOptimize.h"
#include <iostream>
#if ! (defined(__ppc) || defined(__powerpc__) || defined(PPC))
#include <xmmintrin.h>
#endif

namespace partest {

PLLModelOptimize::PLLModelOptimize(ParTestOptions * options) :
		ModelOptimize(options) {
//  tr = options->getAlignment()->getTree();
//  empiricalFrequencies = (double **)malloc(sizeof(double *));
//  empiricalFrequencies[0] = (double *)malloc(4*sizeof(double));
//  empiricalFrequencies[0][0] = empiricalFrequencies[0][1] = empiricalFrequencies[0][2] = empiricalFrequencies[0][3] = 0.25;
//  cout << tr << " " << empiricalFrequencies << endl;
// // read_phylip_msa (tr[0], (char *) dataFileName.c_str(), PHYLIP|  _SEQUENTIAL, 0);
//
//#ifndef _MOCK_COMPUTATION
////  initModel(tr, empiricalFrequencies);
//    makeRandomTree(tr);
////  evaluateGeneric(tr, tr->start, TRUE);
//#endif
//  ModelSet modelset(options->getRateVariation(), options->getDataType(), options->getAlignment()->getNumSeqs());
//  cout << tr << " " << empiricalFrequencies << endl;
//  for (int i=0; i<modelset.getNumberOfModels(); i++) optimizeModel(modelset.getModel(i), i);

}

PLLModelOptimize::~PLLModelOptimize() {
}

int PLLModelOptimize::optimizeModel(Model * model,
		PartitionElement * partitionElement, int index, int groupCount) {
#ifdef _MOCK_COMPUTATION
	model->setLnL(-1 * (rand()%500 + 10000));
	if (model->isPInv()) {
		model->setpInv(0.5d);
	}
	if (model->isGamma()) {
		model->setAlpha(0.5d);
	}
	if (model->isPF()) {
		double *freqs = (double *)malloc(model->getNumberOfFrequencies() * sizeof(double));
		int i;
		for (i=0; i < model->getNumberOfFrequencies(); i++) {
			freqs[i] = 1.0/model->getNumberOfFrequencies();
		}
		model->setFrequencies(freqs);
		free(freqs);
	}
#else
	model->print(cout);
	cout << tr << " " << empiricalFrequencies << endl;

//  treeEvaluate(tr, 32);

	printf("tree evaluated: %f\n", tr->likelihood);
	model->print(cout);
	return 0;
#endif
}

} /* namespace partest */
