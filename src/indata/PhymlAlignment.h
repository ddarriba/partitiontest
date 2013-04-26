/*
 * PhymlAlignment.h
 *
 *  Created on: Jan 21, 2013
 *      Author: diego
 */

#ifndef PHYMLALIGNMENT_H_
#define PHYMLALIGNMENT_H_

#include "Alignment.h"
extern "C" {
typedef double phydbl;
typedef struct __Calign {
	struct __Align **c_seq; /*! compressed sequences      */
	phydbl *b_frq; /*! observed state frequencies */
	short int *invar; /*! 1 -> states are identical, 0 states vary */
	int *wght; /*! # of each site in c_align */
	short int *ambigu; /*! ambigu[i]=1 is one or more of the sequences at site
	 i display an ambiguous character */
	phydbl obs_pinvar;
	int n_otu; /*! number of taxa */
	int clean_len; /*! uncrunched sequences lenghts without gaps */
	int crunch_len; /*! crunched sequences lengths */
	int init_len; /*! length of the uncompressed sequences */
	int *sitepatt; /*! this array maps the position of the patterns in the
	 compressed alignment to the positions in the uncompressed
	 one */
	int format; /*! 0 (default): PHYLIP. 1: NEXUS. */
};

typedef struct __Align {
	char *name; /*! sequence name */
	int len; /*! sequence length */
	char *state; /*! sequence itself */
	short int *is_ambigu; /*! is_ambigu[site] = 1 if state[site] is an ambiguous character. 0 otherwise */
} align;

}
namespace partest {

class PhymlAlignment: public Alignment {
public:
	PhymlAlignment(PhymlAlignment * alignment, int firstPosition,
			int lastPosition);
	PhymlAlignment(PhymlAlignment * alignment, int * firstPosition,
			int * lastPosition, int numberOfSections);
	PhymlAlignment(std::string alignmentFile, DataType dataType);
	virtual ~PhymlAlignment();
	struct __Calign *getCData();
	Alignment * splitAlignment(int firstPosition, int lastPosition);
	Alignment * splitAlignment(int * firstPosition, int * lastPosition,
			int numberOfSections);
private:
	struct __Calign * build_cdata(PhymlAlignment * alignment,
			int * firstPosition, int * lastPosition, int numberOfSections);
	struct __Calign * cdata;
};

} /* namespace partest */
#endif /* PHYMLALIGNMENT_H_ */
