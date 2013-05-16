#ifndef PARTEST_TO_PHYML_H
#define PARTEST_TO_PHYML_H

#define DATATYPE_NT 1
#define DATATYPE_AA 2
#define FREQTYPE_MODEL     1
#define FREQTYPE_EMPIRICAL 2
#define FREQTYPE_CUSTOM    3

typedef struct {
	const char * ioFile;
	const char * treeFile;
	struct __Tree * tree;
	int dataType;
	int freqType;
	double * frequencies;
	int optimizeTopology;
	int optimizeBranchLengths;
	int optimizeRates;
	int optimizePInvar;
	int optimizeGamma;
	const char * model;
} phyml_indata;

typedef struct {
	double lnL;
	double pinv;
	double alpha;
	double *rates;
	double *frequencies;
	char *tree;
} phyml_outdata;

void get_aa_freqs(struct __Calign *data);
void get_nt_freqs(struct __Calign *data);
void free_cdata(struct __Calign * cdata);
struct __Calign *read_data(const char * ioFile, int dataType);
struct __Option * build_options(phyml_indata indata);
double phyml_lk(struct __Option *io, phyml_outdata * outdata);

void compute_distances();

#endif /* PARTEST_TO_PHYML_H */
