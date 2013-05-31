/*
 * ParTestToPhyml.cpp
 *
 *  Created on: Jan 17, 2013
 *      Author: diego
 */

#include "spr.h"
#include "lk.h"
#include "optimiz.h"
#include "bionj.h"
#include "models.h"
#include "free.h"
#include "simu.h"
#include "eigen.h"
#include "pars.h"
#include "alrt.h"
#include "partest_to_phyml.h"

void getRates(int i, phydbl *matrix1, phydbl *freqs) {
	switch (i) {
	case 0:
		Init_Qmat_LG(matrix1, freqs);
		break;
	case 1:
		Init_Qmat_WAG(matrix1, freqs);
		break;
	case 2:
		Init_Qmat_JTT(matrix1, freqs);
		break;
	case 3:
		Init_Qmat_MtREV(matrix1, freqs);
		break;
	case 4:
		Init_Qmat_Dayhoff(matrix1, freqs);
		break;
	case 5:
		Init_Qmat_DCMut(matrix1, freqs);
		break;
	case 6:
		Init_Qmat_RtREV(matrix1, freqs);
		break;
	case 7:
		Init_Qmat_CpREV(matrix1, freqs);
		break;
	case 8:
		Init_Qmat_VT(matrix1, freqs);
		break;
	case 9:
		Init_Qmat_Blosum62(matrix1, freqs);
		break;
	case 10:
		Init_Qmat_MtMam(matrix1, freqs);
		break;
	case 11:
		Init_Qmat_MtArt(matrix1, freqs);
		break;
	case 12:
		Init_Qmat_HIVw(matrix1, freqs);
		break;
	case 13:
		Init_Qmat_HIVb(matrix1, freqs);
		break;
	}
}

void compute_distances() {
	phydbl * matrix1 = alloca(400 * sizeof(phydbl));
	phydbl * matrix2 = alloca(400 * sizeof(phydbl));
	phydbl * freqs = alloca(20 * sizeof(phydbl));
	int i, j, k;
	for (i = 1; i < 14; i++) {
		for (j = 0; j < i; j++) {
			getRates(i, matrix1, freqs);
			getRates(j, matrix2, freqs);

			double distance = 0.0;
			for (k = 0; k < 20; k++) {
				distance += (matrix2[k] - matrix1[k])
						* (matrix2[k] - matrix1[k]);
			}
			distance = sqrt(distance);
			/* Compute euclidean distances */
			printf("%6.2f\t", distance);
		}
		printf("\n");
	}
	exit(0);
}

void get_aa_freqs(calign *data) {
	int i, j, k;
	phydbl A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y;
	phydbl fA, fC, fD, fE, fF, fG, fH, fI, fK, fL, fM, fN, fP, fQ, fR, fS, fT,
			fV, fW, fY;
	int w;
	phydbl sum;

	fA = fC = fD = fE = fF = fG = fH = fI = fK = fL = fM = fN = fP = fQ = fR =
			fS = fT = fV = fW = fY = 1. / 20.;

	For(k,8)
	{
		A = C = D = E = F = G = H = I = K = L = M = N = P = Q = R = S = T = V =
				W = Y = .0;

		For(i,data->n_otu)
		{
			For(j,data->crunch_len)
			{
				w = data->wght[j];
				if (w) {
					switch (data->c_seq[i]->state[j]) {
					case 'A':
						A += w;
						break;
					case 'C':
						C += w;
						break;
					case 'D':
						D += w;
						break;
					case 'E':
						E += w;
						break;
					case 'F':
						F += w;
						break;
					case 'G':
						G += w;
						break;
					case 'H':
						H += w;
						break;
					case 'I':
						I += w;
						break;
					case 'K':
						K += w;
						break;
					case 'L':
						L += w;
						break;
					case 'M':
						M += w;
						break;
					case 'N':
						N += w;
						break;
					case 'P':
						P += w;
						break;
					case 'Q':
						Q += w;
						break;
					case 'R':
						R += w;
						break;
					case 'S':
						S += w;
						break;
					case 'T':
						T += w;
						break;
					case 'V':
						V += w;
						break;
					case 'W':
						W += w;
						break;
					case 'Y':
						Y += w;
						break;
					case 'Z':
						Q += w;
						break;
					case 'X':
					case '?':
					case 'O':
					case '-':
						A += w * fA;
						C += w * fC;
						D += w * fD;
						E += w * fE;
						F += w * fF;
						G += w * fG;
						H += w * fH;
						I += w * fI;
						K += w * fK;
						L += w * fL;
						M += w * fM;
						N += w * fN;
						P += w * fP;
						Q += w * fQ;
						R += w * fR;
						S += w * fS;
						T += w * fT;
						V += w * fV;
						W += w * fW;
						Y += w * fY;
						break;
					default:
						break;
					}
				}
			}
		}
		sum = (A + C + D + E + F + G + H + I + K + L + M + N + P + Q + R + S + T
				+ V + W + Y);
		fA = A / sum;
		fC = C / sum;
		fD = D / sum;
		fE = E / sum;
		fF = F / sum;
		fG = G / sum;
		fH = H / sum;
		fI = I / sum;
		fK = K / sum;
		fL = L / sum;
		fM = M / sum;
		fN = N / sum;
		fP = P / sum;
		fQ = Q / sum;
		fR = R / sum;
		fS = S / sum;
		fT = T / sum;
		fV = V / sum;
		fW = W / sum;
		fY = Y / sum;
	}

	data->b_frq[0] = fA;
	data->b_frq[1] = fR;
	data->b_frq[2] = fN;
	data->b_frq[3] = fD;
	data->b_frq[4] = fC;
	data->b_frq[5] = fQ;
	data->b_frq[6] = fE;
	data->b_frq[7] = fG;
	data->b_frq[8] = fH;
	data->b_frq[9] = fI;
	data->b_frq[10] = fL;
	data->b_frq[11] = fK;
	data->b_frq[12] = fM;
	data->b_frq[13] = fF;
	data->b_frq[14] = fP;
	data->b_frq[15] = fS;
	data->b_frq[16] = fT;
	data->b_frq[17] = fW;
	data->b_frq[18] = fY;
	data->b_frq[19] = fV;
}

void get_nt_freqs(calign *data) {
	int i, j, k;
	phydbl A, C, G, T;
	phydbl fA, fC, fG, fT;
	int w;

	fA = fC = fG = fT = .25;

	For(k,8)
	{
		A = C = G = T = .0;
		For(i,data->n_otu)
		{
			For(j,data->crunch_len)
			{
				w = data->wght[j];
				if (w) {
					switch (data->c_seq[i]->state[j]) {
					case 'A':
						A += w;
						break;
					case 'C':
						C += w;
						break;
					case 'G':
						G += w;
						break;
					case 'T':
						T += w;
						break;
					case 'U':
						T += w;
						break;
					case 'M':
						C += w * fC / (fC + fA);
						A += w * fA / (fA + fC);
						break;
					case 'R':
						G += w * fG / (fA + fG);
						A += w * fA / (fA + fG);
						break;
					case 'W':
						T += w * fT / (fA + fT);
						A += w * fA / (fA + fT);
						break;
					case 'S':
						C += w * fC / (fC + fG);
						G += w * fG / (fC + fG);
						break;
					case 'Y':
						C += w * fC / (fC + fT);
						T += w * fT / (fT + fC);
						break;
					case 'K':
						G += w * fG / (fG + fT);
						T += w * fT / (fT + fG);
						break;
					case 'B':
						C += w * fC / (fC + fG + fT);
						G += w * fG / (fC + fG + fT);
						T += w * fT / (fC + fG + fT);
						break;
					case 'D':
						A += w * fA / (fA + fG + fT);
						G += w * fG / (fA + fG + fT);
						T += w * fT / (fA + fG + fT);
						break;
					case 'H':
						A += w * fA / (fA + fC + fT);
						C += w * fC / (fA + fC + fT);
						T += w * fT / (fA + fC + fT);
						break;
					case 'V':
						A += w * fA / (fA + fC + fG);
						C += w * fC / (fA + fC + fG);
						G += w * fG / (fA + fC + fG);
						break;
					case 'N':
					case 'X':
					case '?':
					case 'O':
					case '-':
						A += w * fA;
						C += w * fC;
						G += w * fG;
						T += w * fT;
						break;
					default:
						break;
					}
				}
			}
		}
		fA = A / (A + C + G + T);
		fC = C / (A + C + G + T);
		fG = G / (A + C + G + T);
		fT = T / (A + C + G + T);
	}

	data->b_frq[0] = fA;
	data->b_frq[1] = fC;
	data->b_frq[2] = fG;
	data->b_frq[3] = fT;
}

void print_options(option * io) {

	printf("**********BEGIN IO*************\n\n");
	printf(" in_tree:       %d\n", io->in_tree);
	printf(" in_align_file: %s\n", io->in_align_file);
	printf(" in_tree_file:  %s\n", io->in_tree_file);
	printf(" in_constraint_tree_file:  %s\n", io->in_constraint_tree_file);
	printf(" datatype:      %d\n", io->datatype);
	printf(" n_data_sets:   %d\n", io->n_data_sets);
	printf(" n_trees:       %d\n", io->n_trees);
	printf(" init_len:      %d\n", io->init_len);
	printf(" n_otu:         %d\n", io->n_otu);
	printf(" multigene:         %d\n", io->multigene);
	printf(" config_multigene:  %d\n", io->config_multigene);
	printf(" n_part:        %d\n", io->n_part);
	printf(" curr_gt:       %d\n", io->curr_gt);
	printf(" ratio_test:    %d\n", io->ratio_test);
	printf(" ready_to_go:   %d\n", io->ready_to_go);
	printf(" r_seed:        %d\n", io->r_seed);

	printf("\n**********BEGIN MODEL*************\n\n");
	struct __Model * mod = io->mod;

	printf(" modelname:   %s\n", mod->modelname);
	printf(" mod_num:       %d\n", mod->mod_num);
	printf(" n_diff_rr:     %d\n", mod->n_diff_rr);
	printf(" update_eigen:  %d\n", mod->update_eigen);
	printf(" whichmodel:    %d\n", mod->whichmodel);
	printf(" n_catg:        %d\n", mod->n_catg);
	printf(" invar:        %d\n", mod->invar);
	printf(" ns:           %d\n", mod->ns);
	printf(" use_m4mod:    %d\n", mod->use_m4mod);
	printf(" gamma_median: %d\n", mod->gamma_median);

	printf("\n**********BEGIN OPT*************\n\n");
	struct __Optimiz * opt = mod->s_opt;
	printf(" a: %d k: %d l: %d i: %d f: %d r: %d, p: %d\n", opt->opt_alpha,
			opt->opt_kappa, opt->opt_lambda, opt->opt_pinvar,
			opt->opt_state_freq, opt->opt_rr, opt->opt_subst_param);
	printf(" COV alpha %d delta %d free_rates %d\n", opt->opt_cov_alpha,
			opt->opt_cov_delta, opt->opt_cov_free_rates);
	printf(" Opt BL %d TOPO %d TOPO_SEARCH %d\n", opt->opt_bl, opt->opt_topo,
			opt->topo_search);
	printf(" INIT_LK %f n_it_max %d last_opt %d\n", opt->init_lk, opt->n_it_max,
			opt->last_opt);
	printf(" random_input_tree %d\n", opt->random_input_tree);
	printf(" n_rand_starts     %d\n", opt->n_rand_starts);
	printf(" brent_it_max      %d\n", opt->brent_it_max);
	printf(" steph_spr         %d\n", opt->steph_spr);
	printf(" user_state_freq   %d\n", opt->user_state_freq);
	printf(" opt_five_branch   %d\n", opt->opt_five_branch);

	printf(" min_diff_lk_local  %f\n", opt->min_diff_lk_local);
	printf(" min_diff_lk_global %f\n", opt->min_diff_lk_global);
	printf(" min_diff_lk_move   %f\n", opt->min_diff_lk_move);
	printf(" p_moves_to_examine %f\n", opt->p_moves_to_examine);

	printf(" fast_nni   %d\n", opt->fast_nni);
	printf(" greedy     %d\n", opt->greedy);
	printf(" spr_pars   %d\n", opt->spr_pars);
	printf(" spr_lnL    %d\n", opt->spr_lnL);
	printf(" constrained_br_len   %d\n", opt->constrained_br_len);
	printf(" opt_gamma_br_len     %d\n", opt->opt_gamma_br_len);
	printf(" opt_free_mixt_rates  %d\n", opt->opt_free_mixt_rates);
}

void print_cdata(calign * cdata) {
	int i, j, k;
	printf("**********BEGIN CDATA*************\n\n");
	printf(" n_otu:      %d\n", cdata->n_otu);
	printf(" clean_len:  %d\n", cdata->clean_len);
	printf(" crunch_len: %d\n", cdata->crunch_len);
	printf(" init_len:   %d\n", cdata->init_len);
	printf(" format:     %d\n", cdata->format);
	printf(" obs_pinvar: %f\n", cdata->obs_pinvar);
	int wght = 0;
	int ambigu = 0;
	int invar = 0;
	for (i = 0; i < cdata->crunch_len; i++) {
		wght += cdata->wght[i];
		ambigu += cdata->ambigu[i];
		invar += cdata->invar[i];
	}
	printf(" wght:   %d\n", wght);
	printf(" ambigu: %d\n", ambigu);
	printf(" invar:  %d\n", invar);
}

void free_cdata(calign * cdata) {
	int i;

	Free(cdata->invar);
	Free(cdata->wght);
	Free(cdata->ambigu);
	Free(cdata->b_frq);
	Free(cdata->sitepatt);
	For(i,cdata->n_otu)
	{
		Free(cdata->c_seq[i]->name);
		if (cdata->c_seq[i]->state) {
			Free(cdata->c_seq[i]->state);
			if (cdata->c_seq[i]->is_ambigu)
				Free(cdata->c_seq[i]->is_ambigu);
		}
		Free(cdata->c_seq[i]);
	}
	Free(cdata->c_seq);
	Free(cdata);
}

calign *clone_cdata(calign * indata, int numberOfStates) {
	calign *cdata;

	int i, j, k, site;
	int n_patt, which_patt;
	cdata = (calign *) mCalloc(1, sizeof(calign));
	cdata->n_otu = indata->n_otu;
	//cdata->c_seq = (align **) mCalloc(indata->n_otu, sizeof(align *));
	cdata->wght = (int *) mCalloc(indata->crunch_len, sizeof(int));
	cdata->b_frq = (phydbl *) mCalloc(numberOfStates, sizeof(phydbl));
	cdata->ambigu = (short int *) mCalloc(indata->crunch_len,
			sizeof(short int));
	cdata->invar = (short int *) mCalloc(indata->crunch_len, sizeof(short int));

	memcpy(cdata->wght, indata->wght, indata->crunch_len * sizeof(int));
	memcpy(cdata->b_frq, indata->b_frq, numberOfStates * sizeof(phydbl));
	memcpy(cdata->ambigu, indata->ambigu,
			indata->crunch_len * sizeof(short int));
	memcpy(cdata->invar, indata->invar, indata->crunch_len * sizeof(short int));

	cdata->crunch_len = indata->crunch_len;
	cdata->c_seq = indata->c_seq;
	cdata->init_len = indata->init_len;

	if (numberOfStates == 20) {
		get_aa_freqs(cdata);
	} else {
		get_nt_freqs(cdata);
	}
	//Get_AA_Freqs(cdata);

	return cdata;
}

void ParTest_Printf(char *format, ...) {
	return;
}

void ParTest_Fprintf(FILE *fp, char *format, ...) {
	return;
}

struct __Calign * read_data(const char * ioFile, int dataType) {
	option *io;
	model * mod;
	optimiz *s_opt;
	m4 *m4mod;
	io = (option *) Make_Input();
	mod = (model *) Make_Model_Basic();
//	s_opt = (optimiz *) Make_Optimiz();
	m4mod = (m4 *) M4_Make_Light();
	Set_Defaults_Input(io);
	Set_Defaults_Model(mod);
//	Set_Defaults_Optimiz(s_opt);

	io->ratio_test = 0;

	strcpy(io->in_align_file, ioFile);
	io->fp_in_align = Openfile(io->in_align_file, 0);

	io->mod = mod;
	io->mod->m4mod = m4mod;

	mod->io = io;
//	mod->s_opt = s_opt;

//	io->fp_in_align = Openfile(io->in_align_file, 0);

	/* DATATYPE */
	if (dataType == DATATYPE_NT) {
		io->datatype = NT;
	} else if (dataType == DATATYPE_AA) {
		io->datatype = AA;
	}

	struct __Calign *cdata = 0;

	Get_Seq(io);

	if (io->data) {
		cdata = Compact_Data(io->data, io);

		Free_Seq(io->data, cdata->n_otu);

		if (cdata)
			Check_Ambiguities(cdata, io->datatype, io->mod->state_len);
		else {
			printf("\n. Err in file %s at line %d\n", __FILE__, __LINE__);
			return NULL ;
		}
	}

	Free_Optimiz(mod->s_opt);
	Free_Custom_Model(mod);
	Free_Model_Basic(mod);
	M4_Free_M4_Model(mod->m4mod);
	Free(mod);

	if (dataType == DATATYPE_NT) {
		get_nt_freqs(cdata);
	} else if (dataType == DATATYPE_AA) {
		get_aa_freqs(cdata);
	}

	if (io->fp_in_align)
		fclose(io->fp_in_align);

	Free_Input(io);

	return cdata;
}

double memory_required(t_tree *tree) {
	/* Rough estimate of the amount of memory that has to be used */

	long int nbytes;
	int n_otu;
	model *mod;

	mod = tree->mod;
	n_otu = tree->io->n_otu;
	nbytes = 0;

	/* Partial Pars */
	nbytes += (2 * n_otu - 3) * 2 * tree->data->crunch_len * sizeof(int);
	nbytes += (2 * n_otu - 3) * 2 * tree->data->crunch_len
			* sizeof(unsigned int);
	nbytes += (2 * n_otu - 3) * 2 * tree->data->crunch_len * mod->ns
			* sizeof(int);

	/* Pmat */
	nbytes += (2 * n_otu - 3) * mod->n_catg * mod->ns * mod->ns
			* sizeof(phydbl);

	/* Partial Lk */
	nbytes += ((2 * n_otu - 3) * 2 - tree->n_otu) * tree->data->crunch_len
			* mod->n_catg * mod->ns * sizeof(phydbl);

	/* Scaling factors */
	nbytes += ((2 * n_otu - 3) * 2 - tree->n_otu) * tree->data->crunch_len
			* sizeof(int);

	return (double) nbytes / (1.E+06);

}

option * build_options(phyml_indata indata) {

	const char * ioFile = indata.ioFile;
	const char * treeFile = indata.treeFile;
	int dataType = indata.dataType;
	int freqType = indata.freqType;
	double * frequencies = indata.frequencies;
	int optimizeTopology = indata.optimizeTopology;
	int optimizeBranchLengths = indata.optimizeBranchLengths;
	int optimizeRates = indata.optimizeRates;
	int optimizePInvar = indata.optimizePInvar;
	int optimizeGamma = indata.optimizeGamma;
	char * pModel = malloc(sizeof(char) * (strlen(indata.model) + 2));
	strcpy(pModel, indata.model);

	option *io;
	model * mod;
	optimiz *s_opt;
	m4 *m4mod;
	io = (option *) Make_Input();
	mod = (model *) Make_Model_Basic();
	s_opt = (optimiz *) Make_Optimiz();
	m4mod = (m4 *) M4_Make_Light();
	Set_Defaults_Input(io);
	Set_Defaults_Model(mod);
	Set_Defaults_Optimiz(s_opt);

	io->mod = mod;
	io->mod->m4mod = m4mod;

	mod->io = io;
	mod->s_opt = s_opt;

	int i;

	io->ratio_test = 0;
	io->print_site_lnl = NO;

	strcpy(io->in_align_file, ioFile);
	io->fp_in_align = Openfile(io->in_align_file, 0);

	if (treeFile) {
		strcpy(io->in_tree_file, treeFile);
		io->in_tree = 2;
		io->fp_in_tree = Openfile(treeFile, 0);
	}

	/* DATATYPE */
	if (dataType == DATATYPE_NT) {
		io->datatype = NT;
		io->mod->ns = 4;
		io->mod->state_len = 1;
		io->mod->m4mod->n_o = 4;

		strcpy(io->mod->custom_mod_string, pModel);

		Make_Custom_Model(io->mod);
		Translate_Custom_Mod_String(io->mod);

		io->mod->whichmodel = CUSTOM;
		strcpy(io->mod->modelname, "custom");
		io->mod->s_opt->opt_kappa = 0;
		io->mod->s_opt->opt_rr = 1;

	} else if (dataType == DATATYPE_AA) {

		io->datatype = AA;
		io->mod->state_len = 1;
		io->mod->s_opt->opt_kappa = 0;
		io->mod->ns = 20;
		io->mod->m4mod->n_o = 20;

		if ((io->mod->whichmodel == JC69) || (io->mod->whichmodel == K80)
				|| (io->mod->whichmodel == F81)
				|| (io->mod->whichmodel == HKY85)
				|| (io->mod->whichmodel == F84) || (io->mod->whichmodel == TN93)
				|| (io->mod->whichmodel == GTR)
				|| (io->mod->whichmodel == CUSTOM)) {
			io->mod->whichmodel = LG;
			strcpy(io->mod->modelname, "LG\0");
		}
	}

	/* FREQUENCIES */
	if (freqType == FREQTYPE_EMPIRICAL) {
		if (io->datatype == NT)
			io->mod->s_opt->opt_state_freq = NO;
		else if (io->datatype == AA)
			io->mod->s_opt->opt_state_freq = YES;
	} else if (freqType == FREQTYPE_MODEL) {
		if (io->datatype == NT)
			io->mod->s_opt->opt_state_freq = YES;
		else if (io->datatype == AA)
			io->mod->s_opt->opt_state_freq = NO;
	} else if (freqType == FREQTYPE_CUSTOM) {
		phydbl sum;

		io->mod->s_opt->opt_state_freq = 0;
		io->mod->s_opt->user_state_freq = 1;

		io->mod->user_b_freq->v[0] = (phydbl) frequencies[0];
		io->mod->user_b_freq->v[1] = (phydbl) frequencies[1];
		io->mod->user_b_freq->v[2] = (phydbl) frequencies[2];
		io->mod->user_b_freq->v[3] = (phydbl) frequencies[3];

		sum = (io->mod->user_b_freq->v[0] + io->mod->user_b_freq->v[1]
				+ io->mod->user_b_freq->v[2] + io->mod->user_b_freq->v[3]);

		io->mod->user_b_freq->v[0] /= sum;
		io->mod->user_b_freq->v[1] /= sum;
		io->mod->user_b_freq->v[2] /= sum;
		io->mod->user_b_freq->v[3] /= sum;

		if (io->mod->user_b_freq->v[0] < .0 || io->mod->user_b_freq->v[1] < .0
				|| io->mod->user_b_freq->v[2] < .0
				|| io->mod->user_b_freq->v[3] < .0
				|| io->mod->user_b_freq->v[0] > 1.
				|| io->mod->user_b_freq->v[1] > 1.
				|| io->mod->user_b_freq->v[2] > 1.
				|| io->mod->user_b_freq->v[3] > 1.) {
			Warn_And_Exit("\n. Invalid base frequencies.\n");
		}
	}

	/* TREE OPTIMIZATION */
	io->mod->s_opt->opt_topo = optimizeTopology;
	io->mod->s_opt->opt_bl = optimizeBranchLengths;
	io->mod->s_opt->opt_subst_param = optimizeRates;

	/* PROPORTION OF INVARIABLE SITES */
	if (optimizePInvar) {
		io->mod->s_opt->opt_pinvar = 1;
		io->mod->invar = 1;
	} else {
		io->mod->pinvar->v = (phydbl) 0.0;
		io->mod->invar = 0;
		io->mod->s_opt->opt_pinvar = 0;
	}

	/* GAMMA RATES */
	if (optimizeGamma) {
		io->mod->s_opt->opt_alpha = 1;
		io->mod->n_catg = 4;
	} else {
		io->mod->s_opt->opt_alpha = 0;
		io->mod->alpha->v = 100.0;
		io->mod->n_catg = 1;
	}

	/* MODEL */
	for (i = 0; i < strlen(indata.model); i++)
		Uppercase(pModel + i);
	if (!isalpha(pModel[0])) {

		strcpy(io->mod->custom_mod_string, pModel);

		if (strlen(io->mod->custom_mod_string) != 6) {
			Warn_And_Exit("\n. The string should be of length 6.\n");
		} else {
			Make_Custom_Model(io->mod);
			Translate_Custom_Mod_String(io->mod);
		}

		io->datatype = NT;
		io->mod->whichmodel = CUSTOM;
		strcpy(io->mod->modelname, "custom");
		io->mod->s_opt->opt_kappa = 0;
		io->mod->s_opt->opt_rr = 1;

	} else if (strcmp(pModel, "JC69") == 0) {
		io->datatype = NT;
		io->mod->whichmodel = JC69;
	} else if (strcmp(pModel, "K80") == 0) {
		io->datatype = NT;
		io->mod->whichmodel = K80;
	} else if (strcmp(pModel, "F81") == 0) {
		io->datatype = NT;
		io->mod->whichmodel = F81;
	} else if (strcmp(pModel, "HKY85") == 0) {
		io->datatype = NT;
		io->mod->whichmodel = HKY85;
	} else if (strcmp(pModel, "F84") == 0) {
		io->datatype = NT;
		io->mod->whichmodel = F84;
	} else if (strcmp(pModel, "TN93") == 0) {
		io->datatype = NT;
		io->mod->whichmodel = TN93;
	} else if (strcmp(pModel, "GTR") == 0) {
		io->datatype = NT;
		io->mod->whichmodel = GTR;
	} else if (strcmp(pModel, "DAYHOFF") == 0) {
		io->datatype = AA;
		io->mod->whichmodel = DAYHOFF;
	} else if (strcmp(pModel, "JTT") == 0) {
		io->datatype = AA;
		io->mod->whichmodel = JTT;
	} else if (strcmp(pModel, "MTREV") == 0) {
		io->datatype = AA;
		io->mod->whichmodel = MTREV;
	} else if (strcmp(pModel, "LG") == 0) {
		io->datatype = AA;
		io->mod->whichmodel = LG;
	} else if (strcmp(pModel, "WAG") == 0) {
		io->datatype = AA;
		io->mod->whichmodel = WAG;
	} else if (strcmp(pModel, "DCMUT") == 0) {
		io->datatype = AA;
		io->mod->whichmodel = DCMUT;
	} else if (strcmp(pModel, "RTREV") == 0) {
		io->datatype = AA;
		io->mod->whichmodel = RTREV;
	} else if (strcmp(pModel, "CPREV") == 0) {
		io->datatype = AA;
		io->mod->whichmodel = CPREV;
	} else if (strcmp(pModel, "VT") == 0) {
		io->datatype = AA;
		io->mod->whichmodel = VT;
	} else if (strcmp(pModel, "BLOSUM62") == 0) {
		io->datatype = AA;
		io->mod->whichmodel = BLOSUM62;
	} else if (strcmp(pModel, "MTMAM") == 0) {
		io->datatype = AA;
		io->mod->whichmodel = MTMAM;
	} else if (strcmp(pModel, "MTART") == 0) {
		io->datatype = AA;
		io->mod->whichmodel = MTART;
	} else if (strcmp(pModel, "HIVW") == 0) {
		io->datatype = AA;
		io->mod->whichmodel = HIVW;
	} else if (strcmp(pModel, "HIVB") == 0) {
		io->datatype = AA;
		io->mod->whichmodel = HIVB;
	} else if (strcmp(pModel, "CUSTOM") == 0) {
		io->datatype = AA;
		io->mod->whichmodel = CUSTOMAA;
	} else {
		PhyML_Printf(
				"\n. The model name is incorrect. Please see the documentation.\n");
		Exit("\n");
	}

	Set_Model_Name(io->mod);

	//TODO: SOME TEMPORARY THINGS
	//io->mod->s_opt->constrained_br_len = YES;
	io->mod->s_opt->min_diff_lk_local = 0.01;
	io->mod->s_opt->min_diff_lk_global = 0.01;
	//TODO: *********************

	if (io->mod->s_opt->constrained_br_len == YES) {
		io->mod->s_opt->opt_topo = NO;
	}

	if (io->datatype == UNDEFINED)
		io->datatype = NT;

	if ((io->mod->s_opt->n_rand_starts)
			&& (io->mod->s_opt->topo_search == NNI_MOVE)
			&& (io->mod->s_opt->random_input_tree)) {
		Warn_And_Exit(
				"\n. The random starting tree option is only compatible with SPR based search options.\n");
	}

	io->mod->s_opt->opt_cov_delta = 0;
	io->mod->s_opt->opt_cov_alpha = 0;
	io->mod->s_opt->opt_cov_free_rates = 0;

	io->print_trace = 0;
	io->print_boot_trees = 0;

	if (io->appebr_run_ID) {
		strcat(io->out_tree_file, "_");
		strcat(io->out_stats_file, "_");
		strcat(io->out_tree_file, io->run_id_string);
		strcat(io->out_stats_file, io->run_id_string);
	}
	strcat(io->out_tree_file, ".txt");
	strcat(io->out_stats_file, ".txt");

	if (io->mod->n_catg == 1)
		io->mod->s_opt->opt_alpha = 0;

	if (!io->mod->s_opt->opt_subst_param) {
		io->mod->s_opt->opt_alpha = 0;
		io->mod->s_opt->opt_kappa = 0;
		io->mod->s_opt->opt_lambda = 0;
		io->mod->s_opt->opt_pinvar = 0;
		io->mod->s_opt->opt_rr = 0;
	}

	if (io->mod->whichmodel != K80 && io->mod->whichmodel != HKY85
			&& io->mod->whichmodel != F84 && io->mod->whichmodel != TN93) {
		io->mod->s_opt->opt_kappa = 0;
	}

	if (io->datatype == AA && io->mod->whichmodel == CUSTOMAA
			&& !io->fp_aa_rate_mat) {
		PhyML_Printf(
				"\n. Custom model option with amino-acid requires you to specify a rate matrix file through the '--aa_rate_file' option.\n");
		Exit("\n");
	}

	// TODO: THIS OUTFILES !!!
	io->fp_out_tree = NULL;	//Openfile(io->out_tree_file, 1);
	io->fp_out_stats = NULL;	//Openfile(io->out_stats_file, 1);

	if (io->mod->whichmodel == GTR) {
		Make_Custom_Model(io->mod);
		io->mod->s_opt->opt_rr = 1;
	}

	free(pModel);
	return io;
}

double phyml_lk(option *io, phyml_outdata * outdata) {

	calign *cdata = clone_cdata(io->cdata, io->mod->ns);
	t_tree *tree;
	int n_otu, num_data_set;
	int num_tree, tree_line_number, num_rand_tree;
	matrix *mat;
	model * mod;
	time_t t_beg, t_end;
	phydbl best_lnL, most_likely_size, tree_size;
	int r_seed;
	char *most_likely_tree = NULL;
#ifdef DEBUG
	printf("[TRACE PHYML] SITES %d %d\n", io->cdata->init_len, io->cdata->crunch_len);
#endif
#ifdef QUIET
	//setvbuf(stdout,NULL,_IOFBF,2048);
#endif

	tree = NULL;
	mod = NULL;
	best_lnL = UNLIKELY;
	most_likely_size = -1.0;
	tree_size = -1.0;

	r_seed = (io->r_seed < 0) ? (time(NULL )) : (io->r_seed);
	srand(r_seed);
	io->r_seed = r_seed;
	io->quiet = 1;
	if (io->in_tree == 2)
		Test_Multiple_Data_Set_Format(io);
	else
		io->n_trees = 1;

	if (io->n_trees == 0 && io->in_tree == 2) {
		ParTest_Printf(
				"\n. The input tree file does not provide a tree in valid format.");
		return 0.0;
	}
	mat = NULL;
	tree_line_number = 0;

	n_otu = 0;
	best_lnL = UNLIKELY;
	//Get_Seq(io);
	Make_Model_Complete(io->mod);
	Set_Model_Name(io->mod);

#ifdef PHYML_DEBUG
	Print_Settings(io);
#endif

	mod = io->mod;

	if (!io->mod->s_opt->random_input_tree)
		io->mod->s_opt->n_rand_starts = 1;

	if ((io->mod->s_opt->random_input_tree)
			&& (io->mod->s_opt->topo_search != NNI_MOVE))
		if (!io->quiet)
			ParTest_Printf("\n. [Random start %3d/%3d]\n", num_rand_tree + 1,
					io->mod->s_opt->n_rand_starts);

	Init_Model(cdata, mod, io);

	if (io->mod->use_m4mod)
		M4_Init_Model(mod->m4mod, cdata, mod);

	switch (io->in_tree) {
	case 0:
	case 1: {
		tree = Dist_And_BioNJ(cdata, mod, io);
		break;
	}
	case 2: {
		tree = Read_User_Tree(cdata, mod, io);
		break;
	}
	}

	time(&t_beg);
	time(&(tree->t_beg));

	tree->mod = mod;
	tree->io = io;
	tree->data = cdata;
	tree->both_sides = YES;
	tree->n_pattern = tree->data->crunch_len;

	Prepare_Tree_For_Lk(tree);

	if (io->in_tree == 1)
		Spr_Pars(tree);

	if (io->do_alias_subpatt) {
		tree->update_alias_subpatt = YES;
		Lk(tree);
		tree->update_alias_subpatt = NO;
	}

#ifdef _MOCK_COMPUTATION
	tree->c_lnL = 100.0;
#else
	if (tree->mod->s_opt->opt_topo) {
		if (tree->mod->s_opt->topo_search == NNI_MOVE)
			Simu_Loop(tree);
		else if (tree->mod->s_opt->topo_search == SPR_MOVE)
			Speed_Spr_Loop(tree);
		else
			Best_Of_NNI_And_SPR(tree);

		if (tree->n_root)
			Add_Root(tree->t_edges[0], tree);
	} else {
		if (tree->mod->s_opt->opt_subst_param || tree->mod->s_opt->opt_bl)
			Round_Optimize(tree, tree->data, ROUND_MAX);
		else
			Lk(tree);
	}
#endif

#ifdef PHYML_DEBUG
	print_cdata(cdata);
	print_options(io);
#endif

	tree->both_sides = 1;
	Lk(tree);
	Pars(tree);
	Get_Tree_Size(tree);
	ParTest_Printf("\n\n. Log likelihood of the current tree: %f.\n",
			tree->c_lnL);

	Br_Len_Involving_Invar(tree);
	Rescale_Br_Len_Multiplier_Tree(tree);

	if (!tree->n_root)
		Get_Best_Root_Position(tree);

	/* Record the most likely tree in a string of characters */
	if (tree->c_lnL > best_lnL) {
		best_lnL = tree->c_lnL;
		if (most_likely_tree)
			Free(most_likely_tree);
		most_likely_tree = Write_Tree(tree, NO);
		most_likely_size = Get_Tree_Size(tree);
	}

	/* 		  JF(tree); */

	time(&t_end);

	/* Start from BioNJ tree */
	if ((num_rand_tree == io->mod->s_opt->n_rand_starts - 1)
			&& (tree->mod->s_opt->random_input_tree)) {
		/* Do one more iteration in the loop, but don't randomize the tree */
		num_rand_tree--;
		tree->mod->s_opt->random_input_tree = 0;
	}

	outdata->lnL = best_lnL;
	outdata->pinv = io->mod->pinvar->v;
	outdata->alpha = io->mod->alpha->v;
	outdata->tree = (char *) malloc(strlen(most_likely_tree) + 1);
	strcpy(outdata->tree, most_likely_tree);

	if (io->datatype == NT) {
		memcpy(outdata->rates,io->mod->rr->v,6 * sizeof(double));
		memcpy(outdata->frequencies,tree->mod->pi->v, 4*sizeof(double));
	} else {
		memcpy(outdata->frequencies,tree->mod->pi->v, 20*sizeof(double));
	}

	Free_Spr_List(tree);
	Free_One_Spr(tree->best_spr);
	if (tree->mat)
		Free_Mat(tree->mat);
	Free_Triplet(tree->triplet_struct);
	Free_Tree_Pars(tree);
	Free_Tree_Lk(tree);
	Free_Tree(tree);

	if (io->fp_out_tree) {
		if (io->n_data_sets == 1)
			rewind(io->fp_out_tree);
		ParTest_Fprintf(io->fp_out_tree, "%s\n", most_likely_tree);
	}

	free(cdata->wght);
	free(cdata->b_frq);
	free(cdata->ambigu);
	free(cdata->invar);
	free(cdata);

	Free_Model_Complete(mod);

	if (most_likely_tree)
		Free(most_likely_tree);

	if (mod->s_opt->n_rand_starts > 1)
		ParTest_Printf("\n. Best log likelihood: %f\n", best_lnL);

	Free_Optimiz(mod->s_opt);
	Free_Custom_Model(mod);
	Free_Model_Basic(mod);
	M4_Free_M4_Model(mod->m4mod);
	Free(mod);

	if (io->fp_in_align)
		fclose(io->fp_in_align);
	if (io->fp_in_tree)
		fclose(io->fp_in_tree);
	if (io->fp_out_lk)
		fclose(io->fp_out_lk);
	if (io->fp_out_tree)
		fclose(io->fp_out_tree);
	if (io->fp_out_trees)
		fclose(io->fp_out_trees);
	if (io->fp_out_stats)
		fclose(io->fp_out_stats);

	Free_Input(io);

	time(&t_end);

	return best_lnL;
}

