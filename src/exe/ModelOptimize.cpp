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

#include "ModelOptimize.h"

#include "exe/ModelSelector.h"
#include "util/Utilities.h"
#include "indata/PartitionMap.h"
#include "indata/TreeManager.h"

#include <pll/parsePartition.h>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <cstring>
#include <cassert>

using namespace std;

namespace partest
{

  ModelOptimize::ModelOptimize ()
  {

  }

  ModelOptimize::~ModelOptimize ()
  {

  }

  string ModelOptimize::buildStartingFixedTree (int do_ml)
  {
    if (I_AM_ROOT)
    {
      bool loadedTree = false;

      if (ckpAvailable)
      {
        fstream ofs ((ckpPath + os_separator + ckpStartingTree).c_str (),
                     ios::in);

        if (ofs)
        {
          cout << timestamp () << " Loading topology from checkpoint..."
              << endl;
          ofs.seekg (0);
          size_t treeLen;
          ofs.read ((char *) &(treeLen), sizeof(size_t));
          starting_tree = (char *) malloc (treeLen + 1);
          ofs.read ((char *) starting_tree, (long) treeLen);

          pllNewickTree * nt;
          nt = pllNewickParseString (starting_tree);
          pllTreeInitTopologyNewick (tree, nt, PLL_FALSE);
          pllNewickParseDestroy (&nt);
          strcpy (tree->tree_string, starting_tree);
          free (starting_tree);

          size_t npergenetrees;
          ofs.read ((char *) &(npergenetrees), sizeof(size_t));
          if (npergenetrees)
          {
            pergene_starting_tree = (char **) malloc (
                (size_t) number_of_genes * sizeof(char *));
          }
          for (size_t i = 0; i < npergenetrees; i++)
          {
            ofs.read ((char *) &(treeLen), sizeof(size_t));
            pergene_starting_tree[i] = (char *) malloc (
                treeLen * sizeof(char) + 1);
            ofs.read ((char *) pergene_starting_tree[i], (long) treeLen);
          }

          ofs.close ();
          loadedTree = true;
        }
      }

      if (!loadedTree)
      {
        cout << timestamp () << " Computing fixed topology..." << endl;
        pllAlignmentData * alignData = 0;
        Utilities::duplicateAlignmentData (&alignData, phylip);

        partitionList * compParts = pllPartitionsCommit (pllPartsQueue,
                                                         alignData);

        if (pergene_branch_lengths)
        {
          compParts->perGeneBranchLengths = PLL_TRUE;
        }

        pllAlignmentRemoveDups (alignData, compParts);

        pllTreeInitTopologyForAlignment (tree, alignData);
        pllLoadAlignment (tree, alignData, compParts);

        pllComputeRandomizedStepwiseAdditionParsimonyTree (tree, compParts);

        if (do_ml)
        {
          tree->start = tree->nodep[1];

          switch (data_type)
            {
            case DT_PROTEIC:
              for (int cur_part = 0; cur_part < compParts->numberOfPartitions;
                  cur_part++)
              {
                pInfo * current_part = compParts->partitionData[cur_part];
                current_part->dataType = PLL_AA_DATA;
                current_part->states = 20;
                current_part->protUseEmpiricalFreqs = PLL_FALSE;
                current_part->optimizeBaseFrequencies = PLL_FALSE;
                current_part->optimizeAlphaParameter = PLL_TRUE;
                current_part->optimizeSubstitutionRates = PLL_FALSE;
                current_part->protModels = PLL_AUTO;
                current_part->alpha = 0.0;
              }
              cout << timestamp () << " Loading AUTO models" << endl;
              break;
            case DT_NUCLEIC:
              for (int cur_part = 0; cur_part < compParts->numberOfPartitions;
                  cur_part++)
              {
                pInfo * current_part = compParts->partitionData[cur_part];
                current_part->dataType = PLL_DNA_DATA;
                current_part->states = 4;
                current_part->optimizeBaseFrequencies = PLL_TRUE;
                current_part->optimizeAlphaParameter = PLL_TRUE;
                current_part->optimizeSubstitutionRates = PLL_TRUE;
              }
              cout << timestamp () << " Loading GTR models" << endl;
              break;
            default:
              assert(0);
            }

          pllInitModel (tree, compParts);

          tree->doCutoff = ML_PARAM_CUTOFF;
          if (epsilon == AUTO_EPSILON)
          {
            tree->likelihoodEpsilon = -0.001 * tree->likelihood;
          }
          else
          {
            tree->likelihoodEpsilon = epsilon;
          }
          tree->stepwidth = ML_PARAM_STEPWIDTH;
          tree->max_rearrange = ML_PARAM_MAXREARRANGE;
          tree->initial = tree->bestTrav = ML_PARAM_BESTTRAV;
          tree->initialSet = ML_PARAM_INITIALSET;

          cout << timestamp ()
              << " Building ML topology (this might take a while...)" << endl;
          pllRaxmlSearchAlgorithm (tree, compParts, PLL_TRUE);
        }
        else
        {
          cout << timestamp () << " Updating branch lengths..." << endl;
          pllInitModel (tree, compParts);
          double lk = 0;
          double optEpsilon = 1;
          pllEvaluateLikelihood (tree, compParts, tree->start,
          PLL_TRUE,
                                 PLL_FALSE);
          if (!reoptimize_branch_lengths)
          {
            while (fabs (lk - tree->likelihood) > optEpsilon)
            {
              lk = tree->likelihood;
              pllOptRatesGeneric (tree, compParts, optEpsilon,
                                  compParts->rateList);
              pllOptBaseFreqs (tree, compParts, optEpsilon,
                               compParts->freqList);
              pllEvaluateLikelihood (tree, compParts, tree->start,
              PLL_TRUE,
                                     PLL_FALSE);
              pllOptAlphasGeneric (tree, compParts, optEpsilon,
                                   compParts->alphaList);
              pllEvaluateLikelihood (tree, compParts, tree->start,
              PLL_TRUE, PLL_FALSE);
            }
          }
        }

        pllTreeToNewick (tree->tree_string, tree, compParts, tree->start->back,
        PLL_TRUE,
                         PLL_TRUE,
                         PLL_FALSE,
                         PLL_FALSE, PLL_FALSE, PLL_SUMMARIZE_LH, PLL_FALSE,
                         PLL_FALSE);
        tree->tree_string[tree->treeStringLength - 1] = '\0';
        if (pergene_branch_lengths)
        {
          pergene_starting_tree = (char **) malloc (
              (size_t) number_of_genes * sizeof(char *));
          for (size_t i = 0; i < number_of_genes; i++)
          {
            /* create one set of branch lengths for each gene */
            pergene_starting_tree[i] = (char *) malloc (
                (size_t) tree->treeStringLength * sizeof(char) + 1);
            pllTreeToNewick (pergene_starting_tree[i], tree, compParts,
                             tree->start->back,
                             PLL_TRUE,
                             PLL_TRUE,
                             PLL_FALSE,
                             PLL_FALSE, PLL_FALSE, (int) i,
                             PLL_FALSE,
                             PLL_FALSE);
          }
        }
        else
        {
          pergene_starting_tree = 0;
        }

        pllPartitionsDestroy (tree, &compParts);
        pllAlignmentDataDestroy (alignData);

        if (ckpAvailable)
        {
          /* store tree */
          fstream ofs ((ckpPath + os_separator + ckpStartingTree).c_str (),
                       ios::out);
          ofs.seekg (0);
          size_t treeLen = strlen (tree->tree_string) + 1;
          ofs.write ((char *) &treeLen, sizeof(size_t));
          ofs.write ((char *) tree->tree_string, (long) treeLen);
          if (pergene_starting_tree)
          {
            ofs.write ((char *) &number_of_genes, sizeof(size_t));
            for (size_t i = 0; i < number_of_genes; i++)
            {
              treeLen = strlen (pergene_starting_tree[i]) + 1;
              ofs.write ((char *) &treeLen, sizeof(size_t));
              ofs.write ((char *) pergene_starting_tree[i], (long) treeLen);
            }
          }
          else
          {
            size_t zero = 0;
            ofs.write ((char *) &zero, sizeof(size_t));
          }
          ofs.close ();
        }
      }

      starting_tree = tree->tree_string;
      cout << timestamp () << " Starting tree loaded " << starting_tree << endl;
    }
#ifdef HAVE_MPI
    unsigned long treelen;
    if (I_AM_ROOT)
      treelen = strlen (starting_tree);
    MPI_Bcast (&treelen, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
    if (!I_AM_ROOT)
      starting_tree = (char *) malloc (sizeof(char) * treelen);
    MPI_Bcast (starting_tree, treelen, MPI_CHAR, 0, MPI_COMM_WORLD);
    starting_tree[treelen] = '\0';
    MPI_Barrier (MPI_COMM_WORLD);
#endif
    return string (starting_tree);
  }

  string ModelOptimize::buildFinalTreeLinking (PartitioningScheme * finalScheme,
                                               bool reoptimizeParameters)
  {
    cout << timestamp () << " Computing fixed topology..." << endl;
    pllAlignmentData * alignData = 0;
    Utilities::duplicateAlignmentData (&alignData, phylip);

    partitionList * compParts = pllPartitionsCommit (pllPartsQueue, alignData);

    pllAlignmentRemoveDups (alignData, compParts);

    pllTreeInitTopologyForAlignment (tree, alignData);
    pllLoadAlignment (tree, alignData, compParts);

    pllComputeRandomizedStepwiseAdditionParsimonyTree (tree, compParts);
    tree->start = tree->nodep[1];

    for (size_t i = 0; i < finalScheme->getNumberOfElements (); i++)
    {
      PartitionElement * pe = finalScheme->getElement (i);
      switch (data_type)
        {
        case DT_PROTEIC:
          for (size_t cur_part = 0;
              cur_part < (size_t) compParts->numberOfPartitions; cur_part++)
          {
            pe = finalScheme->getElement (cur_part);
            ProteicModel * pModel =
                static_cast<ProteicModel *> (pe->getBestModel ()->getModel ());
            int matrix = pModel->getMatrix ();
            pInfo * current_part = compParts->partitionData[cur_part];
            current_part->dataType = PLL_AA_DATA;
            current_part->states = 20;
            current_part->protUseEmpiricalFreqs = pModel->isPF ();
            current_part->optimizeBaseFrequencies = PLL_FALSE;
            current_part->optimizeAlphaParameter =
                reoptimizeParameters ? PLL_TRUE : PLL_FALSE;
            current_part->optimizeSubstitutionRates = PLL_FALSE;
            current_part->protModels = matrix;
            current_part->alpha = pModel->getAlpha ();
          }
          break;
        case DT_NUCLEIC:
          {
            for (size_t cur_part = 0;
                cur_part < (size_t) compParts->numberOfPartitions; cur_part++)
            {
              pe = finalScheme->getElement (cur_part);
              NucleicModel * nModel =
                  static_cast<NucleicModel *> (pe->getBestModel ()->getModel ());
              //int matrix = nModel->getMatrix();
              pInfo * current_part = compParts->partitionData[cur_part];
              current_part->dataType = PLL_DNA_DATA;
              current_part->states = 4;
              current_part->optimizeBaseFrequencies = nModel->isPF ();
              current_part->optimizeAlphaParameter =
                  reoptimizeParameters ? PLL_TRUE : PLL_FALSE;
              current_part->optimizeSubstitutionRates =
                  reoptimizeParameters ? PLL_TRUE : PLL_FALSE;
            }
            break;
          }
        default:
          assert(0);
        }
    }

    pllInitModel (tree, compParts);

    pllRaxmlSearchAlgorithm (tree, compParts, PLL_FALSE);

    pllTreeToNewick (tree->tree_string, tree, compParts, tree->start->back,
    PLL_TRUE,
                     PLL_TRUE,
                     PLL_FALSE,
                     PLL_FALSE, PLL_FALSE, PLL_SUMMARIZE_LH, PLL_FALSE,
                     PLL_FALSE);
    tree->tree_string[tree->treeStringLength - 1] = '\0';

    pllPartitionsDestroy (tree, &compParts);
    pllAlignmentDataDestroy (alignData);

    starting_tree = tree->tree_string;
    cout << timestamp () << " Starting tree loaded" << endl << endl;
    return string (tree->tree_string);
  }

  string ModelOptimize::buildFinalTree (PartitioningScheme * finalScheme,
                                        bool reoptimizeParameters)
  {
    bool loadedTree = false;

    char * final_tree = 0;

    if (ckpAvailable)
    {
      fstream ofs ((ckpPath + os_separator + ckpFinalTree).c_str (), ios::in);

      if (ofs)
      {
        cout << endl << timestamp () << " Loading topology from checkpoint..."
            << endl;
        ofs.seekg (0);
        int treeLen;
        ofs.read ((char *) &(treeLen), sizeof(int));
        final_tree = (char *) malloc ((size_t) treeLen + 1);
        ofs.read ((char *) final_tree, treeLen);
        ofs.close ();

        loadedTree = true;
      }
    }

    if (!loadedTree)
    {

      pllInstanceAttr attr;
      attr.fastScaling = PLL_FALSE;
      attr.randomNumberSeed = 0x54321;
      attr.rateHetModel = PLL_GAMMA;
      attr.saveMemory = PLL_FALSE;
      attr.useRecom = PLL_FALSE;
      attr.numberOfThreads = number_of_threads;

      pllInstance * fTree = pllCreateInstance (&attr);

      pllQueue * parts;
      pllPartitionRegion * pregion;
      pllPartitionInfo * pinfo;

      pllQueueInit (&parts);
      for (size_t i = 0; i < finalScheme->getNumberOfElements (); i++)
      {
        PartitionElement * pe = finalScheme->getElement (i);
        pinfo = (pllPartitionInfo *) malloc (sizeof(pllPartitionInfo));
        pinfo->ascBias = PLL_FALSE;
        pllQueueInit (&(pinfo->regionList));
        pllQueueAppend (parts, (void *) pinfo);

        pinfo->partitionName = (char *) malloc (5);
        strcpy (pinfo->partitionName, "PART");
        pinfo->partitionModel = (char *) malloc (1);
        pinfo->protModels = -1;
        pinfo->protUseEmpiricalFreqs = -1;
        pinfo->dataType =
            (data_type == DT_NUCLEIC) ? PLL_DNA_DATA : PLL_AA_DATA;
        pinfo->optimizeBaseFrequencies = PLL_TRUE;
        for (size_t j = 0; j < pe->getNumberOfSections (); j++)
        {
          pregion = (pllPartitionRegion *) malloc (sizeof(pllPartitionRegion));
          pregion->start = (int) pe->getSection (j).start;
          pregion->end = (int) pe->getSection (j).end;
          pregion->stride = 1;
          pllQueueAppend (pinfo->regionList, (void *) pregion);
        }
      }
      partitionList * compParts = pllPartitionsCommit (parts, phylip);
      pllQueuePartitionsDestroy (&parts);

      cout << endl << timestamp ()
          << " Conducting final topology optimization... " << endl;

      pllAlignmentRemoveDups (phylip, compParts);

      pllTreeInitTopologyForAlignment (fTree, phylip);
      pllLoadAlignment (fTree, phylip, compParts);

      switch (starting_topology)
        {
        case StartTopoMP:
        case StartTopoML:
          pllComputeRandomizedStepwiseAdditionParsimonyTree (fTree, compParts);
          fTree->start = fTree->nodep[1];
          break;
        case StartTopoFIXED:
        case StartTopoFIXEDML:
          {
            pllNewickTree * nt;
            nt = pllNewickParseString (starting_tree);
            pllTreeInitTopologyNewick (fTree, nt, PLL_FALSE);
            pllNewickParseDestroy (&nt);
            break;
          }
        case StartTopoUSER:
          pllNewickTree * nt;
          nt = pllNewickParseFile (user_tree->c_str ());
          pllTreeInitTopologyNewick (fTree, nt, PLL_FALSE);
          pllNewickParseDestroy (&nt);
          break;
        default:
          assert(0);
          break;
        }

      switch (data_type)
        {
        case DT_PROTEIC:
          {
            for (size_t cur_part = 0;
                cur_part < (size_t) compParts->numberOfPartitions; cur_part++)
            {
              PartitionElement * pe = finalScheme->getElement (cur_part);
              ProteicModel * pModel =
                  static_cast<ProteicModel *> (pe->getBestModel ()->getModel ());
              int matrix = pModel->getMatrix ();
              pInfo * current_part = compParts->partitionData[cur_part];
              current_part->dataType = PLL_AA_DATA;
              current_part->states = 20;
              current_part->protUseEmpiricalFreqs = pModel->isPF ();
              current_part->optimizeBaseFrequencies = PLL_FALSE;
              current_part->optimizeAlphaParameter =
                  reoptimizeParameters ? PLL_TRUE : PLL_FALSE;
              current_part->optimizeSubstitutionRates = PLL_FALSE;
              current_part->protModels = matrix;
              current_part->alpha = pModel->getAlpha ();
              if (!reoptimizeParameters)
              {
                current_part->alpha = pModel->getAlpha ();
                memcpy (current_part->frequencies, pModel->getFrequencies (),
                NUM_PROT_FREQS);
                memcpy (current_part->substRates, pModel->getRates (),
                NUM_AA_RATES);
              }
            }
            break;
          }
        case DT_NUCLEIC:
          {
            for (size_t cur_part = 0;
                cur_part < (size_t) compParts->numberOfPartitions; cur_part++)
            {
              PartitionElement * pe = finalScheme->getElement (cur_part);
              NucleicModel * nModel =
                  static_cast<NucleicModel *> (pe->getBestModel ()->getModel ());
              pInfo * current_part = compParts->partitionData[cur_part];
              current_part->dataType = PLL_DNA_DATA;
              current_part->states = 4;
              current_part->optimizeBaseFrequencies = nModel->isPF ();
              current_part->optimizeAlphaParameter =
                  reoptimizeParameters ? PLL_TRUE : PLL_FALSE;
              current_part->optimizeSubstitutionRates =
                  reoptimizeParameters ? PLL_TRUE : PLL_FALSE;
            }
            break;
          }
        default:
          assert(0);
        }

      pllInitModel (fTree, compParts);

      if (!reoptimizeParameters)
      {
        for (size_t cur_part = 0;
            cur_part < (size_t) compParts->numberOfPartitions; cur_part++)
        {
          pInfo * current_part = compParts->partitionData[cur_part];
          PartitionElement * pe = finalScheme->getElement (cur_part);
          current_part->alpha = pe->getBestModel ()->getModel ()->getAlpha ();
          memcpy (current_part->frequencies,
                  pe->getBestModel ()->getModel ()->getFrequencies (),
                  NUM_NUC_FREQS);
          memcpy (current_part->substRates,
                  pe->getBestModel ()->getModel ()->getRates (),
                  NUM_DNA_RATES);
          if (data_type == DT_NUCLEIC)
          {
            if (!pe->getBestModel ()->getModel ()->isPF ())
            {
              for (int i = 0; i < 4; i++)
              {
                current_part->frequencies[i] = 0.25;
              }
            }
          }
        }
      }

      if (data_type == DT_NUCLEIC)
      {
        for (size_t cur_part = 0;
            cur_part < (size_t) compParts->numberOfPartitions; cur_part++)
        {
          PartitionElement * pe = finalScheme->getElement (cur_part);
          NucleicModel * nModel =
              static_cast<NucleicModel *> (pe->getBestModel ()->getModel ());
          const char * m = nModel->getMatrixName ().c_str ();
          char * symmetryPar = (char *) malloc (12 * sizeof(char));
          symmetryPar[0] = m[0];
          symmetryPar[11] = '\0';
          for (int j = 1; j < 6; j++)
          {
            symmetryPar[(j - 1) * 2 + 1] = ',';
            symmetryPar[j * 2] = m[j];
          }
          pllSetSubstitutionRateMatrixSymmetries (symmetryPar, compParts,
                                                  (int) cur_part);
          free (symmetryPar);
        }
      }
      pllRaxmlSearchAlgorithm (fTree, compParts, PLL_FALSE);

      pllTreeToNewick (fTree->tree_string, fTree, compParts, fTree->start->back,
      PLL_TRUE,
                       PLL_TRUE,
                       PLL_FALSE,
                       PLL_FALSE, PLL_FALSE, PLL_SUMMARIZE_LH, PLL_FALSE,
                       PLL_FALSE);
      fTree->tree_string[strlen (fTree->tree_string) - 1] = '\0';

      pllPartitionsDestroy (fTree, &compParts);

      int treeLen = (int) strlen (fTree->tree_string) + 1;
      final_tree = (char *) malloc ((size_t) treeLen);
      strcpy (final_tree, fTree->tree_string);

      if (ckpAvailable)
      {
        /* store tree */
        fstream ofs ((ckpPath + os_separator + ckpFinalTree).c_str (),
                     ios::out);
        ofs.seekg (0);

        ofs.write ((char *) &treeLen, sizeof(int));
        ofs.write ((char *) final_tree, treeLen);
        ofs.close ();
      }

      cout << timestamp () << " Final tree lnL: " << fixed << setprecision (4)
          << fTree->likelihood << endl;
      pllDestroyInstance (fTree);
    }

    cout << timestamp () << " Final tree: " << final_tree << endl;

    string finalTreeStr (final_tree);
    free (final_tree);

    return finalTreeStr;
  }

  int ModelOptimize::optimizePartitioningScheme (PartitioningScheme * scheme,
                                                 int index, int limit)
  {

    cout << timestamp () << " - scheme "
        << setw (Utilities::iDecLog (limit) + 1) << setfill ('0') << right
        << index + 1 << "/" << limit << setfill (' ') << endl;

    /* check number of elements to optimize */
    int elementsToOptimize = 0;
    for (size_t cur_element = 0; cur_element < scheme->getNumberOfElements ();
        cur_element++)
    {
      if (!scheme->getElement (cur_element)->isOptimized ())
      {
        elementsToOptimize++;
      }
    }

    int currentElementToOptimize = 0;
    for (size_t cur_element = 0; cur_element < scheme->getNumberOfElements ();
        cur_element++)
    {
      PartitionElement * element = scheme->getElement (cur_element);
      if (!element->isOptimized ())
      {
        optimizePartitionElement (element, currentElementToOptimize,
                                  elementsToOptimize);
        currentElementToOptimize++;
      }
    }

    return EX_OK;
  }

  int ModelOptimize::optimizePartitionElement (PartitionElement * element,
                                               int index, int limit)
  {

    if (element->isOptimized ())
    {
      return EX_OK;
    }
    element->setupStructures ();

    cout << timestamp () << " - - -";
#ifdef HAVE_MPI
    cout << " [" << myRank << "]";
#endif
    cout << " element " << setw (Utilities::iDecLog (limit) + 1)
        << setfill ('0') << right << index + 1 << "/" << limit << setfill (' ');
#ifdef DEBUG
    if (epsilon == AUTO_EPSILON)
      cout << " [" << element->getEpsilon () << "]";
#endif
    cout << endl;

    for (size_t modelIndex = 0; modelIndex < element->getNumberOfModels ();
        modelIndex++)
    {
      optimizeModel (element, modelIndex, (int) element->getNumberOfModels ());
    }

    ModelSelector ms (element, ic_type, element->getSampleSize ());

    element->destroyStructures ();

    return EX_OK;
  }

  void ModelOptimize::optimizeModel (PartitionElement * element,
                                     size_t modelIndex, int limit)
  {

    TreeManager * treeManager = element->getTreeManager ();
    Model * model = element->getModel (modelIndex);
    double lk;

    /* set parameters for single partition element */
    treeManager->setModelParameters (model, 0, false);

    if (starting_topology == StartTopoML)
    {
      treeManager->searchMlTopology (true);
    }
    else
    {
      /* main optimization loop */
      double cur_epsilon = epsilon;
      if (epsilon == AUTO_EPSILON)
      {
        cur_epsilon = element->getEpsilon ();
      }
      int smoothIterations = 32;
      int iters = 5;
      do
      {
        lk = treeManager->getLikelihood ();
        treeManager->optimizeBranchLengths (smoothIterations);
        treeManager->optimizeModelParameters (cur_epsilon);
        iters--;
      }
      while (fabs (lk - treeManager->getLikelihood ()) > cur_epsilon && iters>0);
    }

    if (!isfinite (treeManager->getLikelihood ()))
    {
      cerr << "[ERROR] Likelihood score for partition " << element->getName ()
          << " is " << treeManager->getLikelihood () << endl;
      exit_partest (EX_DATAERR);
    }

    /* set newick tree for optimized model */
    model->setLnL (treeManager->getLikelihood ());
    model->setTree (treeManager->getNewickTree ());
    model->setBranchLengthsScaler (treeManager->getBranchLengthMultiplier ());

    model->setFrequencies (treeManager->getFrequencies ());
    if (model->isGamma ())
      model->setAlpha (treeManager->getAlpha ());
    model->setRates (treeManager->getRates ());

    if (data_type == DT_PROTEIC && optimize_mode == OPT_GTR)
    {
      /* set chosen model */
      ProteicModel * pModel = static_cast<ProteicModel *> (model);
      pModel->setMatrix (
          static_cast<ProtMatrix> (treeManager->getAutoProtModel ()));
      model->setName (
          Utilities::getProtMatrixName (
              static_cast<ProtMatrix> (treeManager->getAutoProtModel ())));
    }

    if (verbosity)
    {
      cout << timestamp () << " - - - - -";
#ifdef HAVE_MPI
      cout << " [" << myRank << "]";
#endif
      cout << " " << setw (Utilities::iDecLog (limit) + 1) << setfill ('0')
          << right << modelIndex + 1 << "/" << limit << " " << model->getName ()
          << " (" << fixed << setprecision (4) << treeManager->getLikelihood ()
          << ")" << setfill (' ') << endl;
    }
  }

} /* namespace partest */
