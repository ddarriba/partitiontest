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
 * @file Model.h
 *
 * @brief Generic substitution model.
 */

#ifndef MODEL_H_
#define MODEL_H_

#include "util/GlobalDefs.h"
#include <string>

namespace partest
{

  /**
   * @brief Generic substitution model.
   */
  class Model
  {
  public:
    /**
     * Generic constructor of a model.
     *
     * @param rateVariation The rate variation and frequencies parameters (+I, +G, +F).
     * @param numberOfTaxa Number of taxa (required for computing
     *        the free parameters.
     */
    Model (bitMask rateVariation, int numberOfTaxa);

    /**
     * @brief Gets the name of the model.
     *
     * @return The name of the model.
     */
    std::string getName (void) const;

    /**
     * @brief Sets the name of the model.
     *
     * @param name The name of the model.
     */
    void setName (std::string name);

    /**
     * @brief Gets the name of the substitution scheme.
     *
     * @return The name of the substitution scheme.
     */
    std::string getMatrixName (void) const;

    /**
     * @brief Sets the name of the substitution scheme.
     *
     * @param matrixName The name of the substitution scheme.
     */
    void setMatrixName (std::string matrixName);

    /**
     * @brief Gets the rate variation and frequencies parameters as a bitmask.
     *
     * @return The rate variation and frequencies parameters as a bitmask.
     */
    bitMask getRateVariation (void) const;

    /**
     * @brief Gets the model tree in Newick format.
     *
     * @return The model tree in Newick format.
     */
    std::string getTree (void) const;

#ifdef _IG_MODELS
    /**
     * @brief Gets whether the model considers a proportion of invariable sites.
     *
     * @return True, if it is a +I model.
     */
    bool isPInv(void);
#endif

    /**
     * @brief Gets whether the model considers rate variation among sites.
     *
     * @return True, if it is a +G model.
     */
    bool isGamma (void) const;

    /**
     * @brief Gets whether the model considers model-driven/empirical/equal frequencies.
     *
     * For nucleotide models, this member function returns True if the model uses equal base
     * frequencies, or False for empirical frequencies.
     * For protein models, this member function returns True for empirical amino-acid
     * frequencies, or False for model-driven frequencies.
     *
     * @return True, if it is a +F model.
     */
    bool isPF (void) const;

    /**
     * @brief Gets whether the model was already optimized.
     *
     * Gets whether the model was already optimized. Thus, if True, the model has a lnL score
     * and gamma shape / pInv / frequencies.
     *
     * @return True, if the model was optimized.
     */
    bool isOptimized (void) const;

    /**
     * @brief Gets the positive log likelihood score.
     *
     * @return -lnL.
     */
    double getLnL (void) const
    {
      return lnL;
    }

    /**
     * @brief Gets the gamma shape (alpha parameter).
     *
     * This value makes sense only for +G models. Otherwise, it is 100.0
     *
     * @return The alpha parameter for the gamma distribution.
     */
    double getAlpha (void) const
    {
      return alpha;
    }

#ifdef _IG_MODELS
    /**
     * @brief Gets the proportion of invariable sites.
     *
     * @return The proportion of invariable sites.
     */
    double getpInv(void) const
    {
      return pInv;
    }
#endif

    /**
     * @brief Gets the number of state frequencies.
     *
     * @return 4 for nucleotide and 20 for amino-acids.
     */
    int getNumberOfFrequencies (void) const
    {
      return numberOfFrequencies;
    }

    /**
     * @brief Gets the total number of free parameters.
     *
     * @return The number of free parameters.
     */
    int getNumberOfFreeParameters (void) const
    {
      return (modelFreeParameters + treeFreeParameters);
    }

    /**
     * @brief Gets the branch lengths scaler
     *
     * If branch lengths are re-optimized, the scaler is 1.0.
     *
     * @return The branch lengths scaler.
     */
    double getBranchLengthsScaler (void) const
    {
      return (branchLengthsScaler);
    }

    void setBranchLengthsScaler (double value)
    {
      branchLengthsScaler = value;
    }

    /**
     * @brief Gets the number of free parameters of the model (whithout the branch lengths parameters).
     *
     * @return The number of free parameters of the model (without the tree).
     */
    int getModelFreeParameters (void) const
    {
      return modelFreeParameters;
    }

    /**
     * @brief Gets the number of free parameters of the tree (branch lengths).
     *
     * @return The number of free parameters of the tree (branch lengths).
     */
    int getTreeFreeParameters (void) const
    {
      return treeFreeParameters;
    }

    /**
     * @brief Sets the positive log likelihood score.
     *
     * @param _lnL Log likelihood in absolute value.
     */
    void setLnL (double _lnL)
    {
      this->lnL = _lnL;
    }

    /**
     * @brief Sets the gamma shape (alpha parameter).
     *
     * This value makes sense only for +G models. Otherwise, it is 100.0
     *
     * @param alpha The alpha parameter for the gamma distribution.
     */
    void setAlpha (double alpha);

#ifdef _IG_MODELS
    /**
     * @brief Sets the proportion of invariable sites.
     *
     * @param pInv The proportion of invariable sites.
     */
    void setpInv(double pInv);
#endif

    /**
     * @brief Sets the model tree in Newick format.
     *
     * @param tree The model tree in Newick format.
     */
    void setTree (std::string tree);

    /**
     * @brief Sets the model tree in Newick format.
     *
     * @param tree The model tree in Newick format.
     */
    void setTree (char * tree);

    double * getFrequencies (void) const;

    /**
     * @brief Sets the number of state frequencies.
     *
     * @param frequencies Array with the state frequencies: 4 for nucleotide and 20 for amino-acids.
     */
    virtual void setFrequencies (const double * frequencies) = 0;

    /**
     * Gets the exchangeability rates.
     */
    double * getRates (void) const;

    /**
     * @brief Sets the exchangeability rates.
     *
     * @param rates Array with the exchangeability rates.
     */
    virtual void setRates (const double * rates) = 0;
    /**
     * @brief Computes the distance to other model.
     *
     * @param other The other model.
     */
    virtual double distanceTo (Model * other) const = 0;

    /**
     * @brief Prints the model details and parameters.
     */
    virtual void print (std::ostream& out, const char * prefix = "") const;

    virtual ~Model ();
  protected:
    /** Rate variation (+I, +G, +F) */
    bitMask rateVariation;
    /** Likelihood */
    double lnL;
    /** Alpha value of gamma distribution */
    double alpha;
#ifdef _IG_MODELS
    /** Proportion of invariable sites */
    double pInv;
#endif
    /** Substitution rates */
    double *rates;
    /** State frequencies */
    double *frequencies;
    /** Number of state frequencies */
    int numberOfFrequencies;
    /** Full name of the model. */
    std::string name;
    /** Name of the model matrix. */
    std::string matrixName;
    /** Number of free parameters of the model. */
    int modelFreeParameters;
    /** Number of free parameters of the tree. */
    int treeFreeParameters;
    /** Most likely tree */
    std::string tree;
    /** Branch lengths scaler */
    double branchLengthsScaler;
  };

} /* namespace partest */

#endif /* MODEL_H_ */
