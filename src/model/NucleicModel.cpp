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

#include "NucleicModel.h"
#include "util/Utilities.h"
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <cassert>

using namespace std;

namespace partest
{

  NucleicModel::NucleicModel (NucMatrix _matrix, bitMask rateVariation,
                              int numberOfTaxa) :
      Model (rateVariation, numberOfTaxa), matrix (_matrix)
  {
    /* treeFreeParameters is already initialized to the number of branches */
    this->numberOfFrequencies = NUM_NUC_FREQS;
    this->frequencies = (double *) malloc (
        (size_t) numberOfFrequencies * sizeof(double));
    for (int i = 0; i < numberOfFrequencies; i++)
      this->frequencies[i] = 1.0 / numberOfFrequencies;

    this->rates = (double *) malloc (NUM_DNA_RATES * sizeof(double));
    for (int i = 0; i < NUM_DNA_RATES; i++)
      this->rates[i] = 1.0;

    switch (matrix)
      {
      case NUC_MATRIX_JC:
        name.append ("JC");
        matrixName.append ("000000");
        break;
      case NUC_MATRIX_F81:
        name = "F81";
        matrixName.append ("000000");
        break;
      case NUC_MATRIX_K80:
        name.append ("K80");
        matrixName.append ("010010");
        modelFreeParameters += 1;
        break;
      case NUC_MATRIX_HKY:
        name.append ("HKY");
        matrixName.append ("010010");
        modelFreeParameters += 1;
        break;
      case NUC_MATRIX_TrNef:
        name.append ("TrNef");
        matrixName.append ("010020");
        modelFreeParameters += 2;
        break;
      case NUC_MATRIX_TrN:
        name.append ("TrN");
        matrixName.append ("010020");
        modelFreeParameters += 2;
        break;
      case NUC_MATRIX_TPM1:
        name.append ("TPM1");
        matrixName.append ("012210");
        modelFreeParameters += 2;
        break;
      case NUC_MATRIX_TPM1uf:
        name.append ("TPM1uf");
        matrixName.append ("012210");
        modelFreeParameters += 2;
        break;
      case NUC_MATRIX_TPM2:
        name.append ("TPM2");
        matrixName.append ("010212");
        modelFreeParameters += 2;
        break;
      case NUC_MATRIX_TPM2uf:
        name.append ("TPM2uf");
        matrixName.append ("010212");
        modelFreeParameters += 2;
        break;
      case NUC_MATRIX_TPM3:
        name.append ("TPM3");
        matrixName.append ("012012");
        modelFreeParameters += 2;
        break;
      case NUC_MATRIX_TPM3uf:
        name.append ("TPM3uf");
        matrixName.append ("012012");
        modelFreeParameters += 2;
        break;
      case NUC_MATRIX_TIM1ef:
        name.append ("TIM1ef");
        matrixName.append ("012230");
        modelFreeParameters += 3;
        break;
      case NUC_MATRIX_TIM1:
        name.append ("TIM1");
        matrixName.append ("012230");
        modelFreeParameters += 3;
        break;
      case NUC_MATRIX_TIM2ef:
        name.append ("TIM2ef");
        matrixName.append ("010232");
        modelFreeParameters += 3;
        break;
      case NUC_MATRIX_TIM2:
        name.append ("TIM2");
        matrixName.append ("010232");
        modelFreeParameters += 3;
        break;
      case NUC_MATRIX_TIM3ef:
        name.append ("TIM3ef");
        matrixName.append ("012032");
        modelFreeParameters += 3;
        break;
      case NUC_MATRIX_TIM3:
        name.append ("TIM3");
        matrixName.append ("012032");
        modelFreeParameters += 3;
        break;
      case NUC_MATRIX_TVMef:
        name.append ("TVMef");
        matrixName.append ("012314");
        modelFreeParameters += 4;
        break;
      case NUC_MATRIX_TVM:
        name.append ("TVM");
        matrixName.append ("012314");
        modelFreeParameters += 4;
        break;
      case NUC_MATRIX_SYM:
        name.append ("SYM");
        matrixName.append ("012345");
        modelFreeParameters += 5;
        break;
      case NUC_MATRIX_GTR:
        name.append ("GTR");
        matrixName.append ("012345");
        modelFreeParameters += 5;
        break;
      default:
        assert(0);
      }

    if (rateVariation & RateVarF)
    {
      modelFreeParameters += 3;
    }

    if (rateVariation & RateVarI)
    {
      name.append ("+I");
      /* proportion of invariable sites free parameter */
      modelFreeParameters++;
    }

    if (rateVariation & RateVarG)
    {
      name.append ("+G");
      /* alpha free parameter */
      modelFreeParameters++;
    }
  }

  NucleicModel::~NucleicModel ()
  {
    // NOTHING
  }

  void NucleicModel::setFrequencies (const double * _frequencies)
  {
    for (int i = 0; i < 4; i++)
    {
      this->frequencies[i] = _frequencies[i];
    }
  }

  NucMatrix NucleicModel::getMatrix (void) const
  {
    return matrix;
  }

  void NucleicModel::allocateRates (void)
  {
    this->rates = (double *) malloc (NUM_DNA_RATES * sizeof(double));
  }

  void NucleicModel::setRates (const double * _rates)
  {
    for (int i = 0; i < NUM_DNA_RATES; i++)
      this->rates[i] = _rates[i];
  }

  double NucleicModel::distanceTo (Model * otherModel) const
  {
    NucleicModel * other = static_cast<NucleicModel *> (otherModel);

    /* we calculate a factor K such that mse(Ra*K, Rb) is minimized */
    double numFactor = 0.0, denFactor = 0.0;
    for (int i = 0; i < 6; i++)
    {
      numFactor += rates[i] * other->rates[i];
      denFactor += rates[i] * rates[i];
    }
    double kFactor = numFactor / denFactor;

//	double rDistance = Utilities::normalizedEuclideanDistance(rates,
//			other->rates, 6);
    double rDistance = Utilities::euclideanDistance (rates, other->rates, 6,
                                                     kFactor);
    double fDistance = Utilities::euclideanDistance (frequencies,
                                                     other->frequencies, 4);
    double aDistance = fabs(alpha - other->alpha);

    double distance = wgt_r * rDistance +
                      wgt_f * fDistance +
                      wgt_a * aDistance;

    return distance;
  }

  void NucleicModel::print (ostream& cout, const char * prefix) const
  {
    cout << prefix << "Name:  " << name << endl;
    if (isOptimized ())
    {
      cout << prefix << "lnL:   " << lnL << endl;
#ifdef _IG_MODELS
      if (isPInv())
      {
        cout << prefix << "pInv:  " << pInv << endl;
      }
#endif
      if (isGamma ())
      {
        cout << prefix << "alpha: " << alpha << endl;
      }
      cout << prefix << "Frequencies:" << endl;
      cout << prefix << "  f(A) : " << frequencies[0] << endl;
      cout << prefix << "  f(C) : " << frequencies[1] << endl;
      cout << prefix << "  f(G) : " << frequencies[2] << endl;
      cout << prefix << "  f(T) : " << frequencies[3] << endl;
      cout << prefix << "Rates:" << endl;
      cout << prefix << "  R(a) : " << rates[0] << endl;
      cout << prefix << "  R(b) : " << rates[1] << endl;
      cout << prefix << "  R(c) : " << rates[2] << endl;
      cout << prefix << "  R(d) : " << rates[3] << endl;
      cout << prefix << "  R(e) : " << rates[4] << endl;
      cout << prefix << "  R(f) : " << rates[5] << endl;
      cout << prefix << "Most Likely Tree: " << tree;
    }
    else
    {
      cout << prefix << "State: Unoptimized" << endl;
    }
  }

} /* namespace partest */
