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
 * @file PartitionElement.h
 *
 * @brief Description of a single partition
 */

#ifndef PARTITIONELEMENT_H_
#define PARTITIONELEMENT_H_

#include "util/GlobalDefs.h"
#include "indata/PllTreeManager.h"
#include "model/SelectionModel.h"
#include "model/NucleicModel.h"
#include "model/ProteicModel.h"

#include <string>
#include <vector>

namespace partest
{

  class PartitionElement
  {

  public:
    PartitionElement (t_partitionElementId id);
    t_partitionElementId getId (void)
    {
      return id;
    }
    std::string & getName (void)
    {
      return name;
    }

    int setupStructures (void);
    int destroyStructures (void);
    virtual ~PartitionElement ();

    size_t getNumberOfSites (void) const;
    size_t getNumberOfPatterns (void) const;
    size_t getNumberOfSections (void) const;
    PEsection getSection (size_t i);

    size_t getNumberOfModels (void) const;
    Model * getModel (size_t index);
    std::vector<Model *> getModels (void) const;

    double getLnL (void) const;
    SelectionModel * getBestModel (void);
    void setBestModel (SelectionModel * model);

    double * getBranchLengths (void);

    void setTagged (bool tag_status)
    {
      tag = tag_status;
    }
    bool isTagged ()
    {
      return tag;
    }

    bool isReady (void);
    bool isOptimized (void);
    double getEpsilon (void);

    TreeManager * getTreeManager (void);

    double getSampleSize (void);

    int loadData (void);
    int storeData (void);

    void print (std::ostream & out);
  private:

    bool ready;

    t_partitionElementId id;
    size_t numberOfSections;

    size_t numberOfSites;
    size_t numberOfPatterns;

    std::vector<Model *> models;
    SelectionModel * bestModel;

    std::string name, ckpname, ckphash;
    double sampleSize;

    PllTreeManager * treeManager;

    std::vector<PEsection> sections;

    bool ckpLoaded;
    bool tag;

    double * branchLengths;
  };

}
#endif /* PARTITIONELEMENT_H_ */
