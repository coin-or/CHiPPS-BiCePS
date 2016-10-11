/*===========================================================================*
 * This file is part of the Branch, Constrain and Price Software (BiCePS)    *
 *                                                                           *
 * BiCePS is distributed under the Eclipse Public License as part of the     *
 * COIN-OR repository (http://www.coin-or.org).                              *
 *                                                                           *
 * Authors:                                                                  *
 *                                                                           *
 *          Yan Xu, Lehigh University                                        *
 *          Ted Ralphs, Lehigh University                                    *
 *                                                                           *
 * Conceptual Design:                                                        *
 *                                                                           *
 *          Yan Xu, Lehigh University                                        *
 *          Ted Ralphs, Lehigh University                                    *
 *          Laszlo Ladanyi, IBM T.J. Watson Research Center                  *
 *          Matthew Saltzman, Clemson University                             *
 *                                                                           *
 * Copyright (C) 2001-2017, Lehigh University, Yan Xu, and Ted Ralphs.       *
 * All Rights Reserved.                                                      *
 *===========================================================================*/


#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include <iostream>

#include "BcpsBranchStrategy.h"
#include "BcpsModel.h"

BcpsBranchStrategy::BcpsBranchStrategy(BcpsModel * model)
  : model_(model), numBranchObjects_(0), branchObjects_(NULL), bestIndex_(-1) {
}

BcpsBranchStrategy::~BcpsBranchStrategy() {
  for (int k=0; k<numBranchObjects_; ++k) {
    delete branchObjects_[k];
  }
  delete[] branchObjects_;
}

void BcpsBranchStrategy::setBranchObjects(int num, BcpsBranchObject **& obj) {
  // clear members
  clearBranchObjects();
  // take ownership of objects
  branchObjects_ = obj;
  obj = NULL;
  numBranchObjects_ = num;
  // reset index of the best branch object
  bestIndex_ = -1;
  // set best branch index
  bestBranchObject();
}

void
BcpsBranchStrategy::setBranchObjects(std::vector<BcpsBranchObject*> & obj) {
  // clear members
  clearBranchObjects();
  // take ownership of objects
  numBranchObjects_ = obj.size();
  branchObjects_ = new BcpsBranchObject*[numBranchObjects_];
  std::copy(obj.begin(), obj.end(), branchObjects_);
  for(int i=0; i<numBranchObjects_; ++i) {
    obj[i] = NULL;
  }
  // reset index of the best branch object
  bestIndex_ = -1;
  // set best branch index
  bestBranchObject();
}

/* Compare N branch objects and identify bestObject_. Return index
   of best and sets way of branch bestObject_. */
BcpsBranchObject * BcpsBranchStrategy::bestBranchObject() {
  if (numBranchObjects_==0) {
    // there are no branch objects
    std::cerr << "No branch objects in the branch strategy!" << std::endl;
    std::cerr << "This might mean all columns are feasible!" << std::endl;
    throw std::exception();
  }
  if (bestIndex_!=-1) {
    return branchObjects_[bestIndex_];
  }
  bestIndex_ = 0;
  for (int i=1; i<numBranchObjects_; ++i) {
    BcpsBranchObject * current = branchObjects_[i];
    int curr_better = betterBranchObject(current, branchObjects_[bestIndex_]);
    if (curr_better) {
      bestIndex_ = i;
    }
  }
  return branchObjects_[bestIndex_];
}

void BcpsBranchStrategy::clearBranchObjects() {
  for (int i=0; i<numBranchObjects_; ++i) {
    delete branchObjects_[i];
  }
  if (branchObjects_) {
    delete[] branchObjects_;
    branchObjects_ = NULL;
  }
  numBranchObjects_ = 0;
  // reset index of the best branch object.
  bestIndex_ = -1;
}
