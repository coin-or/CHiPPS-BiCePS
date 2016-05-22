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

#ifndef BcpsBranchStrategy_h_
#define BcpsBranchStrategy_h_

#include "BcpsBranchObject.h"

class BcpsModel;
class BcpsTreeNode;

//#############################################################################
// NOTE: Borrow ideas from COIN/Cbc.
//#############################################################################


/*!
  This class is an Abstract Base Class (ABC) for a branching strategy. A
  branching strategy specifies:
  <ul>
  <li> how to select a candidate set of branching objects
  <li> how to compare two branching objects
  </ul>

*/

class BcpsBranchStrategy {
  /// Type of branching strategy.
  int type_;
  /// Pointer to model.
  BcpsModel * model_;

  /// Following members are used to store candidate branching objects. NOTE:
  /// They are required to be cleared before starting another round of
  /// selecting.
  //@{ Number of candidate branching objects.
  int numBranchObjects_;
  /// The set of candiate branching objects.
  BcpsBranchObject ** branchObjects_;
  //@}

  /** Following members are used to store information about best
      branching object found so far.
      NOTE: They are required to be cleared before starting another
      round of selecting.*/
  //@{
  /// Index of the best branching object found so far.
  int bestIndex_;
  //@}

public:
  ///@name Constructors and destructors.
  //@{
  /// Useful Constructor.
  BcpsBranchStrategy(BcpsModel *m);
  /// Destructor.
  virtual ~BcpsBranchStrategy();
  //@}

  ///@name Get data fields
  //@{
  /// Get type.
  int type() { return type_; }
  BcpsModel * model() const { return model_; }
  int numBranchObjects() { return numBranchObjects_; }
  int bestIndex() { return bestIndex_; }
  // BcpsBranchObject * bestBranchObject() { return bestBranchObject_; }
  // BcpsBranchObject ** branchObjects() { return branchObjects_; }
  //@}

  ///@name Set data fields
  //@{
  /// Set type.
  void setType(int type) { type_ = type; }
  /// Set model.
  void setModel(BcpsModel * model) { model_ = model; }
  // set branch objects, takes ownership of input objects.
  void setBranchObjects(int num, BcpsBranchObject **obj);
  void setBranchObjects(std::vector<BcpsBranchObject*> & obj);
  void setBestIndex(int index) { bestIndex_ = index; }
  //@}

  ///@name Selecting and Creating branches.
  //@{
  /// Create a set of candidate branching objects from the given node.
  virtual int createCandBranchObjects(BcpsTreeNode * node) = 0;
  /// Compare current to other, return 1 if current is better, 0 otherwise
  virtual int betterBranchObject(BcpsBranchObject const * current,
                                 BcpsBranchObject const * other) = 0;
  /// Compares branching objects stored. Returns the best branching object
  /// with the smallest index.
  virtual BcpsBranchObject * bestBranchObject();
  /// Clear branch objects stored.
  void clearBranchObjects();
  //@}

private:
  /// Disable default constructor.
  BcpsBranchStrategy();
  /// Disable copy constructor.
  BcpsBranchStrategy(BcpsBranchStrategy const & other);
  /// Disable copy assignment operator.
  BcpsBranchStrategy & operator=(const BcpsBranchStrategy & rhs);

};

#endif
