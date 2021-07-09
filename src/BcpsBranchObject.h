/*===========================================================================*
 * This file is part of the Branch, Constrain and Price Software (BiCePS)    *
 *                                                                           *
 * BiCePS is distributed under the Eclipse Public License as part of the     *
 * COIN-OR repository (http://www.coin-or.org).                              *
 *                                                                           *
 * Authors:                                                                  *
 *                                                                           *
 *          Yan Xu, Lehigh University                                        *
 *          Aykut Bulut, Lehigh University                                   *
 *          Ted Ralphs, Lehigh University                                    *
 *                                                                           *
 * Conceptual Design:                                                        *
 *                                                                           *
 *          Yan Xu, Lehigh University                                        *
 *          Ted Ralphs, Lehigh University                                    *
 *          Laszlo Ladanyi, IBM T.J. Watson Research Center                  *
 *          Matthew Saltzman, Clemson University                             *
 *                                                                           *
 * Copyright (C) 2001-2019, Lehigh University, Yan Xu, and Ted Ralphs.       *
 * All Rights Reserved.                                                      *
 *===========================================================================*/


//#############################################################################
// Borrow ideas from COIN/Cbc
//#############################################################################

#ifndef BcpsBranchObject_h_
#define BcpsBranchObject_h_

#include "BcpsConfig.h"
#include "BcpsModel.h"

#include "Alps.h"
#include "AlpsEncoded.h"

/*!

  # BcpsBranchObject

  BcpsBranchObject contains the member data required when branching a node.
  Branching is creating new nodes from current one. This is an abstract class
  that represents data needed for the branching process in the most general
  form. Branching objects can be simple integer variables or more
  complicated objects like SOS.

  Branching in general form is not necessarily binary. Branching process may
  create more than two nodes.

  This interface lets creating arbitrary number of nodes. score() is updated
  while the new branches are created. next() returns an integer that
  identifies the new branch that will be created. This can be thought as
  direction of the branch in the binary case. The indexing starts at 0. if
  next() is 5, this indicates the branching object already created 5 nodes.

  BcpsBranchObject has no idea what next() returns or means. The super class
  will put a meaning to it. For example, in a binary branching object
  implementing this ABC, 0 might mean down node and 1 might mean up node.
  When next() is 2, it means the branch object created both down and up nodes.

  score() is the quality of this object. It is used to compare this branching
  object to others.

  # Notes(aykut):

  I removed up and down score data fields. I added a single field to store
  the score. Up and down assumes binary branching and this is a general ABC.
  Up and down should be implemented in the super class in case needed.

*/

class BCPSLIB_EXPORT BcpsBranchObject: virtual public AlpsKnowledge {
  /// Type of branching. This will be set by the application built on top of
  /// Bcps.
  int type_;
  /// Branch object index. The index is not necessarily the same as variable
  /// index. It will be set by the user to identify this object.
  int index_;
  /// Quality/Goodness of this object. It is  used when comparing two branching
  /// enities. Derived class can add more metrics like this.
  double score_;
  /// Current branching value. When branching on integer variables, it can be
  /// the fractional value branched. Its meaning will be defined by the super
  /// class.
  double value_;

public:
  ///@name Constructors and Destructor.
  //@{
  /// Constructor.
  BcpsBranchObject(int type, int index, int score);
  /// Constructor.
  BcpsBranchObject(int type, int index, double score, double value);
  /// Copy constructor.
  BcpsBranchObject(BcpsBranchObject const & other);
  /// Copy assignment operator
  BcpsBranchObject & operator=(BcpsBranchObject const & rhs);
  /// Destructor.
  virtual ~BcpsBranchObject() { /* Do nothing */}
  //@}

  ///@name Get functions.
  //@{
  /// Get type.
  int type() const { return type_; }
  /// Get index.
  int index() const { return index_; }
  /// Return score.
  double score() const { return score_; }
  /// Return object branching value.
  double value() const { return value_; }
  //@}

  ///@name Set functions.
  //@{
  /// Set score.
  void setScore(double score) { score_ = score; }
  //@}

  ///@name Pure virtual functions.
  /// The number of branch arms created for this branch object.
  virtual int numBranches() const = 0;
  /// The number of branch arms left to be evaluated.
  virtual int numBranchesLeft() const = 0;
  /// Spit out a branch and, update this or superclass fields if necessary.
  virtual double branch(bool normalBranch = false) = 0;

  /// Encode the content of this into the given AlpsEncoded object.
  virtual AlpsReturnStatus encode(AlpsEncoded * encoded) const;
  /// Decode the given AlpsEncoded object into this.
  virtual AlpsReturnStatus decodeToSelf(AlpsEncoded & encoded);
  //@}

private:
  /// Disable default constructor.
  BcpsBranchObject();
};

#endif
