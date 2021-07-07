/*===========================================================================*
 * This file is part of the Bcps Linear Solver (BLIS).                       *
 *                                                                           *
 * ALPS is distributed under the Eclipse Public License as part of the       *
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
 *                                                                           *
 * Copyright (C) 2001-2017, Lehigh University, Yan Xu, and Ted Ralphs.       *
 * All Rights Reserved.                                                      *
 *===========================================================================*/


//#############################################################################
// NOTE: Borrow ideas from COIN/Cbc
//#############################################################################


#ifndef BlisBranchStrategyPseudo_h_
#define BlisBranchStrategyPseudo_h_

#include "BcpsBranchObject.h"
#include "BcpsBranchStrategy.h"
#include "BlisModel.h"

class BlisTreeNode;


/** Blis branching strategy default class
    This class implements a simple default algorithm, betterBranchObject(),
    for choosing a branching variable. */
class BlisBranchStrategyPseudo : public BcpsBranchStrategy {
    /// score factor used. See class documentation.
    double score_factor_;
    ///@name Statistics
    //@{
    /// number of observations for each integer variable
    int * down_num_;
    int * up_num_;
    /// estimated improvement in the objective value per change in each variable
    /// derivative_[i] is the average of all observations for variable i.
    /// these are \f$ \varphi \f$ variables in the documentation
    double * down_derivative_;
    double * up_derivative_;
    /// reverse map of relaxed columns, rev_relaxed_[index] gives the index of
    /// the varaible in relaxed columns array.
    std::map<int,int> rev_relaxed_;
    /// update scores of the stored branch objects.
    void update_statistics(BlisTreeNode * node);

 public:
  BlisBranchStrategyPseudo(BlisModel * model);
  virtual ~BlisBranchStrategyPseudo();
  virtual int createCandBranchObjects(BcpsTreeNode * node);
  /// Compare current to other, return 1 if current is better, 0 otherwise
  virtual int betterBranchObject(BcpsBranchObject const * current,
                                 BcpsBranchObject const * other);
private:
  /// Disable default constructor.
  BlisBranchStrategyPseudo();
  /// Disable copy constructor.
  BlisBranchStrategyPseudo(BlisBranchStrategyPseudo const & other);
  /// Disable copy assignment operator.
  BlisBranchStrategyPseudo & operator=(BlisBranchStrategyPseudo const & rhs);
};

#endif
