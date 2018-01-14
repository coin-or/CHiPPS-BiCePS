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
 * Copyright (C) 2001-2018, Lehigh University, Yan Xu, and Ted Ralphs.       *
 * All Rights Reserved.                                                      *
 *===========================================================================*/

#ifndef Bcps_h_
#define Bcps_h_

//! \page handle BcpsMainpage

/*! \mainpage

  BiCePS (Branch, Constrain, and Price Software) is an abstract library build
  on top of Alps. Alps does not have any assumption about the tree search
  algorithm. Well BiCePS does. It assumes we are building Branch and Cut type
  algorithm. It provides the means for this.

  BiCePS assuming branch and bound is also very flexible. It only assumes
  this, it does not assume much, hence a flexible library. It does not have an
  assumption about the underlying optimization problem.

  BLIS (BiCePS Linear Integer Solver) is built on top of BiCePS to solve MILP
  problems. Similarly DisCO (Discrete Conic Optimization) solver is build
  on top of BiCePS to solve second order conic optimization problems.

 */

enum BcpsReturnStatus {
   BcpsReturnStatusOk = 0,
   BcpsReturnStatusErr
};

enum BcpsKnowledgeType{
   BcpsKnowledgeTypeConstraint  = 11,
   BcpsKnowledgeTypeVariable  = 12
};

enum BcpsValidRegion{
    BcpsValidLocal = 0,
    BcpsValidGlobal
};

enum BcpsSubproblemStatus {
  BcpsSubproblemStatusOptimal,
  BcpsSubproblemStatusAbandoned,
  BcpsSubproblemStatusPrimalInfeasible,
  BcpsSubproblemStatusDualInfeasible,
  BcpsSubproblemStatusPrimalObjLim,
  BcpsSubproblemStatusDualObjLim,
  BcpsSubproblemStatusIterLim,
  BcpsSubproblemStatusUnknown
};

/** Node branch direction, is it a left node or right */
enum BcpsNodeBranchDir {
    BcpsNodeBranchDirectionDown = 0,
    BcpsNodeBranchDirectionUp
};

#endif
