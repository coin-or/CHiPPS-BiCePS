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


/** Blis branching strategy default class
    This class implements a simple default algorithm, betterBranchObject(),
    for choosing a branching variable. */
class BlisBranchStrategyPseudo : public BcpsBranchStrategy {

 private:
    /** Illegal Assignment operator.*/
    BlisBranchStrategyPseudo& operator=(const BlisBranchStrategyPseudo& rhs);
    BlisBranchStrategyPseudo(const BlisBranchStrategyPseudo& rhs);
    BlisBranchStrategyPseudo();

    int reliability_;

 public:

    /** Useful Constructor. */
    BlisBranchStrategyPseudo(BlisModel *model, int rel)
        :
        BcpsBranchStrategy(model),
        reliability_(rel)
        {}

    /** Destructor. */
    virtual ~BlisBranchStrategyPseudo() {}

    /** Set reliability. */
    void setReliability(int rel) { reliability_ = rel; }

    /** Clone a brancing strategy. */
    /* virtual BcpsBranchStrategy * clone() const { */
    /*     return new BlisBranchStrategyPseudo(*this); */
    /* } */

    /** Compare branching object thisOne to bestSoFar. If thisOne is better
	than bestObject, return branching direction(1 or -1), otherwise
	return 0.
	If bestSorFar is NULL, then always return branching direction(1 or -1).
    */
    virtual int betterBranchObject(BcpsBranchObject const * current,
				   BcpsBranchObject const * other);

    /** Create a set of candidate branching objects. */
    virtual int createCandBranchObjects(BcpsTreeNode * node);
};

#endif
