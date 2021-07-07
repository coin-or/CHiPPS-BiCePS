/*===========================================================================*
 * This file is part of the Bcps Linear Solver (BLIS).                       *
 *                                                                           *
 * BLIS is distributed under the Eclipse Public License as part of the       *
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
// Borrow ideas from COIN/Cbc
//#############################################################################

#ifndef BlisBranchObjectInt_h_
#define BlisBranchObjectInt_h_


#include "BcpsBranchObject.h"

#include "BlisModel.h"


//#############################################################################


class BlisBranchObjectInt : public BcpsBranchObject {

    /// upper bound of the down branch
    double ubDownBranch_;
    /// lower bound of the up branch
    double lbUpBranch_;

 public:
    /** Construct a branching object, which branching on variable varInd.
        \param index     the index of integer variable in object set
        \param score   the double score/goodness
        \param direction  the direction of first branching: 1(up), -1(down)
        \param value      the fractional solution value of variable varInd
    */
    BlisBranchObjectInt(int index,
                        double score,
                        double value)
        :
        BcpsBranchObject(BLIS_BO_INT, index, score, value)
        {
          ubDownBranch_ = floor(value);
          lbUpBranch_ = ceil(value);
        }

    /** Copy constructor. */
    BlisBranchObjectInt(const BlisBranchObjectInt &);
    BlisBranchObjectInt(const BcpsBranchObject *);
    /** Assignment operator. */
    BlisBranchObjectInt & operator = (const BlisBranchObjectInt& rhs);
    /** Destructor. */
    virtual ~BlisBranchObjectInt() {}

    ///@name Virtual functions inherited from BcpsBranchObject
    /// The number of branch arms created for this branch object.
    virtual int numBranches() const;
    /// The number of branch arms left to be evaluated.
    virtual int numBranchesLeft() const;
    /// Spit out a branch and, update this or superclass fields if necessary.
    virtual double branch(bool normalBranch = false);
    //@}

    ///@name Bound getting functions.
    //@{
    /// Get upper bound of the down branch.
    double ubDownBranch() const { return ubDownBranch_; }
    /// Get lower bound of the up branch.
    double lbUpBranch() const { return lbUpBranch_; }
    //@}

    ///@name Encode and Decode functions
    ///@{
    using AlpsKnowledge::encode;
    /// Encode the content of this into the given AlpsEncoded object.
    virtual AlpsReturnStatus encode(AlpsEncoded * encoded) const;
    /// Decode the given AlpsEncoded object into a new AlpsKnowledge object and
    /// return a pointer to it.
    virtual AlpsKnowledge * decode(AlpsEncoded & encoded) const;
    /// Decode the given AlpsEncoded object into this.
    virtual AlpsReturnStatus decodeToSelf(AlpsEncoded & encoded);
    ///@}

 private:
  /// Disable default constructor.
  BlisBranchObjectInt();
};

#endif
