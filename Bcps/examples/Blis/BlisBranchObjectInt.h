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


#include "BcpsBranchObject.h"

#include "BlisModel.h"


//#############################################################################


class BlisBranchObjectInt : public BcpsBranchObject {

 protected:

    int direction_;

    /** Down_[0]: the lower bound of down arm;
        Down_[1]: the upper bound of down arm; */
    double down_[2];

    /** Up_[0]: the lower bound of upper arm;
        Up_[1]: the upper bound of upper arm; */
    double up_[2];

 public:
    /** Construct a branching object, which branching on variable varInd.
        \param index     the index of integer variable in object set
        \param score   the double score/goodness
        \param direction  the direction of first branching: 1(up), -1(down)
        \param value      the fractional solution value of variable varInd
    */
    BlisBranchObjectInt(int index,
                        double score,
                        int direction,
                        double value)
        :
        BcpsBranchObject(BLIS_BO_INT, index, score, value)
        {
            direction_ = direction;
            BlisModel * model = dynamic_cast <BlisModel*> (broker()->getModel());
            //double value = model->solver()->getColSolution()[index];
            down_[0] = model->solver()->getColLower()[index];
            down_[1] = floor(value);
            up_[0] = ceil(value);
            up_[1] = model->getColUpper()[index];
        }

    /** Create a degenerate branching object.
        Specifies a `one-direction branch'. Calling branch() for this
        object will always result in lowerValue <= x <= upperValue.
        Used to fix a variable when lowerValue = upperValue.
    */
    BlisBranchObjectInt(int index,
                        double score,
                        int direction,
                        double lowerValue,
                        double upperValue)
        :
        BcpsBranchObject(BLIS_BO_INT, index, score, lowerValue)
        {
            direction_ = direction;
            down_[0] = lowerValue;
            down_[1] = upperValue;
            up_[0] = lowerValue;
            up_[1] = upperValue;
        }

    /** Copy constructor. */
    BlisBranchObjectInt(const BlisBranchObjectInt &);

    /** Assignment operator. */
    BlisBranchObjectInt & operator = (const BlisBranchObjectInt& rhs);

    ///@name Pure virtual functions inherited.
    //@{
    /// The number of branch arms created for this branch object.
    virtual int numBranches() const;
    /// The number of branch arms left to be evaluated.
    virtual int numBranchesLeft() const;
    //@}


    /** Clone. */
    virtual BcpsBranchObject * clone() const {
        return (new BlisBranchObjectInt(*this));
    }

    /** Destructor. */
    virtual ~BlisBranchObjectInt() {}

    /** Set the bounds for the variable according to the current arm
	of the branch and advances the object state to the next arm.
	Returns change in guessed objective on next branch. */
    virtual double branch(bool normalBranch = false);

    /** \brief Print something about branch - only if log level high. */
    virtual void print(bool normalBranch);

    /** Get down arm bounds. */
    const double *getDown() const { return down_; }

    /** Get upper arm bounds. */
    const double *getUp() const { return up_; }

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

};
