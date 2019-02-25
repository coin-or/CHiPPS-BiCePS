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

#ifndef BlisNodeDesc_h_
#define BlisNodeDesc_h_

//#############################################################################

#include "CoinWarmStartBasis.hpp"

#include "AlpsNodeDesc.h"
#include "BcpsNodeDesc.h"

#include "BlisHelp.h"
#include "BlisModel.h"

//#############################################################################


class BlisNodeDesc : public BcpsNodeDesc {

 private:

    /** Branched direction to create it. For updating pseudocost. */
    int branchedDir_;

    /** Branched object index to create it. For updating pseudocost. */
    int branchedInd_;

    /** Branched value to create it. For updating pseudocost. */
    double branchedVal_;

    /** Warm start. */
    CoinWarmStartBasis *basis_;

 public:

    /** Default constructor. */
    BlisNodeDesc() :
        BcpsNodeDesc(),
        branchedDir_(0),
        branchedInd_(-1),
        branchedVal_(0.0),
	basis_(NULL)
        {}

    /** Useful constructor. */
    BlisNodeDesc(BlisModel* m)
	:
	BcpsNodeDesc(m),
        branchedDir_(0),
        branchedInd_(-1),
        branchedVal_(0.0),
	basis_(NULL)
	{}

    /** Destructor. */
    virtual ~BlisNodeDesc() { delete basis_; }

    /** Set basis. */
    void setBasis(CoinWarmStartBasis *&ws) {
        if (basis_) { delete basis_; }
        basis_= ws;
        ws = NULL;
    }

    /** Get warm start basis. */
    CoinWarmStartBasis * getBasis() const { return basis_; }

    /** Set branching direction. */
    void setBranchedDir(int d) { branchedDir_ = d; }

    /** Get branching direction. */
    int getBranchedDir() const { return branchedDir_; }

    /** Set branching object index. */
    void setBranchedInd(int d) { branchedInd_ = d; }

    /** Get branching object index. */
    int getBranchedInd() const { return branchedInd_; }

    /** Set branching value. */
    void setBranchedVal(double d) { branchedVal_ = d; }

    /** Get branching direction. */
    double getBranchedVal() const { return branchedVal_; }

    ///@name Encode and Decode functions
    ///@{
    using AlpsKnowledge::encode;
    /// Encode the content of this into the given AlpsEncoded object.
    virtual AlpsReturnStatus encode(AlpsEncoded * encoded) const {
        AlpsReturnStatus status;
        status = AlpsNodeDesc::encode(encoded);
        status = BcpsNodeDesc::encodeBcps(encoded);
        assert(status==AlpsReturnStatusOk);
        encoded->writeRep(branchedDir_);
        encoded->writeRep(branchedInd_);
        encoded->writeRep(branchedVal_);
        // Encode basis if available
        int available = 0;
        if (basis_) {
            available = 1;
            encoded->writeRep(available);
            int numCols = basis_->getNumStructural();
            int numRows = basis_->getNumArtificial();
            encoded->writeRep(numCols);
            encoded->writeRep(numRows);
            // Pack structural.
            int nint = (basis_->getNumStructural() + 15) >> 4;
            encoded->writeRep(basis_->getStructuralStatus(), nint * 4);
            // Pack artificial.
            nint = (basis_->getNumArtificial() + 15) >> 4;
            encoded->writeRep(basis_->getArtificialStatus(), nint * 4);
        }
        else {
            encoded->writeRep(available);
        }
        return status;
    }
    /// Decode the given AlpsEncoded object into a new AlpsKnowledge object and
    /// return a pointer to it.
    virtual AlpsKnowledge * decode(AlpsEncoded & encoded) const {
        // get pointers for message logging
        AlpsReturnStatus status;
        BlisNodeDesc * new_desc = new BlisNodeDesc();
        //status = new_desc->decodeToSelf(encoded);
        //assert(status==AlpsReturnStatusOk);
        new_desc->decodeToSelf(encoded);
        return new_desc;
    }

    /// Decode the given AlpsEncoded object into this.
    virtual AlpsReturnStatus decodeToSelf(AlpsEncoded & encoded) {
        AlpsReturnStatus status = AlpsReturnStatusOk;
        status = AlpsNodeDesc::decodeToSelf(encoded);
        // todo(aykut) rename this function in Bcps level?
        status = BcpsNodeDesc::decodeBcps(encoded);
        assert(status==AlpsReturnStatusOk);
        encoded.readRep(branchedDir_);
        encoded.readRep(branchedInd_);
        encoded.readRep(branchedVal_);
        // decode basis if available
        int available;
        encoded.readRep(available);
        if (available==1) {
            if (basis_) {
                delete basis_;
            }
            int numCols;
            int numRows;
            encoded.readRep(numCols);
            encoded.readRep(numRows);
            int tempInt;
            // Structural
            int nint = (numCols + 15) >> 4;
            char * structuralStatus = new char[4 * nint];
            encoded.readRep(structuralStatus, tempInt);
            assert(tempInt == nint*4);
            // Artificial
            nint = (numRows + 15) >> 4;
            char * artificialStatus = new char[4 * nint];
            encoded.readRep(artificialStatus, tempInt);
            assert(tempInt == nint*4);
            basis_ = new CoinWarmStartBasis();
            if (!basis_) {
                throw CoinError("Out of memory", "BlisDecodeWarmStart", "HELP");
            }
            basis_->assignBasisStatus(numCols, numRows,
                              structuralStatus, artificialStatus);
            assert(!structuralStatus);
            assert(!artificialStatus);
        }
        else {
            basis_ = NULL;
        }
        return status;
    }
    ///@}
};
#endif
