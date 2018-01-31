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


#include "CoinHelperFunctions.hpp"

#include "BlisBranchObjectInt.h"


//#############################################################################

// Copy constructor
BlisBranchObjectInt::
BlisBranchObjectInt(const BlisBranchObjectInt & rhs)
    :
    BcpsBranchObject(rhs)
{
    ubDownBranch_ = rhs.ubDownBranch();
    lbUpBranch_ = rhs.lbUpBranch();
}

BlisBranchObjectInt::
BlisBranchObjectInt(const BcpsBranchObject * rhs)
  : AlpsKnowledge(rhs->getType(), rhs->broker_),
    BcpsBranchObject(*rhs)
{
    BlisBranchObjectInt const * blis_other =
        dynamic_cast<BlisBranchObjectInt const *>(rhs);
    if (blis_other==NULL) {
        std::cerr << "Fatal error!" << std::endl;
        throw std::exception();
    }
    ubDownBranch_ = blis_other->ubDownBranch();
    lbUpBranch_ = blis_other->lbUpBranch();
}

// Assignment operator
BlisBranchObjectInt &
BlisBranchObjectInt::operator=(const BlisBranchObjectInt& rhs)
{
    if (this != &rhs) {
	BcpsBranchObject::operator=(rhs);
        ubDownBranch_ = rhs.ubDownBranch();
        lbUpBranch_ = rhs.lbUpBranch();
    }
    return *this;
}

/// The number of branch arms created for this branch object.
int BlisBranchObjectInt::numBranches() const {
  std::cerr << "Not implemented yet!" << std::endl;
  throw std::exception();
  return 0;
}


/// The number of branch arms left to be evaluated.
int BlisBranchObjectInt::numBranchesLeft() const {
  std::cerr << "Not implemented yet!" << std::endl;
  throw std::exception();
  return 0;
}


/// Spit out a branch and, update this or superclass fields if necessary.
double BlisBranchObjectInt::branch(bool normalBranch) {
  std::cerr << "Not implemented yet!" << std::endl;
  throw std::exception();
  return 0.0;
}


AlpsReturnStatus BlisBranchObjectInt::encode(AlpsEncoded * encoded) const {
  assert(encoded);
  AlpsReturnStatus status;
  status = BcpsBranchObject::encode(encoded);
  assert(status==AlpsReturnStatusOk);
  encoded->writeRep(ubDownBranch_);
  return status;
}

AlpsKnowledge * BlisBranchObjectInt::decode(AlpsEncoded &encoded) const {
  // create a new object with default values,
  // Bcps decode will decode right values into them.
  AlpsKnowledge * new_bo = new BlisBranchObjectInt(-1, 0.0, 0.0);
  new_bo->decodeToSelf(encoded);
  return new_bo;
}

/// Decode the given AlpsEncoded object into this.
AlpsReturnStatus BlisBranchObjectInt::decodeToSelf(AlpsEncoded & encoded) {
  AlpsReturnStatus status;
  // decode Bcps part.
  status = BcpsBranchObject::decodeToSelf(encoded);
  assert(status=AlpsReturnStatusOk);
  // decode fields of BlisBranchObjectInt
  encoded.readRep(ubDownBranch_);
  encoded.readRep(lbUpBranch_);
  return status;
}
