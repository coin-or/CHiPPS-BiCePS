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
    direction_ = rhs.direction_;
    down_[0] = rhs.down_[0];
    down_[1] = rhs.down_[1];
    up_[0] = rhs.up_[0];
    up_[1] = rhs.up_[1];
}

//#############################################################################

// Assignment operator
BlisBranchObjectInt &
BlisBranchObjectInt::operator=(const BlisBranchObjectInt& rhs)
{
    if (this != &rhs) {
	BcpsBranchObject::operator=(rhs);
        direction_ = rhs.direction_;
	down_[0] = rhs.down_[0];
	down_[1] = rhs.down_[1];
	up_[0] = rhs.up_[0];
	up_[1] = rhs.up_[1];
    }
    return *this;
}

//#############################################################################

/*
  Perform a branch by adjusting the bounds of the specified variable.
  Note that each arm of the branch advances the object to the next arm by
  advancing the value of direction_.

  Providing new values for the variable's lower and upper bounds for each
  branching direction gives a little bit of additional flexibility and will
  be easily extensible to multi-direction branching.
  Returns change in guessed objective on next branch
*/
double
BlisBranchObjectInt::branch(bool normalBranch)
{
    BlisModel * model = dynamic_cast<BlisModel *>(broker()->getModel());

    int iColumn = model->getIntVars()[index()];

    if (direction_<0) {
#ifdef BLIS_DEBUG_MORE
	{
	    double olb, oub ;
	    olb = model->solver()->getColLower()[iColumn];
	    oub = model->solver()->getColUpper()[iColumn];
	    printf("branching down on var %d: [%g,%g] => [%g,%g]\n",
		   iColumn,olb,oub,down_[0],down_[1]);
	}
#endif
	model->solver()->setColLower(iColumn, down_[0]);
	model->solver()->setColUpper(iColumn, down_[1]);
	direction_ = 1;
    }
    else {
#ifdef BLIS_DEBUG_MORE
	{
	    double olb, oub ;
	    olb = model->solver()->getColLower()[iColumn];
	    oub = model->solver()->getColUpper()[iColumn];
	    printf("branching up on var %d: [%g,%g] => [%g,%g]\n",
		   iColumn,olb,oub,up_[0],up_[1]);
	}
#endif
	model->solver()->setColLower(iColumn, up_[0]);
	model->solver()->setColUpper(iColumn, up_[1]);
	direction_ = -1;	  // Swap direction
    }

    return 0.0;
}

//#############################################################################

// Print what would happen
void
BlisBranchObjectInt::print(bool normalBranch)
{
    BlisModel * model = dynamic_cast<BlisModel*>(broker()->getModel());
    int iColumn = model->getIntVars()[index()];
    double olb, oub ;

    if (direction_ < 0) {
        olb = model->solver()->getColLower()[iColumn] ;
        oub = model->solver()->getColUpper()[iColumn] ;
        printf("BlisInteger would branch down on var %d: [%g,%g] => [%g,%g]\n",
               iColumn,olb,oub,down_[0],down_[1]);

    }
    else {
        olb = model->solver()->getColLower()[iColumn] ;
        oub = model->solver()->getColUpper()[iColumn] ;
        printf("BlisInteger would branch up on var %d: [%g,%g] => [%g,%g]\n",
               iColumn,olb,oub,up_[0],up_[1]);
    }
}

//#############################################################################



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
// double BlisBranchObjectInt::branch(bool normalBranch = false) {
//   std::cerr << "Not implemented yet!" << std::endl;
//   throw std::exception();
//   return 0.0;
// }


AlpsReturnStatus BlisBranchObjectInt::encode(AlpsEncoded * encoded) const {
  assert(encoded);
  AlpsReturnStatus status;
  status = BcpsBranchObject::encode(encoded);
  assert(status==AlpsReturnStatusOk);
  encoded->writeRep(down_[0]);
  encoded->writeRep(down_[1]);
  encoded->writeRep(up_[0]);
  encoded->writeRep(up_[1]);
  return status;
}

AlpsKnowledge * BlisBranchObjectInt::decode(AlpsEncoded &encoded) const {
  // create a new object with default values,
  // Bcps decode will decode right values into them.
  AlpsKnowledge * new_bo = new BlisBranchObjectInt(-1, 0.0, -1, 0.0);
  new_bo->decodeToSelf(encoded);
  return new_bo;
}

/// Decode the given AlpsEncoded object into this.
AlpsReturnStatus BlisBranchObjectInt::decodeToSelf(AlpsEncoded & encoded) {
  AlpsReturnStatus status;
  // decode Bcps part.
  status = BcpsBranchObject::decodeToSelf(encoded);
  assert(status=AlpsReturnStatusOk);
  // decode fields of DcoBranchObject
  encoded.readRep(down_[0]);
  encoded.readRep(down_[1]);
  encoded.readRep(up_[0]);
  encoded.readRep(up_[1]);
  return status;
}
