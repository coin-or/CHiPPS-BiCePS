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
 * Copyright (C) 2001-2013, Lehigh University, Yan Xu, and Ted Ralphs.       *
 * All Rights Reserved.                                                      *
 *===========================================================================*/


#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include "BcpsModel.h"
#include "BcpsBranchObject.h"


//#############################################################################

// Copy constructor 
BcpsBranchObject::BcpsBranchObject(const BcpsBranchObject &rhs)
{
    model_ = rhs.model_;
    objectIndex_ = rhs.objectIndex_;
    direction_ = rhs.direction_;
    value_ = rhs.value_;
    numBranchesLeft_ = rhs.numBranchesLeft_;
}

//#############################################################################

// Assignment operator 
BcpsBranchObject & 
BcpsBranchObject::operator = ( const BcpsBranchObject& rhs)
{
    if (this != &rhs) {
	model_ = rhs.model_;
	objectIndex_ = rhs.objectIndex_;
	direction_ = rhs.direction_;
	value_ = rhs.value_;
	numBranchesLeft_ = rhs.numBranchesLeft_;
    }

    return *this;
}

//#############################################################################

