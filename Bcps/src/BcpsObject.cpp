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
 * Copyright (C) 2001-2015, Lehigh University, Yan Xu, and Ted Ralphs.       *
 * All Rights Reserved.                                                      *
 *===========================================================================*/


#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif


#include "BcpsModel.h"
#include "BcpsObject.h"


//#############################################################################


// Assignment operator 
BcpsObject & 
BcpsObject::operator = ( const BcpsObject& rhs)
{
    if (this!=&rhs) {
        objectIndex_ = rhs.objectIndex_;
	repType_ = rhs.repType_;
	intType_ = rhs.intType_;
	status_ = rhs.status_;
	lbHard_ = rhs.lbHard_;
	ubHard_ = rhs.ubHard_;
	lbSoft_ = rhs.lbSoft_;
	ubSoft_ = rhs.ubSoft_;
        hashValue_ = rhs.hashValue_;
	numInactive_ = rhs.numInactive_;
    }
    
    return *this;
}

//#############################################################################

// Returns floor and ceiling i.e. closest valid points
void 
BcpsObject::floorCeiling(double & floorValue, 
			 double & ceilingValue, 
			 double value,
			 double tolerance) const
{
    if (fabs(floor(value + 0.5) - value) > tolerance) {
	floorValue = floor(value);
    } 
    else {
	floorValue = floor(value + 0.5);
    }
    ceilingValue = floorValue + 1.0;
}

//#############################################################################

