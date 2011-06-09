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
 * Copyright (C) 2001-2011, Lehigh University, Yan Xu, and Ted Ralphs.       *
 * All Rights Reserved.                                                      *
 *===========================================================================*/

//#############################################################################

#include "BcpsSolution.h"

//#############################################################################

BcpsSolution *
BcpsSolution::selectNonzeros(const double etol) const
{
    BcpsSolution *sol = NULL;
    
    return sol;
}

//#############################################################################

BcpsSolution* 
BcpsSolution::selectFractional(const double etol) const
{
    BcpsSolution *sol = NULL;
    
    return sol;
}

//#############################################################################

/** Pack Bcps part of solution into an encoded objects. */
AlpsReturnStatus 
BcpsSolution::encodeBcps(AlpsEncoded *encoded) const
{
    AlpsReturnStatus status = AlpsReturnStatusOk;
    encoded->writeRep(size_);
    encoded->writeRep(values_, size_);
    encoded->writeRep(quality_);

    return status;
}

//#############################################################################

/** Unpack Bcps part of solution from an encoded objects. */
AlpsReturnStatus 
BcpsSolution::decodeBcps(AlpsEncoded & encoded)
{
    AlpsReturnStatus status = AlpsReturnStatusOk;
    
    encoded.readRep(size_);
    encoded.readRep(values_, size_);
    encoded.readRep(quality_);

    return status;
}



//#############################################################################


//#############################################################################
