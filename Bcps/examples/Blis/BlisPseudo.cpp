/*===========================================================================*
 * This file is part of the Bcps Linear Solver (BLIS).                       *
 *                                                                           *
 * ALPS is distributed under the Common Public License as part of the        *
 * COIN-OR repository (http://www.coin-or.org).                              *
 *                                                                           *
 * Authors: Yan Xu, Lehigh University                                       *
 *          Ted Ralphs, Lehigh University                                    *
 *          Laszlo Ladanyi, IBM T.J. Watson Research Center                  *
 *          Matthew Saltzman, Clemson University                             *
 *                                                                           * 
 *                                                                           *
 * Copyright (C) 2001-2005, International Business Machines                  *
 * Corporation, Lehigh University, Yan Xu, Ted Ralphs, Matthew Salzman and   *
 * others. All Rights Reserved.                                              *
 *===========================================================================*/

#include <cmath>
#include <cassert>

#include "Alps.h"

#include "BlisPseudo.h"

//#############################################################################

void 
BlisPseudocost::update(const int dir,
                       const double parentObjValue,
                       const double objValue,
                       const double solValue)
{
    
    double objDiff = objValue - parentObjValue;    

#ifdef BLIS_DEBUG
    assert(objDiff/(1.0+objValue) >= -1.0e-5);
#endif

    update(dir, objDiff, solValue);    
}

//#############################################################################

void
BlisPseudocost::update(int dir,
                       double objDiff,
                       double solValue)
{
    double fraction;
    double cost;
    
#ifdef BLIS_DEBUG
    if (objDiff < -1.0e-1) {
        std::cout << "objDiff=" << objDiff
                  << ", dir="<< dir << ", x=" << solValue <<std::endl;
        assert(0);
    }
#endif

    if (dir == 1) {
        fraction = ceil(solValue) - solValue;
        if (fraction >= 1.0e-5) {
            cost = objDiff / (fraction + 1.0e-9);
            upCost_ = (upCost_ * upCount_ + cost) / (1 + upCount_);
            ++upCount_;
        }
        else {
#ifdef BLIS_DEBUG
            ALPS_PRINTF("WARNING: small fraction %.10g, shouldn't happen.\n",
                        fraction);
            assert(0);
#endif
        }
    }
    else if (dir == -1) {
        fraction = solValue - floor(solValue);
        if (fraction >= 1.0e-5) {
            cost = objDiff / (fraction + 1.0e-9);
            downCost_ = (downCost_ * downCount_ + cost) / (1 + downCount_);
            ++downCount_;
        }
        else {
#ifdef BLIS_DEBUG
            ALPS_PRINTF("WARNING: small fraction %.10g, shouldn't happen.\n",
                        fraction);
            assert(0);
#endif
        }
    }
    else {
        ALPS_PRINTF("ERROR: wrong direction %d.\n", dir);
        assert(0);
    }
    
    score_ = weight_* ALPS_MIN(upCost_, downCost_) + 
        (1.0 - weight_) * ALPS_MAX(upCost_, downCost_);
    
}

//#############################################################################
