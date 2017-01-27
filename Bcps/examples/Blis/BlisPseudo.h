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

#ifndef BlisPseudo_h_
#define BlisPseudo_h_

#include "CoinError.hpp"

//#############################################################################

class BlisPseudocost 
{
 private:
    /** Use to calculate score. */
    double weight_;

    /** How many times being branched up. */
    int upCount_;

    /** Average object change when branching up. */
    double upCost_;

    /** How many times being branched down. */
    int downCount_;

    /** Average object change when branching down. */
    double downCost_;

    /** The estimated importance. 
        Score = weight * MIN(downCost_, upCost_) +
        (1.0 - weight) * MAX(downCost_, upCost_)
    */
    double score_;
    
 public:
    /** Default constructor. */
    BlisPseudocost() : 
        weight_(1.0),
        upCount_(0),
        upCost_(0.0), 
        downCount_(0), 
        downCost_(0.0),
        score_(0.0)
        {}
    
    /** Useful constructor. */
    BlisPseudocost(double uc, 
                   int un,
                   double dc, 
                   int dn,
                   double s)
        :
        weight_(1.0),
        upCount_(un),
        upCost_(uc), 
        downCount_(dn),
        downCost_(dc),
        score_(s) 
        {}
    
    /** Set weigth. */
    void setWeight(double w) { 
        if (w < 0.0 || w > 1.0) {   
            throw CoinError("weight is not in range [0,1]", "setWeight", 
                            "BlisPseudo");
        }
        weight_= w; 
    }

    /** Update pseudocost. */
    void update(const int dir,
                const double parentObjValue,
                const double objValue,
                const double solValue);

    /** Update pseudocost. */
    void update(const int dir,
                const double objDiff,
                const double solValue);
    
    /** Get up branching count. */
    int getUpCount() { return upCount_; }

    /** Get up branching cost. */
    double getUpCost() { return upCost_; }
    
    /** Get down branching count. */
    int getDownCount() { return downCount_; }
    
    /** Get down branching cost. */
    double getDownCost() { return downCost_; } 

   /** Get importance. */
    double getScore() { return score_; } 
};

#endif
