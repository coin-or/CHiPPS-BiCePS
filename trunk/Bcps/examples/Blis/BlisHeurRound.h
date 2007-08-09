/*===========================================================================*
 * This file is part of the Bcps Linear Solver (BLIS).                       *
 *                                                                           *
 * BLIS is distributed under the Common Public License as part of the        *
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
 * Copyright (C) 2001-2005, International Business Machines                  *
 * Corporation, Lehigh University, Yan Xu, Ted Ralphs, Matthew Salzman and   *
 * others. All Rights Reserved.                                              *
 *===========================================================================*/

//#############################################################################
// This file is modified from CbcHeuristic.hpp
//#############################################################################

#ifndef BlisHeurRound_h_
#define BlisHeurRound_h_

#include <string>
#include <vector>

#include "CoinPackedMatrix.hpp"
#include "OsiCuts.hpp"

#include "BlisHeuristic.h"

class BlisModel;

//#############################################################################

/** Rounding Heuristic.  */
class BlisHeurRound : public BlisHeuristic {
 private:
    /** Illegal Assignment operator. */ 
    BlisHeurRound & operator=(const BlisHeurRound& rhs);

 protected:
    /** Column majored matrix. */
    CoinPackedMatrix matrix_;
    
    /** Row majored matrix. */
    CoinPackedMatrix matrixByRow_;
    
    /** Seed for random stuff. */
    int seed_;
    
 public:
    /** Default Constructor. */
    BlisHeurRound() {}
    
    /** Constructor with model - assumed before cuts. */
    BlisHeurRound(BlisModel * model, const char *name, int strategy)
        :
        BlisHeuristic(model, name, strategy)
        {
            assert(model->solver());
        }

    /** Destructor. */
    ~BlisHeurRound() {}
  
    /** Copy constructor. */
    BlisHeurRound( const BlisHeurRound &);
    
    /** Clone a rounding heuristic */
    virtual BlisHeuristic * clone() const;
    
    /** update model (This is needed if cliques update matrix etc). */
    virtual void setModel(BlisModel * model);
    
    /** returns 0 if no solution, 1 if valid solution
        with better objective value than one passed in
        Sets solution values if good, sets objective value (only if good)
        This is called after cuts have been added - so can not add cuts
    */
    virtual bool searchSolution(double & objectiveValue,
                               double * newSolution);
    
    /** Set seed */
    void setSeed(int value) { seed_ = value; }

};
#endif
