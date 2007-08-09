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
// This file is modified from COIN/Cbc/CbcHeuristic.hpp
//#############################################################################


#ifndef BlisHeuristic_h_
#define BlisHeuristic_h_

#include <string>
#include <vector>

#include "CoinPackedMatrix.hpp"
#include "OsiCuts.hpp"

class BlisModel;


//#############################################################################


/** Heuristic base class */
class BlisHeuristic {

 private:

    /** Illegal Assignment operator */
    BlisHeuristic & operator = (const BlisHeuristic& rhs);
    
 protected:

    /** Pointer to the model */
    BlisModel *model_;
    
    /** Heuristics name */
    char *name_;
    
    /** When to call findSolution() routine. 
        -2:                    disable
        -1:                    just root
         0:                    automatically decided by BLIS
         any positive integer: the node interval between the call
    */
    int strategy_;

    /** Number of solutions found. */
    int numSolutions_;

    /** Used CPU/User time. */
    double time_;

    /** The times of calling this heuristic. */
    int calls_;
    
 public:

    /** Default Constructor. */
    BlisHeuristic() {
        model_ = NULL;
        name_ = strdup("Unknown");
        strategy_ = 0;
        numSolutions_ = 0;
        time_ = 0.0;
        calls_ = 0;
    }
    
    /** Useful constructor. */
    BlisHeuristic(BlisModel *model, const char *name, int strategy) {
        model_ = model;
        if (name) {
            name_ = strdup(name);
        }
        else {
            name_ = strdup("Unknown");
        }
        strategy_ = strategy;
        numSolutions_ = 0;
        time_ = 0.0;
        calls_ = 0;
    }

    /** Distructor. */
    virtual ~BlisHeuristic() { if (name_) free(name_); }

    /** Copy constructor. */
    BlisHeuristic(const BlisHeuristic & rhs) {        
        model_ = rhs.model_;
        name_ = strdup(rhs.name_);
        strategy_ = rhs.strategy_; // What if disabled?
        numSolutions_ = 0;
        time_ = 0.0;
        calls_ = 0;
    }
 
    /** update model (This is needed if cliques update matrix etc). */
    virtual void setModel(BlisModel * model) { model_ = model ;}
    
    /** Get/set strategy. */
    //@{
    virtual void setStrategy(int strategy) { strategy_ = strategy; }
    virtual int strategy() { return strategy_; }
    //@]

    /** Clone a heuristic. */
    virtual BlisHeuristic * clone() const = 0;
    
    /** returns 0 if no solution, 1 if valid solution
        with better objective value than one passed in
        Sets solution values if good, sets objective value 
        This is called after cuts have been added - so can not add cuts
    */
    virtual bool searchSolution(double & objectiveValue, 
                                double * newSolution)=0;
    
    /** returns 0 if no solution, 1 if valid solution, -1 if just
        returning an estimate of best possible solution
        with better objective value than one passed in
        Sets solution values if good, sets objective value 
        (only if nonzero code) 
        This is called at same time as cut generators - so can add cuts
        Default is do nothing
    */
    virtual int searchSolution(double & objectiveValue,
                               double * newSolution,
                               OsiCuts & cs) { return 0; }
    
    /** Number of solutions found. */
    int numSolutions() { return numSolutions_; }
    
    /** Cpu time used. */
    inline double time() { return time_; }
    
    /** Cpu time used. */
    inline int calls() { return calls_; }
};

#endif

