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
 * Copyright (C) 2001-2013, Lehigh University, Yan Xu, and Ted Ralphs.       *
 * All Rights Reserved.                                                      *
 *===========================================================================*/

#include "CoinHelperFunctions.hpp"
#include "CoinTime.hpp"
#include "CoinWarmStartBasis.hpp"

#include "Alps.h"
#include "AlpsKnowledgeBroker.h"

#include "BlisBranchStrategyStrong.h"
#include "BlisSolution.h"
#include "BlisObjectInt.h"

//#############################################################################

// Copy constructor 
BlisBranchStrategyStrong::BlisBranchStrategyStrong (
    const BlisBranchStrategyStrong & rhs
    )
    :
    BcpsBranchStrategy()
{
    bestChangeUp_ = rhs.bestChangeUp_;
    bestNumberUp_ = rhs.bestNumberUp_;
    bestChangeDown_ = rhs.bestChangeDown_;
    bestNumberDown_ = rhs.bestNumberDown_;
}

//#############################################################################

/** Create a set of candidate branching objects. */
int 
BlisBranchStrategyStrong::createCandBranchObjects(int numPassesLeft)
{
    int bStatus = 0;
    int i, j, pass;

    int col, ind;
    int numInfs = 0;
    int numIntegerInfs = 0;  // For integer objects.
    int numObjectInfs = 0;   // For non-integer objects.
    double sumDeg = 0.0;     // For solution estimate.
    double lpX;
    
    BlisObjectInt *intObject = NULL;
    
    BlisModel *model = dynamic_cast<BlisModel *>(model_);
    OsiSolverInterface * solver = model->solver();
    
    int numCols = model->getNumCols();
    int numObjects = model->numObjects();
    bool beforeSolution = (model->getNumSolutions() == 0);
    
    int givenStrongLen = dynamic_cast<BlisParams*>
        (model->BlisPar())->entry(BlisParams::strongCandSize);
    int strongLen = givenStrongLen;
    int maxStrongLen = CoinMax(CoinMin(givenStrongLen, numObjects), 1);
    
    BlisStrong * candStrongs = new BlisStrong [maxStrongLen];

    //------------------------------------------------------
    // Backup solver status
    //------------------------------------------------------

    // Objective value. 
    double objValue = solver->getObjSense() * solver->getObjValue();

    // Column bounds.
    const double * lower = solver->getColLower();
    const double * upper = solver->getColUpper();
    double * saveUpper = new double[numCols];
    double * saveLower = new double[numCols];
    for (i = 0; i < numCols; ++i) {
	saveLower[i] = lower[i];
	saveUpper[i] = upper[i];
    }
    
    // Primal solution.
    double * saveSolution = new double[numCols];
    memcpy(saveSolution, solver->getColSolution(),numCols * sizeof(double));

    //------------------------------------------------------
    // Select a set of objects based on feasibility. 
    // NOTE: we might go round this loop twice if we are feed in
    //       a "feasible" solution.
    //------------------------------------------------------
    
    for (pass = 0; pass < 2; ++pass) {
	
	// Compute how many infeasible objects. 
        // NOTE: it set model->savedLpSolution
	model->feasibleSolution(numIntegerInfs, numObjectInfs);

	sumDeg = 0.0; 
	numInfs = 0;

        int preferDir;        
	int leastFrac = 0;
        double infeasibility;
        double minInf = ALPS_ZERO;
        BlisObjectInt * object = NULL;
        
	for (i = 0; i < maxStrongLen; ++i) {
	    candStrongs[i].bObject = NULL;
	}
	
	strongLen = 0;
        
	for (i = 0; i < numObjects; ++i) {
            
            // TODO: currently all integer object.
	    object = dynamic_cast<BlisObjectInt *>(model->objects(i));
	    infeasibility = object->infeasibility(model, preferDir);
            
	    if (infeasibility) {
		++numInfs;
                
		// Increase estimated degradation to solution
		sumDeg += object->pseudocost().getScore();
		
		// Check for suitability based on infeasibility.
                if (infeasibility > minInf) {
		    
                    if (candStrongs[leastFrac].bObject) {
                        // The slot already has one, free it.
                        delete candStrongs[leastFrac].bObject;
                    }
                    
                    // Create new branching object.
		    candStrongs[leastFrac].bObject = 
			object->createBranchObject(model, preferDir);

		    candStrongs[leastFrac].bObject->setUpScore(infeasibility);
                    candStrongs[leastFrac].bObject->setObjectIndex(i);
                    candStrongs[leastFrac].objectIndex = i;
                    
		    strongLen = CoinMax(strongLen, leastFrac + 1);
                    
                    // Find an empty or the worst slot.
		    leastFrac = -1;
		    minInf = ALPS_INFINITY;
                    
		    for(j = 0; j < maxStrongLen; ++j) {
                        if (!candStrongs[j].bObject) {
                            // j is an empty slots.
                            minInf = ALPS_ZERO;
                            leastFrac = j;
                            break;
                        }
			else if(candStrongs[j].bObject->getUpScore() < minInf){
			    minInf = candStrongs[j].bObject->getUpScore();
			    leastFrac = j;
			}
		    }
		}
	    }
	}

	if (numInfs) {
#ifdef BLIS_DEBUG_MORE
            std::cout << "STRONG: numInfs = " << numInfs
                      << std::endl;
#endif
            model->setSolEstimate(objValue + sumDeg);
            
	    break;
	}
	else if (pass == 0) {
	    // The first pass and ip feasible
            
#ifdef BLIS_DEBUG_MORE
	    std::cout << "STRONG: given a feasible sol" << std::endl;
#endif
            
	    bool roundAgain = false;
            
	    CoinWarmStartBasis * ws = 
		dynamic_cast<CoinWarmStartBasis*>(solver->getWarmStart());
	    if (!ws) break;

	    // Force solution values within bounds
	    for (i = 0; i < numCols; ++i) {
		double lpX = saveSolution[i];
		if (lpX < lower[i]) {
		    saveSolution[i] = lower[i];
		    roundAgain = true;
		    ws->setStructStatus(i, CoinWarmStartBasis::atLowerBound);
		} 
		else if (lpX > upper[i]) {
		    saveSolution[i] = upper[i];
		    roundAgain = true;
		    ws->setStructStatus(i, CoinWarmStartBasis::atUpperBound);
		} 
	    }
	    
	    if (roundAgain) {
		// Need resolve, then do the second round selection.
		solver->setWarmStart(ws);
		delete ws;

		// Resolve.
		solver->resolve();
		
		// Save new lp solution.
		memcpy(saveSolution, 
		       solver->getColSolution(),
		       numCols * sizeof(double));

		if (!solver->isProvenOptimal()) {
		    // Become infeasible, can do nothing. 
		    bStatus = -2;
		    break;
		}
	    } 
	    else {
		delete ws;
		break;
	    }
	}
    }

    //------------------------------------------------------
    // Now we have a set of candidate branching object, 
    // evaluate them.
    //------------------------------------------------------
    
    int bestBO = 0;
    double bestScore = 0.0;
    
    for (i = 0; i < strongLen; ++i) {

	candStrongs[i].numIntInfUp = numInfs;
	candStrongs[i].numIntInfDown = numInfs;

	if ( !model->objects(candStrongs[i].bObject->getObjectIndex())->
             boundBranch(model) ) {
            // Weild branching: not by modifying variable bounds.
            // Exit selection routine.
	    strongLen = 0;
        }
        
	// Find the most fractional object in case of doing simple branch
	if (candStrongs[i].bObject->getUpScore() > bestScore) {
	    bestScore = candStrongs[i].bObject->getUpScore();
	    bestBO = i;
	}
    }
    
    // If we have hit max time don't do strong branching
    double timeLimit = model->AlpsPar()->entry(AlpsParams::timeLimit);
    bool maxTimeReached = (CoinCpuTime() - model->startTime_  > timeLimit);
    
#ifdef BLIS_DEBUG_MORE
    printf("1. strongLen = %d, maxTimeReached %d, numPassesLeft %d\n", 
	   strongLen, maxTimeReached, numPassesLeft);
#endif
    
    if (strongLen <= 0 || maxTimeReached) {
        
        //--------------------------------------------------
        // Simple max infeasibility branching.
        //--------------------------------------------------
        
#ifdef BLIS_DEBUG
        std::cout << "NOT STRONG: maxTimeReached=" << maxTimeReached
                  << "; numPassesLeft=" << numPassesLeft 
                  << std::endl;
#endif

        // Create candidate set.
        numBranchObjects_ = 1;
        branchObjects_ = new BcpsBranchObject * [1];
        branchObjects_[0] = candStrongs[bestBO].bObject;
        candStrongs[bestBO].bObject = NULL;
        strongLen = 0;
    }
    else {
        
        //--------------------------------------------------------
        /* Strong branching */
        //--------------------------------------------------------
        
        // set true to say look at all even if some fixed (experiment)
	bool solveAll = false;
	int saveLimit;
        
	CoinWarmStart * ws = solver->getWarmStart();
	solver->getIntParam(OsiMaxNumIterationHotStart, saveLimit);
	if (beforeSolution) {
	    solver->setIntParam(OsiMaxNumIterationHotStart, 10000);
        }
        
        // Mark hot start
        solver->markHotStart();
        
#ifdef BLIS_DEBUG_MORE
	printf("BEFORE LOOP: strongLen = %d\n",strongLen);
#endif
        
	for (i = 0; i < strongLen; ++i) {
	    double objChange;
	    double newObjValue = ALPS_DBL_MAX;
            
	    // status is 0 finished, 1 infeasible and other
	    int lpStatus;
            
            //----------------------------------------------
	    // Branching down.
            //----------------------------------------------
            
            candStrongs[i].bObject->setDirection(-1);
            candStrongs[i].bObject->branch();
            solver->solveFromHotStart();

            if (solver->isProvenOptimal()) {
                lpStatus = 0; // optimal
            }
            else if (solver->isIterationLimitReached()
                     &&!solver->isDualObjectiveLimitReached()) {
                lpStatus = 2; // unknown 
            }
            else {
                lpStatus = 1; // infeasible
            }
            
            newObjValue = solver->getObjSense() * solver->getObjValue();
            objChange = newObjValue-objValue ;
            
#ifdef BLIS_DEBUG_MORE
	    std::cout << "Down: lpStatus = " << lpStatus << std::endl;
#endif

	    if (lpStatus == 0) {

                // Update pseudocost
                ind = candStrongs[i].objectIndex;
                intObject = dynamic_cast<BlisObjectInt *>(model->objects(ind));
                col = intObject->columnIndex();
                lpX = saveSolution[col];
                intObject->pseudocost().update(-1, objChange, lpX);

		candStrongs[i].finishedDown = true ;
		if (newObjValue >= model->getCutoff()) {
		    objChange = ALPS_DBL_MAX; // say infeasible
		} 
                else {    
		    if(model->feasibleSolution(candStrongs[i].numIntInfDown)){
			printf("STRONG:easy:down:found a feasible solution\n");
		    }
		    
		    // See if integer solution
		    if (model->feasibleSolution(candStrongs[i].numIntInfDown,
						candStrongs[i].numObjInfDown)) {
#ifdef BLIS_DEBUG_MORE
			printf("STRONG:down:found a feasible solution\n");
#endif
			
			model->setBestSolution(BLIS_SOL_STRONG,
					       newObjValue,
					       solver->getColSolution());
			(model->getKnowledgeBroker())->
                            getNodeSelection()->setWeight(0.0);
			BlisSolution* ksol = 
			    new BlisSolution(solver->getNumCols(), 
					       solver->getColSolution(), 
					       newObjValue);
			model->getKnowledgeBroker()->
                            addKnowledge(AlpsKnowledgeTypeSolution, 
                                         ksol, 
                                         newObjValue);

			objChange = ALPS_DBL_MAX ;
		    }
		}
	    }
	    else if (lpStatus == 1) {
	      objChange = ALPS_DBL_MAX ;
	    } 
	    else {
		// Can't say much as we did not finish
		candStrongs[i].finishedDown = false ;
	    }
            
	    candStrongs[i].bObject->setDownScore(objChange);
	    
	    // restore bounds
            int numDiff = 0;
            for (j = 0; j < numCols; ++j) {
                if (saveLower[j] != lower[j]) {
                    solver->setColLower(j, saveLower[j]);
                    ++numDiff;
                }
                if (saveUpper[j] != upper[j]) {
                    solver->setColUpper(j, saveUpper[j]);
                    ++numDiff;
                }
            }

#ifdef BLIS_DEBUG_MORE
            std::cout << "numDiff = " << numDiff << std::endl;
#endif	    
            
            //----------------------------------------------
	    // Branching up.
            //----------------------------------------------

            candStrongs[i].bObject->branch();
            solver->solveFromHotStart();

            if (solver->isProvenOptimal()) {
                lpStatus = 0; // optimal
            }
            else if (solver->isIterationLimitReached()
                     &&!solver->isDualObjectiveLimitReached()) {
                lpStatus = 2; // unknown 
            }
            else {
                lpStatus = 1; // infeasible
            }
            
            newObjValue = solver->getObjSense() * solver->getObjValue();
            objChange = newObjValue - objValue ;
            
#ifdef BLIS_DEBUG_MORE
	    std::cout << "STRONG: Up: lpStatus = " << lpStatus << std::endl;
#endif      

	    if (lpStatus == 0) {
                // Update pseudocost
                ind = candStrongs[i].objectIndex;
                intObject = dynamic_cast<BlisObjectInt *>(model->objects(ind));
                col = intObject->columnIndex();
                lpX = saveSolution[col];
                intObject->pseudocost().update(1, objChange, lpX);

		candStrongs[i].finishedUp = true ;
		if (newObjValue >= model->getCutoff()) {
		    objChange = ALPS_DBL_MAX; // Cutoff
		} 
		else {
		    // See if integer solution
		    if (model->feasibleSolution(candStrongs[i].numIntInfUp,
						candStrongs[i].numObjInfUp)) {
                        
#ifdef BLIS_DEBUG_MORE
			printf("STRONG:up:found a feasible solution\n");
#endif
                        
			model->setBestSolution(BLIS_SOL_STRONG,
					       newObjValue,
					       solver->getColSolution());
                        
			model->getKnowledgeBroker()->
                            getNodeSelection()->setWeight(0.0);
                        
			BlisSolution* ksol =
			    new BlisSolution(solver->getNumCols(), 
					     solver->getColSolution(), 
					     newObjValue);
                        
			model->getKnowledgeBroker()->
                            addKnowledge(AlpsKnowledgeTypeSolution, 
                                         ksol, 
                                         newObjValue);
                        
                        objChange = ALPS_DBL_MAX;
		    }
		}
	    } 
	    else if (lpStatus == 1) {
		objChange = ALPS_DBL_MAX;
	    } 
            else {
		// Can't say much as we did not finish
		candStrongs[i].finishedUp = false ;
	    }
	    candStrongs[i].bObject->setUpScore(objChange);
	    
	    // restore bounds
            for (j = 0; j < numCols; ++j) {
                if (saveLower[j] != lower[j]) {
                    solver->setColLower(j,saveLower[j]);
                }
                if (saveUpper[j] != upper[j]) {
                    solver->setColUpper(j,saveUpper[j]);
                }
	    }
            
            //----------------------------------------------
            // End of evaluation for this branching object.
            // Possibilities are:
	    // 1) Both sides below cutoff; this variable is a 
            //    candidate for branching.
	    // 2) Both sides infeasible or above the obj cutoff: 
	    //    no further action here. Break from the evaluation loop and 
	    //    assume the node will be purged by the caller.
	    // 3) One side below cutoff: Install the branch (i.e., fix the 
	    //    variable). Break from the evaluation loop and assume the 
	    //    node will be reoptimised by the caller.
            //----------------------------------------------
            
	    if (candStrongs[i].bObject->getUpScore() < ALPS_INFINITY) {
		if(candStrongs[i].bObject->getDownScore() < ALPS_INFINITY) {
		    // feasible - no action
		} 
                else {
		    // up feasible, down infeasible
		    bStatus = -1;
		    if (!solveAll) {
			candStrongs[i].bObject->setDirection(1);
			candStrongs[i].bObject->branch();
			break;
		    }
		}
	    } 
            else {
		if(candStrongs[i].bObject->getDownScore() < ALPS_INFINITY) {
		    // down feasible, up infeasible
		    bStatus = -1;
		    if (!solveAll) {
			candStrongs[i].bObject->setDirection(-1);
			candStrongs[i].bObject->branch();
			break;
		    }
		} 
                else {
		    // neither side feasible
		    bStatus = -2;
		    break;
		}
	    }
	}// EOF the loop for checking each candiate
        
        
        //--------------------------------------------------
        // Unmark hotstart and reset lp solver.
        //--------------------------------------------------
        
        solver->unmarkHotStart();
        solver->setIntParam(OsiMaxNumIterationHotStart, saveLimit);
        solver->setWarmStart(ws);
        delete ws;
    }
     

    if (bStatus >= 0) {
        // Store the set of candidate branching objects. 
        numBranchObjects_ = strongLen;
        
        branchObjects_ = new BcpsBranchObject* [strongLen];        
        for (i = 0; i < strongLen; ++i) {
            branchObjects_[i] = candStrongs[i].bObject;
            candStrongs[i].bObject = NULL;
        }
    }
    
    // Restore solution
    solver->setColSolution(saveSolution);
    
    // Cleanup.
    delete [] saveSolution;
    delete [] saveLower;
    delete [] saveUpper;

    for (i = 0; i < maxStrongLen; ++i) {
        if (candStrongs[i].bObject) {
            delete candStrongs[i].bObject;
        }
    }
    delete [] candStrongs;
    

    return bStatus;
}

//#############################################################################

/** Strong method to compare object thisOne to the bestObject_. 
    Compare based on number of infeasibility objects (numInfUp, numInfDn)
    until a solution is found by search, then swith to change 
    (changeUp, changeDn) in objective. 
    The lesser of the MIN(numInfUp, numInfDown), the better.
    The larger of the MIN(changeUp, changeDown), the better. */
int
BlisBranchStrategyStrong::betterBranchObject(BcpsBranchObject * thisOne,
					     BcpsBranchObject * bestSoFar)
{
    int betterDirection = 0;
    double bestChange;
    
    //BlisModel *model = dynamic_cast<BlisModel *>(model_);    

    if (!bestSoFar) {
	bestChange = -1.0;
    }
    else {
	bestChange = ALPS_MIN(bestChangeUp_, bestChangeDown_);
    }
    
    double upCost = thisOne->getUpScore();
    double downCost = thisOne->getDownScore();
    
    if (upCost <= downCost) {
	if (upCost > bestChange) {
	    betterDirection = 1;    // Branching up.
	}
    } 
    else {
	if (downCost > bestChange) {
	    betterDirection = -1;   // Branching down
	}
    }
    
    if (betterDirection) {
	bestChangeUp_ = upCost;
	bestChangeDown_ = downCost;
    }

    return betterDirection;
}

//#############################################################################
