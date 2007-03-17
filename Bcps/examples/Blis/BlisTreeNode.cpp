/*===========================================================================*
 * This file is part of the Bcps Linear Solver (BLIS).                       *
 *                                                                           *
 * BLIS is distributed under the Common Public License as part of the        *
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

#include <cassert>
#include <iostream>
#include <utility>
#include <cmath>
#include <vector>

#include "CoinUtility.hpp"
#include "OsiRowCut.hpp"
#include "OsiColCut.hpp"
#include "OsiRowCutDebugger.hpp"
#include "OsiCuts.hpp"

#include "AlpsKnowledge.h"
#include "AlpsEnumProcessT.h"
#include "AlpsKnowledgeBroker.h"
#include "AlpsTreeNode.h"

#include "BcpsBranchStrategy.h"

#include "BlisBranchObjectInt.h"
#include "BlisConstraint.h"
#include "BlisHelp.h"
#include "BlisTreeNode.h"
#include "BlisModel.h"
#include "BlisNodeDesc.h"
#include "BlisModel.h"
#include "BlisObjectInt.h"
#include "BlisParams.h"
#include "BlisSolution.h"
//#include "BlisVariable.h"

#define REMOVE_SLACK 1

//#############################################################################

AlpsTreeNode*
BlisTreeNode::createNewTreeNode(AlpsNodeDesc *&desc) const
{
    // Create a new tree node
    BlisTreeNode *node = new BlisTreeNode(desc);

    // Set solution estimate for this nodes.
    // solEstimate = quality_ + sum_i{min{up_i, down_i}}
    
    node->setSolEstimate(solEstimate_);

#ifdef BLIS_DEBUG_MORE
    printf("BLIS:createNewTreeNode: quality=%g, solEstimate=%g\n",
           quality_, solEstimate_);
#endif
    
    desc = NULL;

    return node;
}

//#############################################################################

// NOTE: if rampup,
// - parent must be explicit if not NULL,
// - this node is explicit.

int
BlisTreeNode::process(bool isRoot, bool rampUp)
{
    int status = BLIS_OK;
    int j, k = -1;
    int numCols, numRows, numCoreCols, numCoreRows;
    int numStartRows, origNumStartRows;
    int maxNumCons, tempNum;

    int numIntInfs = 0;
    int numObjInfs = 0;

    int origNumOldCons = 0;
    int currNumOldCons = 0;
    int numNewCons = 0;
    int maxNewCons = 0;

    int numAppliedCons = 0;
    int pass = 0;
    int maxPass = 20;

    // Only autmatic stategy has depth limit.
    int maxConstraintDepth = 20;

    int numPassesLeft = 20;      
    int bStatus = -1;
 
    double cutoff = getKnowledgeBroker()->getIncumbentValue();
    double parentObjValue = getQuality();
    double primalTolerance = 1.0e-6;
    double tailOffTol = 1.0e-7;
    double preObjValue = -ALPS_OBJ_MAX;
    double improvement = 100.0;

    double  heurObjValue;
    double *heurSolution = NULL;
    double *currLpSolution = NULL;

    bool keepOn = true;
    bool feasibleIP = false;
    bool needBranch = false;
    bool betterSolution = false;
    bool lpFeasible = false;
    bool foundSolution = false;
    bool genConsHere = false;

    CoinWarmStartBasis::Status rowStatus;
    BlisConstraint *aCon = NULL;
    BcpsObject **newConstraints = NULL;

    int numDelRows = 0;
    int *delRow = NULL;
    int *oldConsPos = NULL;
    
    std::vector<int> delIndices;
    
    BlisModel* model = dynamic_cast<BlisModel*>(desc_->getModel());

    AlpsPhase phase = knowledgeBroker_->getPhase();
    
    //------------------------------------------------------
    // Check if this can be fathomed by objective cutoff.
    //------------------------------------------------------
    
    //std::cout << "parentObjValue = " << parentObjValue << std::endl;
    if (parentObjValue - primalTolerance > cutoff) {
	setStatus(AlpsNodeStatusFathomed);
        //std::cout << "fathom!" <<std::endl;
	goto TERM_PROCESS;
    }
    else {
        model->setActiveNode(this);
        model->addNumNodes();
    }
    
    //------------------------------------------------------
    // Get model information and parameters.
    //------------------------------------------------------

    numCols = model->solver()->getNumCols();
    numRows = model->solver()->getNumRows();
    numCoreCols = model->getNumCoreVariables();
    numCoreRows = model->getNumCoreConstraints();
    numAppliedCons =  numRows - numCoreRows;
    
    maxNumCons = model->getMaxNumCons();

    heurSolution = new double [numCols];
    currLpSolution = new double [numCols];

    //------------------------------------------------------
    // Decides if can generate constraints.
    //------------------------------------------------------

    // Mark if this node is root or not.
    model->isRoot_ = isRoot;

    genConsHere = false;
    if (model->useCons() == -2) {
	genConsHere = false;
    }
    else if (model->useCons() == -1) {
	// The original root only
	if (isRoot && index_ == 0) genConsHere = true;
    }
    else if (model->useCons() == 0) {
	// NOTE: Only automatic has depth limit.
	if (depth_ < maxConstraintDepth) {
            if (!diving_ || isRoot) genConsHere = true;
	}
    }
    else {
	genConsHere = true;
    }
    
    //======================================================
    // Restore, load and solve the subproblem.
    // (1) LP infeasible
    //     a. set status to be fathom.
    // (2) LP feasible
    //     a. MILP feasible. Check whether need update incumbent.
    //     b. LP feasible but not MIP feasible. Check whether can be 
    //        fathomed, if not, choose a branch variable.
    //======================================================

    //------------------------------------------------------
    // Extract info from this node and load subproblem into lp solver.
    //------------------------------------------------------
    
    installSubProblem(model);

    //------------------------------------------------------
    // Bounding, heuristic searching, constraint generating.
    //------------------------------------------------------

    numStartRows = model->solver()->getNumRows();
    origNumStartRows = numStartRows;

    origNumOldCons = numStartRows - numCoreRows;
    currNumOldCons = origNumOldCons;
    
#ifdef BLIS_DEBUG
    std::cout << "PROCESS: genConsHere=" << genConsHere
	      << ", useCons=" << model->useCons()
	      << ", numCoreRows=" << numCoreRows
	      << ", numStartRows=" << numStartRows
	      << ", currNumOldCons=" << currNumOldCons << std::endl;
#endif

    if (currNumOldCons > 0) {
	oldConsPos = new int [currNumOldCons];
	for (k = 0; k < currNumOldCons; ++k) {
	    oldConsPos[k] = k;
	}
    }

    if (genConsHere) {
        maxNewCons = maxNumCons;
        newConstraints = new BcpsObject* [maxNewCons];    
    }

    while (keepOn && (pass < maxPass)) {
        ++pass;
        keepOn = false;
        
        //--------------------------------------------------
        // Bounding to get the quality of this node.
        //--------------------------------------------------
        
        status = bound(model);
	if (pass == 1) {
	    int iter = model->solver()->getIterationCount();
	    model->addNumIterations(iter);
	}
        
        switch(status) {
        case BLIS_LP_OPTIMAL:
            // Check if IP feasible 
            feasibleIP = model->feasibleSolution(numIntInfs, numObjInfs);
            
            if (feasibleIP) {         
                // IP feasible 
		
		if (quality_ < cutoff) {  
                    // Better than incumbent
                    betterSolution = true;
                    model->setBestSolution(BLIS_SOL_BOUNDING,
                                           quality_, 
                                           model->getLpSolution());
                    getKnowledgeBroker()->getNodeSelection()->setWeight(0.0);
                    BlisSolution* ksol = 
                        new BlisSolution(numCols, 
                                           model->getLpSolution(), 
                                           quality_);
                    getKnowledgeBroker()->addKnowledge(ALPS_SOLUTION, 
                                                       ksol, 
                                                       quality_); 
                    // Update cutoff
                    cutoff = getKnowledgeBroker()->getIncumbentValue();
                }
                setStatus(AlpsNodeStatusFathomed);
		goto TERM_PROCESS;
            }
            else {
                cutoff = getKnowledgeBroker()->getIncumbentValue();
                if (quality_ > cutoff) {
                    setStatus(AlpsNodeStatusFathomed);
                    goto TERM_PROCESS;
                }
                needBranch = true;
                reducedCostFix(model);
                
                //------------------------------------------
                // Check if tailoff
                //------------------------------------------

                if (pass > 1) {
                    improvement = quality_ - preObjValue;
                    if (improvement > tailOffTol) {
                        // NOTE: still need remove slacks, although
                        //       tailoff.
                        keepOn = true;
                    }
                    
#ifdef BLIS_DEBUG_MORE
                    std::cout << "PROCESS: pass["<< pass << "], improvement=" 
                              << improvement << ", tailOffTol=" << tailOffTol
                              << std::endl;
#endif
                }
                else {
                    keepOn = true;
                }
                // Update previous objective value.
                preObjValue = quality_;

                //------------------------------------------
                // Remove non-core slack constraints. 
                //------------------------------------------

                numRows = model->getNumRows();
                
                if ( genConsHere &&
                     //(improvement > tailOffTol) && 
                     //(numRows > numStartRows) ) {
                     (numRows > numCoreRows) ) {   
                 
#ifdef BLIS_DEBUG                  
                    if ( (numStartRows + numNewCons != numRows) ||
                         (numCoreRows + currNumOldCons + numNewCons != numRows) ) {
                        
                        std::cout << "ERROR: numRows=" << numRows
                                  << "; numCoreRows=" << numCoreRows
                                  << "; numStartRows=" << numStartRows
                                  << "; numNewCons=" << numNewCons
                                  << "; currNumOldCons=" << currNumOldCons
                                  << std::endl;
                        
                        assert(numRows - numStartRows == numNewCons);
                    }
#endif
                    
                    int *oldDelMark = NULL;
                    if (currNumOldCons > 0) {
                        oldDelMark = new int [currNumOldCons];
                        CoinZeroN(oldDelMark, currNumOldCons);
                    }
                    int *newDelMark = NULL;
                    if (numNewCons > 0) {
                        newDelMark = new int [numNewCons];
                        CoinZeroN(newDelMark, numNewCons);
                    }
                    
                    const CoinWarmStartBasis* ws= 
                        dynamic_cast<CoinWarmStartBasis*>
                        (model->solver()->getWarmStart());
                    
                    // Make sure delIndices is empty.
                    assert(delIndices.size()==0);

#if REMOVE_SLACK
                    for (k = numCoreRows; k < numRows; ++k) {
                        rowStatus = ws->getArtifStatus(k);

                        if (rowStatus == CoinWarmStartBasis::basic) {
                            delIndices.push_back(k);

                            if (k < numStartRows) {
                                oldDelMark[(k-numCoreRows)] = 1;
                            }
                            else {
				newDelMark[(k-numStartRows)] = 1;
                            }
                        }
                    }
#endif
                    numDelRows = delIndices.size();
		    
#ifdef BLIS_DEBUG
                    std::cout << "PROCESS: new cuts=" << numNewCons 
                              << ", slack cuts=" << numDelRows 
                              << ", numRows=" << numRows 
                              << ", numStartRows=" <<numStartRows << std::endl;
#endif
                    
                    if (numDelRows > 0) {
                        delRow = new int [numDelRows];
                        for (k = 0; k < numDelRows; ++k) {
                            delRow[k] = delIndices[k];
#ifdef BLIS_DEBUG
			    std::cout << "REMOVE: slack row " << delRow[k] 
				      << std::endl;
#endif
                        }
                        
                        //----------------------------------
                        // Delete from lp solver.
                        //----------------------------------
                        
                        model->solver()->deleteRows(numDelRows, delRow);
                        
                        delete [] delRow;
                        delRow = NULL;
                        delIndices.clear();

                        //----------------------------------
                        // Delete from the old cut position array.
                        //----------------------------------

                        tempNum = 0;
                        for (k = 0; k < currNumOldCons; ++k) {
                            if (oldDelMark[k] != 1) {
                                // Survived
                                oldConsPos[tempNum++] = oldConsPos[k];    
                            }
                        }
                        currNumOldCons = tempNum;
                        numStartRows = numCoreRows + currNumOldCons;
                        
                        //----------------------------------
                        // Delete from new cut vector.
                        //----------------------------------

                        //std::cout << std::endl;
                        tempNum = 0;
                        for (k = 0; k < numNewCons; ++k) {
                            if (newDelMark[k] == 1) {
                                // Deleted
#ifdef BLIS_DEBUG_MORE
                                std::cout << "delete cut " << k 
                                          << ", size=" 
                                          << dynamic_cast<BlisConstraint*>(newConstraints[k])->getSize()
                                          << std::endl;
#endif
                                
                                delete newConstraints[k];
                                newConstraints[k] = NULL;
                            }
                            else {
                                // Survived
                                newConstraints[tempNum++] = newConstraints[k];
                            }
                        }
                        //assert(tempNum + numDelRows == numNewCons);
                        numAppliedCons -= numNewCons;
                        numAppliedCons += tempNum;
                        numNewCons = tempNum;                        
                        
                        //----------------------------------
                        // Resolve to update solution info in lp solver.
                        //----------------------------------
                        
                        int easy = 2;
                        model->solver()->setHintParam(OsiDoInBranchAndCut,
                                                      true, OsiHintDo, &easy);
                        model->solver()->resolve();
                        model->solver()->setHintParam(OsiDoInBranchAndCut,
                                                      true, OsiHintDo, NULL) ;
#ifdef BLIS_DEBUG
                        if (model->solver()->getIterationCount() != 0) {
                            // TODO: maybe some cuts become slack again
#ifdef BLIS_DEBUG
                            std::cout << "SLACK: resolve changed solution!"
                                      << ", iter=" 
				      << model->solver()->getIterationCount()
				      << std::endl;
#endif
                        }
			else {
#ifdef BLIS_DEBUG
			    std::cout<<"SLACK: resolve don't changed solution!"
                                     << std::endl;
#endif
			}
			
                        double tempOV = model->solver()->getObjValue();
                        double ovDiff = fabs(quality_ - tempOV);
                        
                        if (ovDiff /(1.0 + tempOV) > 1.0e-3) {
                            std::cout << "ERROR: SLACK: quality_("<<quality_ 
                                      << ") != tempOV(" << tempOV
                                      << ")" << std::endl;
                            assert(0);
                        }
                        else {
                            std::cout << "AFTER SLACK: quality_("<<quality_ 
                                      << ") == tempOV(" << tempOV
                                      << ")" << std::endl;
                        }
#endif
                    }
                    
                    delete ws;
                    delete [] newDelMark;
                    delete [] oldDelMark;
                }
            }
            
            break;
        case BLIS_LP_ABANDONED:
#ifdef BLIS_DEBUG
            assert(0);
#endif
            status = BLIS_ERR_LP;
            goto TERM_PROCESS;
        case BLIS_LP_DUAL_INF:
            // FIXME: maybe also primal infeasible
#ifdef BLIS_DEBUG
	    assert(0);
#endif
            status = BLIS_UNBOUND;
            goto TERM_PROCESS;
        case BLIS_LP_PRIMAL_INF:
            setStatus(AlpsNodeStatusFathomed);
            quality_ = -ALPS_OBJ_MAX;       // Remove it as soon as possilbe
            goto TERM_PROCESS;
        case BLIS_LP_DUAL_LIM:
            setStatus(AlpsNodeStatusFathomed);
            quality_ = -ALPS_OBJ_MAX;       // Remove it as soon as possilbe
            goto TERM_PROCESS;
        case BLIS_LP_PRIMAL_LIM:
        case BLIS_LP_ITER_LIM:
            /* Can say much, need branch */
            needBranch = true;
#ifdef BLIS_DEBUG
            assert(0);
#endif
            goto TERM_BRANCH;
            break;
        default:
#ifdef BLIS_DEBUG
            std::cout << "PROCESS: unknown status "  <<  status << std::endl;
            assert(0);
#endif
            break;
        }

        //--------------------------------------------------
        // Apply heuristics.
        //--------------------------------------------------
        
        if (keepOn && model->useHeuristics_) {
            heurObjValue = getKnowledgeBroker()->getIncumbentValue();
            for (k = 0; k < model->numHeuristics(); ++k) {
                foundSolution = false;
                foundSolution = 
                    model->heuristics(k)->searchSolution(heurObjValue,
                                                         heurSolution);
		if (foundSolution) {
                    model->setBestSolution(BLIS_SOL_ROUNDING,
                                           heurObjValue,
                                           heurSolution);
                    BlisSolution* ksol = new BlisSolution(numCols, 
							  heurSolution, 
							  heurObjValue);
                    getKnowledgeBroker()->addKnowledge(ALPS_SOLUTION,
                                                       ksol,
                                                       heurObjValue);
                    // Update cutoff
                    cutoff = getKnowledgeBroker()->getIncumbentValue();
                    if (quality_ > cutoff) {
                        setStatus(AlpsNodeStatusFathomed);
                        goto TERM_PROCESS;
                    }
		}
            }
        }
        
        //--------------------------------------------------
        // Generate constraints.
        //--------------------------------------------------
        
        if ( keepOn && genConsHere && (numAppliedCons < maxNumCons) ) {
            
            OsiCuts newCutPool;
            
            memcpy(currLpSolution, 
                   model->getLpSolution(),
                   numCols * sizeof(double));
            
            status = generateConstraints(model, newCutPool);
            
            if (status != BLIS_LP_OPTIMAL) {
                setStatus(AlpsNodeStatusFathomed);
                quality_ = -ALPS_OBJ_MAX; // Remove it as soon as possilbe
                goto TERM_PROCESS;
            }
            
            // TODO: column cuts, which are already installed by probing.
            tempNum = newCutPool.sizeRowCuts();
            
            if (tempNum > 0) {
                applyConstraints(model, newCutPool, currLpSolution);
                keepOn = true;
                tempNum = newCutPool.sizeRowCuts();

                // Move cuts from OsiCuts to vector newConstraints.
                for (k = 0; k < tempNum; ++k) {
                    aCon = BlisOsiCutToConstraint(&(newCutPool.rowCut(k)));

		    //aCon->hashing(model);

                    newConstraints[numNewCons++] = aCon;
                    if (numNewCons >= maxNewCons) {
                        // No space, need resize
#ifdef BLIS_DEBUG
                        std::cout << "NEWCUT: resize, maxNewCons = " 
                                  << maxNewCons << std::endl;
#endif
                        maxNewCons *= 2;
                        BcpsObject **tempNews = new BcpsObject* [maxNewCons];
                        memcpy(tempNews, 
                               newConstraints,
                               numNewCons * sizeof(BcpsObject *));
                        delete [] newConstraints;
                        newConstraints = tempNews;
                    }
                }
                numAppliedCons += tempNum;
            }
        }
        else {
            keepOn = false;
        }
    }
    
#ifdef BLIS_DEBUG_MORE
    printf("needBranch = %d\n", needBranch);
#endif    
    
    //------------------------------------------------------
    // Select branching object
    //------------------------------------------------------
    
 TERM_BRANCH:
    
    if (needBranch) { 
	
        // Find the best branching object
        numPassesLeft = 20;      
        bStatus = -1;
        
        while (bStatus == -1) { 
            foundSolution = false;
            if(getKnowledgeBroker()->getProcRank() == -1) {
                std::cout << "*** I AM RANK ONE: before choose:bStatus = " 
                          << bStatus << std::endl;
            }
            bStatus = selectBranchObject(model, 
                                         foundSolution, 
                                         numPassesLeft);
            --numPassesLeft;
            
            if (bStatus == -1) { 
                lpFeasible = model->resolve();
                
                //resolved = true ;
#ifdef BLIS_DEBUG_MORE
                printf("Resolve since some col fixed, Obj value %g, numRows %d, cutoff %g\n",
                       model->solver()->getObjValue(),
                       model->solver()->getNumRows(),
                       model->getCutoff());
#endif
                if (lpFeasible) {
		    // Update new quality.
                    setQuality(model->solver()->getObjValue() *
			       model->solver()->getObjSense());
                    if (getQuality() > model->getCutoff()) {
                        bStatus = -2;
                    }
                    feasibleIP = model->feasibleSolution(numIntInfs, 
                                                         numObjInfs);
                    if (feasibleIP) {
                        bStatus = -2;
                        setStatus(AlpsNodeStatusFathomed);
                        //FIXME: is it true? I think this must be true
                        // - 11/16/05, Probably not really. It's possible
                        //   some lp solution are integral, some are not.
                        //   If lp solver has some randomness of producing
                        //   optimal solution, then assertion fails. (p0033)
                        // assert(!feasibleIP);
                    }
                }
                else {
                    // Should not happen. No, it will happen when other
                    // branch is ip feasible, and cause this branch to fathom
                    // when resolving. Test enigma.
                    // assert(0);
                    bStatus = -2;
                    setStatus(AlpsNodeStatusFathomed);
                }
            }
            
            if(getKnowledgeBroker()->getProcRank() == -1) {
                std::cout << "*** I AM RANK ONE: bStatus = " << bStatus
                          << std::endl;
            }
        }
        
        assert(bStatus != -1);
        
        //----------------------------------------------------
        // If found a branching object:
        // 1. Record basis
        // 2. Record soft var bound difference. 
        // 3. Record add/del constraints.
        // NOTE: Hard var bound differences have been recorded when branch().
        //       startXXXXX have branching bounds for this node
        //----------------------------------------------------
        
        if (bStatus >= 0) {
            
#ifdef BLIS_DEBUG_MORE
            //FIXME: Can we do this check? I think so.
            //bool fea = false;
            //int numInfObj;
            //fea = model->feasibleSolution(numUnsatisfied_, numInfObj);
            //assert(!fea);
            
            BlisBranchObjectInt *branchObject =
                dynamic_cast<BlisIntegerBranchObject *>(branchObject_);
            std::cout << "SetPregnant: branchedOn = " 
                      << model->getIntVars()[branchObject->variable()]
                      << std::endl;
#endif
            //--------------------------------------------------
            // Mark as pregnant.
            //--------------------------------------------------

            setStatus(AlpsNodeStatusPregnant);
	    
            //--------------------------------------------------
            // Save basis.
            //--------------------------------------------------

            CoinWarmStartBasis *ws = dynamic_cast<CoinWarmStartBasis*>
                (model->solver()->getWarmStart());
            BlisNodeDesc *desc = dynamic_cast<BlisNodeDesc *>(desc_);
            desc->setBasis(ws);

            //----------------------------------------------
            // Save variable/constraint bound, non-core constraints
	    // and non-core variable.
	    // The sizes of hard variable bound vectors are numCols.
	    // The sizes of soft variable bound vectors are number
	    // of modified.
            //----------------------------------------------

	    int *tempVarLBPos = model->tempVarLBPos();
	    int *tempVarUBPos = model->tempVarUBPos();
	    //int *tempConLBPos = model->tempConLBPos();
	    //int *tempConUBPos = model->tempConUBPos();
	    
	    int numModSoftColLB = 0;
	    int numModSoftColUB = 0;
	    const double *currColLB = model->solver()->getColLower();
	    const double *currColUB = model->solver()->getColUpper();
	    //const double *currRowLB = model->solver()->getRowLower();
	    //const double *currRowUB = model->solver()->getRowUpper();
	    
	    double *startColLB = model->startVarLB();
	    double *startColUB = model->startVarUB();
	    //double *startRowLB = model->startConLB();
	    //double *startRowUB = model->startConUB();

	    //BlisConstraint **tempCons = NULL;


#ifdef BLIS_DEBUG_MORE
	    // Debug survived old constraints.
	    for (k = 0; k < currNumOldCons; ++k) {
		int oldPos = oldConsPos[k];
		BlisConstraint *aCon = model->oldConstraints()[oldPos];
		assert(aCon);
		std::cout << "SAVE: DBG: oldPos=" << oldPos
			  << ", k=" << k << ", len=" << aCon->getSize()
			  << ", node=" << index_ << std::endl;
	    }
#endif

	    //----------------------------------------------
	    // Decide if save explicit decription.
	    //----------------------------------------------
	    
	    if (depth_ % 30 == 0 || isRoot || (phase == ALPS_PHASE_RAMPUP)) {
		explicit_ = 1;
		std::cout << "SAVE: node "<< index_ <<" explicitly, "
			  << "depth=" << depth_ << std::endl;
	    }
	    else {
		explicit_ = 0;
		//std::cout << "SAVE: node "<< index_ <<" relatively, "
		//  << "depth=" << depth_ << std::endl;
	    }
	    
	    //explicit_ = 1;

	    if (explicit_ || (phase == ALPS_PHASE_RAMPUP) ) {
		// NOTE: full hard bound has been stored. 
		
		int index;
		int numModify = 0;
		int numSoftVarLowers = 0;
		int numSoftVarUppers = 0;
		double value;
		
		double *fVarHardLB = new double [numCols];
		double *fVarHardUB = new double [numCols];
		int *fVarHardLBInd = new int [numCols];
		int *fVarHardUBInd = new int [numCols];
		
		double *fVarSoftLB = new double [numCols];
		double *fVarSoftUB = new double [numCols];
		int *fVarSoftLBInd = new int [numCols];
		int *fVarSoftUBInd = new int [numCols];
		
		for (k = 0; k < numCols; ++k) {
		    fVarSoftLB[k] = ALPS_DBL_MAX;
		    fVarSoftUB[k] = -ALPS_DBL_MAX;
		    fVarHardLB[k] = ALPS_DBL_MAX;
		    fVarHardUB[k] = -ALPS_DBL_MAX;
		    fVarHardLBInd[k] = k;
		    fVarHardUBInd[k] = k;
		}

		//------------------------------------------
		// Build the path to an explicit node.
		//------------------------------------------
		
		AlpsTreeNode *parent = parent_;

		// First push this node since it has branching hard bounds.
		model->leafToRootPath.push_back(this);
		BlisNodeDesc* pathDesc = NULL;
		
		if (phase != ALPS_PHASE_RAMPUP) {
		    while(parent) {
			model->leafToRootPath.push_back(parent);
			if (parent->getExplicit()) {
			    // Reach an explicit node, then stop.
			    break;
			}
			else {
			    parent = parent->getParent();
			}
		    }
		}

#ifdef BLIS_DEBUG_MORE
		std::cout << "SAVE: EXP: path len = "<<model->leafToRootPath.size()
			  << std::endl;
#endif
		//------------------------------------------
		// Summarize bounds.
		//------------------------------------------
		
		for(j = model->leafToRootPath.size() - 1; j > -1; --j) {
		    
		    pathDesc = dynamic_cast<BlisNodeDesc*>
			((model->leafToRootPath.at(j))->getDesc());
		    
		    //--------------------------------------
		    // Full variable hard bounds.
		    //--------------------------------------
		    
		    numModify = pathDesc->getVars()->lbHard.numModify;
		    for (k = 0; k < numModify; ++k) {
			index = pathDesc->getVars()->lbHard.posModify[k];
			value = pathDesc->getVars()->lbHard.entries[k];
			fVarHardLB[index] = value;
		    }
		    
		    numModify = pathDesc->getVars()->ubHard.numModify;
		    for (k = 0; k < numModify; ++k) {
			index = pathDesc->getVars()->ubHard.posModify[k];
			value = pathDesc->getVars()->ubHard.entries[k];
			fVarHardUB[index] = value;
		    }
		    
		    //--------------------------------------
		    // Full variable soft bounds.
		    //--------------------------------------
		    
		    numModify = pathDesc->getVars()->lbSoft.numModify;
#ifdef BLIS_DEBUG_MORE
		    std::cout << "SAVE: EXP: j=" << j << ", numModify soft lb="
			      << numModify << std::endl;
#endif
		    for (k = 0; k < numModify; ++k) {
			index = pathDesc->getVars()->lbSoft.posModify[k];
			value = pathDesc->getVars()->lbSoft.entries[k];
			fVarSoftLB[index] = value;
		    }
		    
		    numModify = pathDesc->getVars()->ubSoft.numModify;
#ifdef BLIS_DEBUG_MORE
		    std::cout << "SAVE: EXP: j=" << j << ", numModify soft ub="
			      << numModify << std::endl;
#endif
		    for (k = 0; k < numModify; ++k) {
			index = pathDesc->getVars()->ubSoft.posModify[k];
			value = pathDesc->getVars()->ubSoft.entries[k];
			fVarSoftUB[index] = value;
		    }

		} // EOF of for(path)

		//------------------------------------------
		// Collect modified soft bounds at this node.
		// NOTE: Do this after collecting previous soft bounds.
		//------------------------------------------
		
		numModSoftColLB = 0;
		numModSoftColUB = 0;
		for (k = 0; k < numCoreCols; ++k) {
		    if (currColLB[k] != startColLB[k]) {
			fVarSoftLB[k] = currColLB[k];
			++numModSoftColLB;
#ifdef BLIS_DEBUG_MORE
			printf("Col %d, soft lb change, start %g, curr %g\n",
			       k, startColLB[k], currColLB[k]);
#endif
			
		    }
		    if (currColUB[k] != startColUB[k]) {
			fVarSoftUB[k] = currColUB[k];
			++numModSoftColUB;
		    }
		}
		
#ifdef BLIS_DEBUG_MORE
		std::cout << "SAVE: EXP: THIS: numModSoftColLB = "<<numModSoftColLB 
			  << ", numModSoftColUB = " << numModSoftColUB << std::endl;
#endif
		
		//--------------------------------------
		// Debug if bounds are consistant.
		//--------------------------------------

#ifdef BLIS_DEBUG
		for (k = 0; k < numCols; ++k) {
		    
		    //std::cout << "EXP: COL[" << k <<"]: " 
		    //      <<"hardLB=" << fVarHardLB[k] 
		    //      <<", hardUB=" << fVarHardUB[k] << std::endl;
		    
		    // Hard lower bound should not greater than
		    // hard upper bound.
		    if (fVarHardLB[k] > fVarHardUB[k] + ALPS_GEN_TOL) {
			printf("ERROR: Col %d, \tlb %g,  \tub %g\n",
			       k, fVarHardLB[k], fVarHardUB[k]);
			assert(0);
		    }
		    
		    if (fVarSoftLB[k] < ALPS_BND_MAX) {
			// Soft lower changed, and should not greater 
			// than its hard upper bound.
			if (fVarSoftLB[k] > fVarHardUB[k] + ALPS_GEN_TOL) {
			    printf("ERROR: Col %d, \tlb %g,  \tub %g\n",
				   k, fVarSoftLB[k], fVarHardUB[k]);
			    assert(0);	
			}
		    }
		    
		    if (fVarSoftUB[k] > -ALPS_BND_MAX) {
			// Soft upper changed, and should not less 
			// than its hard lower bound.
			if (fVarSoftUB[k] < fVarHardLB[k] - ALPS_GEN_TOL) {
			    printf("ERROR: Col %d, \tlb %g,  \tub %g\n",
				   k, fVarHardLB[k], fVarSoftUB[k]);
			    assert(0);
			}
		    }
		}
#endif
		
		//------------------------------------------
		// Record hard variable bounds. FULL set.
		//------------------------------------------
		
		desc->assignVarHardBound(numCols,
					 fVarHardLBInd,
					 fVarHardLB,
					 numCols,
					 fVarHardUBInd,
					 fVarHardUB);
		
		//------------------------------------------
		// Recode soft variable bound. Modified.
		//------------------------------------------
		
		for (k = 0; k < numCols; ++k) {
		    if (fVarSoftLB[k] < ALPS_BND_MAX) {
			fVarSoftLBInd[numSoftVarLowers] = k;
			fVarSoftLB[numSoftVarLowers++] = fVarSoftLB[k];
		    }
		    if (fVarSoftUB[k] > -ALPS_BND_MAX) {
			fVarSoftUBInd[numSoftVarUppers] = k;
			fVarSoftUB[numSoftVarUppers++] = fVarSoftUB[k];
		    }
		}

		
#ifdef BLIS_DEBUG_MORE
		// Print soft bounds.
		std::cout << "SAVE: EXP: numSoftVarLowers=" << numSoftVarLowers
			  << ", numSoftVarUppers=" << numSoftVarUppers
			  << std::endl;
		for (k = 0; k < numSoftVarLowers; ++k) {
		    std::cout << "Col[" << fVarSoftLBInd[k] << "]: soft lb="
			      << fVarSoftLB[k] << std::endl;		    
		}
		std::cout << "------------------" << std::endl;
		for (k = 0; k < numSoftVarUppers; ++k) {
		    std::cout << "Col[" << fVarSoftUBInd[k] << "]: soft ub="
			      << fVarSoftUB[k] << std::endl;
		}
		std::cout << "------------------" << std::endl << std::endl;
#endif

		//if ( (numSoftVarUppers > 0) || (numSoftVarLowers > 0) ) {
                
                // Assign it anyway so to delete memory(fVarSoftLBInd,etc.)
                desc->assignVarSoftBound(numSoftVarLowers, 
                                         fVarSoftLBInd,
                                         fVarSoftLB,
                                         numSoftVarUppers, 
                                         fVarSoftUBInd,
                                         fVarSoftUB);

                //------------------------------------------
                // Full set of active non-core constraints.
                //------------------------------------------
                // old constraint: model->oldConstraints_, currNumOldCons.
                // new constraint: newConstraints, numNewCons.

                BcpsObject **toAddCons = new BcpsObject * [currNumOldCons + 
                                                           numNewCons];
                if (currNumOldCons > 0) {
                    // Hard copy of the survived old constraints.
                    for (k = 0; k < currNumOldCons; ++k) {
			int oldPos = oldConsPos[k];
                        BlisConstraint *aCon = model->oldConstraints()[oldPos];
                        assert(aCon);
#ifdef BLIS_DEBUG
			std::cout << "SAVE: EXP: currNumOldCons=" << currNumOldCons 
				  << ", k=" << k << ", len=" << aCon->getSize()
                                  << ", node=" << index_ << std::endl;
#endif
			
                        BlisConstraint *newCon = new BlisConstraint(*aCon);
                        toAddCons[k] = newCon;
                    }
                }
                
		if (numNewCons > 0) {
		    // Include new constraints.
		    memcpy(toAddCons + currNumOldCons, 
			   newConstraints,
			   numNewCons * sizeof(BcpsObject *));
		}

		//------------------------------------------
                // Save in description. Add first delete exiting, then add.
		// Pointers in model->oldConstraints_ can be dangling.
		// It is not safe to use model->oldConstraints_ after adding.
		
		// If this node is the root of a subtree, and before processing
		// it has a list of cuts, then model->oldConstraints_ 
		// stores pointer to the cuts when installing.
		
		// Need update model->constraints_ here OR do not be smart!
		// 1/6/06: Choose to be dumn.
		//------------------------------------------

		//------------------------------------------
		// Generating constraints,
		// also means that slack ones might been removed.
		//------------------------------------------

		int numTotal = currNumOldCons + numNewCons;
                desc->setAddedConstraints(numTotal, toAddCons);
                
#ifdef BLIS_DEBUG
		std::cout << "SAVE: EXP: currNumOldCons=" << currNumOldCons
			  << ", numNewCons=" << numNewCons
			  << std::endl;
#endif	
		
                //------------------------------------------
                // Full set of active non-core variables.
                //------------------------------------------
                // Todo.
                
		//------------------------------------------
		// Clear path vector.
		//------------------------------------------
		
		model->leafToRootPath.clear();
		assert(model->leafToRootPath.size() == 0);
	    }
	    else { // Relative.
		//------------------------------------------
		// Record soft bound changes for core vars.
		//------------------------------------------

		// Variable bound change 
		numModSoftColLB = 0;
		numModSoftColUB = 0;
		for (k = 0; k < numCoreCols; ++k) {
		    if (currColLB[k] != startColLB[k]) {
			tempVarLBPos[numModSoftColLB] = k;
			/* startColLB as a temporary storage vector */
			startColLB[numModSoftColLB] = currColLB[k];
			
#ifdef BLIS_DEBUG_MORE
			printf("Col %d, soft lb change, start %g, curr %g\n",
			       k, startColLB[k], currColLB[k]);
#endif
			
			++numModSoftColLB;
		    }
		    if (currColUB[k] != startColUB[k]) {
			tempVarUBPos[numModSoftColUB] = k;
			startColUB[numModSoftColUB] = currColUB[k];
			++numModSoftColUB;
		    }
		}

#ifdef BLIS_DEBUG_MORE
		std::cout << "SAVE: REL: numModSoftColLB = " << numModSoftColLB 
			  << ", numModSoftColUB = " << numModSoftColUB 
			  << std::endl;
#endif
            
		if (numModSoftColLB > 0 || numModSoftColUB > 0) {
#ifdef BLIS_DEBUG
		    //assert(0);
#endif
		    desc->setVarSoftBound(numModSoftColLB, 
					  tempVarLBPos,
					  startColLB,
					  numModSoftColUB, 
					  tempVarUBPos,
					  startColUB);
		}
            
		//------------------------------------------
		// TODO: Constraint bounds change.
		//------------------------------------------

#if 0
		for (k = 0; k < numCoreRows; ++k) {
		    if (currRowLB[k] != startRowLB[k]) {
			tempConLBPos[numModSoftRowLB] = k;
			startRowLB[numModSoftRowLB] = currRowLB[k];
			++numModSoftRowLB;
		    }
		    if (currRowUB[k] != startRowUB[k]) {
			tempConUBPos[numModSoftRowUB] = k;
			startRowUB[numModSoftRowUB] = currRowUB[k];
			++numModSoftRowUB;
		    }
		}
		if (numModSoftRowLB > 0 || numModSoftRowUB > 0) {
		    desc->setConSoftBound(numModSoftRowLB, 
					  tempConLBPos,
					  startRowLB,
					  numModSoftRowUB, 
					  tempConUBPos,
					  startRowUB);
		}
#endif	
            
		if (genConsHere) {
		    // NOTE: explicit_ can NOT do this if, since genConsHere maybe
		    //       false here, but there are maybe cons from parents.

		    //--------------------------------------
		    // Record add constraints.
		    //--------------------------------------

		    // Apend will copy old, then add new. 
		    // If this node has a list of cuts before pointers in 
		    // model->oldConstraints() will be kept. Safe!
		    if (numNewCons > 0) {
			desc->appendAddedConstraints(numNewCons, newConstraints);
                    }
		    
		    //--------------------------------------
		    // Record deleted constraint positions.
		    //--------------------------------------
		    
		    int *oldLeft = new int [origNumOldCons];
		    int leftCon;
		    CoinZeroN(oldLeft, origNumOldCons);
		    
		    for (k = 0; k < currNumOldCons; ++k) {
			leftCon = oldConsPos[k];
			assert(leftCon >= 0 && leftCon < origNumOldCons);
			oldLeft[leftCon] = 1;
		    }
		    // 
		    leftCon = 0;
		    for (k = 0; k < origNumOldCons; ++k) {
			if (oldLeft[k] == 0) {
			    // Deleted. Now oldLeft stores delete position.
			    oldLeft[leftCon++] = k;
			}
                        //FIXME: clean
                        //assert(k < 15196924);
		    }
		    desc->delConstraints(leftCon, oldLeft);
		    
#ifdef BLIS_DEBUG
		    std::cout << "PROCESS: ADD: new cuts=" << numNewCons 
			      << ", numRows=" << model->solver()->getNumRows() 
			      << ", numStartRows="<< numStartRows 
                              << ", origNumStartRows="<< origNumStartRows 
                              << ", num removed=" << leftCon << std::endl;                    
#endif
		    
		}// EOF of if(genConsHere)
	    } // EOF of relative
        }
        else if (bStatus == -2) {
            
#ifdef BLIS_DEBUG_MORE
            std::cout << "bStatus = -2" << std::endl;
#endif
            //branchObject->getDown()[0], branchObject->getDown()[1]);
	    
            setStatus(AlpsNodeStatusFathomed);
        }
        else {
            throw CoinError("No branch object found", "process", 
                            "BlisTreeNode");
        }
    }
    
    //------------------------------------------------------
    // End of process()
    //------------------------------------------------------
    
 TERM_PROCESS:

    delete [] heurSolution;
    delete [] currLpSolution;

    if (status_ == AlpsNodeStatusFathomed) {
	// Delete new cuts since no use anymore.
        for (k = 0; k < numNewCons; ++k) {
            delete newConstraints[k];
        }
    }
    delete [] newConstraints;
    delete [] oldConsPos;
    
    model->isRoot_ = false;

#ifdef BLIS_DEBUG_MORE
    // Debug survived old constraints.
    //int currNumOldCons = model->getNumOldConstraints();
    for (k = 0; k < currNumOldCons; ++k) {
	BlisConstraint *aCon = model->oldConstraints()[k];
	assert(aCon);
	std::cout << "SAVE: DBG: TERM: "
		  << "k=" << k << ", len=" << aCon->getSize()
		  << ", node=" << index_ << std::endl;
    }
#endif
    

    return status;
}

//#############################################################################

std::vector< CoinTriple<AlpsNodeDesc*, AlpsNodeStatus, double> >
BlisTreeNode::branch()
{

    //------------------------------------------------------
    // Change one var hard bound and record the change in nodedesc:
    // THINK: how about constraint bounds? When to update?
    // TODO: how about other SOS object, etc.?
    //------------------------------------------------------

    AlpsPhase phase = knowledgeBroker_->getPhase();

    double objVal = getQuality();

    BlisNodeDesc* childDesc = NULL;  
      
    std::vector< CoinTriple<AlpsNodeDesc*, AlpsNodeStatus, double> > 
	childNodeDescs;
    
    BlisModel* model = dynamic_cast<BlisModel*>(desc_->getModel());    

    int numCols = model->getNumCols();


#ifdef BLIS_DEBUG_MORE
    // Debug survived old constraints.
    int currNumOldCons = model->getNumOldConstraints();
    for (int k = 0; k < currNumOldCons; ++k) {
	BlisConstraint *aCon = model->oldConstraints()[k];
	assert(aCon);
	std::cout << "BRANCH: DBG: "
		  << "k=" << k << ", len=" << aCon->getSize()
		  << ", node=" << index_ << std::endl;
    }
#endif

    //------------------------------------------------------
    // Get branching object. TODO: Assume integer branching object. 
    //------------------------------------------------------
    
    BlisBranchObjectInt *branchObject =
	dynamic_cast<BlisBranchObjectInt *>(branchObject_);
    int objInd = branchObject->getObjectIndex();

    double bValue = branchObject->getValue();

    BlisObjectInt *obj = dynamic_cast<BlisObjectInt *>(model->objects(objInd));
    int branchVar = obj->columnIndex();
    
#ifdef BLIS_DEBUG
    if ( (branchVar < 0) || (branchVar >= numCols) ) {
	std::cout << "ERROR: BRANCH(): branchVar = " << branchVar 
		  << "; numCols = " << numCols  << std::endl;
	throw CoinError("branch index is out of range", 
			"branch", "BlisTreeNode");
    }
#endif

#ifdef BLIS_DEBUG
    printf("BRANCH(): on %d, phase %d\n", branchVar, phase);
    printf("DOWN: lb %g, up %g\n",
	   branchObject->getDown()[0], branchObject->getDown()[1]);
    printf("UP  : lb %g, up %g\n",
	   branchObject->getUp()[0], branchObject->getUp()[1]);
#endif

    BlisNodeDesc* thisDesc = dynamic_cast<BlisNodeDesc*>(desc_);

        
    //======================================================
    //------------------------------------------------------
    // Create down-branch node description.
    //------------------------------------------------------
    //======================================================
    
    childDesc = new BlisNodeDesc(model);
    
    if (phase == ALPS_PHASE_RAMPUP) {
	
	//--------------------------------------------------
	// Store a full description since each node will be the root of
	// a subtree.
	// NOTE: this desc must be explicit during rampup.
	//--------------------------------------------------	

	int index, k;
	int numModify = -1;
	double value;
	
	double *fVarHardLB = new double [numCols];
	double *fVarHardUB = new double [numCols];
	int *fVarHardLBInd = new int [numCols];
	int *fVarHardUBInd = new int [numCols];
		
	double *fVarSoftLB = NULL;
	double *fVarSoftUB = NULL;
	int *fVarSoftLBInd = NULL;
	int *fVarSoftUBInd = NULL;

	//--------------------------------------------------
	// Full hard variable bounds.
	//--------------------------------------------------

	numModify = thisDesc->getVars()->lbHard.numModify;
	assert(numModify == numCols);
	for (k = 0; k < numModify; ++k) {
	    index = thisDesc->getVars()->lbHard.posModify[k];
	    assert(index == k);
	    value = thisDesc->getVars()->lbHard.entries[k];
	    fVarHardLB[k] = value;
	    fVarHardLBInd[k] = index;
	}
	
	numModify = thisDesc->getVars()->ubHard.numModify;
	assert(numModify == numCols);
	for (k = 0; k < numModify; ++k) {
	    index = thisDesc->getVars()->ubHard.posModify[k];
	    assert(index == k);
	    value = thisDesc->getVars()->ubHard.entries[k];
	    fVarHardUB[k] = value;
	    fVarHardUBInd[k] = index;
	}
	
	// Branching bounds.
	fVarHardLB[branchVar] = branchObject->getDown()[0];
	fVarHardUB[branchVar] = branchObject->getDown()[1];


	childDesc->assignVarHardBound(numCols,
				      fVarHardLBInd,
				      fVarHardLB,
				      numCols,
				      fVarHardUBInd,
				      fVarHardUB);

	//--------------------------------------------------
	// Soft variable bounds.
	//--------------------------------------------------

	int numSoftVarLowers = thisDesc->getVars()->lbSoft.numModify;
	assert(numSoftVarLowers >= 0 && numSoftVarLowers <= numCols);
	if (numSoftVarLowers > 0) {
	    fVarSoftLB = new double [numSoftVarLowers];
	    fVarSoftLBInd = new int [numSoftVarLowers];
	    for (k = 0; k < numSoftVarLowers; ++k) {
		index = thisDesc->getVars()->lbSoft.posModify[k];
		value = thisDesc->getVars()->lbSoft.entries[k];
		fVarSoftLB[k] = value;
		fVarSoftLBInd[k] = index;
	    }
	}
		    
	int numSoftVarUppers = thisDesc->getVars()->ubSoft.numModify;
	assert(numSoftVarUppers >= 0 && numSoftVarUppers <= numCols);
	if (numSoftVarUppers > 0) {
	    fVarSoftUB = new double [numSoftVarUppers];
	    fVarSoftUBInd = new int [numSoftVarUppers];
	    for (k = 0; k < numSoftVarUppers; ++k) {
		index = thisDesc->getVars()->ubSoft.posModify[k];
		value = thisDesc->getVars()->ubSoft.entries[k];
		fVarSoftUB[k] = value;
		fVarSoftUBInd[k] = index;
	    }
	}

#ifdef BLIS_DEBUG_MORE
	// Print soft bounds.
	std::cout << "\nBRANCH: numSoftVarLowers=" << numSoftVarLowers
		  << ", numSoftVarUppers=" << numSoftVarUppers
		  << std::endl;
	for (k = 0; k < numSoftVarLowers; ++k) {
	    std::cout << "Col[" << fVarSoftLBInd[k] << "]: soft lb="
		      << fVarSoftLB[k] << std::endl;		    
	}
	std::cout << "------------------" << std::endl;
	for (k = 0; k < numSoftVarUppers; ++k) {
	    std::cout << "Col[" << fVarSoftUBInd[k] << "]: soft ub="
		      << fVarSoftUB[k] << std::endl;
	}
	std::cout << "------------------" << std::endl << std::endl;
#endif

	// Assign it anyway so to transfer ownership of memory(fVarSoftLBInd,etc.)
	childDesc->assignVarSoftBound(numSoftVarLowers, 
				      fVarSoftLBInd,
				      fVarSoftLB,
				      numSoftVarUppers, 
				      fVarSoftUBInd,
				      fVarSoftUB);

	//--------------------------------------------------
	// Full set of non-core constraints.
	// NOTE: non-core constraints have been saved in description
	//       when process() during ramp-up.
	//--------------------------------------------------

	BcpsObject **tempCons = NULL;
	int tempInt = 0;
	
	tempInt = thisDesc->getCons()->numAdd;
	if (tempInt > 0) {
	    tempCons = new BcpsObject* [tempInt];
	    for (k = 0; k < tempInt; ++k) {
		BlisConstraint *aCon = dynamic_cast<BlisConstraint *>
		    (thisDesc->getCons()->objects[k]);
                
		assert(aCon);
		assert(aCon->getSize() > 0);
		assert(aCon->getSize() < 100000);
		BlisConstraint *newCon = new BlisConstraint(*aCon);
		tempCons[k] = newCon;
	    }
	}
	    

#if 0
	else {
	    // No cons or only root cons.
	    tempInt = model->getNumOldConstraints();

	    if (tempInt > 0) {
		tempCons = new BcpsObject* [tempInt];
	    }
	    for (k = 0; k < tempInt; ++k) {
		BlisConstraint *aCon = model->oldConstraints()[k];                
		assert(aCon);
		assert(aCon->getSize() > 0);
		assert(aCon->getSize() < 100000);
		BlisConstraint *newCon = new BlisConstraint(*aCon);
		tempCons[k] = newCon;
	    }
	}
#endif

#ifdef BLIS_DEBUG_MORE
	std::cout << "BRANCH: down: tempInt=" << tempInt <<std::endl;
#endif
	// Fresh desc, safely add.
	childDesc->setAddedConstraints(tempInt, tempCons);
    }
    else {
	
	//--------------------------------------------------
	// Relative: Only need to record hard var bound change. 
	// NOTE: soft var bound changes are Record after selectBranchObject.
	//--------------------------------------------------
	
	childDesc->setVarHardBound(1,
				   &branchVar,
				   &(branchObject->getDown()[0]),
				   1,
				   &branchVar,
				   &(branchObject->getDown()[1]));
    }

    childDesc->setBranchedDir(-1);
    childDesc->setBranchedInd(objInd);
    childDesc->setBranchedVal(bValue);

    // Copy warm start.
    CoinWarmStartBasis *ws = thisDesc->getBasis();
    CoinWarmStartBasis *newWs = new CoinWarmStartBasis(*ws);
    childDesc->setBasis(newWs);
    
    childNodeDescs.push_back(CoinMakeTriple(static_cast<AlpsNodeDesc *>
					    (childDesc),
					    AlpsNodeStatusCandidate,
					    objVal));

    //======================================================
    //------------------------------------------------------
    // Create up-branch node description.
    //------------------------------------------------------
    //======================================================
    
    childDesc = new BlisNodeDesc(model);

    if (phase == ALPS_PHASE_RAMPUP) {

	//--------------------------------------------------
	// Store a full description since each node will be the root of
	// a subtree.
	// NOTE: parent must be explicit during rampup.
	//--------------------------------------------------	

	int index, k;
	int numModify = -1;
	double value;
	
	double *fVarHardLB = new double [numCols];
	double *fVarHardUB = new double [numCols];
	int *fVarHardLBInd = new int [numCols];
	int *fVarHardUBInd = new int [numCols];
		
	double *fVarSoftLB = NULL;
	double *fVarSoftUB = NULL;
	int *fVarSoftLBInd = NULL;
	int *fVarSoftUBInd = NULL;

	//--------------------------------------------------
	// Full hard variable bounds.
	//--------------------------------------------------
	
	numModify = thisDesc->getVars()->lbHard.numModify;
	assert(numModify == numCols);
	for (k = 0; k < numModify; ++k) {
	    index = thisDesc->getVars()->lbHard.posModify[k];
	    assert(index == k);
	    value = thisDesc->getVars()->lbHard.entries[k];
	    fVarHardLB[k] = value;
	    fVarHardLBInd[k] = index;
	}
	
	numModify = thisDesc->getVars()->ubHard.numModify;
	assert(numModify == numCols);
	for (k = 0; k < numModify; ++k) {
	    index = thisDesc->getVars()->ubHard.posModify[k];
	    assert(index == k);
	    value = thisDesc->getVars()->ubHard.entries[k];
	    fVarHardUB[k] = value;
	    fVarHardUBInd[k] = index;
	}
	
	// Branching bounds.
	fVarHardLB[branchVar] = branchObject->getUp()[0];
	fVarHardUB[branchVar] = branchObject->getUp()[1];

	childDesc->assignVarHardBound(numCols,
				      fVarHardLBInd,
				      fVarHardLB,
				      numCols,
				      fVarHardUBInd,
				      fVarHardUB);

	//--------------------------------------------------
	// Soft variable bounds.
	//--------------------------------------------------

	int numSoftVarLowers = thisDesc->getVars()->lbSoft.numModify;
	assert(numSoftVarLowers >= 0 && numSoftVarLowers <= numCols);
	if (numSoftVarLowers > 0) {
	    fVarSoftLB = new double [numSoftVarLowers];
	    fVarSoftLBInd = new int [numSoftVarLowers];
	    for (k = 0; k < numSoftVarLowers; ++k) {
		index = thisDesc->getVars()->lbSoft.posModify[k];
		value = thisDesc->getVars()->lbSoft.entries[k];
		fVarSoftLB[k] = value;
		fVarSoftLBInd[k] = index;
	    }
	}
		    
	int numSoftVarUppers = thisDesc->getVars()->ubSoft.numModify;
	assert(numSoftVarUppers >= 0 && numSoftVarUppers <= numCols);
	if (numSoftVarUppers > 0) {
	    fVarSoftUB = new double [numSoftVarUppers];
	    fVarSoftUBInd = new int [numSoftVarUppers];
	    for (k = 0; k < numSoftVarUppers; ++k) {
		index = thisDesc->getVars()->ubSoft.posModify[k];
		value = thisDesc->getVars()->ubSoft.entries[k];
		fVarSoftUB[k] = value;
		fVarSoftUBInd[k] = index;
	    }
	}

	// Assign it anyway so to transfer ownership of memory(fVarSoftLBInd,etc.)
	childDesc->assignVarSoftBound(numSoftVarLowers, 
				      fVarSoftLBInd,
				      fVarSoftLB,
				      numSoftVarUppers, 
				      fVarSoftUBInd,
				      fVarSoftUB);

	//--------------------------------------------------
	// Full set of non-core constraints.
	// NOTE: non-core constraints have been saved in description
	//       when process() during ramp-up.
	//--------------------------------------------------
	
	BcpsObject **tempCons = NULL;
	int tempInt = 0;
	
	tempInt = thisDesc->getCons()->numAdd;
	if (tempInt > 0) {
	    tempCons = new BcpsObject* [tempInt];
	
	    for (k = 0; k < tempInt; ++k) {
		BlisConstraint *aCon = dynamic_cast<BlisConstraint *>
		    (thisDesc->getCons()->objects[k]);
	    
		assert(aCon);
		assert(aCon->getSize() > 0);
		assert(aCon->getSize() < 1000);
		BlisConstraint *newCon = new BlisConstraint(*aCon);
		tempCons[k] = newCon;
	    }
	}
	
#if 0
	else {
	    // No cons or only root cons.
	    tempInt = model->getNumOldConstraints();
	    if (tempInt > 0) {
		tempCons = new BcpsObject* [tempInt];
	    }
	    for (k = 0; k < tempInt; ++k) {
		BlisConstraint *aCon = model->oldConstraints()[k];                
		assert(aCon);
		assert(aCon->getSize() > 0);
		assert(aCon->getSize() < 1000);
		BlisConstraint *newCon = new BlisConstraint(*aCon);
		tempCons[k] = newCon;
	    }
	}
#endif

#ifdef BLIS_DEBUG_MORE
	std::cout << "BRANCH: up: tempInt=" << tempInt <<std::endl;
#endif
	// Fresh desc, safely add.
	childDesc->setAddedConstraints(tempInt, tempCons);
    }
    else {

	//--------------------------------------------------
	// Relative: Only need to record hard var bound change. 
	// NOTE: soft var bound changes are Record after selectBranchObject.
	//--------------------------------------------------

	childDesc->setVarHardBound(1,
				   &branchVar,
				   &(branchObject->getUp()[0]),
				   1,
				   &branchVar,
				   &(branchObject->getUp()[1]));
    }

    childDesc->setBranchedDir(1);
    childDesc->setBranchedInd(objInd);
    childDesc->setBranchedVal(bValue);
    
    // Copy warm start.
    CoinWarmStartBasis *newWs2 = new CoinWarmStartBasis(*ws);
    childDesc->setBasis(newWs2);
    
    childNodeDescs.push_back(CoinMakeTriple(static_cast<AlpsNodeDesc *>
					    (childDesc),
					    AlpsNodeStatusCandidate,
					    objVal));  

    // Change node status to branched.
    status_ = AlpsNodeStatusBranched;
    
    return childNodeDescs;
}

//#############################################################################

/* FIXME: need rewrite from scratch */
/* 0: find a branch var, -1 no branch var (should not happen) */

int BlisTreeNode::selectBranchObject(BlisModel *model, 
                                     bool& foundSol, 
                                     int numPassesLeft) 
{
    int bStatus = 0;

    if(branchObject_) {
        delete branchObject_;
        branchObject_ = NULL;
    }
    
    //------------------------------------------------------
    // Get branching strategy.
    //------------------------------------------------------

    BcpsBranchStrategy *strategy = model->branchStrategy();
    if (!strategy) {
        throw CoinError("No branch strategy.", "process()","BlisTreeNode");
    }

    //------------------------------------------------------
    // Create branching object candidates.
    //-----------------------------------------------------
    
    bStatus = strategy->createCandBranchObjects(numPassesLeft);   
    
    //------------------------------------------------------
    // Select the best branching objects.
    //-----------------------------------------------------
    
    if (bStatus >= 0) {
        
        int bestIndex = strategy->bestBranchObject();
        
        if (bestIndex >= 0) {
            // Move best branching object to node.
            branchObject_ = strategy->getBestBranchObject();

#ifdef BLIS_DEBUG_MORE
            std::cout << "SELECTBEST: Set branching obj" << std::endl;
#endif
        }
        else {
#ifdef BLIS_DEBUG
            std::cout << "WARNING: Can't find best branching obj" << std::endl;
            assert(0);
#endif      
        }
        
        // Set guessed solution value
        // solEstimate_ = quality_ + sumDeg;
    }
    
    if (!model->branchStrategy()) {
        delete strategy;
    }

    return bStatus;
}

//#############################################################################

int BlisTreeNode::bound(BcpsModel *model) 
{
    int status = BLIS_OK;
    BlisModel *m = dynamic_cast<BlisModel *>(model);
    
#ifdef BLIS_DEBUG_MORE
    int j;
    int numCols = m->solver()->getNumCols();
    const double * clb = m->solver()->getColLower();
    const double * cub = m->solver()->getColUpper();

    for (j = 0; j < numCols; ++j) {
	std::cout << "c"<< j <<"["<<clb[j]<<", "<< cub[j] << "]" << std::endl;
    }
#endif

    m->solver()->resolve();

    if (m->solver()->isAbandoned()) {
#ifdef BLIS_DEBUG
	std::cout << "BOUND: is abandoned" << std::endl;
#endif
	status = BLIS_LP_ABANDONED;
    }
    else if (m->solver()->isProvenOptimal()) {
#ifdef BLIS_DEBUG
	std::cout << "BOUND: is lp optimal" << std::endl;
#endif
	status = BLIS_LP_OPTIMAL;
        BlisNodeDesc *desc = dynamic_cast<BlisNodeDesc*>(desc_);

        double objValue = m->solver()->getObjValue() *
            m->solver()->getObjSense();
        
        int dir = desc->getBranchedDir();
        if (dir != 0) {
            double objDeg = objValue - quality_;
            int objInd = desc->getBranchedInd();
            double lpX = desc->getBranchedVal();
            BlisObjectInt *intObject = 
                dynamic_cast<BlisObjectInt *>(m->objects(objInd));            
#ifdef BLIS_DEBUG_MORE
            std::cout << "BOUND: col[" << intObject->columnIndex() 
                      << "], dir=" << dir << ", objDeg=" << objDeg
                      << ", x=" << lpX
                      << ", up=" << intObject->pseudocost().getUpCost()
                      << ", down=" << intObject->pseudocost().getDownCost()
                      << ", pre quality=" << quality_
                      << ", objValue=" << objValue
                      << std::endl;
#endif

            intObject->pseudocost().update(dir, objDeg, lpX);
        }

        // Update quality of this nodes.
        quality_ = objValue;
    }
    else if (m->solver()->isProvenPrimalInfeasible()) {
#ifdef BLIS_DEBUG
	std::cout << "BOUND: is primal inf" << std::endl;
#endif
	status = BLIS_LP_PRIMAL_INF;
    }
    else if (m->solver()->isProvenDualInfeasible()) {
#ifdef BLIS_DEBUG
	std::cout << "BOUND: is dual inf" << std::endl;
#endif
	status = BLIS_LP_DUAL_INF;
    }
    else if (m->solver()->isPrimalObjectiveLimitReached()) {
#ifdef BLIS_DEBUG
	std::cout << "BOUND: is primal limit" << std::endl;
#endif
	status = BLIS_LP_PRIMAL_LIM;
    }
    else if (m->solver()->isDualObjectiveLimitReached()) {
#ifdef BLIS_DEBUG
	std::cout << "BOUND: is dual limit" << std::endl;
#endif
	status = BLIS_LP_DUAL_LIM;
    }
    else if (m->solver()->isIterationLimitReached()) {
#ifdef BLIS_DEBUG
	std::cout << "BOUND: is iter limit" << std::endl;
#endif
	status = BLIS_LP_ITER_LIM;
    }
    else {
	std::cout << "UNKNOWN LP STATUS" << std::endl;
	assert(0);
    }
    
    return status;
}

//#############################################################################

int BlisTreeNode::installSubProblem(BcpsModel *m)
{
    AlpsReturnCode status = ALPS_OK;

    int i, k;
    int index;
    double value;

    BlisModel *model = dynamic_cast<BlisModel *>(m);
    assert(model);
    
    BlisNodeDesc *desc = dynamic_cast<BlisNodeDesc*>(desc_);

    int numModify = 0;

    int numCoreVars = model->getNumCoreVariables();
    int numCoreCons = model->getNumCoreConstraints();

    int numCols = model->solver()->getNumCols();
    int numRows = model->solver()->getNumRows();

    //double *varSoftLB = NULL;
    //double *varSoftUB = NULL;
    double *varHardLB = NULL;
    double *varHardUB = NULL;

    //double *conSoftLB = NULL;
    //double *conSoftUB = NULL;
    //double *conHardLB = NULL;
    //double *conHardUB = NULL;
    
    double *startColLB = model->startVarLB();
    double *startColUB = model->startVarUB();
    double *startRowLB = model->startConLB();
    double *startRowUB = model->startConUB();
    
    CoinFillN(startColLB, numCoreVars, -ALPS_DBL_MAX);
    CoinFillN(startColUB, numCoreVars, ALPS_DBL_MAX);
    CoinFillN(startRowLB, numCoreCons, -ALPS_DBL_MAX);
    CoinFillN(startRowUB, numCoreCons, ALPS_DBL_MAX);

#if 0
    //memcpy(startColLB, model->origVarLB(), sizeof(double) * numCoreVars);
    //memcpy(startColUB, model->origVarUB(), sizeof(double) * numCoreVars);
    //memcpy(startRowLB, model->origConLB(), sizeof(double) * numCoreCons);
    //memcpy(startRowUB, model->origConUB(), sizeof(double) * numCoreCons);
#endif

    int numOldCons = 0;
    int tempInt = 0;
    BlisConstraint *aCon = NULL;

    int nodeID = -1;
    nodeID = getIndex();
    //std::cout << "nodeID=" << nodeID << std::endl;

    AlpsPhase phase = knowledgeBroker_->getPhase();

    //======================================================
    // Restore subproblem: 
    //  1. Remove noncore var/con
    //  2. Travel back to root and correct differencing to
    //     full var/con bounds into model->startXXX
    //  3. Set col bounds
    //  4. Set row bounds (is this necessary?)
    //  5. Add contraints except cores
    //  6. Add variables except cores
    //  7. Set basis (should not need modify)
    //======================================================


    //------------------------------------------------------
    // Remove old constraints from lp solver.
    //------------------------------------------------------

    int numDelCons = numRows - numCoreCons;
    
#ifdef BLIS_DEBUG
    std::cout << "INSTALL: numDelCons = " << numDelCons << std::endl;
#endif
	
    if (numDelCons > 0) {
	int *indices = new int [numDelCons];
	if (indices == NULL) {
	    throw CoinError("Out of memory", "installSubProblem", "BlisTreeNode");
	}
	
	for (i = 0; i < numDelCons; ++i) {
	    indices[i] = numCoreCons + i;
	}
	
	model->solver()->deleteRows(numDelCons, indices);
	delete [] indices;
	indices = NULL;
    }
    
    //--------------------------------------------------------
    // Travel back to a full node, then collect diff (add/rem col/row,
    // hard/soft col/row bounds) from the node full to this node.
    //----------------------------
    // Note: if we store full set of logic/agorithm col/row, then
    //       no diff are needed for col/row
    //--------------------------------------------------------
    
    //--------------------------------------------------------
    // Collect differencing bounds. Branching bounds of this node
    // are ALSO collected.
    //--------------------------------------------------------

    BlisNodeDesc* pathDesc = NULL;
    AlpsTreeNode *parent = parent_;    
    
    /* First push this node since it has branching hard bounds. 
       NOTE: during rampup, this desc has full description when branch(). */
    model->leafToRootPath.push_back(this);
    
    if (phase != ALPS_PHASE_RAMPUP) {
	while(parent) {
#ifdef BLIS_DEBUG_MORE
	    std::cout << "Parent id = " << parent->getIndex() << std::endl;
#endif     
	    model->leafToRootPath.push_back(parent);
	    if (parent->getExplicit()) {
		// Reach an explicit node, then stop.
		break;
	    }
	    else {
		parent = parent->getParent();
	    }
	}
    }
    
#ifdef BLIS_DEBUG_MORE
    std::cout << "INSTALL: path len = " << model->leafToRootPath.size()
	      << std::endl;
#endif
    
    //------------------------------------------------------
    // Travel back from this node to the explicit node to
    // collect full description.
    //------------------------------------------------------
    
    for(i = model->leafToRootPath.size() - 1; i > -1; --i) {

#ifdef BLIS_DEBUG_MORE
        if (index_ == 3487) {
            std::cout << "\n----------- NODE ------------" 
                      << model->leafToRootPath.at(i)->getIndex() << std::endl;
        }
#endif
	
	//--------------------------------------------------
	// NOTE: As away from explicit node, bounds become 
	//       tighter and tighter.
	//--------------------------------------------------
        
        pathDesc = dynamic_cast<BlisNodeDesc*>((model->leafToRootPath.at(i))->
                                               getDesc());
        
        varHardLB = pathDesc->getVars()->lbHard.entries;
        varHardUB = pathDesc->getVars()->ubHard.entries;
        
	//--------------------------------------------------      
        // Adjust bounds according to hard var lb/ub.
	// If rampup or explicit, collect hard bounds so far.
	//--------------------------------------------------
	
        numModify = pathDesc->getVars()->lbHard.numModify;
	
#ifdef BLIS_DEBUG_MORE
	std::cout << "INSTALL: numModify lb hard = " << numModify << std::endl;
#endif
	
        for (k = 0; k < numModify; ++k) {
            index = pathDesc->getVars()->lbHard.posModify[k];
            value = pathDesc->getVars()->lbHard.entries[k];
      
#ifdef BLIS_DEBUG_MORE
	    printf("INSTALL: 1, col %d, value %g, startColLB %x\n", 
		   index, value, startColLB);
#endif      
	    // Hard bounds do NOT change according to soft bounds, so
	    // here need std::max.
            startColLB[index] = std::max(startColLB[index], value);
            
#ifdef BLIS_DEBUG_MORE
            if (index_ == 3487) {
                printf("INSTALL: 1, col %d, hard lb %g, ub %g\n", 
                       index, startColLB[index], startColUB[index]);
            }
#endif
        }

#ifdef BLIS_DEBUG_MORE
	std::cout << "INSTALL: numModify ub hard = " << numModify<<std::endl;
#endif

        numModify = pathDesc->getVars()->ubHard.numModify;
        for (k = 0; k < numModify; ++k) {
            index = pathDesc->getVars()->ubHard.posModify[k];
            value = pathDesc->getVars()->ubHard.entries[k];
            startColUB[index] = std::min(startColUB[index], value);
	    
#ifdef BLIS_DEBUG_MORE
            if (index_ == 3487) {
                printf("INSTALL: 2, col %d, hard lb %g, ub %g\n", 
                       index, startColLB[index], startColUB[index]);
            }
            if (startColLB[index] > startColUB[index]) {
                    //assert(0);
            }
#endif
        }
        
        //--------------------------------------------------
        // Adjust bounds according to soft var lb/ub.
	// If rampup or explicit, collect soft bounds so far.
        //--------------------------------------------------

        numModify = pathDesc->getVars()->lbSoft.numModify;
#ifdef BLIS_DEBUG_MORE
	std::cout << "INSTALL: i=" << i << ", numModify soft lb="
		  << numModify << std::endl;
#endif
        for (k = 0; k < numModify; ++k) {
            index = pathDesc->getVars()->lbSoft.posModify[k];
            value = pathDesc->getVars()->lbSoft.entries[k];
            startColLB[index] = std::max(startColLB[index], value);
            
#ifdef BLIS_DEBUG_MORE
            if (index_ == 3487) {
                printf("INSTALL: 3, col %d, soft lb %g, ub %g\n", 
                       index, startColLB[index], startColUB[index]);
            }
            
            if (startColLB[index] > startColUB[index]) {
		//assert(0);
            }
#endif
        }
        numModify = pathDesc->getVars()->ubSoft.numModify;
        
#ifdef BLIS_DEBUG_MORE
	std::cout << "INSTALL: i=" << i << ", numModify soft ub="
		  << numModify << std::endl;
#endif

        for (k = 0; k < numModify; ++k) {
            index = pathDesc->getVars()->ubSoft.posModify[k];
            value = pathDesc->getVars()->ubSoft.entries[k];
            startColUB[index] = std::min(startColUB[index], value);
            
#ifdef BLIS_DEBUG_MORE
            if (index_ == 3487) {
                printf("INSTALL: 4, col %d, soft lb %g, ub %g\n", 
                       index, startColLB[index], startColUB[index]);
            }
            
            if (startColLB[index] > startColUB[index]) {
		//assert(0);
            }
#endif
        }
        
        //--------------------------------------------------
        // TODO: Modify hard/soft row lb/ub.
        //--------------------------------------------------


        //--------------------------------------------------
        // Collect active non-core constraints at parent.
        //--------------------------------------------------

	//----------------------------------------------
	// First collect all generated cuts, then remove 
	// deleted.
	//----------------------------------------------
	
	tempInt = pathDesc->getCons()->numAdd;
	    
#ifdef BLIS_DEBUG_MORE
	std::cout << "\nINSTALL: numAdd = " << tempInt << std::endl;
#endif

	int maxOld = model->getOldConstraintsSize();
	
	for (k = 0; k < tempInt; ++k) {
	    aCon = dynamic_cast<BlisConstraint *>
		(pathDesc->getCons()->objects[k]);
	    
	    assert(aCon);
	    assert(aCon->getSize() > 0);
	    assert(aCon->getSize() < 100000);
	    
#ifdef BLIS_DEBUG_MORE
	    std::cout << "INSTALL: cut  k=" << k 
		      << ", len=" <<aCon->getSize() 
		      << ", node="<< index_ << std::endl;
#endif
	    (model->oldConstraints())[numOldCons++] = aCon;
	    
	    if (numOldCons >= maxOld) {
		// Need resize
#ifdef BLIS_DEBUG_MORE
		std::cout << "INSTALL: resize, maxOld = " 
			  << maxOld << std::endl;
#endif
		maxOld *= 2;
		BlisConstraint **tempCons = new BlisConstraint* [maxOld];
		
		memcpy(tempCons, 
		       model->oldConstraints(), 
		       numOldCons * sizeof(BlisConstraint *));
		
		model->delOldConstraints();
		model->setOldConstraints(tempCons);
		model->setOldConstraintsSize(maxOld);
	    }
	}
	    
	//----------------------------------------------
	// Remove those deleted. 
	// NOTE: model->oldConstraints_ stores all previously 
	// generated active constraints at parent.
	//----------------------------------------------
	
	tempInt = pathDesc->getCons()->numRemove;
            
	if (tempInt > 0) {
	    int tempPos;
	    int *tempMark = new int [numOldCons];
	    CoinZeroN(tempMark, numOldCons);
	    for (k = 0; k < tempInt; ++k) {
		tempPos = pathDesc->getCons()->posRemove[k];
#ifdef BLIS_DEBUG_MORE
		std::cout << "tempPos=" << tempPos 
			  << ", tempInt=" << tempInt 
			  << ", numOldCons=" << numOldCons << std::endl;
#endif
		tempMark[tempPos] = 1;
		
	    }
	    
	    tempInt = 0;
	    for (k = 0; k < numOldCons; ++k) {
		if (tempMark[k] != 1) {
		    // Survived.
		    (model->oldConstraints())[tempInt++]=
			(model->oldConstraints())[k];
		}
	    }
	    if (tempInt + pathDesc->getCons()->numRemove != numOldCons) {
		std::cout << "INSTALL: tempInt=" << tempInt
			  <<", numRemove="<<pathDesc->getCons()->numRemove
			  << ", numOldCons=" << numOldCons << std::endl;
		
		assert(0);
	    }
	    
	    // Update number of old non-core constraints.
	    numOldCons = tempInt;
	    delete [] tempMark;
	}
    } // EOF leafToRootPath.


    //--------------------------------------------------------
    // Debug variable bounds to be installed.
    //--------------------------------------------------------

#ifdef BLIS_DEBUG_MORE
    for (k = 0; k < numCols; ++k) {
        //if (index_ == -1) {
            printf("INSTALL: Col %d, \tlb %g,  \tub %g\n",
                   k, startColLB[k], startColUB[k]);
	    //}
        
        if (startColLB[k] > startColUB[k] + ALPS_GEN_TOL) {
            printf("INSTALL: Col %d, \tlb %g,  \tub %g\n",
                   k, startColLB[k], startColUB[k]);
            assert(0);
        }
    }
#endif   

    //--------------------------------------------------------
    // Clear path vector.
    //--------------------------------------------------------
    
    model->leafToRootPath.clear();
    assert(model->leafToRootPath.size() == 0);
    
    //--------------------------------------------------------
    // Adjust column bounds in lp solver  
    //--------------------------------------------------------
    
    for(i = 0; i < numCols; ++i) {
	model->solver()->setColBounds(i, startColLB[i], startColUB[i]); 
    }
    
    //--------------------------------------------------------
    // TODO: Set row bounds 
    //--------------------------------------------------------
    

    //--------------------------------------------------------
    // Add old constraints, which are collect from differencing.
    //--------------------------------------------------------
    
    // If removed cuts due to local cuts.

    model->setNumOldConstraints(numOldCons);
    
#ifdef BLIS_DEBUG
    std::cout << "INSTALL: after collecting, numOldCons = " << numOldCons 
	      << std::endl;
#endif

    if (numOldCons > 0) {
	const OsiRowCut ** oldOsiCuts = new const OsiRowCut * [numOldCons];
	for (k = 0; k < numOldCons; ++k) {
	    OsiRowCut * acut = 
		BlisConstraintToOsiCut(model->oldConstraints()[k]);
	    oldOsiCuts[k] = acut;
	}
	model->solver()->applyRowCuts(numOldCons, oldOsiCuts);
	for (k = 0; k < numOldCons; ++k) {
	    delete oldOsiCuts[k];
	}
	delete [] oldOsiCuts;
	oldOsiCuts = NULL;
    }
    
    //--------------------------------------------------------
    // Add parent variables, which are collect from differencing.
    //--------------------------------------------------------
    

    //--------------------------------------------------------
    // Set basis
    //--------------------------------------------------------

    CoinWarmStartBasis *pws = desc->getBasis();

    if (pws != NULL) {
	model->solver()->setWarmStart(pws);

#ifdef BLIS_DEBUG
	printf("NODE %d: set warm start\n", getIndex());
	
	numCols = model->solver()->getNumCols();
	numRows = model->solver()->getNumRows();
	int nStr = pws->getNumStructural();
	int nArt = pws->getNumArtificial();
	
	if (numCols != nStr) {
	    std::cout << "nStr=" << nStr << ", numCols=" << numCols 
		      << std::endl;
	    assert(0);
	}
	std::cout << "nArt=" << nArt << ", numRows=" << numRows 
		  << std::endl;
	if (numRows != nArt) {
	    std::cout << "nArt=" << nArt << ", numRows=" << numRows 
		      << std::endl;
	    assert(0);
	}
#endif

    }  

    return status;
}

//#############################################################################

int 
BlisTreeNode::generateConstraints(BlisModel *model, OsiCuts & cutPool) 
{
    int i, numCGs;
    int status = BLIS_LP_OPTIMAL;
    int preNumRowCons = 0;
    int preNumColCons = 0;
    int newCons = 0;
    int strategy = -2;
    int maxStrategy = -2;
    
    bool mustResolve = false;
    bool fullScan = true;
    
    double useTime;
    
    numCGs = model->numCutGenerators();
  
    
    for (i = 0 ; i < numCGs; ++i) {
	
	//----------------------------------------------------
	// Check if call this generator.
	//----------------------------------------------------
	
	strategy =  model->cutGenerators(i)->strategy();
	maxStrategy = ALPS_MAX(strategy, maxStrategy);
      
	bool useThis = false;
	if (strategy == -2) {
	    useThis = false;
	}
	else if (strategy == -1) {
	    if (model->isRoot_) useThis = true;
	}
	else if (strategy == 0) {
	    if (!diving_ || model->isRoot_) useThis = true;
	}
	else if (strategy > 0) {
	    // Num of nodes is set at the beginning of process().
	    int numNodes = model->getNumNodes();
	    if ((numNodes-1) % strategy == 0) {
		useThis = true;
	    }
	}
	
#ifdef BLIS_DEBUG_MORE
	std::cout<<"CUTGEN: " << model->cutGenerators(i)->name() 
		 <<": useThis="<<useThis<<", diving=" <<diving_
		 << ", strategy=" << strategy 
		 << ", num of nodes=" << model->getNumNodes()
		 <<std::endl;
#endif

	//----------------------------------------------------
	// Generator constraints.
	//----------------------------------------------------

	if (useThis) {
	    newCons = 0;
	    preNumRowCons = cutPool.sizeRowCuts();
	    preNumColCons = cutPool.sizeColCuts();
          
	    useTime = CoinCpuTime();
	    mustResolve = 
		model->cutGenerators(i)->generateCons(cutPool, fullScan);
	    useTime = CoinCpuTime() - useTime;

	    if (mustResolve) {
		// TODO: Only probing will return ture.
		status = bound(model);
		if (status == BLIS_LP_OPTIMAL) {
#ifdef BLIS_DEBUG
		    std::cout << "CUTGEN: after probing, this node survived."
			      << std::endl;
#endif             
		}
		else {
#ifdef BLIS_DEBUG
		    std::cout<<"CUTGEN: after probing, this node can fathomed."
			     << std::endl;
#endif
		    break;
		}
	    }
	    
	    //------------------------------------------------
	    // Modify control. 
	    // NOTE: only modify if user choose automatic.
	    //------------------------------------------------
	    
	    if ( (model->useCons() == 0) && 
		 (model->cutGenerators(i)->noConsCalls() > 30) ) {
		// disable.
		model->cutGenerators(i)->setStrategy(-2);
		maxStrategy = ALPS_MAX(strategy, maxStrategy);
	    }
	}
    }
    
    if ( (model->useCons() == 0) && (maxStrategy == -2) ){
	// Previously automatic, now all cut generators have been disabled.
	model->setUseCons(-2);
    }
    
    return status;
}

//#############################################################################

int 
BlisTreeNode::applyConstraints(BlisModel *model, 
                               OsiCuts & osiCutSet,
                               const double *solution)
{
    int status = BLIS_OK;
    int i, k;
    
    //int numColumnCuts = osiCutSet.sizeColCuts() ;
    int numRowCuts = osiCutSet.sizeRowCuts();
    int numToAdd = numRowCuts;
    int numAdded = 0;
    
    if (numRowCuts > 0) {
	BlisParams * BlisPar = model->BlisPar();
	double scaleConFactor = BlisPar->entry(BlisParams::scaleConFactor);
        
        if (numToAdd > 0) { 
	    
            OsiRowCut *rowCut = NULL;
	    
#ifdef BLIS_DEBUG
            printf("\nAPPLYCUT: before selecting, num of new cuts = %d\n",
                   numRowCuts);
#endif
            int numRowsNow = model->solver()->getNumRows();
            int numCols = model->solver()->getNumCols();
            CoinWarmStartBasis *ws = dynamic_cast<CoinWarmStartBasis*>
                (model->solver()->getWarmStart());
            
            const OsiRowCut ** addCuts = new const OsiRowCut * [numToAdd];

            for (i = 0 ; i < numToAdd ; i++) {
		bool keep = true;
		
                rowCut = &(osiCutSet.rowCut(i));
                
		//------------------------------------------
		// Remove:
		//  - empty cuts
		//  - dense cuts
		//  - bad scaled cuts
                //  - weak cuts
		//  - parallel cuts
		//------------------------------------------

		const CoinPackedVector & rowVector = rowCut->row();
		int length = rowVector.getNumElements();
                bool check = true;
                
                while (check) {
                    //--------------------------------------                   
                    // Empty.
                    //--------------------------------------

                    if (length <= 0) {
                        keep = false;

#ifdef BLIS_DEBUG_MORE
                        std::cout << "APPLYCUT: A empty cut." << std::endl;
#endif
                        break;
                    }

                    //--------------------------------------
                    // Dense.
                    //--------------------------------------

                    if(length > model->getDenseConCutoff()){
                        keep = false;
#ifdef BLIS_DEBUG
                        std::cout << "APPLYCUT: A dense cut. length = " 
                                  << length << ", cutoff = " 
                                  << model->getDenseConCutoff() << std::endl;
#endif  
                        break;
                    }

                    //--------------------------------------
                    // Compuate scale factor.
                    //--------------------------------------

                    const double *elements = rowVector.getElements();
                    const int *indices = rowVector.getIndices();
                    
                    int index;
                    double activity = 0.0;
                    
                    double maxElem = 0.0;
                    double minElem = ALPS_DBL_MAX;
                    double scaleFactor;
                    
                    for (k = 0; k < length; ++k) {
                        if (fabs(elements[k]) > maxElem) {
                            maxElem = fabs(elements[k]);
                        }
                        if (fabs(elements[k]) < minElem) {
                            minElem = fabs(elements[k]);
                        }
                        index = indices[k];
                        activity += elements[k] * solution[k];
                    }
                    if(minElem != 0.0) {
                        scaleFactor = maxElem/minElem;
                    }
                    else {
                        assert(0);
                        scaleFactor = ALPS_DBL_MAX;
                    }
                    
#ifdef BLIS_DEBUG
                    std::cout << "APPLYCUT: scaleFactor=" << scaleFactor
                              << ", maxElem=" << maxElem 
                              << ", minElem=" << minElem << std::endl;
#endif
                    if (scaleFactor > scaleConFactor) {
#ifdef BLIS_DEBUG
                        std::cout<< "APPLYCUT: remove a bad scaled cut"
                                 << std::endl;
#endif
                        keep = false;
                        break;
                    }
                    
                    //--------------------------------------
                    // Weak.
                    //--------------------------------------

                    char cutSense = rowCut->sense();
                    double cutRhs = rowCut->rhs();
                    double violation = -1.0;
                    double range;
                    
                    switch(cutSense) {
                    case 'E':
                        violation = fabs(activity - cutRhs);
                        break;
                    case 'G':
                        violation = cutRhs - activity;
                        break;
                    case 'L':
                        violation = activity - cutRhs;
                    case 'R':
                        range = rowCut->range();
                        violation = ALPS_MAX(violation, activity - cutRhs);
                        violation = ALPS_MAX(violation, cutRhs-range-activity);
                        break;
                    case 'N':
                    default:
                        throw CoinError("Unknown cut sense", 
                                        "applyConstraint", "BlisTreeNode");
                        break;
                    }
                    
                    if (violation < 1.0e-6) {
                        // Found a weak cuts.
#ifdef BLIS_DEBUG
                        std::cout<< "APPLYCUT: remove a weak cut, violation="
                                 << violation << std::endl;
#endif
                        keep = false;
                        break;
                    }
                    
                    //--------------------------------------
                    // Parallel cuts.
                    //--------------------------------------
                    
		    bool paral = parallel(model, 
					  &osiCutSet,
					  i,
					  rowCut);
		    if (paral) {
#ifdef BLIS_DEBUG
                        std::cout<< "APPLYCUT: remove a parallel"<< std::endl;
#endif
			keep = false;
			break;
		    }
                    
                    //--------------------------------------
                    // Check once and stop.
                    //--------------------------------------

                    check = false;
                }//while
                
		if (keep) {
                    addCuts[numAdded++] = &(osiCutSet.rowCut(i));
                }
                else {
                    osiCutSet.eraseRowCut(i);
                    --i;
                    --numToAdd;
                }


            }
	    

#ifdef BLIS_DEBUG
            printf("APPLYCUT: after selecting, num of new cuts = %d\n\n",
                   numAdded);
#endif
            
            //----------------------------------------------
            // Add cuts to lp and adjust basis.
            //----------------------------------------------

            model->solver()->applyRowCuts(numAdded, addCuts);
            delete [] addCuts;
            
            ws->resize(numRowsNow + numToAdd, numCols);
            for (i = 0 ; i < numToAdd; ++i) { 
                ws->setArtifStatus(numRowsNow + i,
                                   CoinWarmStartBasis::basic); 
            }
            if (model->solver()->setWarmStart(ws) == false) { 
                throw CoinError("Fail setWarmStart() after cut installation.",
                                "applyConstraints","BlisTreeNode"); 
            }
            delete ws;
        }   
    }
    
    return status;
}

//#############################################################################

int BlisTreeNode::
reducedCostFix(BlisModel *model)
{ 
    int i, var;
    int status = BLIS_OK;

    int numFixedUp = 0;
    int numFixedDown = 0;
    int numTighten = 0;

    double movement;    
    double newBound;
    double boundDistance;
    double dj;

    const double *lb = model->solver()->getColLower();
    const double *ub = model->solver()->getColUpper();
    const double *solution = model->solver()->getColSolution();
    const double *reducedCost = model->solver()->getReducedCost();
    
    double cutup = getKnowledgeBroker()->getIncumbentValue() *
        model->solver()->getObjSense();
    
    if (cutup >= ALPS_OBJ_MAX) return status;
    
    double lpObjValue = model->solver()->getObjValue() * 
        model->solver()->getObjSense();
    double epInt = 1.0e-5;
    
    int numIntegers = model->getNumIntVars();
    const int *intIndices = model->getIntVars();
    
    for (i = 0; i < numIntegers; ++i) { 
	var = intIndices[i];
        
        dj = reducedCost[var];
        
        if (fabs(dj) < epInt) continue;
        
        boundDistance = ub[var] - lb[var];
	if (boundDistance < epInt) continue;
	
        movement = floor((cutup - lpObjValue) / fabs(dj));
        
        if (solution[var] > ub[var] - epInt) {
            /* At upper bound */
            if (movement < boundDistance) {
                /* new lower bound. If movement is 0, then fix. */
                newBound = ub[var] - movement;
                newBound = std::min(newBound, ub[var]);

#ifdef BLIS_DEBUG_MORE
                printf("RED-FIX: dj %g, lb %.10g, ub %.10g, newBound %.10g, movement %g\n", dj, lb[var], ub[var], newBound, movement);
#endif
                
                if (movement <= ALPS_ZERO) {
                    ++numFixedUp;
                }
                else if (newBound < ub[var]){
                    ++numTighten;
                }
                model->solver()->setColLower(var, newBound);
            }
        }
        else if (solution[var] < lb[var] + epInt) {
            /* At lower bound */
            if (movement < boundDistance) {
                newBound = lb[var] + movement;
                newBound = std::max(newBound, lb[var]);
                
#ifdef BLIS_DEBUG_MORE
                printf("RED-FIX: dj %g, lb %g, ub %g, newBound %g, movement %g\n", dj, lb[var], ub[var], newBound, movement);
#endif
		
                if (movement <= ALPS_ZERO) {
                    ++numFixedDown;
                }
                else if(newBound > lb[var] ){
                    ++numTighten;
                }
                /* new upper bound. If movement is 0, then fix. */
                model->solver()->setColUpper(var, newBound);
            }
        }
    }
    
    //int change = numFixedUp + numFixedDown + numTighten;
    //model->reducedCostFixed_ += change;
#ifdef BLIS_DEBUG_MORE
    if (numFixedUp > 0 || numFixedDown > 0 || numTighten > 0) {
        printf("reducedCostFix: numFixedUp = %d, numFixedDown = %d, numTighten %d\n", numFixedUp, numFixedDown, numTighten);
    }
#endif

    return status; 
}

//#############################################################################

AlpsEncoded*
BlisTreeNode::encode() const 
{
#ifdef BLIS_DEBUG
    std::cout << "BlisTreeNode::encode()--start to encode node "
	      << index_ << std::endl;
#endif

    AlpsReturnCode status = ALPS_OK;
    
    // NOTE: "ALPS_NODE" is used as type name.
    AlpsEncoded* encoded = new AlpsEncoded("ALPS_NODE");
    
    // Encode decription.
    status = desc_->encode(encoded);
    
    // Encode Alps portion.
    status = encodeAlps(encoded);
    
    // Encode Bcps portion.
    status = encodeBcps(encoded);
    
    // Nothing to encode for Blis portion.
    
    return encoded;
}

//#############################################################################

AlpsKnowledge* 
BlisTreeNode::decode(AlpsEncoded& encoded) const 
{
    AlpsReturnCode status = ALPS_OK;
    BlisTreeNode* treeNode = NULL;

    BlisModel *model = dynamic_cast<BlisModel*>(desc_->getModel());
    
    //------------------------------------------------------
    // Unpack decription.
    //------------------------------------------------------

    AlpsNodeDesc* nodeDesc = new BlisNodeDesc(model);
    status = nodeDesc->decode(encoded);
    
    //------------------------------------------------------
    // Unpack node.
    //------------------------------------------------------
    
    // Unpack Alps portion.
    treeNode = new BlisTreeNode(nodeDesc);
    treeNode->decodeAlps(encoded);
    
    // Unpack Bcps portion.
    int type = 0;
    encoded.readRep(type);	
    if (type == BLIS_BO_INT) {
	// branchObject_ is simple integer.
	BlisBranchObjectInt *bo = new BlisBranchObjectInt();
	status = bo->decode(encoded);

	// Set bo in treeNode.
	treeNode->setBranchObject(bo);
    }
    
    // Nothing to unpack for Blis portion.
    
    return treeNode;
}

//#############################################################################

void 
BlisTreeNode::convertToExplicit() 
{
#ifdef BLIS_DEBUG
    std::cout << "BLIS: convertToExplicit(); explicit_="<<explicit_ << std::endl;
#endif

    if(!explicit_) {
	
	// Convert to explicit
	explicit_ = 1;
	
	BlisModel* model = dynamic_cast<BlisModel*>(desc_->getModel());
	BlisNodeDesc *desc = dynamic_cast<BlisNodeDesc *>(desc_);
	BlisConstraint *aCon = NULL;
	
	int numCols = model->solver()->getNumCols();

	int i, k, index;
	int tempInt;

	int numModify = 0;
	int numSoftVarLowers = 0;
	int numSoftVarUppers = 0;

	double value;
		
	double *fVarHardLB = new double [numCols];
	double *fVarHardUB = new double [numCols];
	int *fVarHardLBInd = new int [numCols];
	int *fVarHardUBInd = new int [numCols];
	
	double *fVarSoftLB = new double [numCols];
	double *fVarSoftUB = new double [numCols];
	int *fVarSoftLBInd = new int [numCols];
	int *fVarSoftUBInd = new int [numCols];

	for (k = 0; k < numCols; ++k) {
	    fVarSoftLB[k] = ALPS_DBL_MAX;
	    fVarSoftUB[k] = -ALPS_DBL_MAX;
	    fVarHardLB[k] = ALPS_DBL_MAX;
	    fVarHardUB[k] = -ALPS_DBL_MAX;
	    fVarHardLBInd[k] = k;
	    fVarHardUBInd[k] = k;
	}

	int numOldCons = 0;
	int maxOld = model->getOldConstraintsSize();
	BcpsObject ** oldConstraints = new BcpsObject* [maxOld];
	
	//--------------------------------------------------
	// Travel back to a full node, then collect diff (add/rem col/row,
	// hard/soft col/row bounds) from the node full to this node.
	//--------------------------------------------------------
    
	BlisNodeDesc* pathDesc = NULL;
	AlpsTreeNode *parent = parent_;    
	
	model->leafToRootPath.push_back(this);
	
	while(parent) {
#ifdef BLIS_DEBUG_MORE
	    std::cout << "Parent id = " << parent->getIndex() << std::endl;
#endif     
	    model->leafToRootPath.push_back(parent);
	    if (parent->getExplicit()) {
		// Reach an explicit node, then stop.
		break;
	    }
	    else {
		parent = parent->getParent();
	    }
	}
    
#ifdef BLIS_DEBUG
	std::cout << "CONVERT TO EXP: path len = " << model->leafToRootPath.size()
		  << std::endl;
#endif
	
	//------------------------------------------------------
	// Travel back from this node to the explicit node to
	// collect full description.
	//------------------------------------------------------	


	for(i = model->leafToRootPath.size() - 1; i > -1; --i) {

	    pathDesc = dynamic_cast<BlisNodeDesc*>((model->leafToRootPath.at(i))->
						   getDesc());

	    //--------------------------------------
	    // Full variable hard bounds.
	    //--------------------------------------
	    
	    numModify = pathDesc->getVars()->lbHard.numModify;
	    for (k = 0; k < numModify; ++k) {
		index = pathDesc->getVars()->lbHard.posModify[k];
		value = pathDesc->getVars()->lbHard.entries[k];
		fVarHardLB[index] = value;
	    }
		    
	    numModify = pathDesc->getVars()->ubHard.numModify;
	    for (k = 0; k < numModify; ++k) {
		index = pathDesc->getVars()->ubHard.posModify[k];
		value = pathDesc->getVars()->ubHard.entries[k];
		fVarHardUB[index] = value;
	    }
		    
	    //--------------------------------------
	    // Full variable soft bounds.
	    //--------------------------------------
	    
	    numModify = pathDesc->getVars()->lbSoft.numModify;
#ifdef BLIS_DEBUG
	    std::cout << "CONVERT: EXP: i=" << i << ", numModify soft lb="
		      << numModify << std::endl;
#endif
	    for (k = 0; k < numModify; ++k) {
		index = pathDesc->getVars()->lbSoft.posModify[k];
		value = pathDesc->getVars()->lbSoft.entries[k];
		fVarSoftLB[index] = value;
	    }
		    
	    numModify = pathDesc->getVars()->ubSoft.numModify;
#ifdef BLIS_DEBUG
	    std::cout << "CONVERT: EXP: i=" << i << ", numModify soft ub="
		      << numModify << std::endl;
#endif
	    for (k = 0; k < numModify; ++k) {
		index = pathDesc->getVars()->ubSoft.posModify[k];
		value = pathDesc->getVars()->ubSoft.entries[k];
		fVarSoftUB[index] = value;
	    }


            //----------------------------------------------
            // Collect all generated constraints, then remove deleted.
            //----------------------------------------------
            
            tempInt = pathDesc->getCons()->numAdd;
	    
#ifdef BLIS_DEBUG
            std::cout << "\nCONVERT: EXP: numAdd = " << tempInt << std::endl;
#endif
            
            for (k = 0; k < tempInt; ++k) {
                aCon = dynamic_cast<BlisConstraint *>
                    (pathDesc->getCons()->objects[k]);
                
                assert(aCon);
                assert(aCon->getSize() > 0);
                assert(aCon->getSize() < 100000);
                
#ifdef BLIS_DEBUG
		std::cout << "CONVERT: EXP: k=" << k 
			  << ", len=" <<aCon->getSize() << std::endl;
#endif
                oldConstraints[numOldCons++] = aCon;
                
                if (numOldCons >= maxOld) {
                    // Need resize
#ifdef BLIS_DEBUG
                    std::cout << "CONVERT: EXP: resize, maxOld = " 
                              << maxOld << std::endl;
#endif
                    maxOld *= 2;
                    BcpsObject **tempCons = new BcpsObject* [maxOld];
		    
                    memcpy(tempCons, 
                           oldConstraints, 
                           numOldCons * sizeof(BcpsObject *));
		    
		    delete [] oldConstraints;
		    oldConstraints = tempCons;
		    tempCons = NULL;
                }
            }
	    
            //----------------------------------------------
            // Remove those deleted. 
	    // NOTE: oldConstraints stores all previously 
	    // generated active constraints at parent.
            //----------------------------------------------
	    
            tempInt = pathDesc->getCons()->numRemove;
            
            if (tempInt > 0) {
                int tempPos;
                int *tempMark = new int [numOldCons];
                CoinZeroN(tempMark, numOldCons);
                for (k = 0; k < tempInt; ++k) {
                    tempPos = pathDesc->getCons()->posRemove[k];
#ifdef BLIS_DEBUG_MORE
		    std::cout << "tempPos=" << tempPos 
			      << ", tempInt=" << tempInt 
			      << ", numOldCons=" << numOldCons << std::endl;
#endif
                    tempMark[tempPos] = 1;
                }
		
                tempInt = 0;
                for (k = 0; k < numOldCons; ++k) {
                    if (tempMark[k] != 1) {
                        // Survived.
                        oldConstraints[tempInt++] = oldConstraints[k];
                    }
                }
		if (tempInt + pathDesc->getCons()->numRemove != numOldCons) {
		    std::cout << "INSTALL: tempInt=" << tempInt
			      <<", numRemove="<<pathDesc->getCons()->numRemove
			      << ", numOldCons=" << numOldCons << std::endl;
		    
		    assert(0);
		}

                // Update number of old non-core constraints.
                numOldCons = tempInt;
                delete [] tempMark;
            }
	    
	} // EOF for (path)
        
	//------------------------------------------
	// Record hard variable bounds. FULL set.
	//------------------------------------------
	
	desc->assignVarHardBound(numCols,
				 fVarHardLBInd,
				 fVarHardLB,
				 numCols,
				 fVarHardUBInd,
				 fVarHardUB);

	//------------------------------------------
	// Recode soft variable bound. Modified.
	//------------------------------------------
	
	for (k = 0; k < numCols; ++k) {
	    if (fVarSoftLB[k] < ALPS_BND_MAX) {
		fVarSoftLBInd[numSoftVarLowers] = k;
		fVarSoftLB[numSoftVarLowers++] = fVarSoftLB[k];
	    }
	    if (fVarSoftUB[k] > -ALPS_BND_MAX) {
		fVarSoftUBInd[numSoftVarUppers] = k;
		fVarSoftUB[numSoftVarUppers++] = fVarSoftUB[k];
	    }
	}
	// Assign it anyway so to delete memory(fVarSoftLBInd,etc.)
	desc->assignVarSoftBound(numSoftVarLowers, 
				 fVarSoftLBInd,
				 fVarSoftLB,
				 numSoftVarUppers, 
				 fVarSoftUBInd,
				 fVarSoftUB);
	
	//------------------------------------------
	// Recode added constraints.
	//------------------------------------------

	// First make a hard copy.
	for (k = 0; k < numOldCons; ++k) {
	    aCon = dynamic_cast<BlisConstraint *>(oldConstraints[k]);
	    assert(aCon);
	    
	    BlisConstraint *newCon = new BlisConstraint(*aCon);
	    oldConstraints[k] = newCon;
	}
	// Add will first delete, then add. It is safe to use in parallel.
	desc->setAddedConstraints(numOldCons, oldConstraints);

	//------------------------------------------
	// Recode deleted constraints.
	//------------------------------------------
	
	desc->delConstraints(0, NULL);
	
	//--------------------------------------------------
	// Clear path vector.
	//--------------------------------------------------
	
	model->leafToRootPath.clear();
	assert(model->leafToRootPath.size() == 0);


    } // EOF of if.
}

//#############################################################################

// Not defined yet.
void
BlisTreeNode::convertToRelative()
{
    if(explicit_) {
	

    }
}

//#############################################################################

bool 
BlisTreeNode::parallel(BlisModel *model, 
		       OsiCuts *newCutSet,
		       int lastNew,
		       OsiRowCut *rowCut)
{
    bool parallel = false;
    int k;
    double threshold = 0.999;
    
    //------------------------------------------------------
    // Compare with old cuts
    //------------------------------------------------------

    int numOldCons = model->getNumOldConstraints();
    for (k = 0; k < numOldCons; ++k) {
	BlisConstraint *aCon = model->oldConstraints()[k];
	assert(aCon);
	parallel = BlisParallelCutCon(rowCut,
				      aCon,
				      threshold);
	if (parallel) return parallel;
    }
    

    //------------------------------------------------------
    // Compare with old cuts
    //------------------------------------------------------

    for (k = 0; k < lastNew; ++k) {
	OsiRowCut *rowCut2 = &(newCutSet->rowCut(k));
	parallel = BlisParallelCutCut(rowCut,
				      rowCut2,
				      threshold);
	if (parallel) return parallel;
    }
    
    return parallel;
}

//#############################################################################


