/*===========================================================================*
 * This file is part of the Branch, Constrain and Price Software (BiCePS)    *
 *                                                                           *
 * BiCePS is distributed under the Common Public License as part of the      *
 * COIN-OR repository (http://www.coin-or.org).                              *
 *                                                                           *
 * Authors: Yan Xu, SAS Institute Inc.                                       *
 *          Ted Ralphs, Lehigh University                                    *
 *          Laszlo Ladanyi, IBM T.J. Watson Research Center                  *
 *          Matthew Saltzman, Clemson Unniversity                            *
 *                                                                           *
 * Copyright (C) 2001-2006, Lehigh University, Yan Xu, and Ted Ralphs.       *
 * All Rights Reserved.                                                      *
 *===========================================================================*/

#include "BcpsTreeNode.h"
#include "BcpsNodeDesc.h"

//#############################################################################

int
BcpsTreeNode::process(bool isRoot, bool rampUp)
{
    int status = ALPS_OK;
    bool keepOn = true;
    bool fathomed = false;
    bool genCons = true;
    bool genVars = true;
    
    int maxPass = 20; // Should be a parameter.
    int pass = 0;

    BcpsConstraintPool *newConPool = NULL;
    BcpsVariablePool *newVarPool = NULL;

    BcpsModel *model = dynamic_cast<BcpsModel*>(desc_->getModel());
    
    //------------------------------------------------------
    // Extract node information (bounds, constraints, variables) from 
    // this node and load the information into the relaxation solver,
    // such as linear programming solver.
    //------------------------------------------------------
    
    installSubProblem(model);
    
    //------------------------------------------------------
    // Iteratively:
    //  - bounding,
    //  - heuristic searching,
    //  - constraint generating,
    //  - variable generating.
    //------------------------------------------------------

    while (keepOn && (pass < maxPass)) {
        ++pass;
        keepOn = false;
        
        //--------------------------------------------------
        // Bounding to get the quality of this node.
        //--------------------------------------------------
        
        status = bound(model);

        //--------------------------------------------------
        // Handle bounding status:
	//  - relaxed feasible but not integer feasible,
	//  - integer feasible,
	//  - infeasible,
	//  - unbounded,
	//  - fathomed (for instance, reaching objective limit for LP).
	// Set node status accordingly.
        //--------------------------------------------------

	status = handleBoundingStatus(status, keepOn, fathomed);

	if(fathomed || !keepOn) {
	    // Infeasible, fathomed, tailing off or some user's criteriall.
	    break;
	}

        //--------------------------------------------------
        // Generate constraints.
        //--------------------------------------------------

        if (genCons) {
	    status = generateConstraints(model, newConPool);
	}
	
        //--------------------------------------------------
        // Generate variables.
        //--------------------------------------------------
	
        if (genVars) {
	    status = generateVariables(model, newVarPool);
	}
    }
    
    //------------------------------------------------------
    // Select branching object
    //------------------------------------------------------
    
    if (!fathomed) { 
	status = chooseBranchingObject(model);
    }
    

    return status;    
}

//#############################################################################

