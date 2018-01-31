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

// define structs that will be used in creating nodes.
struct SparseVector {
  int * ind;
  double * val;
};

struct Bound {
  SparseVector lower;
  SparseVector upper;
};


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
  AlpsNodeStatus status = getStatus();
  BlisModel * model = dynamic_cast<BlisModel*>(broker()->getModel());
  double cutoff = model->BlisPar()->entry(BlisParams::cutoff);
  double sense = model->BlisPar()->entry(BlisParams::objSense);
  cutoff = sense*cutoff;
  cutoff = CoinMin(cutoff, broker()->getIncumbentValue());
  // check if this can be fathomed
  double rel_gap_limit = model->BlisPar()->entry(BlisParams::optimalRelGap);
  double abs_gap_limit = model->BlisPar()->entry(BlisParams::optimalAbsGap);
  double abs_gap = cutoff-getQuality();
  double rel_gap = abs_gap/fabs(cutoff);
  if (rel_gap_limit>rel_gap or abs_gap_limit>abs_gap) {
    setStatus(AlpsNodeStatusFathomed);
    return AlpsReturnStatusOk;
  }
  if (status==AlpsNodeStatusCandidate or
      status==AlpsNodeStatusEvaluated) {
    boundingLoop();
  }
  else if (status==AlpsNodeStatusBranched or
           status==AlpsNodeStatusFathomed or
           status==AlpsNodeStatusDiscarded) {
    // this should not happen
    std::cerr << "Unexpected node status!" << std::endl;
    throw std::exception();
  }
  return AlpsReturnStatusOk;
}


void BlisTreeNode::boundingLoop() {
  BlisModel * model = dynamic_cast<BlisModel*>(broker_->getModel());
  bool keepBounding = true;
  bool do_branch = false;
  bool genConstraints = 0;
  bool genVariables = false;
  BcpsConstraintPool * constraintPool = new BcpsConstraintPool();
  BcpsVariablePool * variablePool = new BcpsVariablePool();
  installSubProblem();

  while (keepBounding) {
    keepBounding = false;
    // solve subproblem corresponds to this node
    BcpsSubproblemStatus subproblem_status = bound();
    // update number of iterations statistics
    //model->addNumRelaxIterations();
    if ((subproblem_status==BcpsSubproblemStatusOptimal) &&
        (getStatus()==AlpsNodeStatusCandidate or
         getStatus()==AlpsNodeStatusEvaluated)) {
    }
    // end of grumpy message

    // check if this can be fathomed
    double cutoff = model->BlisPar()->entry(BlisParams::cutoff);
    double sense = model->BlisPar()->entry(BlisParams::objSense);
    cutoff = sense*cutoff;
    double rel_gap_limit = model->BlisPar()->entry(BlisParams::optimalRelGap);
    double abs_gap_limit = model->BlisPar()->entry(BlisParams::optimalAbsGap);
    cutoff = CoinMin(cutoff, broker()->getIncumbentValue());
    double abs_gap = cutoff-getQuality();
    double rel_gap = abs_gap/fabs(cutoff);
    //std::cout << "abs " << abs_gap << " limit " << abs_gap_limit << std::endl;
    //std::cout << "rel " << rel_gap << " limit " << rel_gap_limit << std::endl;
    if (rel_gap_limit>rel_gap or abs_gap_limit>abs_gap) {
      setStatus(AlpsNodeStatusFathomed);
      break;
    }
    // call heuristics to search for a solution
    callHeuristics();

    // decide what to do
    branchConstrainOrPrice(subproblem_status, keepBounding, do_branch,
                           genConstraints,
                           genVariables);
    if (getStatus()==AlpsNodeStatusFathomed) {
      // node is fathomed, nothing to do.
      break;
    }
    else if (keepBounding and genConstraints) {
      generateConstraints(constraintPool);
      // add constraints to the model
      applyConstraints(constraintPool);
      // clear constraint pool
      constraintPool->freeGuts();
      // set status to evaluated
      setStatus(AlpsNodeStatusEvaluated);
    }
    else if (keepBounding and genVariables) {
      generateVariables(variablePool);
      // add variables to the model
      // set status to evaluated
      setStatus(AlpsNodeStatusEvaluated);
    }
    else if (keepBounding==false and do_branch==false) {
      // put node back into the list.
      // this means do not change the node status and end processing the node.
      // set status to evaluated
      setStatus(AlpsNodeStatusEvaluated);
    }
    else if (keepBounding==false and do_branch) {
      // branch
      BcpsBranchStrategy * branchStrategy = model->branchStrategy();
      branchStrategy->createCandBranchObjects(this);
      // prepare this node for branching, bookkeeping for differencing.
      // call pregnant setting routine
      processSetPregnant();
    }
    else {
      std::cerr << "This should not happen!" << std::endl;
      throw std::exception();
    }
  }
  delete constraintPool;
  delete variablePool;
}


void BlisTreeNode::processSetPregnant() {
  // get warm start basis from solver
  BlisModel * model = dynamic_cast<BlisModel*>(broker()->getModel());
  CoinWarmStartBasis * ws = NULL;
  if (model->solver()->getWarmStart()!=NULL) {
    ws = dynamic_cast<CoinWarmStartBasis*>(model->solver()->getWarmStart());
  }
  // store basis in the node desciption.
  getDesc()->setBasis(ws);
  // set status pregnant
  setStatus(AlpsNodeStatusPregnant);
}


//#############################################################################

std::vector< CoinTriple<AlpsNodeDesc*, AlpsNodeStatus, double> >
BlisTreeNode::branch()
{
  // get model the node belongs
  BlisModel * model = dynamic_cast<BlisModel*>(broker()->getModel());

  // check node status, this should be a pregnant node.
  if (getStatus()!=AlpsNodeStatusPregnant) {
    std::cerr << "This should be pregnant node!" << std::endl;
    throw std::exception();
  }

  // create return value and push the down and up nodes.
  std::vector< CoinTriple<AlpsNodeDesc*, AlpsNodeStatus, double> > res;

  // check if this can be fathomed
  double rel_gap_limit = model->BlisPar()->entry(BlisParams::optimalRelGap);
  double abs_gap_limit = model->BlisPar()->entry(BlisParams::optimalAbsGap);
  double abs_gap = broker()->getIncumbentValue()-getQuality();
  double rel_gap = abs_gap/fabs(broker()->getIncumbentValue());
  //std::cout << "abs " << abs_gap << " limit " << abs_gap_limit << std::endl;
  //std::cout << "rel " << rel_gap << " limit " << rel_gap_limit << std::endl;
  if (rel_gap_limit>rel_gap or abs_gap_limit>abs_gap) {
    setStatus(AlpsNodeStatusFathomed);
    return res;
  }

  // get Alps phase
  AlpsPhase phase = broker()->getPhase();

  // get branch object
  BlisBranchObjectInt const * branch_object =
    dynamic_cast<BlisBranchObjectInt const *>(branchObject());

  assert(branch_object);
  // get index and value of branch variable.
  //int branch_var = model->relaxedCols()[branch_object->getObjectIndex()];
  int branch_var = branch_object->index();
  double branch_value = branch_object->value();

  // compute child nodes' warm start basis
  CoinWarmStartBasis * child_ws;
  child_ws = (getDesc()->getBasis()==0) ? 0 :
    new CoinWarmStartBasis(*getDesc()->getBasis());

  // create new node descriptions
  BlisNodeDesc * down_node = new BlisNodeDesc(model);
  down_node->setBroker(broker_);
  BlisNodeDesc * up_node = new BlisNodeDesc(model);
  up_node->setBroker(broker_);
  if (phase == AlpsPhaseRampup) {
    // Down Node
    // == copy description of this node to the down node
    copyFullNode(down_node);
    // == update the branching variable hard bounds for the down node
    down_node->vars()->ubHard.entries[branch_var] =
      branch_object->ubDownBranch();

    // Up Node
    // == copy description of this node to up node
    copyFullNode(up_node);
    // == update the branching variable hard bounds for the up node
    up_node->vars()->lbHard.entries[branch_var] =
      branch_object->lbUpBranch();
  }
  else {
    // Store node description relative to the parent.
    // We need to add a hard bound for the branching variable.

    double ub_down_branch = branch_object->ubDownBranch();
    double lb_up_branch = branch_object->lbUpBranch();
    double lb = model->getVariables()[branch_var]->getLbHard();
    double ub = model->getVariables()[branch_var]->getUbHard();
    down_node->setVarHardBound(1,
                             &branch_var,
                             &lb,
                             1,
                             &branch_var,
                             &ub_down_branch);
    up_node->setVarHardBound(1,
                             &branch_var,
                             &lb_up_branch,
                             1,
                             &branch_var,
                             &ub);
  }

  // Down Node
  // == set other relevant fields of down node
  down_node->setBranchedDir(-1);
  down_node->setBranchedInd(branch_object->index());
  down_node->setBranchedVal(branch_value);
  // == set warm start basis for the down node.
  down_node->setBasis(child_ws);
  // Up Node
  // == set other relevant fields of up node
  up_node->setBranchedDir(1);
  up_node->setBranchedInd(branch_object->index());
  up_node->setBranchedVal(branch_value);
  // == set warm start basis for the up node.
  up_node->setBasis(child_ws);

  // push the down and up nodes.
  res.push_back(CoinMakeTriple(static_cast<AlpsNodeDesc*>(down_node),
                               AlpsNodeStatusCandidate,
                               quality_));
  res.push_back(CoinMakeTriple(static_cast<AlpsNodeDesc*>(up_node),
                               AlpsNodeStatusCandidate,
                               quality_));
  return res;
}


// Copies node description of this node to given child node.
// New node is explicitly stored in the memory (no differencing).
void BlisTreeNode::copyFullNode(BlisNodeDesc * child_node) const {
  // get model the node belongs
  BlisModel * model = dynamic_cast<BlisModel*>(broker_->getModel());
  // get description of this node
  BlisNodeDesc * node_desc = getDesc();

  // solver number of columns
  int num_cols = model->solver()->getNumCols();

  // this will keep bound information, lower and upper
  // it is used for both hard and soft bounds.
  Bound bound;
  bound.lower.ind = new int[num_cols];
  bound.lower.val = new double[num_cols];
  bound.upper.ind = new int[num_cols];
  bound.upper.val = new double[num_cols];

  // Hard bounds
  // == Hard lower bounds
  std::copy(node_desc->getVars()->lbHard.posModify,
            node_desc->getVars()->lbHard.posModify + num_cols,
            bound.lower.ind);
  std::copy(node_desc->getVars()->lbHard.entries,
            node_desc->getVars()->lbHard.entries + num_cols,
            bound.lower.val);
  // == Hard upper bounds
  std::copy(node_desc->getVars()->ubHard.posModify,
            node_desc->getVars()->ubHard.posModify + num_cols,
            bound.upper.ind);
  std::copy(node_desc->getVars()->ubHard.entries,
            node_desc->getVars()->ubHard.entries + num_cols,
            bound.upper.val);

  // assign hard bounds for the child node.
  // this takes ownership of the arrays.
  child_node->assignVarHardBound(num_cols,
                                 bound.lower.ind,
                                 bound.lower.val,
                                 num_cols,
                                 bound.upper.ind,
                                 bound.upper.val);

  // Soft bounds.
  // == Soft lower bounds
  // number of entries modified for soft lower bounds
  int sl_num_modify = node_desc->getVars()->lbSoft.numModify;
  // We can reuse bound here since previos arrays are owned by Bcps.
  bound.lower.ind = new int[sl_num_modify];
  bound.lower.val = new double[sl_num_modify];
  std::copy(node_desc->getVars()->lbSoft.posModify,
            node_desc->getVars()->lbSoft.posModify + sl_num_modify,
            bound.lower.ind);
  std::copy(node_desc->getVars()->lbSoft.entries,
            node_desc->getVars()->lbSoft.entries + sl_num_modify,
            bound.lower.val);

  // == Soft upper bounds
  // number of entries modified for soft upper bounds
  int su_num_modify = node_desc->getVars()->ubSoft.numModify;
  bound.upper.ind = new int[su_num_modify];
  bound.upper.val = new double[su_num_modify];
  std::copy(node_desc->getVars()->ubSoft.posModify,
            node_desc->getVars()->ubSoft.posModify + su_num_modify,
            bound.upper.ind);
  std::copy(node_desc->getVars()->ubSoft.entries,
            node_desc->getVars()->ubSoft.entries + su_num_modify,
            bound.upper.val);

  // Assign soft bounds for the child node.
  child_node->assignVarSoftBound(sl_num_modify,
                                 bound.lower.ind,
                                 bound.lower.val,
                                 su_num_modify,
                                 bound.upper.ind,
                                 bound.upper.val);
}


BcpsSubproblemStatus BlisTreeNode::bound()
{
    BcpsSubproblemStatus status;
    BlisModel * model = dynamic_cast<BlisModel*>(broker()->getModel());

#ifdef BLIS_DEBUG_MORE
    int j;
    int numCols = m->solver()->getNumCols();
    const double * clb = m->solver()->getColLower();
    const double * cub = m->solver()->getColUpper();

    for (j = 0; j < numCols; ++j) {
	std::cout << "c"<< j <<"["<<clb[j]<<", "<< cub[j] << "]" << std::endl;
    }
#endif

    model->solver()->resolve();

    if (model->solver()->isAbandoned()) {
#ifdef BLIS_DEBUG
	std::cout << "BOUND: is abandoned" << std::endl;
#endif
	status = BcpsSubproblemStatusAbandoned;
    }
    else if (model->solver()->isProvenOptimal()) {
#ifdef BLIS_DEBUG
	std::cout << "BOUND: is lp optimal" << std::endl;
#endif
        status = BcpsSubproblemStatusOptimal;
        BlisNodeDesc *desc = dynamic_cast<BlisNodeDesc*>(desc_);

        double objValue = model->solver()->getObjValue() *
            model->solver()->getObjSense();

        int dir = desc->getBranchedDir();
        if (dir != 0) {
            double objDeg = objValue - quality_;
            int objInd = desc->getBranchedInd();
            double lpX = desc->getBranchedVal();
            BlisObjectInt *intObject =
                dynamic_cast<BlisObjectInt *>(model->objects(objInd));
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
    else if (model->solver()->isProvenPrimalInfeasible()) {
#ifdef BLIS_DEBUG
	std::cout << "BOUND: is primal inf" << std::endl;
#endif
        status = BcpsSubproblemStatusPrimalInfeasible;
    }
    else if (model->solver()->isProvenDualInfeasible()) {
#ifdef BLIS_DEBUG
	std::cout << "BOUND: is dual inf" << std::endl;
#endif
        status = BcpsSubproblemStatusDualInfeasible;
    }
    else if (model->solver()->isPrimalObjectiveLimitReached()) {
#ifdef BLIS_DEBUG
	std::cout << "BOUND: is primal limit" << std::endl;
#endif
        status = BcpsSubproblemStatusPrimalObjLim;
    }
    else if (model->solver()->isDualObjectiveLimitReached()) {
#ifdef BLIS_DEBUG
	std::cout << "BOUND: is dual limit" << std::endl;
#endif
        status = BcpsSubproblemStatusDualObjLim;
    }
    else if (model->solver()->isIterationLimitReached()) {
#ifdef BLIS_DEBUG
	std::cout << "BOUND: is iter limit" << std::endl;
#endif
        status = BcpsSubproblemStatusIterLim;
    }
    else {
        status = BcpsSubproblemStatusUnknown;
	std::cerr << "UNKNOWN LP STATUS" << std::endl;
        throw std::exception();
    }

    return status;
}

//#############################################################################

int BlisTreeNode::installSubProblem()
{
    AlpsReturnStatus status = AlpsReturnStatusOk;

    int i, k;
    int index;
    double value;

    BlisModel *model = dynamic_cast<BlisModel *>(broker()->getModel());
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

    AlpsPhase phase = broker_->getPhase();

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

    if (phase != AlpsPhaseRampup) {
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

    for(i = static_cast<int> (model->leafToRootPath.size() - 1); i > -1; --i) {

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
BlisTreeNode::generateConstraints(BcpsConstraintPool *conPool)
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

    BlisModel * model = dynamic_cast<BlisModel*>(broker()->getModel());

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
	    useThis = true;
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
		 <<": useThis="<<useThis
		 << ", strategy=" << strategy
		 << ", num of nodes=" << model->getNumNodes()
		 <<std::endl;
#endif

	//----------------------------------------------------
	// Generator constraints.
	//----------------------------------------------------

	if (useThis) {
	    newCons = 0;
            OsiCuts new_cuts;
	    int preNumCons = 0;

	    useTime = CoinCpuTime();
	    mustResolve =
		model->cutGenerators(i)->generateCons(new_cuts, fullScan);
	    useTime = CoinCpuTime() - useTime;

	    if (mustResolve) {
		// TODO: Only probing will return ture.
		status = bound();
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


int BlisTreeNode::generateVariables(BcpsVariablePool *varPool) {
}

int BlisTreeNode::chooseBranchingObject() {
}

void BlisTreeNode::branchConstrainOrPrice(
                     BcpsSubproblemStatus subproblem_status,
                     bool & keepBounding,
                     bool & branch,
                     bool & generateConstraints,
                     bool & generateVariables) {
  keepBounding = false;
  branch = false;
  generateConstraints = false;
  generateVariables = false;
  if (subproblem_status==BcpsSubproblemStatusPrimalInfeasible) {
    setStatus(AlpsNodeStatusFathomed);
    return;
  }
  if (subproblem_status!=BcpsSubproblemStatusOptimal and
      subproblem_status!=BcpsSubproblemStatusDualInfeasible) {
    std::cerr << "This should not happen!" << std::endl;
    throw std::exception();
  }
  // subproblem is solved to optimality. Check feasibility of the solution.
  int numColsInf;
  double colInf;
  BlisModel * model = dynamic_cast<BlisModel*>(broker_->getModel());
  BlisSolution * sol = model->feasibleSolution(numColsInf, colInf);

  int objSense = model->BlisPar()->entry(BlisParams::objSense);
  if (numColsInf) {
    branch = true;
  }
  else if (sol) {
    sol->setDepth(depth_);
    // Store in Alps pool
    broker()->addKnowledge(AlpsKnowledgeTypeSolution,
                           sol,
                           objSense * sol->getQuality());
    // update solver with the new incumbent if better
    double incum_val = broker()->getIncumbentValue();
    model->solver()->setDblParam(OsiDualObjectiveLimit,
                                 objSense*incum_val);
    // set status to fathom
    setStatus(AlpsNodeStatusFathomed);
  }
  else {
    std::cerr << "This should not happen!" << std::endl;
    throw std::exception();
  }
}

void BlisTreeNode::callHeuristics() {
}

//#############################################################################

void
BlisTreeNode::applyConstraints(BcpsConstraintPool const * conPool)
{
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

    double cutup = broker_->getIncumbentValue() *
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


AlpsReturnStatus BlisTreeNode::encode(AlpsEncoded * encoded) const {
  AlpsReturnStatus status;
  status = AlpsTreeNode::encode(encoded);
  assert(status==AlpsReturnStatusOk);
  status = BcpsTreeNode::encode(encoded);
  assert(status==AlpsReturnStatusOk);
  return status;
}


/// Unpack into a new DcoTreeNode object and return a pointer to it.
AlpsKnowledge * BlisTreeNode::decode(AlpsEncoded & encoded) const {
  AlpsReturnStatus status;
  AlpsNodeDesc * new_node_desc = new BlisNodeDesc();
  BlisTreeNode * new_node = new BlisTreeNode(new_node_desc);
  new_node->setBroker(broker_);
  new_node_desc = NULL;
  status = new_node->decodeToSelf(encoded);
  assert(status==AlpsReturnStatusOk);
  return new_node;
}

/// Unpack into this from an encoded object.
AlpsReturnStatus BlisTreeNode::decodeToSelf(AlpsEncoded & encoded) {
  AlpsReturnStatus status;
  status = AlpsTreeNode::decodeToSelf(encoded);
  assert(status==AlpsReturnStatusOk);
  status = BcpsTreeNode::decodeToSelf(encoded);
  assert(status==AlpsReturnStatusOk);
  if (isPregnant()) {
    // delete branch object stored, we will re-create it. Ideally it should not
    // be sent in the first place.
    clearBranchObject();
    setStatus(AlpsNodeStatusEvaluated);
  }
  return status;
}


void
BlisTreeNode::convertToExplicit()
{
#ifdef BLIS_DEBUG
    std::cout << "BLIS: convertToExplicit(); explicit_="<<explicit_ << std::endl;
#endif

    if(!explicit_) {

	// Convert to explicit
	explicit_ = 1;

	BlisModel* model = dynamic_cast<BlisModel*>(broker_->getModel());
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


	for(i = static_cast<int> (model->leafToRootPath.size() - 1); i > -1; --i) {

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
