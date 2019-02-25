/*===========================================================================*
 * This file is part of the Branch, Constrain and Price Software (BiCePS)    *
 *                                                                           *
 * BiCePS is distributed under the Eclipse Public License as part of the     *
 * COIN-OR repository (http://www.coin-or.org).                              *
 *                                                                           *
 * Authors:                                                                  *
 *                                                                           *
 *          Yan Xu, Lehigh University                                        *
 *          Aykut Bulut, Lehigh University                                   *
 *          Ted Ralphs, Lehigh University                                    *
 *                                                                           *
 * Conceptual Design:                                                        *
 *                                                                           *
 *          Yan Xu, Lehigh University                                        *
 *          Ted Ralphs, Lehigh University                                    *
 *          Laszlo Ladanyi, IBM T.J. Watson Research Center                  *
 *          Matthew Saltzman, Clemson University                             *
 *                                                                           *
 * Copyright (C) 2001-2019, Lehigh University, Yan Xu, and Ted Ralphs.       *
 * All Rights Reserved.                                                      *
 *===========================================================================*/

#include "BcpsTreeNode.h"
#include "BcpsNodeDesc.h"
#include <sstream>

extern std::map<BCPS_Grumpy_Msg_Type, char const *> grumpyMessage;
extern std::map<BcpsNodeBranchDir, char> grumpyDirection;

/// Destructor.
BcpsTreeNode::~BcpsTreeNode() {
  clearBranchObject();
}

//todo(aykut) this method is not used yet. it should be implemented in user
//sub-class.
int BcpsTreeNode::process(bool isRoot, bool rampUp) {
  AlpsNodeStatus status = getStatus();
  BcpsModel * model = dynamic_cast<BcpsModel*>(broker()->getModel());
  CoinMessageHandler * message_handler = model->bcpsMessageHandler_;

  // debug stuff
  std::stringstream debug_msg;
  debug_msg << "Processing node ";
  debug_msg << this;
  debug_msg << " index ";
  debug_msg << getIndex();
  debug_msg << " parent ";
  debug_msg << getParent();
  message_handler->message(0, "Bcps", debug_msg.str().c_str(),
                           'G', BCPS_DLOG_PROCESS)
    << CoinMessageEol;
  // end of debug stuff

  // check if this can be fathomed
  if (getQuality() > broker()->getBestQuality()) {
    // debug message
    message_handler->message(0, "Bcps", "Node fathomed due to parent quality.",
                             'G', BCPS_DLOG_PROCESS);
    // end of debug message
    setStatus(AlpsNodeStatusFathomed);
    return AlpsReturnStatusOk;
  }

  if (status==AlpsNodeStatusCandidate ||
      status==AlpsNodeStatusEvaluated) {
    boundingLoop(isRoot, rampUp);
  }
  else if (status==AlpsNodeStatusBranched ||
           status==AlpsNodeStatusFathomed ||
           status==AlpsNodeStatusDiscarded) {
    // this should not happen
    message_handler->message(BCPS_NODE_UNEXPECTEDSTATUS, model->bcpsMessages_)
      << static_cast<int>(status) << CoinMessageEol;
  }
  return AlpsReturnStatusOk;
}

/// Clear branch object stored.
void BcpsTreeNode::clearBranchObject() {
  if (branchObject_) {
    delete branchObject_;
    branchObject_=NULL;
  }
}

int BcpsTreeNode::boundingLoop(bool isRoot, bool rampUp) {
  AlpsNodeStatus status = getStatus();

  BcpsModel * model = dynamic_cast<BcpsModel*>(broker_->getModel());
  CoinMessageHandler * message_handler = model->bcpsMessageHandler_;

  bool keepBounding = true;
  bool fathomed = false;
  bool do_branch = false;
  bool genConstraints = false;
  bool genVariables = false;
  BcpsConstraintPool * constraintPool = new BcpsConstraintPool();
  BcpsVariablePool * variablePool = new BcpsVariablePool();
  installSubProblem();

  while (keepBounding) {
    keepBounding = false;
    // solve subproblem corresponds to this node
    BcpsSubproblemStatus subproblem_status = bound();

    // debug print objective value after bounding
    std::stringstream debug_msg;
    debug_msg << "Subproblem solved. "
              << "status "
              << subproblem_status
              << " Obj value "
              << quality_
              << " estimate "
              << solEstimate_;
    message_handler->message(0, "Bcps", debug_msg.str().c_str(),
                             'G', BCPS_DLOG_PROCESS);
    // end of debug stuff

    // call heuristics to search for a solution
    callHeuristics();

    // decide what to do
    branchConstrainOrPrice(subproblem_status, keepBounding, do_branch,
                           genConstraints,
                           genVariables);

    // debug message
    // reset debug message
    debug_msg.str(std::string());
    debug_msg << "BCP function decided to"
              << " keep bounding "
              << keepBounding
              << " branch "
              << do_branch
              << " generate cons "
              << genConstraints;
    message_handler->message(0, "Bcps", debug_msg.str().c_str(),
                             'G', BCPS_DLOG_PROCESS);
    // end of debug stuff

    if (getStatus()==AlpsNodeStatusFathomed) {
      // node is fathomed, nothing to do.
      break;
    }
    else if (keepBounding && genConstraints) {
      generateConstraints(constraintPool);
      // add constraints to the model
      applyConstraints(constraintPool);
      // clear constraint pool
      constraintPool->freeGuts();
      // set status to evaluated
      setStatus(AlpsNodeStatusEvaluated);
    }
    else if (keepBounding && genVariables) {
      generateVariables(variablePool);
      // add variables to the model
      // set status to evaluated
      setStatus(AlpsNodeStatusEvaluated);
    }
    else if (keepBounding==false && do_branch==false) {
      // put node back into the list.
      // this means do not change the node status and end processing the node.
      // set status to evaluated
      setStatus(AlpsNodeStatusEvaluated);
    }
    else if (keepBounding==false && do_branch) {
      // // prepare for branch() call
      // BcpsBranchStrategy * branchStrategy = model->branchStrategy();
      // // todo(aykut) following should be a parameter
      // // Maximum number of resolve during branching.
      // int numBranchResolve = 10;
      // // todo(aykut) why ub should be an input?
      // branchStrategy->createCandBranchObjects(this);
      // // prepare this node for branching, bookkeeping for differencing.
      // // call pregnant setting routine
      // processSetPregnant();
    }
    else {
      message_handler->message(9998, "Bcps", "This should not happen. "
                               " branchConstrainOrPrice() is buggy.", 'E', 0)
        << CoinMessageEol;
    }

  }
  delete constraintPool;
  delete variablePool;
  return AlpsReturnStatusOk;
}

//todo(aykut) Warm start stuff should be implemented without assuming much
//about the underlying solver.
void BcpsTreeNode::processSetPregnant() {
  // get warm start basis from solver
  // todo(aykut) This does not help much if the underlying solver is an IPM
  // based solver.
  // BcpsModel * model = dynamic_cast<BcpsModel*>(broker()->getModel());
  // CoinWarmStartBasis * ws = dynamic_cast<CoinWarmStartBasis*>
  //   (model->solver()->getWarmStart());
  // // store basis in the node desciption.
  // getDesc()->setBasis(ws);
  // // set status pregnant
  // setStatus(AlpsNodeStatusPregnant);

}


/// Encode the content of this into the given AlpsEncoded object.
AlpsReturnStatus BcpsTreeNode::encode(AlpsEncoded * encoded) const {
  AlpsReturnStatus status = AlpsReturnStatusOk;
  // int type = 0;
  // bool branch_obj_exists;
  // branch_obj_exists = branchObject_ ? true : false;
  // encoded->writeRep(branch_obj_exists);
  // if (branch_obj_exists) {
  //   status = branchObject_->encode(encoded);
  //   assert(status==AlpsReturnStatusOk);
  // }
  return status;
}

/// Decode the given AlpsEncoded object into this.
AlpsReturnStatus BcpsTreeNode::decodeToSelf(AlpsEncoded & encoded) {
  AlpsReturnStatus status = AlpsReturnStatusOk;
  // bool branch_obj_exists;
  // encoded.readRep(branch_obj_exists);
  // if (branch_obj_exists) {
  //   if (branchObject_) {
  //     status = branchObject_->decodeToSelf(encoded);
  //     assert(status==AlpsReturnStatusOk);
  //   }
  // }
  return status;
}
