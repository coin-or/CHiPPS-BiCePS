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


#include "BlisBranchStrategyPseudo.h"
#include "BlisModel.h"
#include "BlisMessage.h"
#include "BlisTreeNode.h"
#include "BlisBranchObjectInt.h"

BlisBranchStrategyPseudo::BlisBranchStrategyPseudo(BlisModel * model):
  BcpsBranchStrategy(model) {
  setType(BLIS_BS_PSEUDOCOST);
  // todo(aykut) think about parametrizing this.
  score_factor_ = 1.0/6.0;
  // all integer variables are relaed.
  int num_relaxed = model->getNumIntVars();
  down_num_ = new int[num_relaxed]();
  up_num_ = new int[num_relaxed]();
  down_derivative_ = new double[num_relaxed]();
  up_derivative_ = new double[num_relaxed]();
  // fill reverse map
  int const * relaxed_cols = model->getIntVars();
  for (int i=0; i<num_relaxed; ++i) {
    rev_relaxed_[relaxed_cols[i]] = i;
  }
}

BlisBranchStrategyPseudo::~BlisBranchStrategyPseudo() {
  if (down_num_) {
    delete[] down_num_;
    down_num_ = NULL;
  }
  if (up_num_) {
    delete[] up_num_;
    up_num_ = NULL;
  }
  if (down_derivative_) {
    delete[] down_derivative_;
    down_derivative_ = NULL;
  }
  if (up_derivative_) {
    delete[] up_derivative_;
    up_derivative_ = NULL;
  }
}

int BlisBranchStrategyPseudo::createCandBranchObjects(BcpsTreeNode * node) {
  // get node
  BlisTreeNode * blis_node = dynamic_cast<BlisTreeNode*>(node);
  // update statistics
  update_statistics(blis_node);
  // get blis model and message stuff
  BlisModel * blis_model = dynamic_cast<BlisModel*>(model());
  // get number of relaxed columns
  // we assume all relaxed columns are integer variables.
  int num_relaxed = blis_model->getNumIntVars();
  // get indices of relaxed object
  int const * relaxed = blis_model->getIntVars();
  // store branch objects in bobjects
  std::vector<BcpsBranchObject*> bobjects;
  // iterate over relaxed columns and populate bobjects
  for (int i=0; i<num_relaxed; ++i) {
    int preferredDir;
    BcpsObject * curr_object = blis_model->getVariables()[relaxed[i]];
    double infeasibility = curr_object->infeasibility(blis_model, preferredDir);
    // check the amount of infeasibility
    if (infeasibility != 0.0) {
      double min = std::min(down_derivative_[i], up_derivative_[i]);
      double max = std::max(down_derivative_[i], up_derivative_[i]);
      // compute score
      double score = score_factor_*max + (1.0-score_factor_)*min;
      // create a branch object for this
      BcpsBranchObject * cb =
        curr_object->createBranchObject(blis_model, preferredDir);
      // set score
      cb->setScore(score);
      bobjects.push_back(cb);
    }
  }
  // add branch objects to branchObjects_
  setBranchObjects(bobjects);
  // bobjects are now owned by BcpsBranchStrategy, do not free them.
  bobjects.clear();
  // set the branch object member of the node
  blis_node->setBranchObject(new BlisBranchObjectInt(bestBranchObject()));
  // compare branch objects and keep the best one at bestBranchObject_
  return 0;
}

int
BlisBranchStrategyPseudo::betterBranchObject(BcpsBranchObject const * current,
                                            BcpsBranchObject const * other) {
  int res;
  if (current->score()>other->score()) {
    res = 1;
  }
  else {
    res = 0;
  }
  return res;
}

void BlisBranchStrategyPseudo::update_statistics(BlisTreeNode * node) {
  // return if this is the root node
  if (node->getParent()==NULL) {
    return;
  }
  // get quality_ of this node, quality is sense*value
  double quality = node->getQuality();
  // get quality_ of the parent node
  double parent_quality = node->getParent()->getQuality();
  // is this node a down or up branch
  int dir = node->getDesc()->getBranchedDir();
  // index of the branched variable for the current node
  int branched_index = rev_relaxed_[node->getDesc()->getBranchedInd()];
  double branched_value = node->getDesc()->getBranchedVal();

  // update statistics
  double frac;
  if (dir==-1) {
    frac = branched_value-floor(branched_value);
    double deriv = (quality-parent_quality) / frac;
    int n = down_num_[branched_index];
    double old = down_derivative_[branched_index];
    down_derivative_[branched_index] = (old*n + deriv)/(n+1);
    down_num_[branched_index]++;
  }
  else if (dir==1) {
    frac = ceil(branched_value)-branched_value;
    double deriv = (quality-parent_quality) / frac;
    int n = up_num_[branched_index];
    double old = up_derivative_[branched_index];
    up_derivative_[branched_index] = (old*n + deriv)/(n+1);
    up_num_[branched_index]++;
  }
  else {
    std::cerr << "Invalid branching direction!" << std::endl;
    throw std::exception();
  }
}
