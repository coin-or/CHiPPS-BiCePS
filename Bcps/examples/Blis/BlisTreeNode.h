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

#ifndef BlisTreeNode_h_
#define BlisTreeNode_h_

//#############################################################################

#include "AlpsNodeDesc.h"

#include "BcpsObjectPool.h"
#include "BcpsTreeNode.h"

#include "BcpsNodeDesc.h"
#include "BlisNodeDesc.h"

class BcpsModel;
class BlisModel;


//#############################################################################
/** This is the class in which we are finally able to concretely define the
    bounding procedure. Here we can assume that we have an LP solver and that
    the objects are cuts and variables, etc. */
//#############################################################################


class BlisTreeNode : public BcpsTreeNode {

 private:

    /** No copy constructor, assignment operator. */
    BlisTreeNode(const BlisTreeNode&);

    BlisTreeNode& operator=(const BlisTreeNode&);

    /** Constraint pool. */
    //BcpsConstraintPool *constraintPool_;

    /** Variable pool. */
    //BcpsVariablePool *variablePool_;

    /** Save an explicit node description. */
    //void saveExplicit();

 public:

    /** Default constructor. */
    BlisTreeNode()
        :
        BcpsTreeNode()
        { init(); }

    /** Useful constructor. */
    BlisTreeNode(BlisModel* m) {
        init();
        desc_ = new BlisNodeDesc(m);
    }

    /** Useful constructor. */
    BlisTreeNode(AlpsNodeDesc *&desc) {
        init();
        desc_ = desc;
        desc = NULL;
    }

    /** Destructor. */
    virtual ~BlisTreeNode() {
        //delete constraintPool_;
        //delete variablePool_;
    }

    /** Initilize member data when constructing a node. */
    void init() {
        //constraintPool_ = new BcpsConstraintPool;
        //variablePool_ = new BcpsVariablePool;
    }

    /** Create a new node based on given desc. */
    AlpsTreeNode* createNewTreeNode(AlpsNodeDesc *&desc) const;

    /** Convert explicit description to difference, and vise-vesa */
    ///@{
    virtual void convertToExplicit();
    virtual void convertToRelative();
    ///@}


    /** Performing the bounding operation. */
    virtual int process(bool isRoot = false, bool rampUp = false);


    /** Takes the explicit description of the current active node and
        creates the children's descriptions, which contain information
        about how the branching is to be done. The stati of the children
        are AlpsNodeStatusCandidate. */
    virtual std::vector< CoinTriple<AlpsNodeDesc*, AlpsNodeStatus, double> >
	branch();

    BlisNodeDesc * getDesc() const {
      return dynamic_cast<BlisNodeDesc*>(AlpsTreeNode::getDesc());
    }

    /** Select a branching object based on give branching strategy. */
    int selectBranchObject(BlisModel *model,
                           bool& foundSol,
                           int numPassesLeft);

    /** Fix and tighten varaibles based optimality conditions. */
    int reducedCostFix(BlisModel *model);

    ///@name Pure virtual functions inherited from Bcps.
    //@{
    virtual int generateConstraints(BcpsConstraintPool *conPool);
    virtual int generateVariables(BcpsVariablePool *varPool);
    virtual int chooseBranchingObject();
    virtual int installSubProblem();
    virtual void branchConstrainOrPrice(BcpsSubproblemStatus subproblem_status,
                                        bool & keepBounding,
                                        bool & branch,
                                        bool & generateConstraints,
                                        bool & generateVariables);
    virtual BcpsSubproblemStatus bound();
    virtual void callHeuristics();
    virtual void applyConstraints(BcpsConstraintPool const * conPool);
    //@}

    /// Sets node status to pregnant and carries necessary operations.
    void processSetPregnant();
    void boundingLoop();
    void copyFullNode(BlisNodeDesc * child_node) const;

    ///@name Encode and Decode functions for parallel execution
    //@{
    /// Get encode from #AlpsKnowledge
    using AlpsKnowledge::encode;
    /// Pack this into an encoded object.
    virtual AlpsReturnStatus encode(AlpsEncoded * encoded) const;
    /// Unpack into this from an encoded object.
    virtual AlpsReturnStatus decodeToSelf(AlpsEncoded & encoded);
    /// Unpack into a new DcoTreeNode object and return a
    /// pointer to it.
    virtual AlpsKnowledge * decode(AlpsEncoded & encoded) const;
    //@}

};

#endif
