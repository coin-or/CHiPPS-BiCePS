/*===========================================================================*
 * This file is part of the Branch, Constrain and Price Software (BiCePS)    *
 *                                                                           *
 * BiCePS is distributed under the Eclipse Public License as part of the     *
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
 * Copyright (C) 2001-2017, Lehigh University, Yan Xu, and Ted Ralphs.       *
 * All Rights Reserved.                                                      *
 *===========================================================================*/

#ifndef BcpsTreeNode_h_
#define BcpsTreeNode_h_

#include "BcpsNodeDesc.h"

#include <vector>

#include "AlpsTreeNode.h"
#include "AlpsNodeDesc.h"

#include "BcpsBranchObject.h"
#include "BcpsObjectPool.h"

/*!
   #BcpsTreeNode class represents a tree node in a Branch, Constraint and Price
   type of tree. Alps does not assume an optimization model, Bcps
   does. BcpsNodeTree holds auxilary data, like branching object
   etc. #BcpsNodeDesc holds the data corresponding to the subproblem the node
   stands for, ie. variable bounds, etc. A differencing scheme is implemented
   in #BcpsNodeDesc. It stores only the differences with respect to parent.

   #BcpsTreeNode fixes and API that should be implemented by the user sub-class.
   Virtual functions that should be implemented by the user are

   <ul>

   <li> ::generateConstraints(). User sub-class should generate constrints and
        populate the input pool.

   <li> ::generateVariables(). User sub-class should generate new varaibles and
        populate the given pool.

   <li> ::chooseBranchingObject(). User sub-class should implement choosing a
        branch object.

   <li> ::installSubProblem(). User sub-class should implement how to install
        the subproblem to the underlying solver. The subproblem data are stored
        in desc_, which is a pointer to an instance of the user node
        description class.

   <li> ::branchConstrainOrPrice(). During bounding loop, this function decides
        what to do next, i.e., keep bounding, branch, generate constraints,
        generate columns.

   </ul>

*/

class BcpsTreeNode : public AlpsTreeNode {
protected:
  /// Branching object for this node, which has information of how to
  /// execute branching.
  BcpsBranchObject * branchObject_;
  /// Bounding loop, generate cols, vars and keeps bounding the subproblem.
  int boundingLoop(bool isRoot, bool rampUp);
  /// Helper function for process method. Sets status to pregnant in
  /// boundingLoop() called from process().
  void processSetPregnant();

public:
  ///@name Constructors and Destructor.
  //@{
  /// Default constructor.
  BcpsTreeNode(): branchObject_(NULL) { }
  /// Destructor.
  virtual ~BcpsTreeNode();
  //@}

  ///@name Pure virtual functions that should be implemented by user
  /// sub-class
  //@{
  /// Generate constraints. The generated constraints are stored
  /// in constraint pool. The default implementation does nothing.
  virtual int generateConstraints(BcpsConstraintPool *conPool) = 0;
  /// Generate variables. The generated varaibles are stored
  /// in variable pool. The default implementation does nothing.
  virtual int generateVariables(BcpsVariablePool *varPool) = 0;
  /// Choose a branching object.
  virtual int chooseBranchingObject() = 0;
  /** Extract node information (bounds, constraints, variables) from
      this node and load the information into the relaxation solver,
      such as linear programming solver.*/
  virtual int installSubProblem() = 0;
  /** Once the solver is done, decide what to do next:
   *  - relaxed feasible but not integer feasible,
   *  - integer feasible,
   *  - infeasible,
   *  - unbounded,
   *  - fathomed (for instance, reaching objective limit for LP).
   * Set node status accordingly.
   * \param subproblem_status   Input The solution status of bounding.
   * \param keepBounding        Output  Whether to keep on bounding.
   * \param branch              Output        Whether to branch this.
   * \param generateConstraints Output Whether to generate constraints for the
   * subproblem.
   * \param generateVariables   Output Whether to generate variables for the
   * subproblem.
   */
  virtual void branchConstrainOrPrice(BcpsSubproblemStatus subproblem_status,
                                     bool & keepBounding,
                                     bool & branch,
                                     bool & generateConstraints,
                                     bool & generateVariables) = 0;
  /// Bounding procedure to estimate quality of this node. Typically involves
  /// solution of the corresponding subproblem.
  virtual BcpsSubproblemStatus bound() = 0;
  /// Call heuristics to search for a solution in the current subproblem.
  virtual void callHeuristics() = 0;
  virtual void applyConstraints(BcpsConstraintPool const * conPool) = 0;
//@}


  ///@name Other functions
  //@{
  /** This methods performs the processing of the node. For branch and bound,
      this would mean performing the bounding operation. This method updates
      the status of the node. Alps level tree manager does the rest.
  */
  virtual int process(bool isRoot = false, bool rampUp = false);
  /// Clear branch object stored.
  void clearBranchObject();
  //@}


  ///@name Get/Set functions
  //@{
  /// Return the branching object.
  const BcpsBranchObject * branchObject() const { return branchObject_; }
  /// Set the branching object.
  void setBranchObject(BcpsBranchObject * b) { branchObject_ = b; }
  //@}

  ///@name Encode and Decode functions
  //@{
  /// Encode the content of this into the given AlpsEncoded object.
  virtual AlpsReturnStatus encode(AlpsEncoded * encoded) const;
  /// Decode the given AlpsEncoded object into this.
  virtual AlpsReturnStatus decodeToSelf(AlpsEncoded & encoded);
  //@}
};

#endif
