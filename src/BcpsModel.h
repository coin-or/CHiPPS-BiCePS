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
 * Copyright (C) 2001-2023, Lehigh University, Yan Xu, and Ted Ralphs.       *
 * All Rights Reserved.                                                      *
 *===========================================================================*/

#ifndef BcpsModel_h_
#define BcpsModel_h_

//#############################################################################

#include <vector>

#include "CoinMessageHandler.hpp"

#include "AlpsModel.h"
#include "AlpsKnowledgeBroker.h"

#include "BcpsConfig.h"
#include "BcpsMessage.h"
#include "BcpsObject.h"

//#############################################################################

class BCPSLIB_EXPORT BcpsModel : public AlpsModel {
protected:
  /** Constraints input by users (before preprocessing). */
  std::vector<BcpsConstraint *> constraints_;
  /** Variables input by users (before preprocessing). */
  std::vector<BcpsVariable *> variables_;
  /** Number of core constraints. By default, all input constraints are
      core. */
  int numCoreConstraints_;
  /** Number of core variables. By default, all input variables are core.*/
  int numCoreVariables_;

public:
  /** Message handler. */
  CoinMessageHandler * bcpsMessageHandler_;
  /** Bcps messages. */
  CoinMessages bcpsMessages_;

public:
  ///@name Constructor and Destructor
  //@{
  /// Default constructor.
  BcpsModel() : numCoreConstraints_(0), numCoreVariables_(0) {
    bcpsMessageHandler_ = new CoinMessageHandler();
    bcpsMessageHandler_->setLogLevel(2);
    bcpsMessages_ = BcpsMessage();
  }
  /// Destructor
  virtual ~BcpsModel() {
    int i = 0;
    int size = static_cast<int> (constraints_.size());
    for (i = 0; i < size; ++i) {
      delete constraints_[i];
    }
    size =  static_cast<int> (variables_.size());
    for (i = 0; i < size; ++i) {
      delete variables_[i];
    }
    delete bcpsMessageHandler_;
  }

  /** Get variables and constraints */
  /**@{*/
  std::vector<BcpsConstraint *> & getConstraints() { return constraints_; }
  std::vector<BcpsVariable *> & getVariables() { return variables_; }

  int getNumCoreConstraints() const { return numCoreConstraints_; }
  int getNumCoreVariables() const { return numCoreVariables_; }
  /**@}*/

  /** Set variables and constraints */
  /**@{*/
  void setConstraints(BcpsConstraint **con, int size) {
    for (int j = 0; j < size; ++j) {
      constraints_.push_back(con[j]);
    }
    int numCons = getNumCoreConstraints();
    setNumCoreConstraints(numCons+size);
  }
  void setNumCoreConstraints(int num) { numCoreConstraints_ = num; }

  void setVariables(BcpsVariable **var, int size) {
    for (int j = 0; j < size; ++j) {
      variables_.push_back(var[j]);
    }
    int numVars = getNumCoreVariables();
    setNumCoreVariables(numVars+size);
  }
  void setNumCoreVariables(int num) { numCoreVariables_ = num; }
  /**@}*/

  /** Add variables and constraints */
  /**@{*/
  void addConstraint(BcpsConstraint *con, bool isCore = true) {
    constraints_.push_back(con);
    if (isCore) {
      numCoreConstraints_++;
    }
  }
  void addVariable(BcpsVariable *var, bool isCore = true) {
    variables_.push_back(var);
    if (isCore) {
      numCoreVariables_++;
    }
  }
  /**@}*/

  /** Return list of variables. */
  std::vector<BcpsVariable *> getVariables() const { return variables_; }

  /** Return list of constraints. */
  std::vector<BcpsConstraint *> getConstrints() const { return constraints_; }

  /** Get the message handler. */
  CoinMessageHandler * bcpsMessageHandler() const
  { return bcpsMessageHandler_; }

  /** Return messages. */
  CoinMessages bcpsMessages() { return bcpsMessages_; }

  ///@name Encode and Decode functions
  ///@{
  /// Get encode function defined in #AlpsKnowledge.
  using AlpsModel::encode;
  /// Pack Bcps data fields
  virtual AlpsReturnStatus encode(AlpsEncoded * encoded) const;
  //@}

};


#endif
