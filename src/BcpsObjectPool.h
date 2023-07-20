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

#ifndef BcpsObjectPool_h_
#define BcpsObjectPool_h_

#include <vector>

#include "AlpsKnowledgePool.h"

#include "BcpsConfig.h"
#include "BcpsObject.h"

//#############################################################################
/** Object pool is used to store objects */
//#############################################################################

class BCPSLIB_EXPORT BcpsObjectPool : public AlpsKnowledgePool {
protected:
  std::vector<AlpsKnowledge *> objects_;

public:
  ///@name Constructors and Destructor.
  //@{
  /// Default constructor.
  BcpsObjectPool();
  /// Destructor.
  virtual ~BcpsObjectPool();
  //@}

  ///@name Other functions
  //@{
  /// Delete object k from pool
  void deleteObject(int k);
  /** Get all objects. */
  std::vector<AlpsKnowledge *> const & getObjects() const { return objects_; }
  /** Get a object. */
  AlpsKnowledge * getObject(int k) const { return objects_[k]; }
  /// Free stored objects.
  void freeGuts();
  //@}

  ///@name Querry methods, inherited from AlpsKnowledgePool
  //@{
  /// Return size of the pool.
  virtual int getNumKnowledges() const;
  /// Check the first item in the pool.
  virtual std::pair<AlpsKnowledge*, double> getKnowledge() const;
  /// Check whether the pool is empty.
  virtual bool hasKnowledge() const;
  /// Query the quantity limit of knowledges.
  virtual int getMaxNumKnowledges() const;
  /// Query the best knowledge in the pool.
  virtual std::pair<AlpsKnowledge*, double> getBestKnowledge() const;
  /// Get a reference to all the knowledges in the pool.*/
  virtual void getAllKnowledges (std::vector<std::pair<AlpsKnowledge*,
                                 double> >& kls) const;
  //@}

  ///@name Knowledge manipulation, inherited from AlpsKnowledgePool
  //@{
  /// Add a knowledge to pool.
  virtual void addKnowledge(AlpsKnowledge * nk, double priority);
  /// Pop the first knowledge from the pool.
  virtual void popKnowledge();
  //@}

  ///@name Other functions
  //@{
  /// Set the quantity limit of knowledges that can be stored in the pool.
  virtual void setMaxNumKnowledges(int num);
  //@}

private:
  BcpsObjectPool(BcpsObjectPool const & other);
  BcpsObjectPool & operator=(BcpsObjectPool const & rhs);
};

//#############################################################################

class BCPSLIB_EXPORT BcpsConstraintPool : public BcpsObjectPool {
public:
  BcpsConstraintPool(): BcpsObjectPool() { }
    virtual ~BcpsConstraintPool() {}

    /** Add a constraint to pool */
    void addConstraint(BcpsConstraint * con) { objects_.push_back(con); }

    /** Delete constraint k from pool */
    void deleteConstraint(int k) { deleteObject(k); }

    /** Query how many constraints are in the pool.*/
    int getNumConstraints() const { return getNumKnowledges(); }

    /** Get the vector of constraints. */
    const std::vector<AlpsKnowledge *>& getConstraints() const {return objects_;}

    /** Get a constraints. */
    AlpsKnowledge *getConstraint(int k) const {return getObject(k);}
private:
  BcpsConstraintPool(BcpsConstraintPool const & other);
  BcpsConstraintPool & operator=(BcpsConstraintPool const & rhs);
};

//#############################################################################

class BCPSLIB_EXPORT BcpsVariablePool : public BcpsObjectPool {
public:
  BcpsVariablePool(): BcpsObjectPool() { }

  virtual ~BcpsVariablePool() {}

    /** Add a variable to pool */
    void addVariable(BcpsVariable * var) { objects_.push_back(var); }

    /** Delete variable k from pool */
    void deleteVariable(int k) { deleteObject(k); }

    /** Query how many variables are in the pool.*/
    int getNumVariables() const { return getNumKnowledges(); }

    /** Get the vector of variables. */
    const std::vector<AlpsKnowledge *>& getVariables() const {return objects_;}

    /** Get the vector of variables. */
    AlpsKnowledge *getVariable(int k) const {return getObject(k);}
private:
  BcpsVariablePool(BcpsVariablePool const & other);
  BcpsVariablePool & operator=(BcpsVariablePool const & rhs);
};

//#############################################################################

#endif // End of file
