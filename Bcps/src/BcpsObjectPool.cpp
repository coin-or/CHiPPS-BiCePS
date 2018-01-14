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
 * Copyright (C) 2001-2015, Lehigh University, Yan Xu, and Ted Ralphs.       *
 * All Rights Reserved.                                                      *
 *===========================================================================*/

#include "BcpsObjectPool.h"

BcpsObjectPool::BcpsObjectPool()
  :AlpsKnowledgePool(AlpsKnowledgePoolTypeUndefined) {
  objects_.clear();
}

BcpsObjectPool::~BcpsObjectPool() {
  freeGuts();
}

void BcpsObjectPool::freeGuts() {
  int size = static_cast<int>(objects_.size());
  for (int i=0; i<size; ++i) {
    delete objects_[i];
  }
  objects_.clear();
}

int BcpsObjectPool::getNumKnowledges() const {
  return static_cast<int>(objects_.size());
}

/// Query a knowledge, but doesn't remove it from the pool
std::pair<AlpsKnowledge*, double> BcpsObjectPool::getKnowledge() const {
  return std::make_pair(objects_[0], 0.0);
}

/// Check whether the pool has knowledge.
bool BcpsObjectPool::hasKnowledge() const {
  return objects_.empty() ? false : true;
}

int BcpsObjectPool::getMaxNumKnowledges() const {
  // todo(aykut) this should be removed in Alps level.
  std::cerr << "Not implemented." << std::endl;
  throw std::exception();
}

std::pair<AlpsKnowledge*, double> BcpsObjectPool::getBestKnowledge() const {
  std::cerr << "Not implemented." << std::endl;
  throw std::exception();
  std::pair<AlpsKnowledge*, double> ret;
  return ret;
}

/// Get a reference to all the knowledges in the pool.*/
void BcpsObjectPool::getAllKnowledges (std::vector<std::pair<AlpsKnowledge*,
                       double> >& kls) const {
  std::cerr << "Not implemented." << std::endl;
  throw std::exception();
}

/// Add a knowledge to pool
void BcpsObjectPool::addKnowledge(AlpsKnowledge * nk, double priority) {
  objects_.push_back(nk);
}

/// Pop the first knowledge from the pool.
void BcpsObjectPool::popKnowledge() {
  std::cerr << "Not implemented." << std::endl;
  throw std::exception();
}

/// Set the quantity limit of knowledges that can be stored in the pool.
void BcpsObjectPool::setMaxNumKnowledges(int num) {
  std::cerr << "Not implemented." << std::endl;
  throw std::exception();
}

/// Delete object k from pool
void BcpsObjectPool::deleteObject(int k) {
  assert(k > -1 && k < ((int)objects_.size()));
  AlpsKnowledge *objectK = getObject(k);
  std::vector<AlpsKnowledge *>::iterator pos;
  pos = objects_.begin() + k;
  objects_.erase(pos);
  // Free memory of object k.
  delete objectK;
}
