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


#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include "BcpsModel.h"
#include "BcpsBranchObject.h"


/// Constructor.
BcpsBranchObject::BcpsBranchObject(int type, int index, int score)
  : AlpsKnowledge(AlpsKnowledgeTypeUndefined, NULL),
    type_(type), index_(index), score_(score), value_(0.0) {

}

/// Constructor.
BcpsBranchObject::BcpsBranchObject(int type, int index, double score,
                                   double value)
  : AlpsKnowledge(AlpsKnowledgeTypeUndefined, NULL),
    type_(type), index_(index), score_(score), value_(value) {
}

/// Copy constructor.
BcpsBranchObject::BcpsBranchObject(BcpsBranchObject const & other)
  : AlpsKnowledge(AlpsKnowledgeTypeUndefined, other.broker_),
    type_(other.type()), index_(other.index()),
    score_(other.score()), value_(other.value()) {
}

/// Copy assignment operator
BcpsBranchObject & BcpsBranchObject::operator=(BcpsBranchObject const & rhs) {
  type_ = rhs.type();
  index_ = rhs.index();
  score_ = rhs.score();
  value_ = rhs.value();
  broker_ = rhs.broker_;
  return *this;
}

/// Pack Bcps portion to an encoded object.
AlpsReturnStatus BcpsBranchObject::encode(AlpsEncoded * encoded) const {
  AlpsReturnStatus status = AlpsReturnStatusOk;
  assert(encoded);
  encoded->writeRep(type_);
  encoded->writeRep(index_);
  encoded->writeRep(score_);
  encoded->writeRep(value_);
  return status;
}

/// Unpack Bcps portion from an encoded object.
AlpsReturnStatus BcpsBranchObject::decodeToSelf(AlpsEncoded & encoded) {
  AlpsReturnStatus status = AlpsReturnStatusOk;
  encoded.readRep(type_);
  encoded.readRep(index_);
  encoded.readRep(score_);
  encoded.readRep(value_);
  return status;
}
