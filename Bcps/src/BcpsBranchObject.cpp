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
BcpsBranchObject::BcpsBranchObject(BcpsModel * model, int type, int index,
                                   int score)
  : model_(model), type_(type), index_(index), score_(score), value_(0.0) {
}

/// Constructor.
BcpsBranchObject::BcpsBranchObject(BcpsModel * model, int type, int index,
                                   double score, double value)
  : model_(model), type_(type), index_(index), score_(score), value_(value) {
}

/// Copy constructor.
BcpsBranchObject::BcpsBranchObject(BcpsBranchObject const & other)
  : model_(other.model()), type_(other.type()), index_(other.index()),
    score_(other.score()), value_(other.value()) {
}

/// Copy assignment operator
BcpsBranchObject & BcpsBranchObject::operator=(BcpsBranchObject const & rhs) {
  model_ = rhs.model();
  type_ = rhs.type();
  index_ = rhs.index();
  score_ = rhs.score();
  value_ = rhs.value();
  return *this;
}

/// Pack Bcps portion to an encoded object.
AlpsReturnStatus BcpsBranchObject::encodeBcps(AlpsEncoded * encoded) const {
  AlpsReturnStatus status = AlpsReturnStatusOk;
  assert(encoded);
  encoded->writeRep(model_);
  encoded->writeRep(type_);
  encoded->writeRep(index_);
  encoded->writeRep(score_);
  encoded->writeRep(value_);
  return status;
}

/// Unpack Bcps portion from an encoded object.
AlpsReturnStatus BcpsBranchObject::decodeBcps(AlpsEncoded & encoded) {
  AlpsReturnStatus status = AlpsReturnStatusOk;
  encoded.readRep(model_);
  encoded.readRep(type_);
  encoded.readRep(index_);
  encoded.readRep(score_);
  encoded.readRep(value_);
  return status;
}
