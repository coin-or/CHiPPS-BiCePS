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
 * Copyright (C) 2001-2018, Lehigh University, Yan Xu, and Ted Ralphs.       *
 * All Rights Reserved.                                                      *
 *===========================================================================*/

#include "Bcps.h"
#include "BcpsModel.h"

AlpsReturnStatus BcpsModel::encode(AlpsEncoded * encoded) const {
  AlpsReturnStatus status = AlpsReturnStatusOk;
  int numCons = static_cast<int> (constraints_.size());
  encoded->writeRep(numCons);
  for (int i=0; i<numCons; ++i) {
    status = constraints_[i]->encode(encoded);
    assert(status==AlpsReturnStatusOk);
  }
  int numVars =  static_cast<int> (variables_.size());
  encoded->writeRep(numVars);
  for (int i=0; i<numVars; ++i) {
    status = variables_[i]->encode(encoded);
    assert(status==AlpsReturnStatusOk);
  }
  // Core node decription
  encoded->writeRep(numCoreConstraints_);
  encoded->writeRep(numCoreVariables_);
  return status;
}
