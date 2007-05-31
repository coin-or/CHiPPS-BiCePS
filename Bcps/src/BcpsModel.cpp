/*===========================================================================*
 * This file is part of the Branch, Constrain and Price Software (BiCePS)    *
 *                                                                           *
 * BiCePS is distributed under the Common Public License as part of the      *
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
 * Copyright (C) 2001-2007, Lehigh University, Yan Xu, and Ted Ralphs.       *
 * All Rights Reserved.                                                      *
 *===========================================================================*/

#include "Bcps.h"
#include "BcpsModel.h"

//#############################################################################

AlpsReturnStatus 
BcpsModel::encodeBcps(AlpsEncoded *encoded) const
{
    AlpsReturnStatus status  = AlpsReturnStatusOk;
    int i;
    int size = static_cast<int> (constraints_.size());
    encoded->writeRep(size);
    for (i = 0; i < size; ++i) {
        constraints_[i]->encode(encoded); 
    }
    
    size =  static_cast<int> (variables_.size());
    encoded->writeRep(size);
    for (i = 0; i < size; ++i) {
        variables_[i]->encode(encoded);
    }
    
    // Core node decription
    encoded->writeRep(numCoreConstraints_);
    encoded->writeRep(numCoreVariables_);
    
    return status;
}

//#############################################################################

AlpsReturnStatus 
BcpsModel::decodeBcps(AlpsEncoded &encoded)
{
    AlpsReturnStatus status  = AlpsReturnStatusOk;
    int i, size;
    
    encoded.readRep(size);
    for (i = 0; i < size; ++i) {
        const AlpsKnowledge* know = 
            broker_->decoderObject(BcpsKnowledgeTypeConstraint);
        BcpsConstraint *con = NULL;
        con = static_cast<BcpsConstraint *>(know->decode(encoded));
        

        //BcpsConstraint *con = static_cast<BcpsConstraint *>
        //( broker_->decoderObject(BcpsKnowledgeTypeConstraint)->decode(encoded) );

        constraints_.push_back(con);
        con = NULL;
    }

    encoded.readRep(size);
    for (i = 0; i < size; ++i) {
        BcpsVariable *var = static_cast<BcpsVariable *>
            ( broker_->decoderObject(BcpsKnowledgeTypeVariable)->decode(encoded) );

        variables_.push_back(var);
        var = NULL;
    }
    
    // Core node decription
    encoded.readRep(numCoreConstraints_);
    encoded.readRep(numCoreVariables_);

    return status;
}

//#############################################################################

