/*===========================================================================*
 * This file is part of the Branch, Constrain and Price Software (BiCePS)    *
 *                                                                           *
 * BiCePS is distributed under the Common Public License as part of the      *
 * COIN-OR repository (http://www.coin-or.org).                              *
 *                                                                           *
 * Authors: Yan Xu, Lehigh University                                       *
 *          Ted Ralphs, Lehigh University                                    *
 *          Laszlo Ladanyi, IBM T.J. Watson Research Center                  *
 *          Matthew Saltzman, Clemson Unniversity                            *
 *                                                                           *
 * Copyright (C) 2001-2007, Lehigh University, Yan Xu, and Ted Ralphs.       *
 * All Rights Reserved.                                                      *
 *===========================================================================*/

#ifndef BcpsModel_h_
#define BcpsModel_h_

//#############################################################################

#include <vector>

#include "AlpsKnowledgeBroker.h"

#include "CoinMessageHandler.hpp"

#include "AlpsModel.h"

#include "BcpsObject.h"
#include "BcpsMessage.h"

//#############################################################################
//#############################################################################

class BcpsModel : public AlpsModel {

 protected:

    /** Constraints input by users (before preprocessing). */
    std::vector<BcpsConstraint *> constraints_;

    /** Number of core constraints. By default, all input constraints are
        core. 
    */
    int numCoreConstraints_;   
    
    /** Variables input by users (before preprocessing). */
    std::vector<BcpsVariable *> variables_;

    /** Number of core variables. By default, all input variables are core.*/
    int numCoreVariables_;
    
    /** Message handler. */
    CoinMessageHandler * bcpsMessageHandler_;

    /** Bcps messages. */
    CoinMessages bcpsMessages_;
    
 public:

    BcpsModel() 
        : 
	numCoreConstraints_(0), 
	numCoreVariables_(0)
        {
	    bcpsMessageHandler_ = new CoinMessageHandler();
	    bcpsMessageHandler_->setLogLevel(2);
	    bcpsMessages_ = BcpsMessage();
	}

        virtual ~BcpsModel() {

        int i = 0, size  = 0;

        size = constraints_.size();
	if (size > 0) {
	    for (i = 0; i < size; ++i) {
		delete constraints_[i]; 
	    }
	}

        size =  variables_.size();
	if (size > 0) {
	    for (i = 0; i < size; ++i) {
		delete variables_[i];
	    }
	}

	delete bcpsMessageHandler_;
    }
    
    /** Get variables and constraints */
    /**@{*/
    std::vector<BcpsConstraint *> & getConstraints() { return constraints_; }
    int getNumCoreConstraints() const { return numCoreConstraints_; }

    std::vector<BcpsVariable *> & getVariables() { return variables_; }
    int getNumCoreVariables() const { return numCoreVariables_; }
    /**@}*/
    
    /** Set variables and constraints */
    /**@{*/
    void setConstraints(BcpsConstraint **con, int size) {
        for (int j = 0; j < size; ++j) {
            constraints_.push_back(con[j]);
        }
    }
    void setNumCoreConstraints(int num) { numCoreConstraints_ = num; }

    void setVariables(BcpsVariable **var, int size) { 
        for (int j = 0; j < size; ++j) {
            variables_.push_back(var[j]);
        }
    }
    void setNumCoreVariables(int num) { numCoreVariables_ = num; }
    /**@}*/

    /** Get the message handler. */
    CoinMessageHandler * bcpsMessageHandler() const 
    { return bcpsMessageHandler_; }
    
    /** Return messages. */
    CoinMessages bcpsMessages() { return bcpsMessages_; }

};


#endif
