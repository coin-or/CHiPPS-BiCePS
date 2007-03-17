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
 * Copyright (C) 2001-2006, Lehigh University, Yan Xu, and Ted Ralphs.       *
 * All Rights Reserved.                                                      *
 *===========================================================================*/

#ifndef BcpsModel_h_
#define BcpsModel_h_

//#############################################################################

#include "AlpsKnowledgeBroker.h"

#include "CoinMessageHandler.hpp"

#include "AlpsModel.h"

#include "BcpsObject.h"
#include "BcpsMessage.h"

//#############################################################################
//#############################################################################

class BcpsModel : public AlpsModel {

 protected:

    int numCoreConstraints_;

    BcpsConstraint **coreConstraints_;

    int numCoreVariables_;

    BcpsVariable **coreVariables_;

    /** Message handler. */
    CoinMessageHandler * bcpsMessageHandler_;

    /** Bcps messages. */
    CoinMessages bcpsMessages_;
    
 public:

    BcpsModel() 
        : 
	numCoreConstraints_(0), 
        coreConstraints_(NULL),
	numCoreVariables_(0), 
        coreVariables_(NULL)
        {
	    bcpsMessageHandler_ = new CoinMessageHandler();
	    bcpsMessageHandler_->setLogLevel(2);
	    bcpsMessages_ = BcpsMessage();
	}

    virtual ~BcpsModel() {
	int i = 0;
	if (numCoreConstraints_ > 0) {
	    for (i = 0; i < numCoreConstraints_; ++i) {
		delete coreConstraints_[i]; 
		coreConstraints_[i] = NULL;
	    }
	    delete [] coreConstraints_;
	    coreConstraints_ = NULL;
	}
	if (numCoreVariables_ > 0) {
	    for (i = 0; i < numCoreVariables_; ++i) {
		delete coreVariables_[i];
		coreVariables_[i] = NULL;
	    }
	    delete [] coreVariables_;
	    coreVariables_ = NULL;
	}
	delete bcpsMessageHandler_;
    }
    
    /** Get core variables and constraints */
    /**@{*/
    int getNumCoreConstraints() const { return numCoreConstraints_; }
    BcpsConstraint ** getCoreConstraints() { return coreConstraints_; }
    int getNumCoreVariables() const { return numCoreVariables_; }
    BcpsVariable ** getCoreVariables() { return coreVariables_; }
    /**@}*/
    
    /** Set core variables and constraints */
    /**@{*/
    void setNumCoreConstraints(int num) { numCoreConstraints_ = num; }
    void setNumCoreVariables(int num) { numCoreVariables_ = num; }
    void setCoreConstraints(BcpsConstraint **con) { coreConstraints_ = con; }
    void setCoreVariables(BcpsVariable **var) { coreVariables_ = var; }
    /**@}*/

    /** Get the message handler. */
    CoinMessageHandler * bcpsMessageHandler() const 
    { return bcpsMessageHandler_; }
    
    /** Return messages. */
    CoinMessages bcpsMessages() { return bcpsMessages_; }

};


#endif
