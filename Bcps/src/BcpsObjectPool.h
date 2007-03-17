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

#ifndef BcpsObjectPool_h_
#define BcpsObjectPool_h_

#include <vector>

#include "AlpsKnowledgePool.h"

#include "BcpsObject.h"

//#############################################################################
/** Object pool is used to store objects */
//#############################################################################

class BcpsObjectPool : public AlpsKnowledgePool {

 private:

    std::vector<AlpsKnowledge *> objects_;

 public:

    /** Default construct. */
    BcpsObjectPool() {}
    virtual ~BcpsObjectPool() {
	if (! objects_.empty()) {
	    freeGuts();
	}
    }

    /** Free object pointers. */
    inline void freeGuts() {
	for (int i = objects_.size() - 1; i >= 0; --i) {
	    delete objects_[i];
	}
        objects_.clear();
    }
    
    /** Reset to empty. Don't free memory. */
    inline void clear(){ objects_.clear(); }
    
    /** Add a knowledge to pool */
    virtual void addKnowledge(AlpsKnowledge * nk, double priority) { 
	objects_.push_back(nk);
    };
 
    /** Query how many knowledges are in the pool.*/
    virtual int getNumKnowledges() const {
	return objects_.size();
    }

    /** Query a knowledge, but doesn't remove it from the pool*/
    virtual std::pair<AlpsKnowledge*, double> getKnowledge() const {
	return std::make_pair(objects_[0], 0.0);
    };

    /** Remove the queried knowledge from the pool*/
    //virtual void popKnowledge() {
    //throw CoinError("Can not call popKnowledge()",
    //		"popKnowledge()", "AlpsKnowledgePool");
    //}

    /** Check whether the pool has knowledge. */
    virtual bool hasKnowledge() const
        { return objects_.empty() ? false : true; };
	

    /** Set the quantity limit of knowledges that can be stored in the pool. */
    //virtual void setMaxNumKnowledges(int num) {
    //std::cout << "Can not call setMaxNumKnowledges without overriding"
    //	  << std::endl;
    //throw CoinError("Can not call  setMaxNumKnowledges()",
    //		"setMaxNumKnowledges()", "AlpsKnowledgePool");
    //}

    /** Query the quantity limit of knowledges. */
    //virtual int getMaxNumKnowledges() const {
	// throw CoinError("Can not call getMaxNumKnowledges()",
	//		    "getMaxNumKnowledges()", "AlpsKnowledgePool");
	//return INT_MAX;
    //}
    
    /** Query the best knowledge in the pool.*/
    //virtual std::pair<AlpsKnowledge*, double> 
    //getBestKnowledge() const {
    //throw CoinError("Can not call  getBestKnowledge()",
    //		"getBestKnowledge()", "AlpsKnowledgePool");
    //}
 
    /** Get a reference to all the knowledges in the pool.*/
    // virtual void getAllKnowledges (std::vector<std::pair<AlpsKnowledge*, 
    //			   double> >& kls) const {
    //std::cout << "Can not call  getAllKnowledge() without overriding"
    //	  << std::endl;
    //throw CoinError("Can not call  getAllKnowledge()",
    //		"getAllKnowledge()", "AlpsKnowledgePool");
    //}
    
    std::vector<AlpsKnowledge *> getObjects() const 
	{
	    return objects_;
	}
};


//#############################################################################

class BcpsConstraintPool : public BcpsObjectPool {
 public:
    BcpsConstraintPool() {}
    virtual ~BcpsConstraintPool() {}
};

//#############################################################################

class BcpsVariablePool : public BcpsObjectPool {
 public:
    BcpsVariablePool() {}
    virtual ~BcpsVariablePool() {}
};


#endif // End of file
