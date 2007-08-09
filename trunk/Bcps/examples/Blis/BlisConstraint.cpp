/*===========================================================================*
 * This file is part of the Bcps Linear Solver (BLIS).                       *
 *                                                                           *
 * BLIS is distributed under the Common Public License as part of the        *
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
 *                                                                           *
 * Copyright (C) 2001-2007, International Business Machines                  *
 * Corporation, Lehigh University, Yan Xu, Ted Ralphs, Matthew Salzman and   *
 * others. All Rights Reserved.                                              *
 *===========================================================================*/

#include "BlisConstraint.h"
#include "BlisModel.h"

//#############################################################################

/* Default constructor. */
BlisConstraint::BlisConstraint() 
    :size_(0), indices_(NULL), values_(NULL) {}

//#############################################################################

/* Useful constructor. */
BlisConstraint::BlisConstraint(int s, const int *ind, const double *val)
{
    size_ = s;
    indices_ = new int [s];
    values_ = new double [s];
    memcpy(indices_, ind, s * sizeof(int));
    memcpy(values_, val, s * sizeof(double));
}

//#############################################################################

/* Useful constructor. */
BlisConstraint::BlisConstraint(double lbh, double ubh, double lbs, double ubs) 
    :
    BcpsConstraint(lbh, ubh, lbs, ubs),
    size_(0), indices_(NULL), values_(NULL) {}

//#############################################################################

/* Useful constructor. */
BlisConstraint::BlisConstraint(double lbh, double ubh, double lbs, double ubs,
			       int s, const int *ind, const double *val)
    : 
    BcpsConstraint(lbh, ubh, lbs, ubs)
{
    size_ = s;
    indices_ = new int [s];
    values_ = new double [s];
    memcpy(indices_, ind, s * sizeof(int));
    memcpy(values_, val, s * sizeof(double));
}

//#############################################################################

/** Destructor. */
BlisConstraint::~BlisConstraint()
{ 
    delete [] indices_; indices_ = NULL;
    delete [] values_; values_ = NULL;
}

//#############################################################################

/** Copy constructor. */
BlisConstraint::BlisConstraint(const BlisConstraint & rhs) 
    : BcpsConstraint(rhs) 
{
    size_ = rhs.size_;
    
    if (size_ <= 0) {
	std::cout << "ERROR: size_ = " << size_ << std::endl;
	assert(size_);
    }
    if (size_ > 0) {
	indices_ = new int [size_];
	values_ = new double [size_];
	memcpy(indices_, rhs.indices_, size_ * sizeof(int));
	memcpy(values_, rhs.values_, size_ * sizeof(double));
    }
    else {
	indices_ = NULL;
	values_ = NULL;
    }
}

//#############################################################################
   
/** Pack Blis part into an encoded object. */
AlpsReturnStatus 
BlisConstraint::encodeBlis(AlpsEncoded *encoded) 
{
    AlpsReturnStatus status = AlpsReturnStatusOk;
    if (size_ <= 0) {
	std::cout << "ERROR: encodeBlis: size_=" << size_<<std::endl;
	assert(size_ > 0);
    }
    std::cout << "encodeBlis: size_=" << size_<<std::endl;
    encoded->writeRep(indices_, size_);
    encoded->writeRep(values_, size_);
    return status;
}    

//#############################################################################

/** Unpack Blis part from a encode object. */
AlpsReturnStatus 
BlisConstraint::decodeBlis(AlpsEncoded &encoded) 
{
    AlpsReturnStatus status = AlpsReturnStatusOk;
    encoded.readRep(indices_, size_);
    if (size_ <= 0) {
	std::cout << "ERROR: decodeBlis: con1, size_=" << size_<<std::endl;
	assert(size_ > 0);
    }
    encoded.readRep(values_, size_);
    if (size_ <= 0) {
	std::cout << "ERROR: decodeBlis: con2, size_=" << size_<<std::endl;
	assert(size_ > 0);
    }
    return status;
}

//#############################################################################

AlpsReturnStatus 
BlisConstraint::encode(AlpsEncoded *encoded) 
{
    AlpsReturnStatus status = AlpsReturnStatusOk;
    status = encodeBcpsObject(encoded);
    status = encodeBlis(encoded);
    return status;
}

//#############################################################################

/** Decode a constraint from an encoded object. */
AlpsKnowledge* 
BlisConstraint::decode(AlpsEncoded& encoded) const 
{
    AlpsReturnStatus status = AlpsReturnStatusOk;
    BlisConstraint* con = new BlisConstraint();    
    
    // Unpack Bcps object part.
    status = con->decodeBcpsObject(encoded);
    if (status) {
	throw CoinError("Failed to decode Bcps part",
			"decode", 
			"BlisObject");
    }
    
    // Unpack Blis part.
    status = con->decodeBlis(encoded);
    if (status) {
	throw CoinError("Failed to decode Blis part", 
			"decode", 
			"BlisObject");
    }
    
    return con;
}

//#############################################################################

/** Compute hash value. */
void 
BlisConstraint::hashing(BcpsModel *model)
{
    assert(model != NULL);
    BlisModel *m = dynamic_cast<BlisModel *>(model);
    
    int k, ind;
    const double * randoms = m->getConRandoms();

    hashValue_ = 0.0;    
    for (k = 0; k < size_; ++k) {
	ind = indices_[k];
	hashValue_ += randoms[ind] * ind;
    }
#ifdef BLIS_DEBUG_MORE
    std::cout << "hashValue_=" << hashValue_ << std::endl;
#endif
}

//#############################################################################
//#############################################################################



