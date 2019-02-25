/*===========================================================================*
 * This file is part of the Bcps Linear Solver (BLIS).                       *
 *                                                                           *
 * BLIS is distributed under the Eclipse Public License as part of the       *
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
 * Copyright (C) 2001-2017, Lehigh University, Yan Xu, and Ted Ralphs.       *
 * All Rights Reserved.                                                      *
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

double BlisConstraint::infeasibility(BcpsModel * m,
                                     int & preferredWay) const {
  std::cerr << "Not implemented, " << std::endl
            << "file: " <<  __FILE__ << std::endl
            << "line: " << __LINE__ << std::endl;
  throw std::exception();
  return 0.0;
}

/// Encode this to an AlpsEncoded object.
AlpsReturnStatus BlisConstraint::encode(AlpsEncoded * encoded) {
  std::cerr << "Not implemented, " << std::endl
            << "file: " <<  __FILE__ << std::endl
            << "line: " << __LINE__ << std::endl;
  throw std::exception();
  return AlpsReturnStatusOk;
}

/// Decode a given AlpsEncoded object to an AlpsKnowledge object and return a
/// pointer to it.
AlpsKnowledge * BlisConstraint::decode(AlpsEncoded & encoded) const {
  std::cerr << "Not implemented, " << std::endl
            << "file: " <<  __FILE__ << std::endl
            << "line: " << __LINE__ << std::endl;
  throw std::exception();
  BlisConstraint * con = new BlisConstraint();
  return con;
}

/// Decode a given AlpsEncoded object into self.
AlpsReturnStatus BlisConstraint::decodeToSelf(AlpsEncoded & encoded) {
  std::cerr << "Not implemented, " << std::endl
            << "file: " <<  __FILE__ << std::endl
            << "line: " << __LINE__ << std::endl;
  throw std::exception();
  return AlpsReturnStatusOk;
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
