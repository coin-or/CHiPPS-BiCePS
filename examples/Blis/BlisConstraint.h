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

#ifndef BlisConstraint_h_
#define BlisConstraint_h_

#include "BcpsObject.h"

//#############################################################################

class BlisConstraint : public BcpsConstraint {

 private:

    int size_;
    int *indices_;
    double *values_;

 public:

    /** Default constructor. */
    BlisConstraint();

    /** Useful constructor. */
    BlisConstraint(int s, const int *ind, const double *val);

    /** Useful constructor. */
    BlisConstraint(double lbh, double ubh, double lbs, double ubs);

    /** Useful constructor. */
    BlisConstraint(double lbh, double ubh, double lbs, double ubs,
                   int s, const int *ind, const double *val);
    /** Destructor. */
    virtual ~BlisConstraint();

    /** Copy constructor. */
    BlisConstraint(const BlisConstraint & rhs);

    /** Return data  */
    /**@{*/
    int getSize() const       { return size_; }
    int* getIndices() const   { return indices_; }
    double* getValues() const { return values_; }
    /**@}*/

    /** Set data  */
    /**@{*/
    void setData(int s, const int *ind, const double *val) {
	if (size_ < s) {
	    delete [] indices_;
	    delete [] values_;
	    indices_ = new int [s];
	    values_ = new double [s];
	}
	size_ = s;
	memcpy(indices_, ind, sizeof(int) * s);
	memcpy(values_, val, sizeof(double) * s);
    }
    /**@}*/

    virtual double infeasibility(BcpsModel * m, int & preferredWay) const;

    ///@name Encode and Decode functions
    //@{
    /// Encode this to an AlpsEncoded object.
    virtual AlpsReturnStatus encode(AlpsEncoded * encoded);
    /// Decode a given AlpsEncoded object to an AlpsKnowledge object and return a
    /// pointer to it.
    virtual AlpsKnowledge * decode(AlpsEncoded & encoded) const;
    // todo(aykut) this should be a pure virtual function in Alps level
    // we can overload this function here due to cv-qualifier.
    /// Decode a given AlpsEncoded object into self.
    AlpsReturnStatus decodeToSelf(AlpsEncoded & encoded);
    //@}

    /** Compute a hash key. */
    virtual void hashing(BcpsModel *model=NULL);
};

//#############################################################################

#endif
