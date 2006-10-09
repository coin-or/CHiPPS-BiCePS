/*===========================================================================*
 * This file is part of the Bcps Linear Solver (BLIS).                       *
 *                                                                           *
 * BLIS is distributed under the Common Public License as part of the        *
 * COIN-OR repository (http://www.coin-or.org).                              *
 *                                                                           *
 * Authors: Yan Xu, SAS Institute Inc.                                       *
 *          Ted Ralphs, Lehigh University                                    *
 *          Laszlo Ladanyi, IBM T.J. Watson Research Center                  *
 *          Matthew Saltzman, Clemson University                             *
 *                                                                           * 
 *                                                                           *
 * Copyright (C) 2001-2005, International Business Machines                  *
 * Corporation, Lehigh University, Yan Xu, Ted Ralphs, Matthew Salzman and   *
 * others. All Rights Reserved.                                              *
 *===========================================================================*/

#ifndef BlisVariable_h_
#define BlisVariable_h_

#include "BcpsObject.h"

//#############################################################################

class BlisVariable : public BcpsVariable {

 private:

    double objCoef_;
    int size_;
    int *indices_;
    double *values_;

 public:

    BlisVariable() : objCoef_(0.0), size_(0), indices_(NULL), values_(NULL) {}
    
    BlisVariable(double obj, int s, const int *ind, const double *val)
	{
            objCoef_ = obj;
	    size_ = s;
	    indices_ = new int [s];
	    values_ = new double [s];
	    memcpy(indices_, ind, s * sizeof(int));
	    memcpy(values_, val, s * sizeof(double));
	}
    
    BlisVariable(double lbh, double ubh, double lbs, double ubs) 
	:
	BcpsVariable(lbh, ubh, lbs, ubs),
        objCoef_(0.0),
	size_(0), indices_(NULL), values_(NULL)
	{}

    BlisVariable(double lbh, double ubh, double lbs, double ubs,
                 double obj, int s, const int *ind, const double *val)
        : 
        BcpsVariable(lbh, ubh, lbs, ubs)
        {
            objCoef_ = obj;
            size_ = s;
            indices_ = new int [s];
            values_ = new double [s];
            memcpy(indices_, ind, s * sizeof(int));
            memcpy(values_, val, s * sizeof(double));
        }
    
    virtual ~BlisVariable(){ 
	if (size_ > 0) {
	    delete [] indices_; indices_ = NULL;
	    delete [] values_; values_ = NULL;
	}
    }
    
    /** Return data  */
    /**@{*/
    int getSize() const     { return size_; }
    int* getIndices() const { return indices_; }
    double* getValues()     { return values_; }    
    /**@}*/
    
    /** Set data  */
    /**@{*/
    void setData(int s, const int *ind, const double *val) {
	if (size_ < s) {
	    delete [] indices_; indices_ = NULL;
	    delete [] values_; values_ = NULL;
	    indices_ = new int [s];
	    values_ = new double [s];
	}
	size_ = s;
	memcpy(indices_, ind, sizeof(int) * s);
	memcpy(values_, val, sizeof(double) * s);
    }
    /**@}*/

 protected:

   /** Pack Blis part into an encoded object. */
    AlpsReturnCode encodeBlis(AlpsEncoded *encoded) {
	AlpsReturnCode status = ALPS_OK;

	encoded->writeRep(objCoef_);
	encoded->writeRep(indices_, size_);
	encoded->writeRep(values_, size_);

	return status;
    }    

    /** Unpack Blis part from a encode object. */
    AlpsReturnCode decodeBlis(AlpsEncoded &encoded) {
	AlpsReturnCode status = ALPS_OK;

	encoded.readRep(objCoef_);
	encoded.readRep(indices_, size_);
	encoded.readRep(values_, size_);
	
	return status;
    }
    
 public:
    
    /** Pack to a encode object. */
    virtual AlpsReturnCode encode(AlpsEncoded *encoded){
	AlpsReturnCode status;

	status = encodeBcps(encoded);
	status = encodeBlis(encoded);
	
	return status;
    }

    /** Decode a variable from an encoded object. */
    virtual AlpsKnowledge* decode(AlpsEncoded &encoded) const {
	AlpsReturnCode status = ALPS_OK;
	BlisVariable * var = new BlisVariable();    
	
	// Unpack Bcps part.
	status = var->decodeBcps(encoded);
	if (status) {
	    throw CoinError("Failed to decode Bcps part of var",
			    "decode",
			    "BlisObject");
	}
	
	// Unpack Blis part.
	status = var->decodeBlis(encoded);
	if (status) {
	    throw CoinError("Failed to decode Blis part of var", 
			    "decode", 
			    "BlisObject");
	}
	return var;
    }
    
};

//#############################################################################

#endif /* End of file */

