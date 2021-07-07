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

    virtual double infeasibility(BcpsModel * bcps_model, int & preferredDir) const {
        bool integral = true;
        if (!integral) {
          return 0.0;
        }
        BlisModel * model = dynamic_cast<BlisModel*>(bcps_model);
        preferredDir = -1;
        // get integer tolerance parameter
        double tolerance = model->BlisPar()->entry(BlisParams::integerTol);
        double value = model->solver()->getColSolution()[getObjectIndex()];
        double dist_to_upper = ceil(value) - value;
        double dist_to_lower = value - floor(value);
        // return the minimum of distance to upper or lower
        double infeas;
        if (dist_to_upper>dist_to_lower) {
          preferredDir = -1;
          infeas = dist_to_lower;
        }
        else {
          preferredDir = 1;
          infeas = dist_to_upper;
        }
        if (infeas<tolerance) {
          infeas = 0.0;
        }
        return infeas;
    }


    virtual BcpsBranchObject * createBranchObject(BcpsModel *m, int way) const {
        BlisModel * model = dynamic_cast<BlisModel*>(m);
        int var_index = getObjectIndex();
        // get current value from solver
        double value = model->solver()->getColSolution()[var_index];
        // we do not know the score of branch object since we do not know
        // the branching strategy.
        double score = 0.0;
        BcpsBranchObject * bo = new BlisBranchObjectInt(var_index, score, value);
        bo->setBroker(broker_);
        return bo;
    }

  ///@name Encode and Decode functions
  //@{
  /// Get encode from #AlpsKnowledge
  using AlpsKnowledge::encode;
  /// Encode this to an AlpsEncoded object.
  virtual AlpsReturnStatus encode(AlpsEncoded * encoded) const {
      std::cerr << "Not implemented, "
                << "file: " <<  __FILE__
                << "line: " << __LINE__
                << std::endl;
      throw std::exception();
      return AlpsReturnStatusOk;
  }
  /// Decode a given AlpsEncoded object to a new BlisVariable object and return
  /// a pointer to it.
  virtual AlpsKnowledge * decode(AlpsEncoded & encoded) const {
      std::cerr << "Not implemented, "
                << "file: " <<  __FILE__
                << "line: " << __LINE__
                << std::endl;
      throw std::exception();
      return NULL;
  }
  /// Decode a given AlpsEncoded object into self.
  virtual AlpsReturnStatus decodeToSelf(AlpsEncoded & encoded) {
      std::cerr << "Not implemented, "
                << "file: " <<  __FILE__
                << "line: " << __LINE__
                << std::endl;
      throw std::exception();
      return AlpsReturnStatusOk;
  }
  //@}
};

//#############################################################################

#endif /* End of file */
