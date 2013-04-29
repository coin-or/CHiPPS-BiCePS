/*===========================================================================*
 * This file is part of the Bcps Linear Solver (BLIS).                       *
 *                                                                           *
 * ALPS is distributed under the Eclipse Public License as part of the       *
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
 * Copyright (C) 2001-2013, Lehigh University, Yan Xu, and Ted Ralphs.       *
 * All Rights Reserved.                                                      *
 *===========================================================================*/

#include "float.h"

#include "CoinFinite.hpp"
#include "CoinTime.hpp"
#include "OsiClpSolverInterface.hpp"

#include "CglGomory.hpp"
#include "CglProbing.hpp"
#include "CglKnapsackCover.hpp"
#include "CglOddHole.hpp"
#include "CglClique.hpp"
#include "CglFlowCover.hpp"
#include "CglMixedIntegerRounding2.hpp"
#include "CglTwomir.hpp"

#include "BcpsObject.h"

#include "BlisBranchObjectInt.h"
#include "BlisBranchStrategyPseudo.h"
#include "BlisBranchStrategyRel.h"
#include "BlisBranchStrategyStrong.h"

#include "BlisConstraint.h"
#include "BlisHeurRound.h"
#include "BlisModel.h"
#include "BlisObjectInt.h"
#include "BlisSolution.h"
#include "BlisTreeNode.h"
#include "BlisVariable.h"

//#############################################################################

void 
BlisModel::init() 
{
    lpSolver_ = NULL;
    numCols_ = 0;
    numRows_ = 0;
    numElems_ = 0;
    colMatrix_ = 0;

    origVarLB_ = NULL;
    origVarUB_ = NULL;
    origConLB_ = NULL;
    origConUB_ = NULL;

    startVarLB_ = NULL;
    startVarUB_ = NULL;
    startConLB_ = NULL;
    startConUB_ = NULL;

    tempVarLBPos_ = NULL;
    tempVarUBPos_ = NULL;
    tempConLBPos_ = NULL;
    tempConUBPos_ = NULL;
    
    objSense_ = 1.0;
    objCoef_ = NULL;

    colType_ = 0;
    numIntVars_ = 0;
    intVars_ = NULL;
    
    numSolutions_ = 0;
    numHeurSolutions_ = 0;
    
    //savedLpSolution_ = NULL;
    incumbent_ = NULL;
    hotstartStrategy_ = 0;
    
    activeNode_ = NULL;
    numStrong_ = 0;
    
    cutoff_ = COIN_DBL_MAX;
    incObjValue_ = COIN_DBL_MAX;
    handler_ = NULL;
    objects_ = NULL;
    numObjects_ = 0;
    
    numNodes_ = 0;
    numIterations_ = 0;
    aveIterations_ = 0;
    
    branchStrategy_ = NULL;
    priority_ = NULL;
    nodeWeight_ = 1.0;

    isRoot_ = true;
    BlisPar_ = new BlisParams;

    /// Timing
    startTime_ = 0.0;
    /// Max solution time allowed
    timeLimit_ = 1.0e75;

    /// Tolerance 
    integerTol_ = 1.0e-5;
    optimalRelGap_ = 1.0e-4;
    optimalAbsGap_ = 1.0e-6;

    /// Heuristic
    useHeuristics_ = true;
    numHeuristics_ = 0;
    heuristics_ = NULL;

    /// Cons related
    useCons_ = 0;
    numCutGenerators_ = 0;
    generators_ = NULL;
    constraintPool_ = NULL;
    oldConstraints_ = NULL;
    oldConstraintsSize_ = 0;
    numOldConstraints_ = 0;
    conRandoms_ = NULL;
}

//#############################################################################

// Read from file (currently MPS format. TODO: LP format).
// Load to solver.
void
BlisModel::readInstance(const char* dataFile)
{
    int j;
    
    //------------------------------------------------------
    // Read in data from MPS file.
    //------------------------------------------------------
    
    CoinMpsIO *mps = new CoinMpsIO;
    int rc = mps->readMps(dataFile, "");
    if(rc) {
        delete mps;
        throw CoinError("Unable to read in instance",
                        "readInstance",
                        "BlisModel");
    }
    

    //------------------------------------------------------
    // Get problem data.
    //------------------------------------------------------

    numCols_ = mps->getNumCols();
    numRows_ = mps->getNumRows();
    numElems_ = mps->getNumElements();

    colMatrix_ = new CoinPackedMatrix();    
    *colMatrix_ = *(mps->getMatrixByCol());
    
    origVarLB_ = new double [numCols_];
    origVarUB_ = new double [numCols_];

    origConLB_ = new double [numRows_];
    origConUB_ = new double [numRows_];
    
    memcpy(origVarLB_, mps->getColLower(), sizeof(double) * numCols_);
    memcpy(origVarUB_, mps->getColUpper(), sizeof(double) * numCols_);
    
    memcpy(origConLB_, mps->getRowLower(), sizeof(double) * numRows_);
    memcpy(origConUB_, mps->getRowUpper(), sizeof(double) * numRows_);
    
    //memcpy(startVarLB_, mps->getColLower(), sizeof(double) * numCols_);
    //memcpy(startVarUB_, mps->getColUpper(), sizeof(double) * numCols_);
    
    //memcpy(startConLB_, mps->getRowLower(), sizeof(double) * numRows_);
    //memcpy(startConUB_, mps->getRowUpper(), sizeof(double) * numRows_);
    
    objSense_ = 1.0; /* Default from MPS is minimization */
    
    objCoef_ = new double [numCols_];
    memcpy(objCoef_, mps->getObjCoefficients(), sizeof(double) * numCols_);
    
    //------------------------------------------------------
    // Classify variable type.
    //------------------------------------------------------

    colType_ = new char [numCols_];
    intVars_ = new int [numCols_];
    numIntVars_ = 0;
    for(j = 0; j < numCols_; ++j) {
	if (mps->isContinuous(j)) {
	    colType_[j] = 'C';
	}
	else {
	    intVars_[numIntVars_++] = j;
	    if (origVarLB_[j] == 0 && origVarUB_[j] == 1.0) {
		colType_[j] = 'B';
	    }
	    else {
		colType_[j] = 'I';
	    }
	}
    }
    
    //------------------------------------------------------
    // Do root preprocessing.
    //------------------------------------------------------
    

    //------------------------------------------------------
    // load problem to lp solver.
    // NOTE: since no presolve, load original one.
    //------------------------------------------------------
    
    if (!lpSolver_) {
        lpSolver_ = new OsiClpSolverInterface();
    }

    lpSolver_->loadProblem(*colMatrix_,
			   origVarLB_, origVarUB_,   
			   objCoef_,
			   origConLB_, origConUB_);

    lpSolver_->setObjSense(objSense_);
    lpSolver_->setInteger(intVars_, numIntVars_);
    
    delete mps;

    if (numIntVars_ == 0) {
	// solve lp and throw error.
	lpSolver_->initialSolve();
	throw CoinError("Input instance is a LP", 
			"readInstance", "BlisModel");
    }
}

//############################################################################ 

/** Read in Alps parameters. */
void 
BlisModel::readParameters(const int argnum, const char * const * arglist)
{
    std::cout << "Reading in ALPS parameters ..." << std::endl;
    AlpsPar_->readFromArglist(argnum, arglist);
    std::cout << "Reading in BLIS parameters ..." << std::endl;
    BlisPar_->readFromArglist(argnum, arglist);
} 

//##############################################################################

/** Write out parameters. */
void 
BlisModel::writeParameters(std::ostream& outstream) const
{
    outstream << "\n================================================"
              <<std::endl;
    outstream << "ALPS Parameters: " << std::endl;
    AlpsPar_->writeToStream(outstream);
    outstream << "\n================================================"
              <<std::endl;
    outstream << "BLIS Parameters: " << std::endl;
    BlisPar_->writeToStream(outstream);
}

//############################################################################ 

/** Do necessary work to make model usable. Return success or not. 
    This function is called when constructing knowledge broker. */
bool 
BlisModel::setupSelf()
{
    int j;

    //------------------------------------------------------
    // Starting time.
    //------------------------------------------------------

    startTime_ = CoinCpuTime();

    //------------------------------------------------------
    // Allocate memory.
    //------------------------------------------------------

    startVarLB_ = new double [numCols_];
    startVarUB_ = new double [numCols_];

    startConLB_ = new double [numRows_];
    startConUB_ = new double [numRows_];

    tempVarLBPos_ = new int [numCols_];
    tempVarUBPos_ = new int [numCols_];

    tempConLBPos_ = new int [numRows_];
    tempConUBPos_ = new int [numRows_];

    //------------------------------------------------------
    // Get parameters.
    //------------------------------------------------------
    
    timeLimit_ = AlpsPar_->entry(AlpsParams::timeLimit);
    
    integerTol_ = BlisPar_->entry(BlisParams::integerTol);
    optimalRelGap_ = BlisPar_->entry(BlisParams::optimalRelGap);
    optimalAbsGap_ = BlisPar_->entry(BlisParams::optimalAbsGap);

    int relibility = BlisPar_->entry(BlisParams::pseudoRelibility);

    //------------------------------------------------------
    // Create core variables and constraints.
    //------------------------------------------------------

#ifdef BLIS_DEBUG
    std::cout << "setupSelf: numCols_ " << numCols_ 
	      << ", numRows_" << numRows_ 
	      << std::endl;
    std::cout << "Create core ..." << std::endl;
#endif

    numCoreVariables_ = numCols_;
    numCoreConstraints_ = numRows_;
    
    // BlisVariable ** tempVars = new BlisVariable* [numCols_];
    // BlisConstraint ** tempCons = new BlisConstraint* [numRows_];

    //variables_ = new BcpsVariable* [numCols_];
    //constraints_ = new BcpsConstraint* [numRows_];
    
    for (j = 0; j < numCols_; ++j) {
	BlisVariable * var = new BlisVariable(origVarLB_[j],
					      origVarUB_[j], 
					      origVarLB_[j], 
					      origVarUB_[j]);
	variables_.push_back(var);
	var = NULL;
	variables_[j]->setObjectIndex(j);
	variables_[j]->setRepType(BCPS_CORE);
	variables_[j]->setIntType(colType_[j]);
	variables_[j]->setStatus(BCPS_NONREMOVALBE);
    }

    for (j = 0; j < numRows_; ++j) {
        BlisConstraint *con = new BlisConstraint(origConLB_[j], 
                                                 origConUB_[j], 
                                                 origConLB_[j], 
                                                 origConUB_[j]);
        constraints_.push_back(con);
        con = NULL;
        constraints_[j]->setObjectIndex(j);
        constraints_[j]->setRepType(BCPS_CORE);
        //coreContraints_[j]->setIntType(colType_[j]);
        constraints_[j]->setStatus(BCPS_NONREMOVALBE);
    }
    
    //------------------------------------------------------
    // Identify integers.
    //------------------------------------------------------
    
    findIntegers(true);
    
    // lpSolver_->initialSolve();
    
#ifdef BLIS_DEBUG_MORE
    std::string problemName;
    lpSolver_->getStrParam(OsiProbName, problemName);
    printf("BLIS: setupSelf: Problem name - %s\n", problemName.c_str());
    lpSolver_->setHintParam(OsiDoReducePrint, false, OsiHintDo, 0);
#endif

    //------------------------------------------------------
    // Set branch strategy.
    //------------------------------------------------------

    int brStrategy = BlisPar_->entry(BlisParams::branchStrategy);
    if (brStrategy == 0) {
        // Max inf
        branchStrategy_ =  new BlisBranchStrategyStrong(this);
    }
    else if (brStrategy == 1) {
        // Pseudocost
        branchStrategy_ =  new BlisBranchStrategyPseudo(this, 1);
    }
    else if (brStrategy == 2) {
        // Relibility
        branchStrategy_ =  new BlisBranchStrategyRel(this, relibility);
    }
    else if (brStrategy == 3) {
        // Strong
        branchStrategy_ =  new BlisBranchStrategyStrong(this);
    }
    else {
        throw CoinError("Unknown branch strategy.", "setupSelf","BlisModel");
    }

    //------------------------------------------------------
    // Add heuristics.
    //------------------------------------------------------

    useHeuristics_ = BlisPar_->entry(BlisParams::useHeuristics);
    int useRound = BlisPar_->entry(BlisParams::heurRound); 

    if (useHeuristics_) {
        if (useRound > -2) {
            // Add rounding heuristic
            BlisHeurRound *heurRound = new BlisHeurRound(this, 
                                                         "Rounding",
                                                         useRound);
            addHeuristic(heurRound);
        }
    }

    //------------------------------------------------------
    // Cut generators settings.
    //------------------------------------------------------

    // Compute dense cutoff.
    
    const CoinPackedMatrix * rowMatrix = lpSolver_->getMatrixByRow();
    const int * rowLen = rowMatrix->getVectorLengths();
    double maxLen = 0.0, minLen = ALPS_DBL_MAX, sumLen = 0.0;
    double aveLen, diffLen, stdLen;
    double denseConFactor = BlisPar_->entry(BlisParams::denseConFactor);

    for (j = 0; j < numRows_; ++j) {
	if (rowLen[j] > maxLen) maxLen = rowLen[j];
	if (rowLen[j] < minLen) minLen = rowLen[j];
	sumLen += rowLen[j];
    }
    assert(numRows_ > 0);
    aveLen = sumLen / numRows_;
    sumLen = 0.0;

    for (j = 0; j < numRows_; ++j) {
	diffLen = rowLen[j] - aveLen;
	diffLen *= diffLen;
	sumLen += diffLen;
    }
    stdLen = sqrt(sumLen/numRows_);
    denseConCutoff_ = static_cast<int>(aveLen + denseConFactor*stdLen);    
    denseConCutoff_ = ALPS_MIN(numCols_/2, denseConCutoff_);
    denseConCutoff_ = ALPS_MAX(10, denseConCutoff_);
    
#ifdef BLIS_DEBUG
    std::cout << "aveLen=" << aveLen << ", minLen=" << minLen
	      << ", maxLen=" << maxLen << ", stdLen=" << stdLen
	      << ", denseConCutoff_=" << denseConCutoff_ << std::endl;
#endif
    
    // NOTE: maxNumCons is valid only for automatic strategy.
    double cutFactor = BlisPar_->entry(BlisParams::cutFactor);

    maxNumCons_ = (int)((cutFactor - 1.0) * numCoreConstraints_);
    
    constraintPool_ = new BcpsConstraintPool();
    oldConstraints_ = new BlisConstraint* [maxNumCons_];
    oldConstraintsSize_ = maxNumCons_;
    
    useCons_ = BlisPar_->entry(BlisParams::useCons); 

#ifdef BLIS_DEBUG
    std::cout << "useCons_ = " << useCons_ << std::endl;
#endif

    int clique = BlisPar_->entry(BlisParams::cutClique);
    int gomory = BlisPar_->entry(BlisParams::cutGomory); 
    int fCover = BlisPar_->entry(BlisParams::cutFlowCover);
    int knap = BlisPar_->entry(BlisParams::cutKnapsack); 
    int mir = BlisPar_->entry(BlisParams::cutMir); 
    int oddHole = BlisPar_->entry(BlisParams::cutOddHole);
    int probe = BlisPar_->entry(BlisParams::cutProbing);
    int twoMir = BlisPar_->entry(BlisParams::cutTwoMir); 

    //------------------------------------------------------
    // Add cut generators.
    //------------------------------------------------------

    if (!useCons_) {
	useCons_ = -2;
    }
    else {
	// Reset to -2, so that addCutGenerator can adjust it properly.
	useCons_ = -2;
        if (probe > -2) {
            CglProbing *probing = new CglProbing;
            probing->setUsingObjective(true);
            probing->setMaxPass(1);
            probing->setMaxPassRoot(5);
            // Number of unsatisfied variables to look at
            probing->setMaxProbe(10);
            probing->setMaxProbeRoot(1000);
            // How far to follow the consequences
            probing->setMaxLook(50);
            probing->setMaxLookRoot(500);
            // Only look at rows with fewer than this number of elements
            probing->setMaxElements(200);
            probing->setRowCuts(3);
            addCutGenerator(probing, "Probing", probe);
        }
        if (clique > -2) {
            CglClique *cliqueCut = new CglClique ;
            cliqueCut->setStarCliqueReport(false);
            cliqueCut->setRowCliqueReport(false);
            addCutGenerator(cliqueCut, "Clique", clique);
        }
        if (oddHole > -2) {
            CglOddHole *oldHoleCut = new CglOddHole;
            oldHoleCut->setMinimumViolation(0.005);
            oldHoleCut->setMinimumViolationPer(0.00002);
            // try larger limit
            oldHoleCut->setMaximumEntries(200);
            addCutGenerator(oldHoleCut, "OddHole", oddHole);
        }
        if (fCover > -2) {
            CglFlowCover *flowGen = new CglFlowCover;
            addCutGenerator(flowGen, "FlowCover", fCover);
        }
        if (knap > -2) {
            CglKnapsackCover *knapCut = new CglKnapsackCover;
            addCutGenerator(knapCut, "Knapsack", knap);
        }
        if (mir > -2) {
            CglMixedIntegerRounding2 *mixedGen = new CglMixedIntegerRounding2;
            addCutGenerator(mixedGen, "MIR", mir);
        }
        if (gomory > -2) {
            CglGomory *gomoryCut = new CglGomory;
            // try larger limit
            gomoryCut->setLimit(300);
            addCutGenerator(gomoryCut, "Gomory", gomory);
        }
        if (twoMir > -2) {
            CglTwomir *twoMirCut =  new CglTwomir;
            addCutGenerator(twoMirCut, "Two mir", twoMir);
        }


        // Random vector
        // TODO, if generating variable, then need more space.
        conRandoms_ = new double [numCols_];
        double deseed = 12345678.0;
        double DSEED2 = 2147483647.0;
        
        for (j = 0; j < numCols_; ++j) {
            deseed *= 16807.0;
            int jseed = (int) (deseed / DSEED2);
            deseed -= (double) jseed * DSEED2;
            double random = deseed / DSEED2;
            conRandoms_[j] = random;
	    //std::cout << "conRandoms_[" << j << "]="
	    //      <<conRandoms_[j]<< std::endl;
        }
    
    }

#ifdef BLIS_DEBUG_MORE
    std::cout << "AFTER: useCons_ = " << useCons_ << std::endl;
#endif

    return true;
}

//############################################################################ 

bool 
BlisModel::feasibleSolution(int & numIntegerInfs)
{
    bool feasible = true;
    numIntegerInfs = 0;
    int i = -1;
    //const int numCols = lpSolver_->getNumCols();
    const double *savedLpSolution = lpSolver_->getColSolution();

#if 0
    if (savedLpSolution_ != 0) {
      delete [] savedLpSolution_;
      savedLpSolution_ = 0;
    }
    
    savedLpSolution_ = new double [numCols];
    memcpy(savedLpSolution_, 
	   lpSolver_->getColSolution(), 
	   sizeof(double) * numCols);
#endif
    
    for (i = 0; i < numIntVars_; ++i) {
      if ( ! checkInteger(savedLpSolution[intVars_[i]]) ) {
	++numIntegerInfs;
	feasible = false;
      }
    }

    return feasible;
}

//############################################################################ 

bool
BlisModel::setBestSolution(BLIS_SOL_TYPE how,
			   double & objectiveValue, 
			   const double * solution, 
			   bool fixVariables)
{
    double cutoff = getCutoff();
    
    // Double check the solution to catch pretenders.
    if (objectiveValue >= cutoff) {  // Bad news
        return false;
    }
    else {  // Better solution
	incObjValue_ = objectiveValue;

	int numColumns = lpSolver_->getNumCols();
	if (incumbent_ == 0) {
	    incumbent_ = new double[numColumns];
	}
	
	memcpy(incumbent_, solution, numColumns*sizeof(double));

        // Update cutoff value in lp solver.
	setCutoff(incObjValue_);
        ++numSolutions_;

	switch (how) {
	case BLIS_SOL_BOUNDING:
#ifdef BLIS_DEBUG
	  std::cout << "Rounding heuristics found a better solution" 
		    <<", old cutoff = " << cutoff 
		    << ", new cutoff = " << getCutoff()  << std::endl;
#endif
	  break;
	case BLIS_SOL_BRANCHING:
#ifdef BLIS_DEBUG
	  std::cout << "Branching found a better solution" 
		    <<", old cutoff = " << cutoff 
		    << ", new cutoff = " << getCutoff()  << std::endl;
#endif
	  break;
	case BLIS_SOL_DIVING:
            ++numHeurSolutions_;
#ifdef BLIS_DEBUG
	  std::cout << "Branching found a better solution" 
		    <<", old cutoff = " << cutoff 
		    << ", new cutoff = " << getCutoff()  << std::endl;
#endif
	  break;
	case BLIS_SOL_ROUNDING:
            ++numHeurSolutions_;
#ifdef BLIS_DEBUG
	  std::cout << "Rounding heuristics found a better solution" 
		    <<", old cutoff = " << cutoff 
		    << ", new cutoff = " << getCutoff()  << std::endl;
#endif
	  break;
	case BLIS_SOL_STRONG:
#ifdef BLIS_DEBUG
	  std::cout << "Strong branching found a better solution" 
		    <<", old cutoff = " << cutoff 
		    << ", new cutoff = " << getCutoff()  << std::endl;
#endif
	  break;
	default:
#ifdef BLIS_DEBUG
	  std::cout << "Nowhere found a better solution" 
		    <<", old cutoff = " << cutoff 
		    << ", new cutoff = " << getCutoff()  << std::endl;
#endif
	  break;
	}


	//setMinimizationObjValue(objectiveValue * lpSolver_->getObjSense());
	
	return true;
    }
}

//############################################################################ 

void 
BlisModel::findIntegers(bool startAgain)
{
#ifdef BLIS_DEBUG
    assert(lpSolver_);
#endif

    if (numIntVars_ && !startAgain && objects_) return;

    int iCol;    
    int numCols = getNumCols();

    const double *colLB = lpSolver_->getColLower();
    const double *colUB = lpSolver_->getColUpper();
    BlisObjectInt *intObject = NULL;

    if (intVars_) {
        delete [] intVars_;
    }
    numIntVars_ = 0;

    for (iCol = 0; iCol < numCols; ++iCol) {
	if (lpSolver_->isInteger(iCol)) ++numIntVars_;
    }

    double weight = BlisPar_->entry(BlisParams::pseudoWeight);
    
    int numObjects = 0;
    int iObject;
    BcpsObject ** oldObject = objects_;    

    for (iObject = 0; iObject < numObjects_; ++iObject) {
	BlisObjectInt * obj =
	    dynamic_cast <BlisObjectInt *>(oldObject[iObject]) ;
        
	if (obj) {
	    delete oldObject[iObject];
        }
	else {
	    oldObject[numObjects++] = oldObject[iObject];
        }
    }

    objects_ = new BcpsObject * [(numIntVars_ + numObjects)];
    intVars_ = new int [numIntVars_];
    numObjects_ = numIntVars_ + numObjects;

    // Walk the variables again, filling in the indices and creating objects 
    // for the integer variables. Initially, the objects hold the indices,
    // variable bounds and pseudocost.
    numIntVars_ = 0;
    for (iCol = 0; iCol < numCols; ++iCol) {
	if(lpSolver_->isInteger(iCol)) {
            
	    intObject = new BlisObjectInt(numIntVars_,
                                          iCol,
                                          colLB[iCol],
                                          colUB[iCol]);
            intObject->pseudocost().setWeight(weight);
            objects_[numIntVars_] = intObject;
	    intVars_[numIntVars_++] = iCol;
	}
    }
    
    // Now append other objects
    memcpy(objects_ + numIntVars_, oldObject, numObjects*sizeof(BcpsObject *));

    // Delete old array (just array)
    delete [] oldObject;
}

//#############################################################################

bool
BlisModel::resolve()
{
    lpSolver_->resolve();
    numIterations_ += lpSolver_->getIterationCount();
    bool feasible = (lpSolver_->isProvenOptimal() &&
                     !lpSolver_->isDualObjectiveLimitReached());
    
    return feasible; 
}

//#############################################################################

// Delete all object information
void 
BlisModel::deleteObjects()
{
    delete [] priority_;
    priority_ = NULL;
    int i;
    for (i = 0; i < numObjects_; ++i) delete objects_[i];
    delete [] objects_;
    objects_ = NULL;
    numObjects_ = 0;
    findIntegers(true);
}

//#############################################################################

BlisModel::~BlisModel()
{
    gutsOfDestructor();
}

//#############################################################################

void 
BlisModel::gutsOfDestructor()
{
    int i;

//    delete [] savedLpSolution_;
//    savedLpSolution_ = NULL;

    delete [] intVars_;
    intVars_ = NULL;

    for (i = 0; i < numObjects_; ++i) delete objects_[i];
    delete [] objects_;
    objects_ = NULL;

    delete [] priority_;
    priority_ = NULL;

    delete [] colType_;
    delete colMatrix_;

    delete [] origVarLB_;
    delete [] origVarUB_;

    delete [] origConLB_;
    delete [] origConUB_;

    delete [] startVarLB_;
    delete [] startVarUB_;

    delete [] startConLB_;
    delete [] startConUB_;

    delete [] tempVarLBPos_;
    delete [] tempVarUBPos_;

    delete [] tempConLBPos_;
    delete [] tempConUBPos_;

    delete [] objCoef_;
    delete [] incumbent_;
    
    if (numHeuristics_ > 0) {
#ifdef BLIS_DEBUG
        std::cout << "MODEL: distructor: numHeuristics =  " 
                  << numHeuristics_ << std::endl;
#endif
        BlisHeuristic *tempH = NULL;
        
        for (i = 0; i < numHeuristics_; ++i) {
            tempH = heuristics_[i];
            delete tempH;
        }
        
        delete [] heuristics_;
        heuristics_ = NULL;
    }

    if (generators_ != NULL) {
#ifdef BLIS_DEBUG
        std::cout << "MODEL: distructor: numCutGenerators = " 
                  << numCutGenerators_ << std::endl;
#endif
        BlisConGenerator *temp = NULL;
        for (i = 0; i < numCutGenerators_; ++i) {
            temp = generators_[i];
            delete temp;
        }
        delete [] generators_;
        generators_ = NULL;
    }

    delete constraintPool_;
    delete [] oldConstraints_;
    delete branchStrategy_;

    delete BlisPar_;
}

//#############################################################################

bool 
BlisModel::feasibleSolution(int & numIntegerInfs, int & numObjectInfs)
{
    int numUnsatisfied = 0;
    double sumUnsatisfied = 0.0;
    int preferredWay;
    int j;
    
#if 0
    if (savedLpSolution_ != 0) {
        delete [] savedLpSolution_;
        savedLpSolution_ = 0;
    }
    
    savedLpSolution_ = new double [lpSolver_->getNumCols()];
    
    // Put current solution in safe place
    memcpy(savedLpSolution_, 
           lpSolver_->getColSolution(),
	   lpSolver_->getNumCols() * sizeof(double));
#endif

    for (j = 0; j < numIntVars_; ++j) {
	const BcpsObject * object = objects_[j];
	
	double infeasibility = object->infeasibility(this, preferredWay);
	if (infeasibility) {
	    assert (infeasibility>0);
	    numUnsatisfied++;
	    sumUnsatisfied += infeasibility;
	}
    }
    numIntegerInfs = numUnsatisfied;
    for (; j < numObjects_; ++j) {
	const BcpsObject * object = objects_[j];
	double infeasibility = object->infeasibility(this, preferredWay);
	if (infeasibility) {
	    assert (infeasibility > 0);
	    numUnsatisfied++;
	    sumUnsatisfied += infeasibility;
	}
    }

    numObjectInfs = numUnsatisfied - numIntegerInfs;

    //printf("numUnsatisfied = %d\n",numUnsatisfied);
    
    return (!numUnsatisfied);
}

//#############################################################################
/*
  Set branching priorities.

  Setting integer priorities looks pretty robust; the call to findIntegers
  makes sure that integer objects are in place. Setting priorities for
  other objects is entirely dependent on their existence, and the routine may
  quietly fail in several directions.
*/
void 
BlisModel::passInPriorities (const int * priorities,
			     bool ifObject, 
			     int defaultValue)
{

    // FIXME: not completed.
    int i;

    findIntegers(false);

    if (!priority_) {
	priority_ = new int[numObjects_];
	for (i = 0; i < numObjects_; ++i) {
            priority_[i] = defaultValue;
        }
    }

    if (priorities) {
	if (ifObject) {
	    memcpy(priority_ + numIntVars_, 
                   priorities,
		   (numObjects_ - numIntVars_) * sizeof(int));
        }
	else {
	    memcpy(priority_, priorities, numIntVars_ * sizeof(int));
        }
    }
}

//#############################################################################

AlpsTreeNode * 
BlisModel::createRoot() {
  
  //-------------------------------------------------------------
  // NOTE: Root will be deleted by ALPS. Root is an explicit node.
  //-------------------------------------------------------------
    
  BlisTreeNode* root = new BlisTreeNode;
  BlisNodeDesc* desc = new BlisNodeDesc(this);
  root->setDesc(desc);

  //-------------------------------------------------------------
  // NOTE: Although original data are stored in model when reading. 
  //   Root desc still store a full copy of col and row bounds when creating.
  //   It will store soft differences after finding a branching object. 
  //   The soft difference are due to reduced cost fixing and probing.
  //   Also the added cols and rows will be stored.
  //-------------------------------------------------------------
  int k;

  std::vector<BcpsVariable *> vars = getVariables();
  std::vector<BcpsConstraint *> cons = getConstraints();

#ifdef BLIS_DEBUG
  std::cout << "BLIS: createRoot(): numCoreVariables_=" << numCoreVariables_
	    << ", numCoreConstraints_=" << numCoreConstraints_ << std::endl;
#endif  

  int *varIndices1 = new int [numCoreVariables_];
  int *varIndices2 = new int [numCoreVariables_];
  int *varIndices3 = NULL; //new int [numCoreVariables_];
  int *varIndices4 = NULL; //new int [numCoreVariables_];
  double *vlhe = new double [numCoreVariables_];
  double *vuhe = new double [numCoreVariables_];
  double *vlse = NULL; //new double [numCoreVariables_];
  double *vuse = NULL; //new double [numCoreVariables_];

  int *conIndices1 = new int [numCoreConstraints_];
  int *conIndices2 = new int [numCoreConstraints_];
  int *conIndices3 = NULL; //new int [numCoreConstraints_];
  int *conIndices4 = NULL; //new int [numCoreConstraints_];
  double *clhe = new double [numCoreConstraints_];
  double *cuhe = new double [numCoreConstraints_];
  double *clse = NULL; //new double [numCoreConstraints_];
  double *cuse = NULL; //new double [numCoreConstraints_];

  //-------------------------------------------------------------
  // Get var bounds and indices.
  //-------------------------------------------------------------

  for (k = 0; k < numCoreVariables_; ++k) {
    vlhe[k] = vars[k]->getLbHard();
    vuhe[k] = vars[k]->getUbHard();
    //vlse[k] = vars[k]->getLbSoft();
    //vuse[k] = vars[k]->getUbSoft();
    varIndices1[k] = k;
    varIndices2[k] = k;
    
    //varIndices3[k] = k;
    //varIndices4[k] = k;

#ifdef BLIS_DEBUG_MORE
    std::cout << "BLIS: createRoot(): var "<< k << ": hard: lb=" << vlhe[k]
	      << ", ub=" << vuhe[k] << std::endl;
#endif  
    
  }

  //-------------------------------------------------------------  
  // Get con bounds and indices.
  //-------------------------------------------------------------

  for (k = 0; k < numCoreConstraints_; ++k) {
    clhe[k] = cons[k]->getLbHard();
    cuhe[k] = cons[k]->getUbHard();
    //clse[k] = cons[k]->getLbSoft();
    //cuse[k] = cons[k]->getUbSoft();
    conIndices1[k] = k;
    conIndices2[k] = k;
    //conIndices3[k] = k;
    //conIndices4[k] = k;
  }

  int *tempInd = NULL;
  BcpsObject **tempObj = NULL;
  
  desc->assignVars(0 /*numRem*/, tempInd,
		   0 /*numAdd*/, tempObj,
		   false, numCoreVariables_, varIndices1, vlhe, /*Var hard lb*/
		   false, numCoreVariables_, varIndices2, vuhe, /*Var hard ub*/
		   false, 0, varIndices3, vlse, /*Var soft lb*/
		   false, 0, varIndices4, vuse);/*Var soft ub*/
  desc->assignCons(0 /*numRem*/, tempInd,
		   0 /*numAdd*/, tempObj,
		   false, numCoreConstraints_,conIndices1,clhe, /*Con hard lb*/
		   false, numCoreConstraints_,conIndices2,cuhe, /*Con hard ub*/
		   false, 0,conIndices3,clse, /*Con soft lb*/
		   false, 0,conIndices4,cuse);/*Con soft ub*/

  //-------------------------------------------------------------  
  // Mark it as an explicit node.
  //-------------------------------------------------------------

  root->setExplicit(1);
  
  return root;
}

//#############################################################################

void 
BlisModel::addHeuristic(BlisHeuristic * heuristic)
{
    BlisHeuristic ** temp = heuristics_;
    heuristics_ = new BlisHeuristic * [numHeuristics_ + 1];

    memcpy(heuristics_, temp, numHeuristics_ * sizeof(BlisHeuristic *));
    delete [] temp;
    
    heuristics_[numHeuristics_++] = heuristic;
}

//#############################################################################

void 
BlisModel::addCutGenerator(CglCutGenerator * generator,
			   const char * name,
			   int strategy,
			   bool normal, 
			   bool atSolution,
			   bool whenInfeasible)
{
#if 0
    if (!generators_) {        
        generators_ = new BlisCutGenerator* [50];
    }
    generators_[numCutGenerators_++] = new BlisConGenerator(this, generator, 
                                                            strategy, name,
                                                            normal, atSolution,
                                                            whenInfeasible);
#else
    BlisConGenerator ** temp = generators_;

    generators_ = new BlisConGenerator * [(numCutGenerators_ + 1)];
    memcpy(generators_, temp, numCutGenerators_ * sizeof(BlisConGenerator *));
    
    generators_[numCutGenerators_++] = 
        new BlisConGenerator(this, generator, name, strategy,
                             normal, atSolution, whenInfeasible);
    delete [] temp;
    temp = NULL;
#endif

    useCons_ = ALPS_MAX(strategy, useCons_);
}

//#############################################################################

AlpsEncoded* 
BlisModel::encode() const 
{ 
    AlpsReturnStatus status = AlpsReturnStatusOk;

    // NOTE: "AlpsKnowledgeTypeModel" is the type name.
    AlpsEncoded* encoded = new AlpsEncoded(AlpsKnowledgeTypeModel);

    //------------------------------------------------------
    // Encode Alps part. 
    // NOTE: Nothing to do for Bcps part.
    //------------------------------------------------------

    status = encodeAlps(encoded);
    
    //------------------------------------------------------
    // Encode Blis part. 
    //------------------------------------------------------

    //------------------------------------------------------
    // Blis parameter.
    //------------------------------------------------------

    BlisPar_->pack(*encoded);

    //------------------------------------------------------
    // Get a column matrix.
    //------------------------------------------------------
    
    const CoinPackedMatrix *matrixByCol = lpSolver_->getMatrixByCol();
    
    const int numRows = lpSolver_->getNumRows();
    encoded->writeRep(numRows);
    
    const int numCols = lpSolver_->getNumCols();
    encoded->writeRep(numCols);

#ifdef BLIS_DEBUG
    std::cout << "BlisModel::encode()-- numRows="<< numRows << "; numCols=" 
	      << numCols << std::endl;
#endif

    //------------------------------------------------------
    // Variable bounds.
    //------------------------------------------------------

    const double* collb = lpSolver_->getColLower();
    // NOTE: when write a array to buffer, the length is written
    //       before the actual array. So don't need send numCols or
    //       numRows.
    encoded->writeRep(collb, numCols);
    const double* colub = lpSolver_->getColUpper();
    encoded->writeRep(colub, numCols);

    //------------------------------------------------------
    // Objective.
    //------------------------------------------------------

    const double* obj = lpSolver_->getObjCoefficients();
    encoded->writeRep(obj, numCols);
    const double objSense = lpSolver_->getObjSense();
    encoded->writeRep(objSense);

    //------------------------------------------------------
    // Constraint bounds.
    //------------------------------------------------------

    const double* rowlb = lpSolver_->getRowLower();
    encoded->writeRep(rowlb, numRows);
    const double* rowub = lpSolver_->getRowUpper();
    encoded->writeRep(rowub, numRows);

    //------------------------------------------------------
    // Matrix.
    //------------------------------------------------------

    int numElements = lpSolver_->getNumElements();
    encoded->writeRep(numElements);
    const double* elementValue = matrixByCol->getElements();
    encoded->writeRep(elementValue, numElements);

    const CoinBigIndex* colStart = matrixByCol->getVectorStarts();

    int numStart = numCols + 1;
    encoded->writeRep(colStart, numStart);

    const int* index = matrixByCol->getIndices();
    encoded->writeRep(index, numElements);

    //------------------------------------------------------
    // Variable type.
    //------------------------------------------------------

    encoded->writeRep(numIntVars_);
    encoded->writeRep(intVars_, numIntVars_);

    //------------------------------------------------------
    // Debug.
    //------------------------------------------------------

#ifdef BLIS_DEBUG
    std::cout << "BlisModel::encode()-- objSense="<< objSense
	      << "; numElements="<< numElements 
	      << "; numIntVars_=" << numIntVars_ 
	      << "; numStart = " << numStart <<std::endl;
#endif

#ifdef BLIS_DEBUG
    std::cout << "rowub=";
    for (int i = 0; i < numRows; ++i){
	std::cout <<rowub[i]<<" ";
    }
    std::cout << std::endl;
    std::cout << "elementValue=";
    for (int j = 0; j < numElements; ++j) {
	std::cout << elementValue[j] << " ";
    }
    std::cout << std::endl;    
#endif

    return encoded;
}

//#############################################################################

void
BlisModel::decodeToSelf(AlpsEncoded& encoded) 
{
    AlpsReturnStatus status = AlpsReturnStatusOk;

    //------------------------------------------------------
    // Decode Alps part. 
    // NOTE: Nothing to do for Bcps part.
    //------------------------------------------------------

    status = decodeAlps(encoded);

    //------------------------------------------------------
    // Decode Blis part. 
    //------------------------------------------------------


    //------------------------------------------------------
    // Blis Parameters.
    //------------------------------------------------------

    BlisPar_->unpack(encoded);

    encoded.readRep(numRows_);
    encoded.readRep(numCols_);    

#ifdef BLIS_DEBUG
    std::cout << "BlisModel::decode()-- numRows_="<< numRows_ 
	      << "; numCols_=" << numCols_ << std::endl;
#endif

    //------------------------------------------------------
    // Variable bounds.
    //------------------------------------------------------

    int numCols;
    encoded.readRep(origVarLB_, numCols);
    assert(numCols == numCols_);

    encoded.readRep(origVarUB_, numCols);
    assert(numCols == numCols_);
    
    //------------------------------------------------------
    // Objective.
    //------------------------------------------------------

    encoded.readRep(objCoef_, numCols);
    assert(numCols == numCols_);
    
    encoded.readRep(objSense_);
    
    //------------------------------------------------------
    // Constraint bounds.
    //------------------------------------------------------
    
    int numRows;
    
    encoded.readRep(origConLB_, numRows);
    assert(numRows == numRows_);

    encoded.readRep(origConUB_, numRows);
    assert(numRows == numRows_);
    
    //------------------------------------------------------
    // Matrix.
    //------------------------------------------------------

    encoded.readRep(numElems_);

    int numElements;
    double* elementValue;
    encoded.readRep(elementValue, numElements);
    assert(numElements == numElems_);

    CoinBigIndex* colStart;

    int numStart;
    encoded.readRep(colStart, numStart);
    assert(numStart == numCols_ + 1);

    int* index;
    encoded.readRep(index, numElements);
    assert(numElements == numElems_);

    //------------------------------------------------------
    // Variable type.
    //------------------------------------------------------

    encoded.readRep(numIntVars_);
    int numInts;
    encoded.readRep(intVars_, numInts);
    assert(numInts == numIntVars_);

    //------------------------------------------------
    // Classify variable type
    //------------------------------------------------

    colType_ = new char [numCols_];
    int j, ind = -1;
    for(j = 0; j < numCols_; ++j) {
	colType_[j] = 'C';
    }
    for(j = 0; j < numInts; ++j) {
	ind = intVars_[j];
	assert(ind >= 0 && ind < numCols_);
	if (origVarLB_[ind] == 0 && origVarUB_[ind] == 1.0) {
	    colType_[ind] = 'B';
	}
	else {
	    colType_[ind] = 'I';
	}
    }
    
    //------------------------------------------------------
    // Debug.
    //------------------------------------------------------

#ifdef BLIS_DEBUG
    std::cout << "BlisModel::decode()-- objSense="<< objSense_
	      <<  "; numElements="<< numElements 
	      << "; numberIntegers_=" << numIntVars_ 
	      << "; numStart = " << numStart <<std::endl;
#endif

#ifdef BLIS_DEBUG
    int i;
    std::cout << "origConUB_=";
    for (i = 0; i < numRows; ++i){
	std::cout <<origConUB_[i]<<" ";
    }
    std::cout << std::endl;
    std::cout << "elementValue=";
    for (j = 0; j < numElements; ++j) {
	std::cout << elementValue[j] << " ";
    }
    std::cout << std::endl;  
    std::cout << "index=";
    for (j = 0; j < numElements; ++j) {
	std::cout << index[j] << " ";
    }
    std::cout << std::endl;  
    std::cout << "colStart=";
    for (j = 0; j < numCols + 1; ++j) {
	std::cout << colStart[j] << " ";
    }
    std::cout << std::endl;   
#endif

    //------------------------------------------------------
    // Check if lpSolver_ is declared in main.
    //------------------------------------------------------

    assert(lpSolver_);

    //------------------------------------------------------
    // Load data to lp solver.
    //------------------------------------------------------

    lpSolver_->loadProblem(numCols, numRows,
			   colStart, index, elementValue,
			   origVarLB_, origVarUB_, 
			   objCoef_,
			   origConLB_, origConUB_);
    
    lpSolver_->setObjSense(objSense_);
    lpSolver_->setInteger(intVars_, numIntVars_);

    //------------------------------------------------------
    // Clean up.
    //------------------------------------------------------

    delete [] colStart;
    colStart = NULL;

    delete [] index;
    index = NULL;
}

//#############################################################################

/** Register knowledge. */
void 
BlisModel::registerKnowledge() {
    // Register model, solution, and tree node
    assert(broker_);
    broker_->registerClass(AlpsKnowledgeTypeModel, new BlisModel);
    std::cout << "Register Alps model." << std::endl;
    
    broker_->registerClass(AlpsKnowledgeTypeNode, new BlisTreeNode(this));
    std::cout << "Register Alps node." << std::endl;
    
    broker_->registerClass(AlpsKnowledgeTypeSolution, new BlisSolution);
    std::cout << "Register Alps solution." << std::endl;
    
    broker_->registerClass(BcpsKnowledgeTypeConstraint, new BlisConstraint);
    std::cout << "Register Bcps constraint." << std::endl;
    
    broker_->registerClass(BcpsKnowledgeTypeVariable, new BlisVariable);
    std::cout << "Register Bcps variable." << std::endl;
}

//#############################################################################

/** Log of specific models. */
void 
BlisModel::modelLog() 
{

    int logFileLevel = AlpsPar_->entry(AlpsParams::logFileLevel);
    
    if (logFileLevel > -1) {
	std::string logfile = AlpsPar_->entry(AlpsParams::logFile);
	std::ofstream logFout(logfile.c_str(), std::ofstream::app);
	writeParameters(logFout);
    }
}

//#############################################################################
