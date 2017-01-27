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

#include <iostream>

#include "CoinError.hpp"
#include "CoinTime.hpp"

#include "OsiSolverInterface.hpp"
#include "OsiClpSolverInterface.hpp"

#include "CglFlowCover.hpp"
#include "CglGomory.hpp"
#include "CglProbing.hpp"
#include "CglKnapsackCover.hpp"
#include "CglOddHole.hpp"

#ifdef COIN_HAS_MPI
#  include "AlpsKnowledgeBrokerMPI.h"
#else
#  include "AlpsKnowledgeBrokerSerial.h"
#endif

#include "BlisModel.h"

//#############################################################################

int main(int argc, char *argv[]) 
{
    try{
	// Set up lp solver
        OsiClpSolverInterface lpSolver;
	lpSolver.getModelPtr()->setDualBound(1.0e10);
	lpSolver.messageHandler()->setLogLevel(0);
	
	// Create BLIS model 
	BlisModel model;
	model.setSolver(&lpSolver);
	
#ifdef  COIN_HAS_MPI
	AlpsKnowledgeBrokerMPI broker(argc, argv, model);
#else
	AlpsKnowledgeBrokerSerial broker(argc, argv, model); 
#endif

	// Search for best solution
	broker.search(&model);
	
	// Report the best solution found and its ojective value
	// broker.printBestResult();
    }
    catch(CoinError& er) {
	std::cerr << "\nBLIS ERROR: \"" << er.message() 
		  << "\""<< std::endl
		  << "             from function \"" << er.methodName()
		  << "\""<< std::endl
		  << "             from class \"" << er.className()
		  << "\"" << std::endl;
    }
    catch(...) {
	std::cerr << "Something went wrong!" << std::endl;
    }
    
    
    return 0;
}

//#############################################################################

