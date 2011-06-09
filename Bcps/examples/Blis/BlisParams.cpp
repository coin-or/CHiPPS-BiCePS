/*===========================================================================*
 * This file is part of the Abstract Library for Parallel Search (ALPS).     *
 *                                                                           *
 * ALPS is distributed under the Common Public License as part of the        *
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
 * Copyright (C) 2001-2011, Lehigh University, Yan Xu, and Ted Ralphs.       *
 * All Rights Reserved.                                                      *
 *===========================================================================*/

#include "BlisParams.h"

using std::make_pair;

void 
BlisParams::createKeywordList() {

  //--------------------------------------------------------
  // Create the list of keywords for parameter file reading
  //--------------------------------------------------------

  //--------------------------------------------------------
  // CharPar
  //--------------------------------------------------------
  
  keys_.push_back(make_pair(std::string("Blis_useCons"),
			    AlpsParameter(AlpsBoolPar, useCons)));

  keys_.push_back(make_pair(std::string("Blis_useHeuristics"),
			    AlpsParameter(AlpsBoolPar, useHeuristics)));

  keys_.push_back(make_pair(std::string("Blis_cutDuringRampup"),
			    AlpsParameter(AlpsBoolPar, cutDuringRampup)));
    
  //--------------------------------------------------------
  // BoolArrayPar
  //--------------------------------------------------------
  
  //--------------------------------------------------------
  // Int Parameters
  //--------------------------------------------------------

  keys_.push_back(make_pair(std::string("Blis_strongCandSize"),
			    AlpsParameter(AlpsIntPar, strongCandSize)));

  // keys_.push_back(make_pair(std::string("Blis_logLevel"),
  //		    AlpsParameter(AlpsIntPar, logLevel)));

  keys_.push_back(make_pair(std::string("Blis_branchStrategy"),
			    AlpsParameter(AlpsIntPar, branchStrategy)));

  keys_.push_back(make_pair(std::string("Blis_heurRound"),
			    AlpsParameter(AlpsIntPar, heurRound)));

  keys_.push_back(make_pair(std::string("Blis_cutClique"),
			    AlpsParameter(AlpsIntPar, cutClique)));

  keys_.push_back(make_pair(std::string("Blis_cutGomory"),
			    AlpsParameter(AlpsIntPar, cutGomory)));

  keys_.push_back(make_pair(std::string("Blis_cutFlowCover"),
			    AlpsParameter(AlpsIntPar, cutFlowCover)));

  keys_.push_back(make_pair(std::string("Blis_cutKnapsack"),
			    AlpsParameter(AlpsIntPar, cutKnapsack)));
  
  keys_.push_back(make_pair(std::string("Blis_cutMir"),
			    AlpsParameter(AlpsIntPar, cutMir)));
  
  keys_.push_back(make_pair(std::string("Blis_cutOddHole"),
			    AlpsParameter(AlpsIntPar, cutOddHole)));
  
  keys_.push_back(make_pair(std::string("Blis_cutProbing"),
			    AlpsParameter(AlpsIntPar, cutProbing)));

  keys_.push_back(make_pair(std::string("Blis_cutTwoMir"),
			    AlpsParameter(AlpsIntPar, cutTwoMir)));

  keys_.push_back(make_pair(std::string("Blis_pseudoRelibility"),
			    AlpsParameter(AlpsIntPar, pseudoRelibility)));

  keys_.push_back(make_pair(std::string("Blis_lookAhead"),
			    AlpsParameter(AlpsIntPar, lookAhead)));

  
  //--------------------------------------------------------
  // Double Parameters.
  //--------------------------------------------------------
  
  keys_.push_back(make_pair(std::string("Blis_integerTol"),
			    AlpsParameter(AlpsDoublePar, integerTol)));
  
  keys_.push_back(make_pair(std::string("Blis_cutoffInc"),
			    AlpsParameter(AlpsDoublePar, cutoffInc)));
  
  keys_.push_back(make_pair(std::string("Blis_optimalRelGap"),
			    AlpsParameter(AlpsDoublePar, optimalRelGap)));
  
  keys_.push_back(make_pair(std::string("Blis_optimalRelGap"),
			    AlpsParameter(AlpsDoublePar, optimalAbsGap)));
  
  keys_.push_back(make_pair(std::string("Blis_pseudoWeight"),
			    AlpsParameter(AlpsDoublePar, pseudoWeight)));
  
  keys_.push_back(make_pair(std::string("Blis_cutFactor"),
			    AlpsParameter(AlpsDoublePar, cutFactor)));
  
  keys_.push_back(make_pair(std::string("Blis_denseConFactor"),
			    AlpsParameter(AlpsDoublePar, denseConFactor)));

  keys_.push_back(make_pair(std::string("Blis_scaleConFactor"),
			    AlpsParameter(AlpsDoublePar, scaleConFactor)));
  
  //--------------------------------------------------------
  // String Parameters.
  //--------------------------------------------------------

}

//#############################################################################

void 
BlisParams::setDefaultEntries() {

  //-------------------------------------------------------------
  // Char Parameters.
  //-------------------------------------------------------------

  setEntry(useHeuristics, true);
  setEntry(cutDuringRampup, false);
  setEntry(useCons, true);
  
  //-------------------------------------------------------------
  // Int Parameters.
  //-------------------------------------------------------------

  setEntry(strongCandSize, 10);
  // setEntry(logLevel, 0);
  setEntry(branchStrategy, 1);
  setEntry(heurRound, 0);
  setEntry(cutClique, 0);
  setEntry(cutGomory, 0);
  setEntry( cutFlowCover, 0);
  setEntry(cutKnapsack, 0);
  setEntry(cutMir, 0);
  setEntry(cutOddHole, -2);
  setEntry(cutProbing, 0);
  setEntry(cutTwoMir, -2);
  setEntry(pseudoRelibility, 8);
  setEntry(lookAhead, 4);
  
  //-------------------------------------------------------------
  // Double Parameters
  //-------------------------------------------------------------

  setEntry(integerTol, 1.0e-5);
  setEntry(cutoffInc, 1.0e-6);
  setEntry(optimalRelGap, 1.0e-6);
  setEntry(optimalAbsGap, 1.0e-4);
  setEntry(pseudoWeight, 1.0);
  setEntry(cutFactor, 4.0);
  setEntry(denseConFactor, 5.0);
  setEntry(scaleConFactor, 1000000.0);

  //-------------------------------------------------------------
  // String Parameters
  //-------------------------------------------------------------
  
}
