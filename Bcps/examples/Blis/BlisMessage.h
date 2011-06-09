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
 * Copyright (C) 2001-2011, Lehigh University, Yan Xu, and Ted Ralphs.       *
 * All Rights Reserved.                                                      *
 *===========================================================================*/

#ifndef BlisMessage_H_
#define BlisMessage_H_

//#############################################################################
// This file is modified from SbbMessage.hpp
//#############################################################################

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

/** This deals with Blis messages (as against Clp messages etc).
    CoinMessageHandler.hpp is the general part of message handling.
    All it has are enum's for the various messages.
    BlisMessage.cpp has text in various languages.
    
    It is trivial to use the .hpp and .cpp file as a basis for
    messages for other components.
 */

#include "CoinMessageHandler.hpp"

enum BLIS_Message
{
  BLIS_END_GOOD,
  BLIS_MAXNODES,
  BLIS_MAXTIME,
  BLIS_MAXSOLS,
  BLIS_SOLUTION,
  BLIS_END,
  BLIS_INFEAS,
  BLIS_STRONG,
  BLIS_SOLINDIVIDUAL,
  BLIS_INTEGERINCREMENT,
  BLIS_STATUS,
  BLIS_GAP,
  BLIS_ROUNDING,
  BLIS_ROOT,
  BLIS_GENERATOR,
  BLIS_BRANCH,
  BLIS_STRONGSOL,
  BLIS_NOINT,
  BLIS_VUB_PASS,
  BLIS_VUB_END,
  BLIS_NOTFEAS1,
  BLIS_NOTFEAS2,
  BLIS_NOTFEAS3,
  BLIS_CUTOFF_WARNING1,
  BLIS_CUTS,
  BLIS_BRANCHSOL,
  BLIS_DUMMY_END
};

class BlisMessage : public CoinMessages {

public:

  /**@name Constructors etc */
  //@{
  /** Constructor */
  BlisMessage(Language language=us_en);
  //@}

};

#endif
