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
 * Copyright (C) 2001-2005, International Business Machines                  *
 * Corporation, Lehigh University, Yan Xu, Ted Ralphs, Matthew Salzman and   *
 * others. All Rights Reserved.                                              *
 *===========================================================================*/

//#############################################################################

#ifndef BlisHelp_h_
#define BlisHelp_h_

#include "AlpsEncoded.h"

class CoinWarmStartBasis;
class OsiRowCut;
class BlisConstraint;
class BlisModel;

//#############################################################################

/** Convert a OsiRowCut to a Blis Contraint. */
BlisConstraint * BlisOsiCutToConstraint(const OsiRowCut *rowCut);

/** Convert a Blis constraint to a OsiRowCut. */
OsiRowCut * BlisConstraintToOsiCut(const BlisConstraint * con);

/** Strong branching on a variable colInd. */
int BlisStrongBranch(BlisModel *model, double objValue, int colInd, double x,
                     const double *saveLower, const double *saveUpper,
		     bool &downKeep, bool &downFinished, double &downDeg,
		     bool &upKeep, bool &upFinished, double &upDeg);

/** Pack coin warm start into an encoded object. */
int BlisEncodeWarmStart(AlpsEncoded *encoded, const CoinWarmStartBasis *ws);

/** Unpack coin warm start from an encoded object. */
CoinWarmStartBasis *BlisDecodeWarmStart(AlpsEncoded &encoded,
					AlpsReturnStatus *rc);

/** Compute and return a hash value of an Osi row cut. */
double BlisHashingOsiRowCut(const OsiRowCut *rowCut, 
			    const BlisModel *model);

/** Check if a row cut parallel with another row cut. */
bool BlisParallelCutCut(OsiRowCut * rowCut1,
			OsiRowCut * rowCut2,
			double threshold = 1.0);

/** Check if a row cut parallel with a constraint. */
bool BlisParallelCutCon(OsiRowCut * rowCut,
			BlisConstraint * con,
			double threshold = 1.0);


#endif
