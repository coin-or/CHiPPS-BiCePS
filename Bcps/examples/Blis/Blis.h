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

//#############################################################################

#ifndef Blis_h_
#define Blis_h_

//#############################################################################

#define BLIS_OK             0
#define BLIS_LP_OPTIMAL     0
#define BLIS_LP_ABANDONED   1
#define BLIS_LP_PRIMAL_INF  2
#define BLIS_LP_DUAL_INF    3
#define BLIS_LP_PRIMAL_LIM  4
#define BLIS_LP_DUAL_LIM    5
#define BLIS_LP_ITER_LIM    6

#define BLIS_ERR_LP         100
#define BLIS_INF            200
#define BLIS_UNBOUND        201
#define BLIS_OPTIMAL          0
#define BLIS_UNKNOWN        202

//#############################################################################

enum BLIS_SOL_TYPE {
    BLIS_SOL_BOUNDING,
    BLIS_SOL_BRANCHING,
    BLIS_SOL_DIVING,
    BLIS_SOL_ROUNDING,
    BLIS_SOL_STRONG
};

//#############################################################################

/** Branching object type. */
enum BLIS_BO_TYPE {
    BLIS_BO_NONE = 0,
    BLIS_BO_INT,
    BLIS_BO_SOS
};

//#############################################################################

#endif
