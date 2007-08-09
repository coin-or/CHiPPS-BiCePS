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

#ifndef BlisSubTree_h_
#define BlisSubTree_h_

//#############################################################################
/** This is the class in which we are finally able to concretely define the
    bounding procedure. Here we can assume that we have an LP solver and that
    the objects are cuts and variables, etc. */
//#############################################################################

class BlisTreeNode : public BcpsTreeNode {

   virtual void process();

   /** Takes the explicit description of the current active node and creates
       the children, adds them to the priority queue, etc. Most likely, the
       node description will contain information about how the branching is
       to be done. */
   void branch(AlpsTreeNode& explicitNode);

   /** Takes the explicit description of the current active node and decides
       whether or not it should be fathomed. */
   virtual bool fathom();
   
};

//#############################################################################
//#############################################################################

#endif
