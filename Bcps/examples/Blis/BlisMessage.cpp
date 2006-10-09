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

#include "BlisMessage.h"

//#############################################################################

typedef struct {
    BLIS_Message internalNumber;
    int externalNumber;              // or continuation
    char detail;
    const char * message;
} Blis_message;

//#############################################################################

static Blis_message us_english[]=
{
    {BLIS_END_GOOD,1,1,"Search completed - best objective %g, took %d iterations and %d nodes"},
    {BLIS_MAXNODES,3,1,"Exiting on maximum nodes"},
    {BLIS_MAXTIME,20,1,"Exiting on maximum time"},
    {BLIS_MAXSOLS,19,1,"Exiting on maximum solutions"},
    {BLIS_SOLUTION,4,1,"Integer solution of %g found after %d iterations and %d nodes"},
    {BLIS_END,5,1,"Partial search took %d iterations and %d nodes"},
    {BLIS_INFEAS,6,1,"The LP relaxation is infeasible or too expensive"},
    {BLIS_STRONG,7,3,"Strong branching on %d (%d), down %g (%d) up %g (%d) value %d"},
    {BLIS_SOLINDIVIDUAL,8,2,"%d has value %g"},
    {BLIS_INTEGERINCREMENT,9,1,"Objective coefficients multiple of %g"},
    {BLIS_STATUS,10,1,"Process[%d]: after %d nodes, %d on tree, %g best solution, best possible %g"},
    {BLIS_GAP,11,1,"Exiting as integer gap of %g less than %g"},
    {BLIS_ROUNDING,12,1,"Integer solution of %g found by rounding after %d iterations and %d nodes"},
    {BLIS_ROOT,13,3,"At root node, %d cuts changed objective from %g to %g in %d passes"},
    {BLIS_GENERATOR,14,2,"Cut generator %d (%s) - %d row cuts (%d active), %d column cuts - new frequency is %d"},
    {BLIS_BRANCH,15,3,"Node %d Obj %g Unsat %d depth %d"},
    {BLIS_STRONGSOL,16,1,"Integer solution of %g found by strong branching after %d iterations and %d nodes"},
    {BLIS_NOINT,3007,0,"No integer variables - nothing to do"},
    {BLIS_VUB_PASS,17,1,"%d solved, %d variables fixed, %d tightened"},
    {BLIS_VUB_END,18,1,"After tightenVubs, %d variables fixed, %d tightened"},
    {BLIS_NOTFEAS1,21,2,"On closer inspection node is infeasible"},
    {BLIS_NOTFEAS2,22,2,"On closer inspection objective value of %g above cutoff of %g"},
    {BLIS_NOTFEAS3,23,2,"Allowing solution, even though largest row infeasibility is %g"},
    {BLIS_CUTOFF_WARNING1,23,1,"Cutoff set to %g - equivalent to best solution of %g"},
    {BLIS_CUTS,24,1, "At node %d, %d cuts changed objective from %g to %g in %d passes"},
    {BLIS_BRANCHSOL,25,1,"Integer solution of %g found by branching after %d iterations and %d nodes"},
    {BLIS_DUMMY_END, 999999, 0, ""}
};

//#############################################################################

/* Constructor */
BlisMessage::BlisMessage(Language language) 
    :
    CoinMessages(sizeof(us_english) / sizeof(Blis_message))
{
    language_ = language;
    strcpy(source_, "Blis");
    Blis_message * message = us_english;

    while (message->internalNumber != BLIS_DUMMY_END) {
	CoinOneMessage oneMessage(message->externalNumber, message->detail,
				  message->message);
	addMessage(message->internalNumber, oneMessage);
	message++;
    }

    // now override any language ones

    switch (language) {

    default:
	message = NULL;
	break;
    }

    // replace if any found
    if (message) {
	while (message->internalNumber != BLIS_DUMMY_END) {
	    replaceMessage(message->internalNumber, message->message);
	    message++;
	}
    }
}
