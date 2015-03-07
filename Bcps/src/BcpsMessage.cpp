/*===========================================================================*
 * This file is part of the Branch, Constrain and Price Software (BiCePS)    *
 *                                                                           *
 * BiCePS is distributed under the Eclipse Public License as part of the     *
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
 * Copyright (C) 2001-2015, Lehigh University, Yan Xu, and Ted Ralphs.       *
 * All Rights Reserved.                                                      *
 *===========================================================================*/

#include "BcpsMessage.h"
#include <cstring>

//#############################################################################

typedef struct {
    BCPS_Message internalNumber;
    int externalNumber;              // or continuation
    char detail;
    const char * message;
} Bcps_message;

//#############################################################################

static Bcps_message us_english[]=
{
    {BCPS_DUMMY_END, 999999, 0, ""}
};

//#############################################################################

/* Constructor */
BcpsMessage::BcpsMessage(Language language) 
    :
    CoinMessages(sizeof(us_english) / sizeof(Bcps_message))
{
    language_ = language;
    strcpy(source_, "Bcps");
    Bcps_message * message = us_english;

    while (message->internalNumber != BCPS_DUMMY_END) {
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
	while (message->internalNumber != BCPS_DUMMY_END) {
	    replaceMessage(message->internalNumber, message->message);
	    message++;
	}
    }
}
