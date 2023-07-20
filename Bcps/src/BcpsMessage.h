/*===========================================================================*
 * This file is part of the Branch, Constrain and Price Software (BiCePS)    *
 *                                                                           *
 * BiCePS is distributed under the Eclipse Public License as part of the     *
 * COIN-OR repository (http://www.coin-or.org).                              *
 *                                                                           *
 * Authors:                                                                  *
 *                                                                           *
 *          Yan Xu, Lehigh University                                        *
 *          Aykut Bulut, Lehigh University                                   *
 *          Ted Ralphs, Lehigh University                                    *
 *                                                                           *
 * Conceptual Design:                                                        *
 *                                                                           *
 *          Yan Xu, Lehigh University                                        *
 *          Ted Ralphs, Lehigh University                                    *
 *          Laszlo Ladanyi, IBM T.J. Watson Research Center                  *
 *          Matthew Saltzman, Clemson University                             *
 *                                                                           *
 * Copyright (C) 2001-2023, Lehigh University, Yan Xu, and Ted Ralphs.       *
 * All Rights Reserved.                                                      *
 *===========================================================================*/

#ifndef BcpsMessage_H_
#define BcpsMessage_H_

//#############################################################################
// This file is modified from SbbMessage.hpp
//#############################################################################

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

/** This deals with Bcps messages. */
#include "CoinMessageHandler.hpp"

enum BCPS_Grumpy_Msg_Type {
  BCPS_GRUMPY_BRANCHED = 0,
  BCPS_GRUMPY_CANDIDATE,
  BCPS_GRUMPY_HEURISTIC,
  BCPS_GRUMPY_INTEGER,
  BCPS_GRUMPY_FATHOMED,
  BCPS_GRUMPY_PREGNANT,
  BCPS_GRUMPY_INFEASIBLE
};

enum BCPS_Message {
  // tree node
  BCPS_NODE_BRANCHONINT,
  BCPS_NODE_UNEXPECTEDSTATUS,
  // grumpy messages
  BCPS_GRUMPY_MESSAGE_LONG,
  BCPS_GRUMPY_MESSAGE_MED,
  BCPS_GRUMPY_MESSAGE_SHORT,
  /// end of messages
  BCPS_DUMMY_END
};

class BcpsMessage : public CoinMessages {

public:

    /**@name Constructors etc */
    //@{
    /** Constructor */
    BcpsMessage(Language language=us_en);
    //@}

};

enum BCPS_Debug_Level {
  BCPS_DLOG_PROCESS = 8,
  BCPS_DLOG_GRUMPY = 16
};

#endif
