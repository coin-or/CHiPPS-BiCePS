#==========================================================================#
# This file is part of the Abstract Library for Parallel Search (ALPS).    #
#                                                                          #
# ALPS is distributed under the Eclipse Public License as part of the      #
# COIN-OR repository (http://www.coin-or.org).                             #
#                                                                          #
# Authors:                                                                 #
#                                                                          #
#          Yan Xu, Lehigh University                                       #
#          Aykut Bulut, Lehigh University                                  #
#          Ted Ralphs, Lehigh University                                   #
#                                                                          #
# Conceptual Design:                                                       #
#                                                                          #
#          Yan Xu, Lehigh University                                       #
#          Ted Ralphs, Lehigh University                                   #
#          Laszlo Ladanyi, IBM T.J. Watson Research Center                 #
#          Matthew Saltzman, Clemson University                            #
#                                                                          #
#                                                                          #
# Copyright (C) 2001-2019, Lehigh University, Yan Xu, Aykut Bulut, and     #
#                          Ted Ralphs.                                     #
# All Rights Reserved.                                                     #
#==========================================================================#

## $Id: Makefile.am 507 2006-06-03 03:20:16Z andreasw $

########################################################################
#                      unitTest for CoinUtils                          #
########################################################################

noinst_PROGRAMS = unitTest

nodist_unitTest_SOURCES = \
	BlisBranchObjectInt.cpp BlisBranchObjectInt.h \
	BlisBranchStrategyPseudo.cpp BlisBranchStrategyPseudo.h \
	BlisConGenerator.cpp BlisConGenerator.h \
	BlisConstraint.cpp BlisConstraint.h \
	Blis.h \
	BlisHelp.cpp BlisHelp.h \
	BlisHeuristic.h \
	BlisHeurRound.cpp BlisHeurRound.h \
	BlisLicense.h \
	BlisMain.cpp \
	BlisMessage.cpp BlisMessage.h \
	BlisModel.cpp BlisModel.h \
	BlisNodeDesc.h \
	BlisObjectInt.cpp BlisObjectInt.h \
	BlisParams.cpp BlisParams.h \
	BlisPseudo.cpp BlisPseudo.h \
	BlisSolution.h \
	BlisSubTree.h \
	BlisTreeNode.cpp BlisTreeNode.h \
	BlisVariable.h

# List libraries of COIN projects
unitTest_LDADD = ../src/libBcps.la $(UNITTEST_LFLAGS) $(BCPSLIB_LFLAGS)

AM_LDFLAGS = $(LT_LDFLAGS)

AM_CPPFLAGS = -I$(srcdir)/../src $(UNITTEST_CFLAGS) $(BCPSLIB_CFLAGS) 	

all: test

test: unitTest$(EXEEXT)
	./unitTest$(EXEEXT) -Alps_instance ./flugpl.mps

.PHONY: test

########################################################################
#                          Cleaning stuff                              #
########################################################################

# Here we list everything that is not generated by the compiler, e.g.,
# output files of a program

DISTCLEANFILES = flugpl.log
