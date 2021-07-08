#==========================================================================#
# This file is part of the Branch, Constrain and Price Software (BiCePS).  #
#                                                                          #
# BiCePS is distributed under the Eclipse Public License as part of the    #
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

## $Id$

# Author:  Andreas Waechter           IBM    2006-04-13
#          Yan Xu                     Lehigh    2006-08-27

include BuildTools/Makemain.inc

########################################################################
#                          Subdirectories                              #
########################################################################

SUBDIRS = src

pkgconfiglib_DATA = bcps.pc

########################################################################
#             Additional files to be included in tarball               #
########################################################################

# Here we need include all files that are not mentioned in other Makefiles

EXTRA_DIST = \
	examples/Blis/BlisBranchObjectInt.cpp \
	examples/Blis/BlisBranchObjectInt.h \
	examples/Blis/BlisBranchStrategyPseudo.cpp \
	examples/Blis/BlisBranchStrategyPseudo.h \
	examples/Blis/BlisBranchStrategyRel.cpp \
	examples/Blis/BlisBranchStrategyRel.h \
	examples/Blis/BlisBranchStrategyStrong.cpp \
	examples/Blis/BlisBranchStrategyStrong.h \
	examples/Blis/BlisConGenerator.cpp \
	examples/Blis/BlisConGenerator.h \
	examples/Blis/BlisConstraint.cpp \
	examples/Blis/BlisConstraint.h \
	examples/Blis/Blis.h \
	examples/Blis/BlisHelp.cpp \
	examples/Blis/BlisHelp.h \
	examples/Blis/BlisHeuristic.h \
	examples/Blis/BlisHeurRound.cpp \
	examples/Blis/BlisHeurRound.h \
	examples/Blis/BlisLicense.h \
	examples/Blis/BlisMain.cpp \
	examples/Blis/BlisMessage.cpp \
	examples/Blis/BlisMessage.h \
	examples/Blis/BlisModel.cpp \
	examples/Blis/BlisModel.h \
	examples/Blis/BlisNodeDesc.h \
	examples/Blis/BlisObjectInt.cpp \
	examples/Blis/BlisObjectInt.h \
	examples/Blis/blis.par \
	examples/Blis/BlisParams.cpp \
	examples/Blis/BlisParams.h \
	examples/Blis/BlisPseudo.cpp \
	examples/Blis/BlisPseudo.h \
	examples/Blis/BlisSolution.h \
	examples/Blis/BlisSubTree.h \
	examples/Blis/BlisTreeNode.cpp \
	examples/Blis/BlisTreeNode.h \
	examples/Blis/BlisVariable.h \
	examples/Blis/Makefile.in \
	examples/Blis/AUTHORS \
	examples/Blis/INSTALL \
	examples/Blis/LICENSE \
	examples/Blis/Makefile.in \
	examples/Blis/README \
	examples/Blis/data/egout.mps \
	examples/Blis/data/flugpl.mps \
	examples/Blis/data/gen.mps \
	examples/Blis/data/knap1.mps \
	examples/Blis/data/knap2.mps \
	examples/Blis/data/knap3.mps


########################################################################
#                           Extra Targets                              #
########################################################################

test: all
	cd test; $(MAKE) test

unitTest: test

clean-local: clean-doxygen-docs
	if test -r test/Makefile; then cd test; $(MAKE) clean; fi

install-exec-local: install-doc

uninstall-local: uninstall-doc uninstall-doxygen-docs

.PHONY: test unitTest doxydoc

