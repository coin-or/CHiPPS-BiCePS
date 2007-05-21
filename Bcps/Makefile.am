#==========================================================================#
# This file is part of the Branch, Constrain and Price Software (BiCePS).  #
#                                                                          #
# BiCePS is distributed under the Common Public License as part of the     #
# COIN-OR repository (http://www.coin-or.org).                             #
#                                                                          #
# Authors:                                                                 #
#                                                                          #
#          Yan Xu, Lehigh University                                       #
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
# Copyright (C) 2001-2006, Lehigh University, Yan Xu, and Ted Ralphs.      #
# All Rights Reserved.                                                     #
#==========================================================================#

## $Id$

# Author:  Andreas Waechter           IBM    2006-04-13
#          Yan Xu                     SAS    2006-08-27

AUTOMAKE_OPTIONS = foreign

########################################################################
#                          Subdirectories                              #
########################################################################

SUBDIRS = src

# We don't want to compile the test subdirectory, unless the test target is
# specified.  But we need to list it as subdirectory to make sure that it is
# included in the tarball

if ALWAYS_FALSE
  SUBDIRS += test
endif

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

clean-local:
	if test -r test/Makefile; then cd test; $(MAKE) clean; fi

distclean-local:
	if test -r test/Makefile; then cd test; $(MAKE) distclean; fi

.PHONY: test unitTest

########################################################################
#                  Installation of the addlibs file                    #
########################################################################

addlibsfile = biceps_addlibs.txt

install-exec-local:
	$(install_sh_DATA) $(addlibsfile) $(DESTDIR)$(libdir)/$(addlibsfile)

uninstall-local:
	rm -f $(DESTDIR)$(libdir)/$(addlibsfile)

########################################################################
#                         Maintainer Stuff                             #
########################################################################

# Files that are generated and should be cleaned with make distclean
DISTCLEANFILES =

include ../BuildTools/Makemain.inc
