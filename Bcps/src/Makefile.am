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
# Copyright (C) 2001-2009, Lehigh University, Yan Xu, and Ted Ralphs.      #
# All Rights Reserved.                                                     #
#==========================================================================#

## $Id$

# Author:  Andreas Waechter           IBM    2006-04-13
#          Yan Xu                     Lehigh    2006-08-27

AUTOMAKE_OPTIONS = foreign

########################################################################
#                               libBcps                                #
########################################################################

# Name of the library compiled in this directory.  We want it to be installed
# in the 'lib' directory
lib_LTLIBRARIES = libBcps.la

# List all source files for this library, including headers
libBcps_la_SOURCES = \
	Bcps.h \
	BcpsBranchObject.cpp BcpsBranchObject.h \
	BcpsBranchStrategy.cpp BcpsBranchStrategy.h \
	BcpsInterface.h \
	BcpsLicense.h \
	BcpsMessage.h BcpsMessage.cpp \
	BcpsModel.h BcpsModel.cpp \
	BcpsNodeDesc.cpp BcpsNodeDesc.h \
	BcpsObject.cpp BcpsObject.h \
	BcpsObjectPool.h BcpsSolution.cpp \
	BcpsSolution.h BcpsSubTree.cpp \
	BcpsSubTree.h BcpsTreeNode.cpp \
	BcpsTreeNode.h

# This is for libtool (on Windows)
libBcps_la_LDFLAGS = $(LT_LDFLAGS)

########################################################################
#                            Additional flags                          #
########################################################################

# Here list all include flags, relative to this "srcdir" directory.  This
# "cygpath" stuff is necessary to compile with native compilers on Windows.
AM_CPPFLAGS = \
	-I`$(CYGPATH_W) $(COINUTILSSRCDIR)/src` \
	-I`$(CYGPATH_W) $(COINUTILSOBJDIR)/inc` \
	-I`$(CYGPATH_W) $(ALPSSRCDIR)/src` \
	-I`$(CYGPATH_W) $(ALPSOBJDIR)/inc`

# Include MPI header directory if MPI is used
if COIN_HAS_MPI
AM_CPPFLAGS +=  -I`$(CYGPATH_W) $(MPIINCDIR)`
endif

# This line is necessary to allow VPATH compilation with MS compilers
# on Cygwin
DEFAULT_INCLUDES = -I. -I`$(CYGPATH_W) $(srcdir)` -I$(top_builddir)/inc

########################################################################
#                Headers that need to be installed                     #
########################################################################

# Here list all the header files that are required by a user of the library,
# and that therefore should be installed in 'install/coin'
includecoindir = $(includedir)/coin
includecoin_HEADERS = \
	../inc/config_bcps.h \
	Bcps.h \
	BcpsBranchObject.h \
	BcpsBranchStrategy.h \
	BcpsInterface.h \
	BcpsLicense.h \
	BcpsMessage.h \
	BcpsModel.h \
	BcpsNodeDesc.h \
	BcpsObject.h \
	BcpsObjectPool.h \
	BcpsSolution.h \
	BcpsSubTree.h \
	BcpsTreeNode.h


#############################################################################
# Create the Config.h file that always defines HAVE_CONFIG_H and install it #
#############################################################################

# You only need to adapt the following line
ConfigHeader = BcpsConfig.h

install-exec-local:
	echo "#ifndef HAVE_CONFIG_H" >bla
	echo "#define HAVE_CONFIG_H" >>bla
	echo "#endif" >> bla
	cat $(srcdir)/$(ConfigHeader) >> bla
	$(install_sh_DATA) bla $(DESTDIR)$(includecoindir)/$(ConfigHeader)
	rm -f bla

uninstall-local:
	rm -f $(DESTDIR)$(includecoindir)/$(ConfigHeader)
