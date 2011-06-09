#==========================================================================#
# This file is part of the Branch, Constrain and Price Software (BiCePS).  #
#                                                                          #
# BiCePS is distributed under the Eclipse Public License as part of the    #
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
# Copyright (C) 2001-2011, Lehigh University, Yan Xu, and Ted Ralphs.      #
# All Rights Reserved.                                                     #
#==========================================================================#

## $Id$

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
AM_CPPFLAGS = $(BCPSLIB_CFLAGS)

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

install-exec-local:
	$(install_sh_DATA) config_bcps.h $(DESTDIR)$(includecoindir)/BcpsConfig.h

uninstall-local:
	rm -f $(DESTDIR)$(includecoindir)/BcpsConfig.h
