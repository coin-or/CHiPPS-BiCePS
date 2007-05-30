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
# Copyright (C) 2001-2007, Lehigh University, Yan Xu, and Ted Ralphs.      #
# All Rights Reserved.                                                     #
#==========================================================================#

## $Id$

# Author:  Andreas Waechter           IBM    2006-04-13
#          Yan Xu                     Lehigh    2006-08-27

#############################################################################
#                       Names and other basic things                        #
#############################################################################

AC_PREREQ(2.59)

AC_INIT([Bcps],[0.8],[https://projects.coin-or.org/Alps/newticket])

AC_COPYRIGHT([
Copyright 2006 International Business Machines and others.
All Rights Reserved.
This file is part of the open source package Coin which is distributed
under the Common Public License.])

# List one file in the package so that the configure script can test
# whether the package is actually there
AC_CONFIG_SRCDIR(src/BcpsModel.h)

# Where should everything be installed by default?  Here, we want it
# to be installed directly in 'bin', 'lib', 'include' subdirectories
# of the directory where configure is run.  The default would be
# /usr/local.
AC_PREFIX_DEFAULT([`pwd`])

#############################################################################
#                         Standard build tool stuff                         #
#############################################################################

# Get the system type
AC_CANONICAL_BUILD

# If this project depends on external projects, the Externals file in
# the source root directory contains definition of where to find those
# externals.  The following macro ensures that those externals are
# retrieved by svn if they are not there yet.
AC_COIN_PROJECTDIR_INIT

# Check if user wants to produce debugging code
AC_COIN_DEBUG_COMPILE(Bcps)

# Get the name of the C++ compiler and appropriate compiler options
AC_COIN_PROG_CXX

# We filter out the -pedantic-errors flag for C and C++ compiler in
# case the user wants to use MPI, since otherwise the compilation
# fails
if test "${with_mpi_lib+set}" = set; then
  CFLAGS=`echo $CFLAGS | sed -e s/-pedantic-errors//g`
  CXXFLAGS=`echo $CXXFLAGS | sed -e s/-pedantic-errors//g`
fi

# Initialize automake and libtool
AC_COIN_INIT_AUTO_TOOLS

#############################################################################
#                              COIN components                              #
#############################################################################

AC_COIN_HAS_PROJECT(Bcps)
AC_COIN_HAS_PROJECT(CoinUtils)
AC_COIN_HAS_PROJECT(Alps)
AC_COIN_HAS_PROJECT(Clp)
AC_COIN_HAS_PROJECT(Cgl)
AC_COIN_HAS_PROJECT(Osi)

#############################################################################
#                             Other dependencies                            #
#############################################################################

AC_COIN_HAS_USER_LIBRARY([Mpi], [MPI], [mpi.h], [MPI_Irecv])
# ToDo: Automatically choose MPI compiler

##############################################################################
#                   VPATH links for example input files                      #
##############################################################################

# In this case, the source code is taken from the Bcps examples directory

AC_CONFIG_LINKS(test/BlisBranchObjectInt.cpp:examples/Blis/BlisBranchObjectInt.cpp
	test/BlisBranchObjectInt.h:examples/Blis/BlisBranchObjectInt.h
	test/BlisBranchStrategyPseudo.cpp:examples/Blis/BlisBranchStrategyPseudo.cpp
	test/BlisBranchStrategyPseudo.h:examples/Blis/BlisBranchStrategyPseudo.h
	test/BlisBranchStrategyRel.cpp:examples/Blis/BlisBranchStrategyRel.cpp
	test/BlisBranchStrategyRel.h:examples/Blis/BlisBranchStrategyRel.h
	test/BlisBranchStrategyStrong.cpp:examples/Blis/BlisBranchStrategyStrong.cpp
	test/BlisBranchStrategyStrong.h:examples/Blis/BlisBranchStrategyStrong.h
	test/BlisConGenerator.cpp:examples/Blis/BlisConGenerator.cpp
	test/BlisConGenerator.h:examples/Blis/BlisConGenerator.h
	test/BlisConstraint.cpp:examples/Blis/BlisConstraint.cpp
	test/BlisConstraint.h:examples/Blis/BlisConstraint.h
	test/Blis.h:examples/Blis/Blis.h
	test/BlisHelp.cpp:examples/Blis/BlisHelp.cpp
	test/BlisHelp.h:examples/Blis/BlisHelp.h
	test/BlisHeuristic.h:examples/Blis/BlisHeuristic.h
	test/BlisHeurRound.cpp:examples/Blis/BlisHeurRound.cpp
	test/BlisHeurRound.h:examples/Blis/BlisHeurRound.h
	test/BlisLicense.h:examples/Blis/BlisLicense.h
	test/BlisMain.cpp:examples/Blis/BlisMain.cpp
	test/BlisMessage.cpp:examples/Blis/BlisMessage.cpp
	test/BlisMessage.h:examples/Blis/BlisMessage.h
	test/BlisModel.cpp:examples/Blis/BlisModel.cpp
	test/BlisModel.h:examples/Blis/BlisModel.h
	test/BlisNodeDesc.h:examples/Blis/BlisNodeDesc.h
	test/BlisObjectInt.cpp:examples/Blis/BlisObjectInt.cpp
	test/BlisObjectInt.h:examples/Blis/BlisObjectInt.h
	test/BlisParams.cpp:examples/Blis/BlisParams.cpp
	test/BlisParams.h:examples/Blis/BlisParams.h
	test/BlisPseudo.cpp:examples/Blis/BlisPseudo.cpp
	test/BlisPseudo.h:examples/Blis/BlisPseudo.h
	test/BlisSolution.h:examples/Blis/BlisSolution.h
	test/BlisSubTree.h:examples/Blis/BlisSubTree.h
	test/BlisTreeNode.cpp:examples/Blis/BlisTreeNode.cpp
	test/BlisTreeNode.h:examples/Blis/BlisTreeNode.h
	test/BlisVariable.h:examples/Blis/BlisVariable.h
	test/flugpl.mps:examples/Blis/data/flugpl.mps)

# In case this is a VPATH configuration we need to make sure that the
# input files for the examples are available in the VPATH directory.

# ToDo
#AC_COIN_VPATH_LINK(examples/Blis/blis.par)
#AC_COIN_VPATH_LINK(examples/input.130)

##############################################################################
#                   Finishing up by writing all the output                   #
##############################################################################

# Here list all the files that configure should create (except for the
# configuration header file)
AC_CONFIG_FILES([Makefile
		 examples/Blis/Makefile
		 examples/Blis/blis.par
                 src/Makefile
                 test/Makefile
		 bcps_addlibs.txt])

# Here put the location and name of the configuration header file
AC_CONFIG_HEADER([inc/config_bcps.h])

# Finally, we let configure write all the output...
AC_COIN_FINALIZE
