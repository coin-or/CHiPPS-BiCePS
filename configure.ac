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
# Copyright (C) 2001-2017, Lehigh University, Yan Xu, and Ted Ralphs.      #
# All Rights Reserved.                                                     #
#==========================================================================#

## $Id$

#############################################################################
#                       Names and other basic things                        #
#############################################################################

AC_INIT([Bcps],[0.95.0],[https://github.com/coin-or/CHiPPS-BiCePS/issues/new],[coin-or-bcps],[https://github.com/coin-or/CHiPPS-BiCePS])

AC_COPYRIGHT([
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
])

# List one file in the package so that the configure script can test
# whether the package is actually there
AC_CONFIG_SRCDIR(src/BcpsModel.h)

# Do some initialization work (version numbers, change prefix default, ...)
AC_COIN_INITIALIZE(9:0:8)

#############################################################################
#                         Standard build tool stuff                         #
#############################################################################

# Get the name of the C++ compiler and appropriate compiler options.
AC_COIN_PROG_CXX

# This is a C++ package, set the language accordingly.
AC_LANG_PUSH(C++)

# Initialize libtool
AC_COIN_PROG_LIBTOOL

# set RPATH_FLAGS to the compiler link flags required to hardcode location
# of the shared objects (expanded_libdir is set somewhere in configure before)
# (use in examples Makefile)
AC_COIN_RPATH_FLAGS([$expanded_libdir])

#############################################################################
#                              COIN components                              #
#############################################################################

AC_COIN_CHK_PKG(Alps, [BcpsLib])
if test $coin_has_alps != yes ; then
  AC_MSG_ERROR([Required package Alps not available.])
fi

# These are only required for the unit test and the Abc example.
AC_COIN_CHK_PKG(Osi, [UnitTest])
AC_COIN_CHK_PKG(Clp, [UnitTest], [osi-clp])
AC_COIN_CHK_PKG(Cgl, [UnitTest])

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

#############################################################################
#                  Check for doxygen                                        #
#############################################################################

AC_COIN_DOXYGEN([CoinUtils Alps])

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
                 bcps.pc])

AC_CONFIG_FILES([doxydoc/doxygen.conf])

# Here put the location and name of the configuration header file
AC_CONFIG_HEADERS([src/config.h src/config_bcps.h])

# Finalize libs and cflags.
AC_COIN_FINALIZE_FLAGS([BcpsLib UnitTest])

# Finally, we let configure write all the output...
AC_COIN_FINALIZE

