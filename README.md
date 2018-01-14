BiCePS Version 0.99 README
==========================

Welcome to BiCePS. BiCePS is distributed under the Eclipse Public License and
is freely redistributable. All source code and documentation is Copyright
2001-2018 by Lehigh University, Yan Xu, Aykut Bulut, and Ted Ralphs. This
README may be distributed freely.

TESTING STATUS
==============

[![Build Status](https://travis-ci.org/coin-or/CHiPPS-BiCePS.svg?branch=master)](https://travis-ci.org/coin-or/CHiPPS-BiCePS)

[![Build status](https://ci.appveyor.com/api/projects/status/lqdk3l72xs5j8qxj?svg=true)](https://ci.appveyor.com/project/tkralphs/chipps-biceps)

DOWNLOAD
========

Binaries available as part of the [CHiPPS-BLIS](http://github.com/coin-or/CHiPPS-BLIS) binary distribution on [BinTray](http://bintray.com):

[ ![Download](https://api.bintray.com/packages/coin-or/download/CHiPPS-BLIS/images/download.svg?version=0.94) ](https://bintray.com/coin-or/download/CHiPPS-BLIS/0.94/link)

CITE
====

[![DOI](https://zenodo.org/badge/23726997.svg)](https://zenodo.org/badge/latestdoi/23726997)

WHAT IS BiCePS?
===============

BiCePS is one of the libraries that make up the CHiPPS (COIN High Performance
Parallel Search Framework) library hierarchy. It implements that data-handling
functions needed to support development of many types of relaxation-based
branch-and-bound algorithms, especially for solving mathematical programs. It
is intended to capture the implementation of methods common to all such
algorithms without assuming anything about the structure of the mathematical
program or the bounding method used. BLIS, which is another layer built on top
of BiCePS, is a concretization of Bcps for the case of mixed integer linear
programs. DIP is another implementation being developed using BLIS that
implements a decomposition-based bounding procedure.

PROJECT WEB PAGE
================

https://projects.coin-or.org/CHiPPS

DOCUMENTATION
=============

The html documentation of the classes in BiCePS (See INSTALL file) can be
created. A user's guide that is under preparation that can help with running
or customizing BiCePS.

INSTALLATION
============

Please see the INSTALL file for a guide to install Bcps.

SUPPORT
=======

1. Authors

Source Code:

Yan Xu (yax2@lehigh.edu)
Aykut Bulut (aykutblt@gmail.com)
Ted Ralphs (ted@lehigh.edu), Project Manager

Original Conceptual Design:

Yan Xu (yax2@lehigh.edu)
Ted Ralphs (ted@lehigh.edu), Project Manager
Laci Ladanyi (ladanyi@us.ibm.com)
Matt Saltzman (mjs@clemson.edu)

2. Bug Reports and Other Issues

Bug reports should be reported on the CHiPPS GitHub web site at

https://github.com/coin-or/CHiPPS-BLIS/issues/new
