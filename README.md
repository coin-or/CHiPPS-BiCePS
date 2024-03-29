# BiCePS

[![A COIN-OR Project](https://coin-or.github.io/coin-or-badge.png)](https://www.coin-or.org)

Projects such as this one are maintained by a small group of volunteers under
the auspices of the non-profit [COIN-OR Foundation](https://www.coin-or.org)
and we need your help! Please consider [sponsoring our
activities](https://github.com/sponsors/coin-or) or [volunteering](mailto:volunteer@coin-or.org) to help!

[![Latest Release](https://img.shields.io/github/v/release/coin-or/CHiPPS-BiCePS?sort=semver)](https://github.com/coin-or/CHiPPS-BiCePS/releases)

_This file is auto-generated from [config.yml](.coin-or/config.yml) using the 
[generate_readme](.coin-or/generate_readme) script.
To make changes, please edit [config.yml](.coin-or/config.yml) or the generation scripts
[here](.coin-or/generate_readme) and [here](https://github.com/coin-or/coinbrew/blob/master/scripts/generate_readme)._

Welcome to BiCePS, (Branch, Constrain, and Price Software), a
data-handling layer of the CHiPPS (COIN-OR HIgh Performance Parallel Search) 
framework. CHiPPS is a framework for implementing parallel graph search 
algorithms. Its methodology generalizes many of the notions of an LP-based 
branch-and-bound algorithm, allowing the implementation of a wide range of 
algorithms with a simplified interface. BiCePS implements the data-handling 
methods required for implementing relaxation-based branch-and-bound. It is an 
intermediate layer of the CHiPPS library hierarchy that that includes a 
library for solving mixed integer linear programs (BLIS).

BiCePS is written in C++ and is released as open source under the [Eclipse Public License 2.0](http://www.opensource.org/licenses/eclipse-2.0).

It is distributed under the auspices of the [COIN-OR Foundation](https://www.coin-or.org)

The BiCePS development site is https://github.com/coin-or/CHiPPS-BiCePS.

## CITE

Code: [![DOI](https://zenodo.org/badge/23726997.svg)](https://zenodo.org/badge/latestdoi/23726997)

Paper: http://dx.doi.org/10.1023/B:SUPE.0000020179.55383.ad

## CURRENT BUILD STATUS

[![Windows Builds](https://github.com/coin-or/CHiPPS-BiCePS/actions/workflows/windows-ci.yml/badge.svg?branch=master)](https://github.com/coin-or/CHiPPS-BiCePS/actions/workflows/windows-ci.yml?query=branch%3Amaster)

[![Linux and MacOS Builds](https://github.com/coin-or/CHiPPS-BiCePS/actions/workflows/linux-ci.yml/badge.svg?branch=master)](https://github.com/coin-or/CHiPPS-BiCePS/actions/workflows/linux-ci.yml?query=branch%3Amaster)

## DOWNLOAD

### Docker image

There is a Docker image that provides BiCePS, as well as other projects
in the [COIN-OR Optimization
Suite](https://github.com/coin-or/COIN-OR-OptimizationSuite) [here](https://hub.docker.com/repository/docker/coinor/coin-or-optimization-suite)

### Binaries

For newer releases, binaries will be made available as assets attached to
releases in Github
[here](https://github.com/coin-or/CHiPPS-BiCePS/releases). Older binaries
are archived as part of CHiPPS-BLIS
[here](https://www.coin-or.org/download/binary/CHiPPS-BLIS).

Due to license incompatibilities, pre-compiled binaries may lack some
functionality. If binaries are not available for your platform for the latest
version and you would like to request them to be built and posted, feel free
to let us know in the discussion formum.

### Source

Source code can be obtained either by

 * Downloading a snapshot of the source code for the latest release version of BiCePS from the
 [releases](https://github.com/coin-or/CHiPPS-BiCePS/releases) page.
 * Cloning this repository from [Github](https://github.com/coin-or/CHiPPS-BiCePS) or 
 * Using the [coinbrew](https://github.com/coin-or/coinbrew) script to get the project and all dependencies (recommended, see below).   

Below is a quick start guide for building on common platforms. More detailed
build instructions are
[here](https://coin-or.github.io/user_introduction.html).

### Dependencies

BiCePS has a number of dependencies, which are detailed in
[config.yml](.coin-or/config.yml). Dependencies on other COIN-OR projects are
automatically downloaded when obtaining the source with `coinbrew`. For some
of the remaining third-party dependencies, automatic download scripts and
build wrappers are provided (and will also be automatically run for required
and recommended dependencies), while other libraries that are aeasy to obtain
must be installed using an appropriate package manager (or may come with your
OS by default). 

## BUILDING from source

The quick start assumes you are in a bash shell. 

### Using `coinbrew`

To download and build BiCePS from source, execute the 
following on the command line. 
```
wget https://raw.githubusercontent.com/coin-or/coinbrew/master/coinbrew
chmod u+x coinbrew
./coinbrew fetch Bcps@master
./coinbrew build Bcps
```
For more detailed instructions on coinbrew, see https://coin-or.github.io/coinbrew.
The `coinbrew` script will fetch the additional projects specified in the Dependencies section of [config.yml](.coin-or/config.yml).

### Without `coinbrew` (Expert users)

 * Download the source code, e.g., by cloning the git repo https://github.com/coin-or/CHiPPS-BiCePS
 * Download and install the source code for the dependencies listed in [config.yml](.coin-or/config.yml)
 * Build the code as follows (make sure to set PKG_CONFIG_PTH to install directory for dependencies).

```
./configure -C
make
make test
make install
```

## BUILDING with MPI (PARALLEL VERSION)

For configuration and compilation of the parallel version, the user has
to specify the location of MPI with options `--with-mpi-cflags`,
`--with-mpi-lflags`, `MPICC`, and `MPICXX`. 

```
./coinbrew build Bcps --enable-static --disable-shared --with-mpi-cflags="$\(pkg-config --cflags mpi\)" --with-mpi-lflags="$\(pkg-config --libs mpi\)" MPICC=mpicc MPICXX=mpiCC
```

## BUILDING EXAMPLES

To build the example codes (which is an old version of Blis, just used for 
testing), configure and build as above. Switch into the appropriate
subdirectory in the source distribution and type `make`.
## Doxygen Documentation

If you have `Doxygen` available, you can build a HTML documentation by typing

`make doxygen-docs` 

in the build directory. If BiCePS was built via `coinbrew`, then the build
directory will be `./build/CHiPPS-BiCePS/master` by default. The doxygen documentation main file
is found at `<build-dir>/doxydoc/html/index.html`.

If you don't have `doxygen` installed locally, you can use also find the
documentation [here](http://coin-or.github.io/CHiPPS-BiCePS/Doxygen).


## Project Links

 * [Code of Conduct](https://www.coin-or.org/code-of-conduct/)
 * [COIN-OR Web Site](http://www.coin-or.org/)
 * [COIN-OR general discussion forum](https://github.com/orgs/coin-or/discussions)
 * [BiCePS Discussion forum](https://github.com/coin-or/CHiPPS-BiCePS/discussions)
 * [Report a bug](https://github.com/coin-or/CHiPPS-BiCePS/issues/new)
 * [Doxygen generated documentation](http://coin-or.github.io/CHiPPS-BiCePS/Doxygen)

## CURRENT TESTING STATUS

  1. Configurations
    - Serial: Well tested.
    - LAMMPI: Well tested.
    - MPICH: Well tested.

  2. Applications (See INSTALL)
    - Blis: an older version of the Blis solver: Well tested.

## Authors

Source Code:

Yan Xu (yax2@lehigh.edu)
Aykut Bulut (aykutblt@gmail.com) 
Ted Ralphs (ted@lehigh.edu), Project Manager

Original Conceptual Design:

Yan Xu (yax2@lehigh.edu)
Ted Ralphs (ted@lehigh.edu), Project Manager
Laci Ladanyi (ladanyi@us.ibm.com)
Matt Saltzman (mjs@clemson.edu)
