## BUILDING AND INSTALLING (SERIAL VERSION)

These instructions are for building and installing ALPS from source. For
instructions on how to obtain pre-built binaries, please see the README file.

IMPORTANT: The build instructions have changed significantly. In most case,
you do not need to clone this repository first! Please follow the instructions
for your platform below.

### Building on Linux

Most Linux distributions come with all the required tools installed. To obtain
the source code, the first step is to get the installer that will then
fetch the source for ALPS and all its dependencies. *You do not need to
clone this repository first, just do the following!* Open a terminal and execute

```
git clone https://www.github.com/coin-or/COIN-OR-OptimizationSuite
```

Next, to check out source code for and build all the necessary projects
(including dependencies), execute the script in the `COIN-OR-OptimizationSuite`
subdirectory. To execute the script, do

```
cd COIN-OR-OptimizationSuite
chmod u+x coin.install.sh
./coin.install.sh
```

(Note: The `chmod` command is only needed if the execute permission is not
automatically set by git on cloning). Once you run the script,
you will be prompted interactively to select a project to fetch and build. The
rest should happen automagically. Alternatively, the following command-line
incantation will execute the procedure non-interactively.

```
./coin.install.sh fetch build --no-prompt --main-proj=Bcps
```

Options that would have been passed to the `configure` script under the old
build system can simply be added to the command-line. For example, to build
with debugging symbols, do

```
./coin.install.sh fetch build --no-prompt --main-proj=Bcps --enable-debug
```

To get help with additional options available in running the script, do

```
./coin/install.sh --help
```

The above procedures will build all required dependencies and Alps itself.
Afterwards, the binaries will be installed in the directory `Mibs/build/bin`
and the libraries in the directory `Bcps/build/lib`. If you wish to
install in a different directory, such as `/usr/local`, then run the command

```
./coin.install.sh install --no-prompt --main-proj=Bcps --prefix=/path/to/install/dir
```

After installation, you will also need to add `/path/to/install/dir/bin` to your
`PATH` variable in your `.bashrc` and also add `/path/to/install/dir/lib`
to your `LD_LIBRARY_PATH` if you want to link to COIN libraries. 

### Building on Windows (MSys2/CYGWIN and MinGW/MSVC)

By far, the easiest way to build on Windows is with the GNU autotools and the
GCC compilers. The first step is to install either
   * [Msys2](https://msys2.github.io/)
   * [CYGWIN](http://cygwin.org/)
   * [Windows Subsystem for Linux](https://docs.microsoft.com/en-us/windows/wsl/install-win10)

If you don't already have CYGWIN installed and don't want to fool around with
WSL (which is a great option if you already know your way around Unix), it is
recommended to use MSys2, since it provides a minimal toolset that is easy to
install. To get MSys2, either download the installer
[here](https://msys2.github.io/) or download and unzip MSys2 base from
[here](http://kent.dl.sourceforge.net/project/msys2/Base/x86_64/msys2-base-x86_64-20150512.tar.xz) 
(this is an out-of-date version, there may be a better place to get an archive
version). 

Following any of the above steps, you should have the `bash` command
(with Msys2, be sure to run `msys2_shell.bat` 
or manually add `msys64\usr\bin`, `msys64\mingw32\bin`, and
`msys64\mingw64\bin` to your Windows path).   

Once you have bash installed and in your `PATH`, open a Windows terminal and
type 

```
bash
pacman -S make wget tar patch dos2unix diffutils git svn
```

To obtain the source code, the first step is to get the installer that will then
fetch the source for Alps and all its dependencies. *You do not need to
clone Alps first, just do the following!* Open a terminal and execute

```
git clone https://www.github.com/coin-or/COIN-OR-OptimizationSuite
```

Next, to check out source code for and build all the necessary projects
(including dependencies), execute the script in the `COIN-OR-OptimizationSuite`
subdirectory. To execute the script, do

```
cd COIN-OR-OptimizationSuite
chmod u+x coi.install.sh
./coin.install.sh
```

(Note: The `chmod` command is only needed if the execute permission is not
automatically set by git on cloning). Once you run the script,
you will be prompted interactively to select a project to fetch and build. the
rest should happen automagically. Alternatively, the following command-line
incantation will execute the procedure non-interactively.

```
./coin.install.sh fetch build --no-prompt --main-proj=Bcps
```
Options that would have been passed to the `configure` script under the old
build system can simply be added to the command-line. For example, to build
with debugging symbols, do

```
./coin.install.sh fetch build --no-prompt --main-proj=Bcps --enable-debug
```

To get help with additional options available in running the script, do

```
./coin/install.sh --help
```

To use the resulting binaries and/or libraries, you will need to add the
full path of the directory `build\bin` to your Windows executable
search `PATH`, or, alternatively, copy the conents of the build directory to 
`C:\Program Files (x86)\Bcps` and add the directory
`C:\Program Files (x86)\Bcps\bin` 
to your Windows executable search `PATH`. You may also consider adding
`C:\Program Files (x86)\Bcps\lib` to the `LIB` path and 
`C:\Program Files (x86)\Bcps\include` to the `INCLUDE` path. 

It is possible to use almost the exact same commands to build with the Visual
Studio compilers. Before doing any of the above commands in the Windows
terminal, first run the `vcvarsall.bat` script for your version of Visual
Studio. Note that you will also need a compatible Fortran compiler if you want
to build any projects requiring Fortran (`ifort` is recommended, but not
free). Then follow all the steps above, but replace the `build` command
with

```
./coin.install.sh fetch build --no-prompt --main-proj=Bcps --enable-msvc
```

### Building with the MSVC++ IDE

These instructions are for MSVC++ Version 10. Instructions for other versions
should be similar. '''The MSVC++ are not regularly tested so please let us
know if they are broken.'''

1. Go to `Bcps/MSVisualStudio/v10` directory and open the solution
file `libBcps.sln`.

2. To test the library, You can build one of the examples by navigating to the
Examples/Knap/MSVisualStudio/v10 directory, for example.


### Building on OS X

OS X is a Unix-based OS and ships with many of the basic components needed to
build COIN-OR, but it's missing some things. For examples, the latest versions
of OS X come with the `clang` compiler but no Fortran compiler. You may also
be missing the `wget` utility and `subversion` and `git` clients (needed for
obtaining source code). The easiest way to get these missing utilitites is to
install Homebrew (see http://brew.sh). After installation, open a terminal and
do

```
brew install gcc wget svn git
```

To obtain the source code, the first step is to get the installer that will then
fetch the source for Alps and all its dependencies. *You do not need to
clone Alps first, just do the following!* Open a terminal and execute

```
git clone https://www.github.com/coin-or/COIN-OR-OptimizationSuite
```

Next, to check out source code for and build all the necessary projects
(including dependencies), execute the script in the `COIN-OR-OptimizationSuite`
subdirectory. To execute the script, do

```
cd COIN-OR-OptimizationSuite
chmod u+x coi.install.sh
./coin.install.sh
```

(Note: The `chmod` command is only needed if the execute permission is not
automatically set by git on cloning). Once you run the script,
you will be prompted interactively to select a project to fetch and build. the
rest should happen automagically. Alternatively, the following command-line
incantation will execute the procedure non-interactively.

```
./coin.install.sh fetch build --no-prompt --main-proj=Bcps
```

With this setup, `clang` will be used for compiling C++ by default and
`gfortran` will be used for Fortran. Since `clang` uses the GNU standard
library, `gfortran` is compatible.

If you want to use the `gcc` compiler provided by Homebrew, then replace the
`build` command above with

```
./coin.install.sh build --no-prompt --main-proj=Bcps CC=gcc-5 CXX=g++-5
```

Options that would have been passed to the `configure` script under the old
build system can simply be added to the command-line. For example, to build
with debugging symbols, do

```
./coin.install.sh fetch build --no-prompt --main-proj=Bcps --enable-debug
```

To get help with additional options available in running the script, do

```
./coin/install.sh --help
```

If you wish to install in a different directory, such as `/usr/local`, then run
the command

```
./coin.install.sh install --no-prompt --main-proj=Bcps --prefix=/path/to/install/dir
```

After installation, you will also need to add `/path/to/install/dir/bin` to your
`PATH` variable in your `.bashrc` and also add `/path/to/install/dir/lib`
to your `DYLD_LIBRARY_PATH` if you want to link to COIN libraries. 

## BUILDING with MPI (PARALLEL VERSION)

For configuration and and compilation of the MPI parallel solver, the user has
to specify the location of MPI with options --with-mpi-incdir, --with-mpi-lib,
MPICC, and MPICXX. Switch to the directory that you just downloaded or
extracted and execute the following commands (substituting the
correct paths and command names).

```
./coin.install.sh build --no-prompt --main-proj=Bcps --enable-static \
--disable-shared --with-mpi-incdir=/usr/local/mpich/include \
--with-mpi-lib='-L/usr/local/mpich/lib  -lmpich' MPICC=mpicc MPICXX=mpiCC

You might want to specify other options, in case you don't want to 
use the default choices that configure makes (e.g., compilers). 

## BUILDING EXAMPLES

To build example codes (ABC and KNAP), configure and build as above. Switch
into the appropriate subdirectory in the source distribution and type "make".

