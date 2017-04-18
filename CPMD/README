This directory contains the first working version of the MiMiC
communication library, to be linked to both GROMACS and CPMD to
facilitate a QM/MM calculation. However, no released version of either
code yet supports such linking. The source code is found in
the src directory.

Building the library requires a C++11 compiler, MPI2.0 compliant MPI library
and cmake version 2.6 or higher. To build, run

mkdir build
cd build
cmake ..
make

To use the library one will need to link against the built library
including the MessageApi.h header

To build tests:
mkdir build
cd build
cmake -DINCLUDE_TESTS=ON ..
make
ctest

GTest library will be downloaded and compiled in order to run unit tests