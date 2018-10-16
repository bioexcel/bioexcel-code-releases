#!/bin/bash

# This script is a convenience for making BioExcel releases.
#
# The resulting tarballs are aggregations of materials suitable for
# BioExcel to release. Their naming reflects either an underlying
# version number, or the date from which they were collated from
# respective git repositories. Details are found in the respective
# README files.
#
# The sha1sum and sha256sum for each tarball is computed and reported
# so we can put that on the BioExcel webpage easily.

# Make tarballs for each component
tar cfz BioExcel-GROMACS-2019-beta1-1.tar.gz GROMACS
tar cfz BioExcel-HADDOCK-20181016-1.tar.gz HADDOCK
tar cfz BioExcel-CPMD-MiMiC-1.tar.gz CPMD/MiMiC

# Make one tarball for everything
tar cfz BioExcel-Software-Release-1.tar.gz GROMACS HADDOCK CPMD

echo "sha1sum values"
for tarballs in *.tar.gz
do
    sha1sum $tarballs
done

echo
echo "sha256sum values"
for tarballs in *.tar.gz
do
    sha256sum $tarballs
done
