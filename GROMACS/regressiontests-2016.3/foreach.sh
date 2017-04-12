#!/bin/bash

# This is just a simple wrapper script to make it easy to do the same
# simple operation on many subdirectories, e.g.
#
# foreach.sh "mv mdout.mdp grompp.mdp" kernel/kernel* complex/tip4p* simple/*
#

command=$1
shift

cwd=$(pwd)

for dir in $*
do
  cd $dir
  $command;
  cd $cwd
done
