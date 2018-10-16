#!/bin/bash
# @ job_name = thegmxtest
# @ step_name = mdrun
# @ output = $(job_name).$(step_name).$(cluster).out
# @ error = $(job_name).$(step_name).$(cluster).err
# @ wall_clock_limit = 0:30:00
# @ notification = error
# @ job_type = bluegene
# @ class = bg32_100
# @ bg_size = 32
# @ group = bg01
# @ environment = COPY_ALL;
# @ queue

# Note $PREFIX needs a final '/' because of how -prefix is
# being abused below
PREFIX=/path/to/gromacs/bin/
#VERBOSE=-verbose
#DOUBLE=-double
BLUEGENE=-bluegene
#TIGHT=-tight

./gmxtest.pl -np 1 kernel ${VERBOSE} ${DOUBLE} -bluegene -prefix ${PREFIX} ${TIGHT}
./gmxtest.pl -np 1 simple ${VERBOSE} ${DOUBLE} -bluegene -prefix ${PREFIX} ${TIGHT}
./gmxtest.pl -np 1 complex ${VERBOSE} ${DOUBLE} -bluegene -prefix ${PREFIX} ${TIGHT}
./gmxtest.pl -np 1 freeenergy ${VERBOSE} ${DOUBLE} -bluegene -prefix ${PREFIX} ${TIGHT}
