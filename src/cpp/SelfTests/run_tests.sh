#!/bin/bash

# ================================================================
# 
# Disclaimer:  IMPORTANT:  This software was developed at the National
# Institute of Standards and Technology by employees of the Federal
# Government in the course of their official duties.  Pursuant to
# title 17 Section 105 of the United States Code this software is not
# subject to copyright protection and is in the public domain.  This
# is an experimental system.  NIST assumes no responsibility
# whatsoever for its use by other parties, and makes no guarantees,
# expressed or implied, about its quality, reliability, or any other
# characteristic.  We would appreciate acknowledgement if the software
# is used.  This software can be redistributed and/or modified freely
# provided that any derivative works bear some notice that they are
# derived from it, and any modified versions bear some notice that
# they have been modified. 
# 
# ================================================================

# ================================================================
# 
# Author: Derek Juba <derek.juba@nist.gov>
# Date:   Wed Sep 14 11:30:15 2016 EDT
# 
# Time-stamp: <2016-09-22 11:57:18 dcj>
# 
# ================================================================

function compare {
    float_re="[[:space:]][-+]\?[[:digit:]]\+[.][[:digit:]]\+\([eE][-+]\?[[:digit:]]\+\)\?\([[:space:]]\|\$\)"

    int_re="[[:space:]][-+]\?[[:digit:]]\+\([[:space:]]\|\$\)"

    exclude_re="Time"

    match_tolerance=0.01

    paste <(grep -v $exclude_re $1 | grep -o $int_re) <(grep -v $exclude_re $2 | grep -o $int_re) | awk '{if ($1 != $2) exit(10)}'

    int_code=$?

    # Don't test attributes ground truth == 0
    paste <(grep -v $exclude_re $1 | grep -o $float_re) <(grep -v $exclude_re $2 | grep -o $float_re) | awk '{if (($1 != 0) && ((($1 - $2)/$1)^2 > '$match_tolerance'^2)) exit(10)}'

    float_code=$?

    if [ $int_code == 0 ] && [ $float_code == 0 ]; then
	echo "PASS"
    else
	echo "*** FAIL ***"
    fi
}

common_params="--num-walks 10000 --num-interior-samples 10000 --seed 0 --print-counts"

function run_test {
    test_file=$1_$2_test.txt
    ground_file=$1_$2_ground.txt

    echo "Testing" $test_file $ground_file "..."

    ../zeno -i $1.bod --num-threads $3 $common_params > $test_file

    compare $test_file $ground_file
}

function run_mpi_test {
    test_file=$1_$2_test.txt
    ground_file=$1_$2_ground.txt

    echo "Testing" $test_file $ground_file "..."

    mpirun -np $4 ../zeno-mpi -i $1.bod --num-threads $3 $common_params > $test_file

    compare $test_file $ground_file
}

run_test unit_cm serial 1
run_test two_spheres_1_1 serial 1
run_test two_spheres_1_4 serial 1
run_test torus_1_4 serial 1
run_test polymer serial 1
run_test 1LYD serial 1

run_test unit_cm threads 8
run_test two_spheres_1_1 threads 8
run_test two_spheres_1_4 threads 8
run_test torus_1_4 threads 8
run_test polymer threads 8
run_test 1LYD threads 8

if [ -f ../zeno-mpi ]; then
    run_mpi_test unit_cm mpi 2 4
    run_mpi_test two_spheres_1_1 mpi 2 4
    run_mpi_test two_spheres_1_4 mpi 2 4
    run_mpi_test torus_1_4 mpi 2 4
    run_mpi_test polymer mpi 2 4
    run_mpi_test 1LYD mpi 2 4
fi

# ================================================================

# Local Variables:
# time-stamp-line-limit: 30
# End:
