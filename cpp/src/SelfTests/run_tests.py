#!/usr/bin/python

# ================================================================
#
# This software was developed by employees of the National Institute of
# Standards and Technology (NIST), an agency of the Federal Government.
# Pursuant to title 17 United States Code Section 105, works of NIST employees
# are not subject to copyright protection in the United States and are
# considered to be in the public domain. Permission to freely use, copy,
# modify, and distribute this software and its documentation without fee is
# hereby granted, provided that this notice and disclaimer of warranty appears
# in all copies.
#
# THE SOFTWARE IS PROVIDED 'AS IS' WITHOUT ANY WARRANTY OF ANY KIND, EITHER
# EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT LIMITED TO, ANY
# WARRANTY THAT THE SOFTWARE WILL CONFORM TO SPECIFICATIONS, ANY IMPLIED
# WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, AND FREEDOM
# FROM INFRINGEMENT, AND ANY WARRANTY THAT THE DOCUMENTATION WILL CONFORM TO
# THE SOFTWARE, OR ANY WARRANTY THAT THE SOFTWARE WILL BE ERROR FREE. IN NO
# EVENT SHALL NIST BE LIABLE FOR ANY DAMAGES, INCLUDING, BUT NOT LIMITED TO,
# DIRECT, INDIRECT, SPECIAL OR CONSEQUENTIAL DAMAGES, ARISING OUT OF,
# RESULTING FROM, OR IN ANY WAY CONNECTED WITH THIS SOFTWARE, WHETHER OR NOT
# BASED UPON WARRANTY, CONTRACT, TORT, OR OTHERWISE, WHETHER OR NOT INJURY WAS
# SUSTAINED BY PERSONS OR PROPERTY OR OTHERWISE, AND WHETHER OR NOT LOSS WAS
# SUSTAINED FROM, OR AROSE OUT OF THE RESULTS OF, OR USE OF, THE SOFTWARE OR
# SERVICES PROVIDED HEREUNDER.
#
# Distributions of NIST software should also include copyright and licensing
# statements of any third-party software that are legally bundled with the
# code in compliance with the conditions of those licenses.
# 
# ================================================================

# ================================================================
# 
# Authors: Derek Juba <derek.juba@nist.gov>
# Created: 2017-06-07
# 
# ================================================================

import argparse
import csv
import subprocess
import glob
import os

allowed_relative_error = 0.01

absolute_error_to_ignore = 1E-15

properties_to_ignore = ["input_file",
                        "initialization_ram",
                        "loading_input_data_ram",
                        "building_spatial_data_structure_ram",
                        "walk_on_spheres_ram",
                        "interior_samples_ram",
                        "initialize_time",
                        "read_time",
                        "broadcast_time",
                        "preprocess_time",
                        "exterior_walk_time",
                        "exterior_reduce_time",
                        "volume_sample_time",
                        "volume_reduce_time",
                        "total_time"]

arg_parser = argparse.ArgumentParser()

arg_parser.add_argument("zeno_executable", type=str,
                        help="Zeno executable")

arg_parser.add_argument("input_directory", type=str,
                        help="Input .bod file directory")

arg_parser.add_argument("ground_directory", type=str,
                        help="Ground truth directory")

arg_parser.add_argument("result_directory", type=str,
                        help="Test results directory")

arg_parser.add_argument("--run-mpi-tests", action="store_true", default=False,
                        help="Run MPI tests")

args = arg_parser.parse_args()

zeno_executable = args.zeno_executable

ground_directory = args.ground_directory
result_directory = args.result_directory
bod_directory    = args.input_directory

run_mpi_tests = args.run_mpi_tests

bod_filenames = glob.glob(os.path.join(bod_directory, "*.bod"))

dev_null = open(os.devnull, "w")

def check_result(result_filename):
    ground_basename = os.path.basename(result_filename).replace("result",
                                                                "ground")
    
    ground_filename = os.path.join(ground_directory, ground_basename)
    
    test_passed = True

    with open(ground_filename) as ground_file, \
         open(result_filename) as result_file:

        ground_reader = csv.reader(ground_file)
        result_reader = csv.reader(result_file)

        for ground_line, result_line in map(None, ground_reader, result_reader):

            #Check if one file ended before the other
            if ground_line == None or result_line == None:
                print("Ground truth file and test results file are different "
                      "lengths")
                test_passed = False

            #Check if property name and type are the same for both lines
            elif ground_line[0:2] != result_line[0:2]:
                print(ground_line[0:2])
                print(result_line[0:2])
                test_passed = False

            else:
                if ground_line[0] in properties_to_ignore:
                    continue
                
                if ground_line[1] == "units":
                    if ground_line[2] != result_line[2]:
                        print(ground_line[0:3])
                        print(result_line[0:3])
                        test_passed = False
                else:
                    ground_float = float(ground_line[2])
                    result_float = float(result_line[2])

                    absolute_error = abs(ground_float - result_float)

                    if absolute_error > absolute_error_to_ignore:
                        if ground_float != 0:
                            relative_error = ((ground_float - result_float) /
                                              ground_float)

                            if abs(relative_error) > allowed_relative_error:
                                print(ground_line[0:3])
                                print(result_line[0:3])
                                test_passed = False
                        else:
                            if result_float != 0:
                                print(ground_line[0:3])
                                print(result_line[0:3])
                                test_passed = False

    if test_passed:
        print("PASS\n")

        return True
    else:
        print("*** FAIL ***\n")

        return False

print("Running tests...\n")

all_tests_passed = True

for bod_filename in sorted(bod_filenames):
    test_name = os.path.splitext(os.path.basename(bod_filename))[0]
    
    result_prefix = os.path.join(result_directory, test_name)
    
    serial_result_filename  = result_prefix + "_serial_result.csv"
    threads_result_filename = result_prefix + "_threads_result.csv"
    mpi_result_filename     = result_prefix + "_mpi_result.csv"

    common_params = ["--input-file", bod_filename,
                     "--num-walks", "10000",
                     "--num-interior-samples", "10000",
                     "--max-rsd-capacitance", "1",
                     "--max-rsd-polarizability", "1",
                     "--max-rsd-volume", "1",
                     "--seed", "0",
                     "--print-counts",
                     "--print-benchmarks"]

    serial_params = [zeno_executable,
                     "--csv-output-file", serial_result_filename,
                     "--num-threads", "1"] + common_params

    threads_params = [zeno_executable,
                      "--csv-output-file", threads_result_filename,
                      "--num-threads", "8"] + common_params

    mpi_params = ["mpirun",
                  "-np", "2",
                  zeno_executable,
                  "--csv-output-file", mpi_result_filename,
                  "--num-threads", "4"] + common_params

    print(test_name + " serial")
    
    subprocess.call(serial_params, stdout=dev_null, cwd=bod_directory)

    serial_test_passed = check_result(serial_result_filename)

    all_tests_passed = all_tests_passed and serial_test_passed

    print(test_name + " threads")
    
    subprocess.call(threads_params, stdout=dev_null, cwd=bod_directory)

    threads_test_passed = check_result(threads_result_filename)

    all_tests_passed = all_tests_passed and threads_test_passed
    
    if run_mpi_tests:
        print(test_name + " mpi")

        subprocess.call(mpi_params, stdout=dev_null, cwd=bod_directory)

        mpi_test_passed = check_result(mpi_result_filename)

        all_tests_passed = all_tests_passed and mpi_test_passed

if all_tests_passed:
    print("ALL TESTS PASSED")
else:
    print("*** TESTS FAILED ***")
    
