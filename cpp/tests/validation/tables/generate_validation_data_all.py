#!/usr/bin/python3

import os
import glob
import subprocess
import argparse
import collections

TestCase = collections.namedtuple("TestCase",
                                  "test_results_file "
                                  "ground_truth_file "
                                  "validation_data_file "
                                  "test_type")

test_cases = [
    TestCase("two_spheres_1_1_1e+06.csv",
             "two_spheres_1_1.csv",
             "two_spheres_1_1_1e+06.csv",
             "analytic"),
    TestCase("two_spheres_1_1_1e+07.csv",
             "two_spheres_1_1.csv",
             "two_spheres_1_1_1e+07.csv",
             "analytic"),
    TestCase("two_spheres_1_4_1e+06.csv",
             "two_spheres_1_4.csv",
             "two_spheres_1_4_1e+06.csv",
             "analytic"),
    TestCase("two_spheres_1_4_1e+07.csv",
             "two_spheres_1_4.csv",
             "two_spheres_1_4_1e+07.csv",
             "analytic"),
    TestCase("voxel_cube_1_1e+06.csv",
             "cube_1_1e+11.csv",
             "cube_1_1e+06.csv",
             "surrogate"),
    TestCase("voxel_cube_1_1e+07.csv",
             "cube_1_1e+11.csv",
             "cube_1_1e+07.csv",
             "surrogate"),
    TestCase("voxel_cuboid_1_2_3_1e+06.csv",
             "cuboid_1_2_3_1e+11.csv",
             "cuboid_1_2_3_1e+06.csv",
             "surrogate"),
    TestCase("voxel_cuboid_1_2_3_1e+07.csv",
             "cuboid_1_2_3_1e+11.csv",
             "cuboid_1_2_3_1e+07.csv",
             "surrogate"),
    TestCase("polymer_1e+06.csv",
             "polymer_1e+11.csv",
             "polymer_1e+06.csv",
             "surrogate"),
    TestCase("polymer_1e+07.csv",
             "polymer_1e+11.csv",
             "polymer_1e+07.csv",
             "surrogate"),
    TestCase("1LYD_1e+06.csv",
             "1LYD_1e+11.csv",
             "1LYD_1e+06.csv",
             "surrogate"),
    TestCase("1LYD_1e+07.csv",
             "1LYD_1e+11.csv",
             "1LYD_1e+07.csv",
             "surrogate"),
    TestCase("spheres_torus_1_4_16_1e+06.csv",
             "spheres_torus_1_4_16_1e+11.csv",
             "spheres_torus_1_4_16_1e+06_surrogate.csv",
             "surrogate"),
    TestCase("spheres_torus_1_4_16_1e+07.csv",
             "spheres_torus_1_4_16_1e+11.csv",
             "spheres_torus_1_4_16_1e+07_surrogate.csv",
             "surrogate"),
    TestCase("spheres_torus_1_4_16_1e+06.csv",
             "torus_1_4.csv",
             "spheres_torus_1_4_16_1e+06_analytic.csv",
             "analytic"),
    TestCase("spheres_torus_1_4_16_1e+07.csv",
             "torus_1_4.csv",
             "spheres_torus_1_4_16_1e+07_analytic.csv",
             "analytic"),
    TestCase("voxels_two_spheres_1_4_32_1e+06.csv",
             "voxels_two_spheres_1_4_32_1e+11.csv",
             "voxels_two_spheres_1_4_32_1e+06_surrogate.csv",
             "surrogate"),
    TestCase("voxels_two_spheres_1_4_32_1e+07.csv",
             "voxels_two_spheres_1_4_32_1e+11.csv",
             "voxels_two_spheres_1_4_32_1e+07_surrogate.csv",
             "surrogate"),
    TestCase("voxels_two_spheres_1_4_32_1e+06.csv",
             "two_spheres_1_4_large.csv",
             "voxels_two_spheres_1_4_32_1e+06_analytic.csv",
             "analytic"),
    TestCase("voxels_two_spheres_1_4_32_1e+07.csv",
             "two_spheres_1_4_large.csv",
             "voxels_two_spheres_1_4_32_1e+07_analytic.csv",
             "analytic")
]

arg_parser = argparse.ArgumentParser()

arg_parser.add_argument("-g", "--generate-script-path", required=True,
                        help="Path to generate_validation_data script")

arg_parser.add_argument("-t", "--test-results-dir", required=True,
                        help="Directory of CSV files containing the combined "
                        "results from all test runs")

arg_parser.add_argument("-s", "--surrogate-ground-truth-dir", required=True,
                        help="Directory of CSV files containing the combined "
                        "results from all surrogate ground-truth runs")

arg_parser.add_argument("-a", "--analytic-ground-truth-dir", required=True,
                        help="Directory of CSV files containing analytic "
                        "ground truth")

arg_parser.add_argument("-v", "--validation-data-dir", required=True,
                        help="Directory into which to write the validation "
                        "table data CSV files")

args = arg_parser.parse_args()

generate_validation_data_script = args.generate_script_path
test_results_dir = args.test_results_dir
surrogate_ground_truth_dir = args.surrogate_ground_truth_dir
analytic_ground_truth_dir = args.analytic_ground_truth_dir
validation_data_dir = args.validation_data_dir

for test_case in test_cases:
    test_results_path = os.path.join(test_results_dir,
                                     test_case.test_results_file)

    if test_case.test_type == "surrogate":
        ground_truth_dir = surrogate_ground_truth_dir
    elif test_case.test_type == "analytic":
        ground_truth_dir = analytic_ground_truth_dir
    else:
        raise Exception("Invalid test type " + test_case.test_type)

    ground_truth_path = os.path.join(ground_truth_dir,
                                     test_case.ground_truth_file)

    validation_data_path = os.path.join(validation_data_dir,
                                        test_case.validation_data_file)

    command = [generate_validation_data_script,
               "-t", test_results_path,
               "-g", ground_truth_path,
               "-v", validation_data_path]

    print(command)

    subprocess.call(command)
