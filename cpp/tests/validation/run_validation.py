#!/usr/bin/python3

import os
import subprocess
import shutil
import argparse

sphinx_dir_path = "../../doc/release/sphinx"

arg_parser = argparse.ArgumentParser()

arg_parser.add_argument("-z", "--zeno-exec-path", required=True,
                        help="Path to the zeno executable file")

args = arg_parser.parse_args()

zeno_exec_path = os.path.abspath(args.zeno_exec_path)

def generate_surrogate_ground_truth():
    generate_surrogate_ground_truth_tests_all_command = [
        "./generate_surrogate_ground_truth_tests_all.py",
        "--output-config-dir",
        "./config_files",
        "--input-bod-dir",
        "./bod_files",
        "--raw-output-dir",
        "./raw_output"
    ]

    print(generate_surrogate_ground_truth_tests_all_command)
    subprocess.call(generate_surrogate_ground_truth_tests_all_command,
                    cwd="./surrogate_ground_truth")

    run_surrogate_ground_truth_tests_all_command = [
        "./run_surrogate_ground_truth_tests_all.py",
        "--input-config-dir",
        "./config_files",
        "--zeno-exec-path",
        zeno_exec_path
    ]

    print(run_surrogate_ground_truth_tests_all_command)
    subprocess.call(run_surrogate_ground_truth_tests_all_command,
                    cwd="./surrogate_ground_truth")

    combine_surrogate_ground_truth_output_all_command = [
        "./combine_surrogate_ground_truth_output_all.py",
        "--input-raw-dir",
        "./raw_output",
        "--output-combined-dir",
        "./combined_output"
    ]

    print(combine_surrogate_ground_truth_output_all_command)
    subprocess.call(combine_surrogate_ground_truth_output_all_command,
                    cwd="./surrogate_ground_truth")

    format_surrogate_ground_truth_output_all_command = [
        "./format_surrogate_ground_truth_output_all.py",
        "--input-combined-dir",
        "./combined_output",
        "--output-formatted-dir",
        "./formatted_output"
    ]

    print(format_surrogate_ground_truth_output_all_command)
    subprocess.call(format_surrogate_ground_truth_output_all_command,
                    cwd="./surrogate_ground_truth")

def run_correctness_tests():
    generate_correctness_tests_all_command = [
        "./generate_correctness_tests_all.py",
        "--output-config-dir",
        "./config_files",
        "--input-bod-dir",
        "./bod_files",
        "--raw-output-dir",
        "./raw_output"
    ]
    
    print(generate_correctness_tests_all_command)
    subprocess.call(generate_correctness_tests_all_command,
                    cwd="./correctness")

    run_correctness_tests_all_command = [
        "./run_correctness_tests_all.py",
        "--input-config-dir",
        "./config_files",
        "--zeno-exec-path",
        zeno_exec_path
    ]

    print(run_correctness_tests_all_command)
    subprocess.call(run_correctness_tests_all_command,
                    cwd="./correctness")

    combine_correctness_output_all_command = [
        "./combine_correctness_output_all.py",
        "--input-raw-dir",
        "./raw_output",
        "--output-combined-dir",
        "./combined_output"
    ]

    print(combine_correctness_output_all_command)
    subprocess.call(combine_correctness_output_all_command,
                    cwd="./correctness")

def generate_tables():
    generate_validation_data_all_command = [ 
        "./generate_validation_data_all.py",
        "--generate-script-path",
        "./generate_validation_data.py",
        "--test-results-dir",
        "../correctness/combined_output",
        "--surrogate-ground-truth-dir",
        "../surrogate_ground_truth/formatted_output",
        "--analytic-ground-truth-dir",
        "../ground_truth/analytic",
        "--validation-data-dir",
        "./csv"
    ]

    print(generate_validation_data_all_command)
    subprocess.call(generate_validation_data_all_command,
                    cwd="./tables")

    convert_csv_to_rst_table_all_command = [
        "./convert_csv_to_rst_table_all.py",
        "--script-path",
        "./convert_csv_to_rst_table.py",
        "--input-csv-dir",
        "./csv",
        "--output-rst-dir",
        "./rst"
    ]

    print(convert_csv_to_rst_table_all_command)
    subprocess.call(convert_csv_to_rst_table_all_command,
                    cwd="./tables")

    copy_source_dir_path = "./tables/rst"
    copy_dest_dir_path = os.path.join(sphinx_dir_path, "validation_tables")

    for copy_source_file_name in os.listdir(copy_source_dir_path):
        copy_source_file_path = os.path.join(copy_source_dir_path,
                                             copy_source_file_name)

        copy_dest_file_name = copy_source_file_name

        copy_dest_file_path = os.path.join(copy_dest_dir_path,
                                           copy_dest_file_name)

        print("Copying " + copy_source_file_path + " to " + copy_dest_file_path)
        shutil.copy2(copy_source_file_path, copy_dest_file_path)

    make_command = ["make", "html"]

    print(make_command)
    subprocess.call(make_command, cwd=sphinx_dir_path)

#generate_surrogate_ground_truth()

run_correctness_tests()

generate_tables()
