#!/usr/bin/python3

import os
import argparse
import subprocess

arg_parser = argparse.ArgumentParser()

arg_parser.add_argument("-i", "--input-combined-dir", required=True,
                        help="Directory containing the combined output files "
                        "from the surrogate ground truth computation")

arg_parser.add_argument("-o", "--output-formatted-dir", required=True,
                        help="Directory into which to write the formatted "
                        "surrogate ground truth data files")

args = arg_parser.parse_args()

combined_output_dir = args.input_combined_dir

formatted_output_dir = args.output_formatted_dir

for combined_output_file in sorted(os.listdir(combined_output_dir)):

    combined_output_file_path = os.path.join(combined_output_dir,
                                             combined_output_file)

    formatted_output_file_path = os.path.join(formatted_output_dir,
                                              combined_output_file)

    format_surrogate_ground_truth_output_command = [
        "./format_surrogate_ground_truth_output.py",
        "--input-combined-file",
        combined_output_file_path,
        "--output-formatted-file",
        formatted_output_file_path
    ]

    print(format_surrogate_ground_truth_output_command)
    subprocess.call(format_surrogate_ground_truth_output_command)
