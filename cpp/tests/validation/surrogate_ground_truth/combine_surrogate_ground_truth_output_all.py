#!/usr/bin/python3

import sys
import os
import argparse

sys.path.append('..') # for ../validation_common.py

from validation_common import combine_output

arg_parser = argparse.ArgumentParser()

arg_parser.add_argument("-i", "--input-raw-dir", required=True,
                        help="Directory containing the raw output files to"
                        "combine")

arg_parser.add_argument("-o", "--output-combined-dir", required=True,
                        help="Directory into which to write the combined "
                        "output files")

args = arg_parser.parse_args()

raw_output_dir = args.input_raw_dir

parsed_output_dir = args.output_combined_dir

for output_csv_files_dir in sorted(os.listdir(raw_output_dir)):
    
        output_csv_files_dir_path = os.path.join(raw_output_dir,
                                                 output_csv_files_dir)

        combined_csv_file_name = output_csv_files_dir + ".csv"

        combined_csv_file_path = os.path.join(parsed_output_dir,
                                              combined_csv_file_name)

        if os.path.isdir(output_csv_files_dir_path):
            combine_output(output_csv_files_dir_path,
                           combined_csv_file_path)
