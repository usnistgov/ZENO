#!/usr/bin/python3

import sys
import os
import argparse

sys.path.append('..') # for ../validation_common.py

from validation_common import generate_configs

arg_parser = argparse.ArgumentParser()

arg_parser.add_argument("-o", "--output-config-dir", required=True,
                        help="Directory into which to write the generated "
                        "config files")

arg_parser.add_argument("-i", "--input-bod-dir", required=True,
                        help="Directory containing bod files from which to "
                        "generate the config files")

arg_parser.add_argument("-r", "--raw-output-dir", required=True,
                        help="Directory to use for raw output in the config "
                        "files")

args = arg_parser.parse_args()

config_files_dir = args.output_config_dir

bod_files_dir = args.input_bod_dir

raw_output_dir = args.raw_output_dir

num_tests = 50

num_walks_list = [1E6, 1E7]

os.makedirs(config_files_dir, exist_ok = True)
os.makedirs(bod_files_dir, exist_ok = True)
os.makedirs(raw_output_dir, exist_ok = True)

print("Generating configs...")
generate_configs(config_files_dir,
                 bod_files_dir,
                 raw_output_dir,
                 num_tests,
                 num_walks_list)
