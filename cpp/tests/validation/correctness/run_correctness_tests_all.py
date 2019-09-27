#!/usr/bin/python3

import sys
import os
import argparse

sys.path.append('..') # for ../validation_common.py

from validation_common import run_tests

arg_parser = argparse.ArgumentParser()

arg_parser.add_argument("-i", "--input-config-dir", required=True,
                        help="Directory containing config files of the tests "
                        "to run")

arg_parser.add_argument("-z", "--zeno-exec-path", required=True,
                        help="Path to the zeno executable file")

args = arg_parser.parse_args()

config_files_dir = args.input_config_dir

zeno_exec = args.zeno_exec_path

for test_dir in sorted(os.listdir(config_files_dir)):
    test_dir_path = os.path.join(config_files_dir,
                                 test_dir)
    
    if os.path.isdir(test_dir_path):
        run_tests(test_dir_path,
                  zeno_exec)
