#!/usr/bin/python3

import os

from validation import generate_configs

zeno_dir = '../../..'

config_files_dir = os.path.join(zeno_dir,
                                'cpp/tests/validation/'
                                'surrogate_ground_truth/config_files')

bod_files_dir = os.path.join(zeno_dir,
                             'cpp/tests/validation/'
                             'surrogate_ground_truth/bod_files')

raw_output_dir = os.path.join(zeno_dir,
                              'cpp/tests/validation/'
                              'surrogate_ground_truth/raw_output')

num_tests = 3

num_walks_list = [1E11]

os.makedirs(config_files_dir, exist_ok = True)
os.makedirs(bod_files_dir, exist_ok = True)
os.makedirs(raw_output_dir, exist_ok = True)

print("Generating configs...")
generate_configs(config_files_dir,
                 bod_files_dir,
                 raw_output_dir,
                 num_tests,
                 num_walks_list)
