#!/usr/bin/python3

import os

from validation import run_tests

zeno_dir = '../../..'

config_files_dir = os.path.join(zeno_dir,
                                'cpp/tests/validation/'
                                'surrogate_ground_truth/config_files')

zeno_exec = os.path.join(zeno_dir,
                         'cpp/build/zeno')

for test_dir in sorted(os.listdir(config_files_dir)):
    test_dir_path = os.path.join(config_files_dir,
                                 test_dir)
    
    if os.path.isdir(test_dir_path):
        run_tests(test_dir_path,
                  zeno_exec)
