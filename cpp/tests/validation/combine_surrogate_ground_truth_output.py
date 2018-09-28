#!/usr/bin/python3

import os

from validation import combine_output

zeno_dir = '../../..'

raw_output_dir = os.path.join(zeno_dir,
                                'cpp/tests/validation/'
                                'surrogate_ground_truth/raw_output')

parsed_output_dir = os.path.join(zeno_dir,
                                 'cpp/tests/validation/'
                                 'surrogate_ground_truth/combined_output')

for output_csv_files_dir in sorted(os.listdir(raw_output_dir)):
    
        output_csv_files_dir_path = os.path.join(raw_output_dir,
                                                 output_csv_files_dir)

        combined_csv_file_name = output_csv_files_dir + ".csv"

        combined_csv_file_path = os.path.join(parsed_output_dir,
                                              combined_csv_file_name)

        if os.path.isdir(output_csv_files_dir_path):
            combine_output(output_csv_files_dir_path,
                           combined_csv_file_path)
