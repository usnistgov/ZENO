#!/usr/bin/python3

import os
import subprocess
import argparse

arg_parser = argparse.ArgumentParser()

arg_parser.add_argument("-s", "--script-path", required=True,
                        help="Path to convert_csv_to_rst script")

arg_parser.add_argument("-i", "--input-csv-dir", required=True,
                        help="Directory of CSV files containing the table data")

arg_parser.add_argument("-o", "--output-rst-dir", required=True,
                        help="Directory into which to write the "
                        "table data RST files")

args = arg_parser.parse_args()

convert_csv_to_rst_script = args.script_path
csv_dir_path = args.input_csv_dir
rst_dir_path = args.output_rst_dir

for csv_file_name in sorted(os.listdir(csv_dir_path)):
    rst_file_name = os.path.splitext(csv_file_name)[0] + ".rst"
    
    csv_file_path = os.path.join(csv_dir_path,
                                 csv_file_name)

    rst_file_path = os.path.join(rst_dir_path,
                                 rst_file_name)

    command = [
        convert_csv_to_rst_script,
        "-i", csv_file_path,
        "-o", rst_file_path
    ]

    print(command)

    subprocess.call(command)
