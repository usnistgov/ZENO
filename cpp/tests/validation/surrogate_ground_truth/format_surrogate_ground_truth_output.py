#!/usr/bin/python3

import math
import csv
import argparse
import statistics
import numpy as np

arg_parser = argparse.ArgumentParser()

arg_parser.add_argument("-i", "--input-combined-file", required=True,
                        help="Path to file containing the combined output "
                        "from the surrogate ground truth computation")

arg_parser.add_argument("-o", "--output-formatted-file", required=True,
                        help="Path to file into which to write the formatted "
                        "surrogate ground truth data")

args = arg_parser.parse_args()

input_filename = args.input_combined_file

output_filename = args.output_formatted_file

property_names = ["capacitance",
                  "electric_polarizability_eigenvalues[0]",
                  "electric_polarizability_eigenvalues[1]",
                  "electric_polarizability_eigenvalues[2]",
                  "mean_electric_polarizability",
                  "hydrodynamic_radius",
                  "volume",
                  "capacitance_sphere_same_volume",
                  "gyration_eigenvalues[0]",
                  "gyration_eigenvalues[1]",
                  "gyration_eigenvalues[2]",
                  "intrinsic_conductivity",
                  "intrinsic_viscosity"]

ground_truth_values = {}

ground_truth_reported_std_devs = {}

# The std devs reported by zeno include sources of uncertainty other than
# just the Monte Carlo algorithm, so they are not what we want to use to
# determine the number of decimal digits for the ground truth.  For that,
# we compute the actual std dev and std uncert of the test results.
ground_truth_actual_std_uncerts = {}

ground_truth_units = {}

with open(input_filename, newline='') as input_file:
    input_reader = csv.reader(input_file)

    for row in input_reader:
        property_name = row[0]
        property_type = row[1]

        if property_name in property_names:
            if property_type == "value":
                property_values = [float(x) for x in row[2:]]

                property_values_mean = statistics.mean(property_values)

                ground_truth_values[property_name] = property_values_mean

                property_values_std_dev = statistics.stdev(property_values)

                property_values_std_uncert = \
                property_values_std_dev / math.sqrt(len(property_values))

                ground_truth_actual_std_uncerts[property_name] = \
                property_values_std_uncert

            elif property_type == "std_dev":
                property_std_devs = [float(x) for x in row[2:]]

                property_std_devs_mean = statistics.mean(property_std_devs)

                ground_truth_reported_std_devs[property_name] = \
                property_std_devs_mean

            elif property_type == "units":
                property_vals = row[2:]

                ground_truth_units[property_name] = property_vals[0]

with open(output_filename, 'w', newline='') as output_file:
    output_writer = csv.writer(output_file)

    for property_name in property_names:
        std_uncert = ground_truth_actual_std_uncerts[property_name]

        num_decimal_places = -1*math.floor(math.log10(std_uncert)) - 1

        rounded_value = np.around(ground_truth_values[property_name],
                                  decimals=num_decimal_places)

        formatted_value = "{0:.{1}f}".format(rounded_value,
                                             num_decimal_places)

        value_row = [property_name,
                     "value",
                     formatted_value]

        units_row = [property_name,
                     "units",
                     ground_truth_units[property_name]]

        std_dev_row = [property_name,
                       "std_dev",
                       ground_truth_reported_std_devs[property_name]]

        output_writer.writerow(units_row)
        output_writer.writerow(value_row)
        output_writer.writerow(std_dev_row)
