#!/usr/bin/python3

import math
import csv
import argparse

def get_mean(vals):
    N = len(vals)
    
    return sum(vals) / N

def get_diff(ground_truth, mean):
    return abs(ground_truth - mean)

def get_relative_diff(ground_truth, diff):
    return (diff / ground_truth) * 100

def get_std_dev(mean, vals):
    N = len(vals)
    
    return math.sqrt(sum((x - mean)**2 for x in vals) / (N - 1))

def get_std_uncert(std_dev, N):
    return std_dev / math.sqrt(N)

def get_expanded_uncert(std_uncert, t):
    return std_uncert * t

def get_relative_expanded_uncert(ground_truth, expanded_uncert):
    return (expanded_uncert / ground_truth) * 100

arg_parser = argparse.ArgumentParser()

arg_parser.add_argument("-t", "--test-results", required=True,
                        help="CSV file containing the combined results from "
                        "all test runs")
                        
arg_parser.add_argument("-g", "--ground-truth", required=True,
                        help="CSV file containing the formatted ground-truth "
                        "values")

arg_parser.add_argument("-v", "--validation-data", required=True,
                        help="Name of a file into which to write the "
                        "validation table data in CSV format")

args = arg_parser.parse_args()

test_results_filename = args.test_results

ground_truth_filename = args.ground_truth

validation_data_filename = args.validation_data

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

ground_truth_means = {}

ground_truth_units = {}

with open(ground_truth_filename, newline='') as ground_truth_file:
    ground_truth_reader = csv.reader(ground_truth_file)

    for row in ground_truth_reader:
        property_name = row[0]
        property_type = row[1]
        property_value = row[2]
        
        if property_name in property_names:
            if property_type == "value":
                ground_truth_means[property_name] = property_value

            elif property_type == "units":
                ground_truth_units[property_name] = property_value

with open(validation_data_filename, 'w', newline='') as validation_data_file:
    validation_data_writer = csv.writer(validation_data_file)

    validation_data_header = ["Property",
                              "Ground Truth",
                              "Units",
                              "Mean",
                              "Difference",
                              "Rel. Diff.",
                              "Std. Dev.",
                              "Std. Uncert.",
                              "Expand Uncert.",
                              "Rel. Exp. Uncert."]

    validation_data_writer.writerow(validation_data_header)

    with open(test_results_filename, newline='') as test_results_file:
        test_results_reader = csv.reader(test_results_file)

        for test_results_row in test_results_reader:
            test_property_name = test_results_row[0]
            test_property_type = test_results_row[1]
        
            if test_property_name in property_names:
                if test_property_type == "value":
                    if test_property_name not in ground_truth_means:
                        print("Warning: " + test_property_name +
                              " not found in " + ground_truth_filename)

                        continue
                    
                    test_property_vals = \
                    [float(x) for x in test_results_row[2:]]
            
                    t = 2.009575 # 95% confidence interval

                    ground_truth = ground_truth_means[test_property_name]

                    unit = ground_truth_units[test_property_name]
                    
                    N = len(test_property_vals)
                    
                    mean = get_mean(test_property_vals)

                    diff = get_diff(float(ground_truth), mean)

                    relative_diff = get_relative_diff(float(ground_truth), diff)

                    std_dev = get_std_dev(mean, test_property_vals)

                    std_uncert = get_std_uncert(std_dev, N)

                    expanded_uncert = \
                    get_expanded_uncert(std_uncert, t)

                    relative_expanded_uncert = \
                    get_relative_expanded_uncert(float(ground_truth),
                                                 expanded_uncert)

                    mean_decimal_places = \
                    0 if "." not in ground_truth else len(ground_truth.split(".")[1])

                    def compute_decimal_places(mean_decimal_places, value):
                        decimal_places = \
                        mean_decimal_places - math.floor(math.log10(diff)) - 1
                        
                        decimal_places = max(0, min(2, decimal_places))

                        return decimal_places

                    diff_decimal_places = compute_decimal_places(mean_decimal_places,
                                                                 diff)

                    std_dev_decimal_places = compute_decimal_places(mean_decimal_places,
                                                                    std_dev)

                    std_uncert_decimal_places = compute_decimal_places(mean_decimal_places,
                                                                       std_uncert)

                    expanded_uncert_decimal_places = compute_decimal_places(mean_decimal_places,
                                                                            expanded_uncert)
                    
                    validation_data_row = [test_property_name,
                                           ground_truth,
                                           unit,
                                           "{0:.{1}f}".format(mean, mean_decimal_places),
                                           "{0:.{1}e}".format(diff, diff_decimal_places),
                                           "{0:.3f}%".format(relative_diff),
                                           "{0:.{1}e}".format(std_dev, std_dev_decimal_places),
                                           "{0:.{1}e}".format(std_uncert, std_uncert_decimal_places),
                                           "{0:.{1}e}".format(expanded_uncert, expanded_uncert_decimal_places),
                                           "{0:.3f}%".format(relative_expanded_uncert)]

                    validation_data_writer.writerow(validation_data_row)
