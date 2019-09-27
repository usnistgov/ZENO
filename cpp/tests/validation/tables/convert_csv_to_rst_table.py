#!/usr/bin/python3

import os
import csv
import argparse
import pytablewriter as ptw

arg_parser = argparse.ArgumentParser()

arg_parser.add_argument("-i", "--input-csv-file", required=True,
                        help="CSV file containing the table data")

arg_parser.add_argument("-o", "--output-rst-file", required=True,
                        help="Name of a file into which to write the "
                        "table data in RST format")

args = arg_parser.parse_args()

csv_filename = args.input_csv_file

rst_filename = args.output_rst_file

property_name_csv_to_rst = {
    "Property": "Property",
    "capacitance": ":math:`C`",
    "electric_polarizability_eigenvalues[0]": "Eigenvalue of :math:`\\mathbf{\\alpha}`",
    "electric_polarizability_eigenvalues[1]": "Eigenvalue of :math:`\\mathbf{\\alpha}`",
    "electric_polarizability_eigenvalues[2]": "Eigenvalue of :math:`\\mathbf{\\alpha}`",
    "mean_electric_polarizability": ":math:`\\langle\\mathbf{\\alpha}\\rangle`",
    "hydrodynamic_radius": ":math:`R_{h}`",
    "volume": ":math:`V`",
    "capacitance_sphere_same_volume": ":math:`C_{0}`",
    "gyration_eigenvalues[0]": "Eigenvalue of :math:`\\mathbf{S}`",
    "gyration_eigenvalues[1]": "Eigenvalue of :math:`\\mathbf{S}`",
    "gyration_eigenvalues[2]": "Eigenvalue of :math:`\\mathbf{S}`",
    "intrinsic_conductivity": "[:math:`\\sigma`]",
    "intrinsic_viscosity": "[:math:`\\eta`]"
}

units_csv_to_rst = {
    "Units": "Units",
    "1": "1",
    "L": "L",
    "L^2": "L :math:`^{2}`",
    "L^3": "L :math:`^{3}`",
    "A": "A",
    "A^2": "A :math:`^{2}`",
    "A^3": "A :math:`^{3}`"
}

table_rows = []

with open(csv_filename, newline='') as csv_file:
    csv_reader = csv.reader(csv_file)

    for row in csv_reader:
        property_name = row[0]
        
        if property_name not in property_name_csv_to_rst:
            raise Exception("Unknown property " + property_name + " in CSV file")
        
        row[0] = property_name_csv_to_rst[property_name]

        units = row[2]

        if units not in units_csv_to_rst:
            raise Exception("Unknown units " + units + " in CSV file")

        row[2] = units_csv_to_rst[units]
        
        table_rows += [row]

with open(rst_filename, 'w') as rst_file:
    rst_writer = ptw.RstGridTableWriter()

    #rst_writer.table_name = os.path.basename(rst_filename)
    
    rst_writer.headers = table_rows[0]

    rst_writer.value_matrix = table_rows[1:]

    # Make sure precision of original numbers is not changed
    rst_writer.type_hints = [ptw.String] * len(rst_writer.headers)

    rst_writer.stream = rst_file

    rst_writer.write_table()
