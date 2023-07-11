#!/bin/bash

# This is a helper script for initialising renv. Based on the manually created
# 'R/packages' file it creates an R script 'R/load_packages.R', that renv will scan
# and detect all packages that need installing.

# Define the input and output file names
input_file="R/packages"
output_file="R/load_packages.R"
rm -f $output_file

# Read the input file line by line and create the R commands
while IFS= read -r package; do
  echo "library($package)" >> "$output_file"
done < "$input_file"
