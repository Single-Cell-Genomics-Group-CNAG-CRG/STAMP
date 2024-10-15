#!/bin/bash

# Create the output directory if it doesn't exist
mkdir -p ./html

# Loop through each sample and render the .qmd file
for sample in iPSC_parental ectoderm mesoderm endoderm; do
  # Create a temporary directory for each sample
  temp_dir=$(mktemp -d)

  # Render the file to the temporary directory
  quarto render QC.qmd -P sample="$sample" --output-dir "$temp_dir" -o "QC_${sample}.html"

  # Move the rendered file to the html directory
  mv "$temp_dir/QC_${sample}.html" ./html/

  # Remove the temporary directory
  rm -rf "$temp_dir"
done
