#!/bin/bash

# Create the output directory if it doesn't exist
mkdir -p ./html

# Loop through each sample and render the .qmd file
for sample in iPSC_parental mesoderm endoderm ectoderm; do
  # Create a temporary directory for each sample
  temp_dir=$(mktemp -d)

  # Render the file to the temporary directory
  quarto render Clustering.qmd -P sample="$sample" --output-dir "$temp_dir" -o "Clustering_${sample}.html"

  # Move the rendered file to the html directory
  mv "$temp_dir/Clustering_${sample}.html" ./html/

  # Remove the temporary directory
  rm -rf "$temp_dir"
done
