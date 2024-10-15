#!/bin/bash

# Create the output directory if it doesn't exist
mkdir -p ./html

# Define an array for the samples and their corresponding "npcs" values
samples=("iPSC_parental" "mesoderm" "endoderm" "ectoderm")
npcs_values=(15 15 15 15)

# Loop through each sample and render the .qmd file
for i in "${!samples[@]}"; do
  sample="${samples[i]}"

  # Get the corresponding "npcs" value from the array
  npcs="${npcs_values[i]}"

  # Create a temporary directory for each sample
  temp_dir=$(mktemp -d)

  # Render the file to the temporary directory with the "npcs" argument
  quarto render PreProc.qmd -P sample="$sample" -P npcs="$npcs" --output-dir "$temp_dir" -o "PreProc_${sample}.html"

  # Move the rendered file to the html directory
  mv "$temp_dir/PreProc_${sample}.html" ./html/

  # Remove the temporary directory
  rm -rf "$temp_dir"
done
