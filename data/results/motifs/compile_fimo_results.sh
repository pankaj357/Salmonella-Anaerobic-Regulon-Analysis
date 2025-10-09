#!/bin/bash

# compile_fimo_results.sh
# This script compiles significant hits from fimo_results_new/ directory into one file per regulator.

# Define a list of your regulators
regulators=("ArcA" "Fnr" "NarL" "NarP" "NsrR" "DcuR" "TtrR")

# Create a final output directory
mkdir -p compiled_results

# Loop through each regulator
for reg in "${regulators[@]}"; do
  echo "Processing $reg..."
  
  # FIMO output is now expected in fimo_results_new/Regulator_fimo.tsv
  fimo_file="fimo_results/${reg}_fimo.tsv"
  
  # Check if the FIMO output file exists
  if [[ -f "$fimo_file" ]]; then
    # Extract significant hits (q-value < 0.05), add a regulator column, and save
    awk -F'\t' 'BEGIN {OFS="\t"} \
        NR==1 {print "Regulator", $0; next} \
        $7 < 0.05 {print "'"$reg"'", $0}' "$fimo_file" > "compiled_results/${reg}_hits.tsv"
    
    echo "  Saved significant hits to compiled_results/${reg}_hits.tsv"
  else
    echo "  WARNING: File $fimo_file not found. Skipping $reg."
  fi
done

echo "All done! Check the 'compiled_results' directory."
