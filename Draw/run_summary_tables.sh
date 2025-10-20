#!/bin/bash

# SMOBench Summary Table Generator
# Adapted from scib RNA Single-Task Summary Tables for vertical integration

echo "=== SMOBench Summary Table Generator ==="
echo "Starting summary table generation..."

# Check if Results/summary_table directory exists, create if not
if [ ! -d "Results/summary_table" ]; then
    echo "Creating Results/summary_table directory..."
    mkdir -p Results/summary_table
fi

echo "Executing R script for summary table generation..."

# Run the R script  
Rscript Draw/plotSMOBenchSummaryTable.R

if [ $? -eq 0 ]; then
    echo "Success: Summary tables generated!"
    echo "Output directory: Results/summary_table/"
    echo ""
    echo "Generated files:"
    ls -la Results/summary_table/SMOBench_summary_table_*.png 2>/dev/null || echo "No summary table files found"
else
    echo "Error: R script execution failed!"
    echo "Please check the error messages above"
fi

echo "=== Summary Table Generation Complete ==="