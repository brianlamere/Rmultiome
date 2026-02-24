#!/bin/bash

# CellBender for samples with VARIABLE quality (2.5% - 8.5% MT)
# Normalizes ambient RNA profiles before integration

RAWDIR="/projects/opioid-per/rawdata"

# Sample-specific parameters based on web_summary cell counts
declare -A SAMPLE_CELLS
SAMPLE_CELLS["LG05"]=15000
SAMPLE_CELLS["LG08"]=12000   # Adjust based on your web summaries
SAMPLE_CELLS["LG22"]=18000
SAMPLE_CELLS["LG23"]=16000
SAMPLE_CELLS["LG25"]=14000
SAMPLE_CELLS["LG26"]=13000
SAMPLE_CELLS["LG31"]=15000
SAMPLE_CELLS["LG33"]=19000
SAMPLE_CELLS["LG38"]=17000   # Pristine quality (2.5% MT)
SAMPLE_CELLS["LG300"]=18000
SAMPLE_CELLS["LG301"]=20000

for SAMPLE in "${!SAMPLE_CELLS[@]}"; do
  EXPECTED=${SAMPLE_CELLS[$SAMPLE]}
  TOTAL=$((EXPECTED * 3))
  
  echo "================================================"
  echo "Sample: ${SAMPLE}"
  echo "Expected cells: ${EXPECTED}"
  echo "Total droplets: ${TOTAL}"
  echo "================================================"
  
  cellbender remove-background \
    --input "${RAWDIR}/${SAMPLE}/filtered_feature_bc_matrix.h5" \
    --output "${RAWDIR}/${SAMPLE}/cellbender_output.h5" \
    --expected-cells $EXPECTED \
    --total-droplets-included $TOTAL \
    --fpr 0.01 \
    --epochs 150 \
    --cuda
    
  echo ""
done
