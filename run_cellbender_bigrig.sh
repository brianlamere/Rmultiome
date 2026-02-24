#!/bin/bash

# CellBender parallelized for high-RAM system
# Your specs: 32 cores, 192GB RAM - can handle 8-10 jobs easily

RAWDIR="/projects1/opioid/rawdata"
MAX_PARALLEL=11  # Conservative; can increase to 10 or even 11

# Exact cell counts from web_summary
declare -A SAMPLE_CELLS
SAMPLE_CELLS["LG05"]=15327
SAMPLE_CELLS["LG08"]=20000
SAMPLE_CELLS["LG22"]=10241
SAMPLE_CELLS["LG23"]=14447
SAMPLE_CELLS["LG25"]=20000
SAMPLE_CELLS["LG26"]=18889
SAMPLE_CELLS["LG31"]=20000
SAMPLE_CELLS["LG33"]=11490
SAMPLE_CELLS["LG38"]=12311
SAMPLE_CELLS["LG300"]=12107
SAMPLE_CELLS["LG301"]=15027

# Function to run CellBender on one sample
run_cellbender_sample() {
  local SAMPLE=$1
  local EXPECTED=$2
  local TOTAL=$((EXPECTED * 3))
  
  local SAMPLE_DIR="${RAWDIR}/${SAMPLE}"
  local INPUT_H5="${RAWDIR}/${SAMPLE}/raw_feature_bc_matrix.h5"
  local OUTPUT_H5="${RAWDIR}/${SAMPLE}/cellbender_output.h5"
  local LOG="${SAMPLE_DIR}/cellbender_output_$(date +%Y%m%d_%H%M%S).log"
  
  cd "$SAMPLE_DIR" || return 1

  echo "[$(date)] Starting ${SAMPLE}: ${EXPECTED} cells" | tee -a "$LOG"
  
  if [ ! -f "$INPUT_H5" ]; then
    echo "ERROR: Input not found: $INPUT_H5" | tee -a "$LOG"
    return 1
  fi
  
  if [ -f "${OUTPUT_H5%.*}_filtered.h5" ]; then
    echo "WARNING: Output exists, skipping ${SAMPLE}" | tee -a "$LOG"
    return 0
  fi
  
  # CPU version, no --cuda
  cellbender remove-background \
    --input "$INPUT_H5" \
    --output "$OUTPUT_H5" \
    --expected-cells $EXPECTED \
    --total-droplets-included $TOTAL \
    --fpr 0.01 \
    --epochs 150 \
    --learning-rate 0.0001 \
    --checkpoint-mins 60 \
    >> "$LOG" 2>&1
  
  local STATUS=$?
  if [ $STATUS -eq 0 ]; then
    echo "[$(date)] ✓ Completed ${SAMPLE}" | tee -a "$LOG"
  else
    echo "[$(date)] ✗ FAILED ${SAMPLE}" | tee -a "$LOG"
  fi
  
  return $STATUS
}

export -f run_cellbender_sample
export RAWDIR

echo "================================================"
echo "System Info:"
echo "  Cores: $(nproc)"
echo "  RAM: $(free -h | grep Mem | awk '{print $2}')"
echo "  Available: $(free -h | grep Mem | awk '{print $7}')"
echo "  Samples: ${#SAMPLE_CELLS[@]}"
echo "  Max parallel: ${MAX_PARALLEL}"
echo "  Started: $(date)"
echo "================================================"
echo ""

# Option 1: Using GNU parallel (if installed)
if command -v parallel &> /dev/null; then
  echo "Using GNU parallel for job control"
  
  # Build job list
  JOBS=()
  for SAMPLE in "${!SAMPLE_CELLS[@]}"; do
    JOBS+=("$SAMPLE ${SAMPLE_CELLS[$SAMPLE]}")
  done
  
  printf '%s\n' "${JOBS[@]}" | \
    parallel -j $MAX_PARALLEL --colsep ' ' \
    run_cellbender_sample {1} {2}

# Option 2: Manual backgrounding with job control
else
  echo "Using bash job control (install GNU parallel for better monitoring)"
  
  RUNNING=0
  for SAMPLE in "${!SAMPLE_CELLS[@]}"; do
    CELLS=${SAMPLE_CELLS[$SAMPLE]}
    
    # Wait if at max parallel jobs
    while [ $(jobs -r | wc -l) -ge $MAX_PARALLEL ]; do
      sleep 10
    done
    
    # Start job in background
    run_cellbender_sample "$SAMPLE" "$CELLS" &
    
    RUNNING=$(jobs -r | wc -l)
    echo "Launched ${SAMPLE} (${RUNNING}/${MAX_PARALLEL} running)"
    
    # Brief pause to prevent stampede
    sleep 2
  done
  
  echo ""
  echo "All jobs launched. Waiting for completion..."
  wait
fi

echo ""
echo "================================================"
echo "All samples complete: $(date)"
echo "================================================"

# Summary report
echo ""
echo "Results:"
for SAMPLE in "${!SAMPLE_CELLS[@]}"; do
  OUTPUT="${RAWDIR}/${SAMPLE}/cellbender_output_filtered.h5"
  if [ -f "$OUTPUT" ]; then
    SIZE=$(du -h "$OUTPUT" | cut -f1)
    echo "  ✓ ${SAMPLE}: ${OUTPUT} (${SIZE})"
  else
    echo "  ✗ ${SAMPLE}: FAILED or INCOMPLETE"
  fi
done
