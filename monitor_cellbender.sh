#!/bin/bash

# Watch CellBender progress

echo "Monitoring CellBender jobs..."
echo "Press Ctrl+C to exit"
echo ""

while true; do
  clear
  echo "=== System Resources ==="
  echo "RAM: $(free -h | grep Mem | awk '{print $3 " / " $2 " (" $7 " available)"}')"
  echo "Swap: $(free -h | grep Swap | awk '{print $3 " / " $2}')"
  echo ""
  
  echo "=== Running CellBender Processes ==="
  ps aux | grep -E "cellbender|PID" | grep -v grep | head -15
  echo ""
  
  echo "=== Job Count ==="
  RUNNING=$(ps aux | grep cellbender | grep -v grep | wc -l)
  echo "Active CellBender jobs: $RUNNING"
  echo ""
  
  echo "=== Completed Samples ==="
  find /projects/opioid-per/rawdata/*/cellbender_output_filtered.h5 2>/dev/null | \
    wc -l | awk '{print $1 " / 11 samples complete"}'
  echo ""
  
  echo "Updated: $(date)"
  sleep 10
done
