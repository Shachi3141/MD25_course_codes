#!/bin/bash

# === Configuration ===
tpr_file="topol.tpr"
plumed_input="plumed.dat"
deffnm="meta"
ntmpi=1
ntomp=24

# === Step 1: Check Required Files ===
if [[ ! -f "$tpr_file" ]]; then
    echo "ERROR: $tpr_file not found. Make sure your TPR is ready."
    exit 1
fi
if [[ ! -f "$plumed_input" ]]; then
    echo "ERROR: $plumed_input not found. Create it before running."
    exit 1
fi

# === Step 2: Run Metadynamics ===
echo "Running Metadynamics..."
gmx mdrun -s $tpr_file -deffnm $deffnm -plumed $plumed_input -ntmpi $ntmpi -ntomp $ntomp

# === Step 3: Build Free Energy Surface ===
echo "Generating Free Energy Surface from HILLS..."
plumed sum_hills --hills HILLS --outfile fes.dat --stride 10 --mintozero

# === Step 4: Organize Outputs for OPES Comparison ===
echo "Organizing output for comparison..."
mkdir -p metadynamics_output
cp COLVAR HILLS fes.dat metadynamics_output/

echo "Metadynamics simulation completed successfully!"

