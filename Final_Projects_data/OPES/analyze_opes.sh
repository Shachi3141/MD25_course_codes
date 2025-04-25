#!/bin/bash
# OPES/Metadynamics Analysis Script
# Handles: FES, CV evolution, simulation quality checks

# ==========================================
# Configuration
# ==========================================
OUTPUT_PREFIX="opes"
PLUMED_FILE="plumed_opes.dat"
GRO_FILE="${OUTPUT_PREFIX}.gro"
CPT_FILE="${OUTPUT_PREFIX}.cpt"
XTC_FILE="${OUTPUT_PREFIX}.xtc"
EDR_FILE="${OUTPUT_PREFIX}.edr"

# ==========================================
# 1. Free Energy Surface Calculation
# ==========================================
echo "=== Calculating Free Energy Surface ==="
plumed sum_hills --hills HILLS --mintozero --outfile fes.dat --stride 1000

# Generate 2D FES plot
gnuplot << EOF
set terminal pngcairo enhanced font "Arial,12"
set output 'fes.png'
set view map
set xlabel "Phi (rad)"
set ylabel "Psi (rad)"
set cblabel "Free Energy (kJ/mol)"
splot 'fes.dat' u 1:2:3 w pm3d
EOF

# ==========================================
# 2. Collective Variable Analysis
# ==========================================
echo "=== Analyzing Collective Variables ==="
# Plot phi/psi over time
gnuplot << EOF
set terminal pngcairo enhanced font "Arial,12"
set output 'cv_evolution.png'
set xlabel "Time (ps)"
set ylabel "Angle (rad)"
set yrange [-pi:pi]
plot \
'COLVAR' u 1:2 w l lw 2 title "Phi", \
'COLVAR' u 1:3 w l lw 2 title "Psi"
EOF

# ==========================================
# 3. Simulation Quality Checks
# ==========================================
echo "=== Running Quality Checks ==="
# Energy analysis
echo "Energy Terms:" | tee energy_analysis.txt
gmx energy -f $EDR_FILE -o energy.xvg << EOF
Potential
Temperature
Pressure
EOF

# RMSD analysis
gmx rms -s $GRO_FILE -f $XTC_FILE -o rmsd.xvg << EOF
Protein
Protein
EOF

# ==========================================
# 4. HILLS File Analysis
# ==========================================
echo "=== Analyzing Metadynamics HILLS ==="
# Plot Gaussian heights
gnuplot << EOF
set terminal pngcairo enhanced font "Arial,12"
set output 'hills_analysis.png'
set xlabel "Time (ps)"
set ylabel "Gaussian Height (kJ/mol)"
plot 'HILLS' u 1:8 w l lw 2 title "Hill Height"
EOF

# ==========================================
# Final Output
# ==========================================
echo -e "\n=== Analysis Complete ==="
echo "Generated files:"
ls -lh fes.dat fes.png cv_evolution.png energy.xvg rmsd.xvg hills_analysis.png

echo -e "\nTo view the FES:"
echo "gnuplot -e \"set view map; splot 'fes.dat' u 1:2:3 w pm3d; pause -1\""
