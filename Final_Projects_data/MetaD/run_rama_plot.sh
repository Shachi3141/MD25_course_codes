#!/bin/bash

# Check if PLUMED is installed
if ! command -v plumed &> /dev/null; then
    echo "PLUMED is not installed or not in PATH. Please install PLUMED and try again."
    exit 1
fi

# Check if gnuplot is installed
if ! command -v gnuplot &> /dev/null; then
    echo "gnuplot is not installed or not in PATH. Please install gnuplot and try again."
    exit 1
fi

# Disable X11 display to avoid authorization errors
export DISPLAY=

# Step 1: Define PLUMED input file
cat > plumed.dat << EOL
# Define dihedral angles for alanine dipeptide
phi: TORSION ATOMS=1,2,3,4
psi: TORSION ATOMS=2,3,4,5

# Print dihedrals to a file for scatter plot
PRINT ARG=phi,psi FILE=phi_psi_data.xvg STRIDE=10

# Metadynamics to generate HILLS file (uncomment if not already done)
# METAD ARG=phi,psi SIGMA=0.1,0.1 HEIGHT=1.2 PACE=500 FILE=HILLS
EOL

echo "PLUMED input file (plumed.dat) created."

# Step 2: Run PLUMED driver to process trajectory
echo "Processing trajectory with PLUMED..."
plumed driver --plumed plumed.dat --igro step5_2.gro --ixtc step5_2.xtc

if [ $? -ne 0 ]; then
    echo "Error running plumed driver. Check plumed.dat, step5_2.gro, and step5_2.xtc."
    exit 1
fi

echo "Trajectory processed. Generated phi_psi_data.xvg and HILLS file."

# Step 3: Reconstruct free energy surface
echo "Reconstructing free energy surface..."
plumed sum_hills --hills HILLS --outfile energy_matrix.txt --bin 36,36 --min -180,-180 --max 180,180

if [ $? -ne 0 ]; then
    echo "Error running sum_hills. Check HILLS file existence. If missing, enable METAD in plumed.dat and rerun the simulation."
    exit 1
fi

echo "Free energy surface saved to energy_matrix.txt."

# Step 4: Generate Ramachandran plot with gnuplot
echo "Generating Ramachandran plot..."
cat > ramachandran.plt << EOL
set terminal pngcairo size 800,600
set output "ramachandran_plot.png"
set title "Ramachandran Plot: ϕ vs ψ"
set xlabel "ϕ (degrees)"
set ylabel "ψ (degrees)"
set xrange [-180:180]
set yrange [-180:180]
set grid
set palette rgbformulae 33,13,10
set colorbox vertical
set cbrange [-20:0]
plot "energy_matrix.txt" matrix with image title "", \
     "phi_psi_data.xvg" using 1:2 with points pointtype 7 pointsize 1 linecolor rgb "red" title "Sampled Conformations"
set output
EOL

gnuplot ramachandran.plt

if [ $? -ne 0 ]; then
    echo "Error generating plot. Check gnuplot installation and files."
    exit 1
fi

echo "Ramachandran plot saved as ramachandran_plot.png."
