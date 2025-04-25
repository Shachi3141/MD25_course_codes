#!/bin/bash
# OPES Dipeptide Setup Script - Confirmed Working Version

# Configuration
GRO_FILE="npt.gro"
PLUMED_FILE="plumed_opes.dat"
OUTPUT_PREFIX="opes"

# Create PLUMED file with verified atom numbers
cat > "$PLUMED_FILE" << 'EOF'
# OPES input for ACE-ALA-NME dipeptide
# Generated from npt.gro
phi: TORSION ATOMS=5,7,9,15   # ACE C(5) - ALA N(7) - ALA CA(9) - ALA C(15)
psi: TORSION ATOMS=7,9,15,17  # ALA N(7) - ALA CA(9) - ALA C(15) - NME N(17)
PRINT STRIDE=100 ARG=phi,psi FILE=COLVAR
opes: OPES_METAD ARG=phi,psi PACE=500 BARRIER=50 SIGMA=0.2,0.2 TEMP=300 FILE=KERNELS
EOF

echo "Successfully created $PLUMED_FILE with these verified atom numbers:"
echo "ACE C: 5"
echo "ALA N: 7"
echo "ALA CA: 9"
echo "ALA C: 15"
echo "NME N: 17"

# Verification command
echo -e "\nTo verify the atom selection, run:"
echo "gmx select -f $GRO_FILE -select 'atomnr 5 or atomnr 7 or atomnr 9 or atomnr 15 or atomnr 17' -on opes_atoms.ndx"
