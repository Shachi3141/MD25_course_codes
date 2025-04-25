#!/bin/bash
# OPES Dipeptide Setup Script - Adaptive Version

# Configuration
GRO_FILE="npt.gro"
PLUMED_FILE="plumed_opes.dat"
OUTPUT_PREFIX="opes"

# Get all protein atoms
PROTEIN_ATOMS=$(grep -E 'ACE|ALA|NME' "$GRO_FILE")

# Find specific atoms
ACE_C=$(echo "$PROTEIN_ATOMS" | awk '/1ACE.* C / {print $3; exit}')
ALA_N=$(echo "$PROTEIN_ATOMS" | awk '/2ALA.* N / {print $3; exit}')
ALA_CA=$(echo "$PROTEIN_ATOMS" | awk '/2ALA.* CA / {print $3; exit}')
ALA_C=$(echo "$PROTEIN_ATOMS" | awk '/2ALA.* C / {print $3; exit}')
NME_N=$(echo "$PROTEIN_ATOMS" | awk '/3NME.* N / {print $3; exit}')

# Verify we found all required atoms
REQUIRED_ATOMS=("$ACE_C" "$ALA_N" "$ALA_CA" "$ALA_C" "$NME_N")
ATOM_NAMES=("ACE C" "ALA N" "ALA CA" "ALA C" "NME N")

for i in "${!REQUIRED_ATOMS[@]}"; do
  if [ -z "${REQUIRED_ATOMS[$i]}" ]; then
    echo "ERROR: Could not find ${ATOM_NAMES[$i]} in $GRO_FILE"
    echo "Available protein atoms:"
    echo "$PROTEIN_ATOMS" | awk '{printf "Atom %3s: %-6s %-4s\n", $3, $1, $2}'
    exit 1
  fi
done

# Create PLUMED file
cat > "$PLUMED_FILE" << EOF
# OPES input for ACE-ALA-NME dipeptide
# Generated from $GRO_FILE
phi: TORSION ATOMS=$ACE_C,$ALA_N,$ALA_CA,$ALA_C   # ACE-C, ALA-N, ALA-CA, ALA-C
psi: TORSION ATOMS=$ALA_N,$ALA_CA,$ALA_C,$NME_N  # ALA-N, ALA-CA, ALA-C, NME-N
PRINT STRIDE=100 ARG=phi,psi FILE=COLVAR
opes: OPES_METAD ARG=phi,psi PACE=500 BARRIER=50 SIGMA=0.2,0.2 TEMP=300 FILE=KERNELS
EOF

echo "Successfully created $PLUMED_FILE with these atom numbers:"
echo "ACE C: $ACE_C"
echo "ALA N: $ALA_N"
echo "ALA CA: $ALA_CA"
echo "ALA C: $ALA_C"
echo "NME N: $NME_N"

echo -e "\nTo verify the atom selection, run:"
echo "gmx select -f $GRO_FILE -select 'atomnr $ACE_C or atomnr $ALA_N or atomnr $ALA_CA or atomnr $ALA_C or atomnr $NME_N' -on opes_atoms.ndx"
