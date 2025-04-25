#!/bin/bash
# OPES Dipeptide Setup Validator and Generator

# Configuration
GRO_FILE="npt.gro"
PLUMED_FILE="plumed_opes.dat"
OUTPUT_PREFIX="opes"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m' # No Color

# Function to verify atom existence
verify_atom() {
    local atom_num=$1
    local expected_name=$2
    local atom_info=$(grep "^.\{15\}$atom_num " "$GRO_FILE" 2>/dev/null)
    
    if [ -z "$atom_info" ]; then
        echo -e "${RED}ERROR: Atom $atom_num not found in $GRO_FILE${NC}"
        return 1
    fi
    
    local atom_name=$(echo "$atom_info" | awk '{print $2}')
    if [ "$atom_name" != "$expected_name" ]; then
        echo -e "${RED}WARNING: Atom $atom_num is '$atom_name' (expected '$expected_name')${NC}"
    else
        echo -e "${GREEN}Atom $atom_num verified as $expected_name${NC}"
    fi
}

# Main validation
echo -e "\n=== Verifying Dipeptide Structure ==="

# Check critical atoms
declare -A required_atoms=(
    [5]="C"    # ACE C
    [7]="N"    # ALA N
    [9]="CA"   # ALA CA
    [15]="C"   # ALA C
    [17]="N"   # NME N
)

all_atoms_valid=true
for atom_num in "${!required_atoms[@]}"; do
    if ! verify_atom "$atom_num" "${required_atoms[$atom_num]}"; then
        all_atoms_valid=false
    fi
done

# Create PLUMED file if all atoms are valid
if $all_atoms_valid; then
    echo -e "\n=== Creating PLUMED Input File ==="
    cat > "$PLUMED_FILE" << EOF
# OPES input for ACE-ALA-NME dipeptide
# Generated from $GRO_FILE
phi: TORSION ATOMS=5,7,9,15   # ACE-C, ALA-N, ALA-CA, ALA-C
psi: TORSION ATOMS=7,9,15,17  # ALA-N, ALA-CA, ALA-C, NME-N
PRINT STRIDE=100 ARG=phi,psi FILE=COLVAR
opes: OPES_METAD ARG=phi,psi PACE=500 BARRIER=50 SIGMA=0.2,0.2 TEMP=300 FILE=KERNELS
EOF

    echo -e "${GREEN}Successfully created $PLUMED_FILE with:${NC}"
    echo "Phi: Atoms 5,7,9,15"
    echo "Psi: Atoms 7,9,15,17"
    
    # Verify file creation
    if [ -f "$PLUMED_FILE" ]; then
        echo -e "\n=== File Contents ==="
        cat "$PLUMED_FILE"
    else
        echo -e "${RED}ERROR: Failed to create $PLUMED_FILE${NC}"
        exit 1
    fi
else
    echo -e "\n${RED}ERROR: Cannot create PLUMED file due to missing atoms${NC}"
    echo "Please check your $GRO_FILE and adjust atom numbers accordingly"
    exit 1
fi

# Final validation step
echo -e "\n=== Running Final Validation ==="
if command -v plumed &>/dev/null; then
    if plumed info < "$PLUMED_FILE" >/dev/null 2>&1; then
        echo -e "${GREEN}PLUMED input file is valid${NC}"
    else
        echo -e "${RED}ERROR: PLUMED validation failed${NC}"
        exit 1
    fi
else
    echo -e "${RED}WARNING: plumed command not found - final validation skipped${NC}"
fi

echo -e "\n=== Setup Complete ==="
echo "You can now run the simulation with:"
echo "gmx mdrun -s $OUTPUT_PREFIX.tpr -plumed $PLUMED_FILE -deffnm $OUTPUT_PREFIX"
