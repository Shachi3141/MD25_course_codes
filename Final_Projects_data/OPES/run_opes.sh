#!/bin/bash
# OPES Simulation Setup and Execution Script

# ==========================================
# CONFIGURATION
# ==========================================
TPR_FILE="opes.tpr"
PLUMED_FILE="plumed_opes.dat"
OUTPUT_PREFIX="opes"
N_CORES=30  # Adjust based on your system

# ==========================================
# FUNCTIONS
# ==========================================
validate_plumed() {
    echo "=== Validating PLUMED input ==="
    # Check for hidden characters
    echo -n "Checking for hidden characters... "
    if grep -q -P '\t' "$PLUMED_FILE"; then
        echo "ERROR: Tabs detected. Replace with spaces."
        exit 1
    fi
    echo "OK"
    
    # Check line endings
    echo -n "Checking line endings... "
    if file "$PLUMED_FILE" | grep -q CRLF; then
        echo "Converting Windows line endings to Unix..."
        dos2unix "$PLUMED_FILE"
    fi
    echo "OK"
    
    # Minimal validation
    echo -n "Basic PLUMED validation... "
    plumed info < "$PLUMED_FILE" > /dev/null 2>&1
    if [ $? -ne 0 ]; then
        echo "ERROR: Invalid PLUMED input"
        exit 1
    fi
    echo "OK"
}

check_atom_numbers() {
    echo "=== Verifying atom numbers ==="
    # Create index file if it doesn't exist
    if [ ! -f index.ndx ]; then
        echo "Creating index file..."
        gmx make_ndx -f npt.gro -o index.ndx <<< "q"
    fi
    
    # Check selected atoms
    for atom in 5 7 9 15 17; do
        if ! gmx select -f npt.gro -select "atomnr $atom" -on /dev/null 2>/dev/null; then
            echo "ERROR: Atom number $atom not found in structure"
            exit 1
        fi
    done
    echo "All specified atoms exist"
}

clean_environment() {
    echo "=== Cleaning environment ==="
    rm -f \#* 
    echo "Cleanup complete"
}

generate_tpr() {
    echo "=== Generating TPR file ==="
    gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o "$TPR_FILE"
    if [ $? -ne 0 ]; then
        echo "ERROR: grompp failed"
        exit 1
    fi
    echo "Successfully generated $TPR_FILE"
}

run_simulation() {
    echo "=== Starting OPES Simulation ==="
    echo "Using $N_CORES CPU cores"
    
    PLUMED_LOAD_DEBUG=on PLUMED_DEBUG=on gmx mdrun \
        -s "$TPR_FILE" \
        -plumed "$PLUMED_FILE" \
        -deffnm "$OUTPUT_PREFIX" \
        -ntmpi 1 -ntomp "$N_CORES" \
        -nb cpu -pme cpu \
        -bonded cpu -update cpu \
        -v 2>&1 | tee simulation.log
    
    if [ $? -ne 0 ]; then
        echo "ERROR: Simulation failed"
        grep -i error simulation.log | tail -n 5
        exit 1
    fi
    
    echo "=== Simulation Completed Successfully ==="
}

# ==========================================
# MAIN EXECUTION
# ==========================================
validate_plumed
check_atom_numbers
clean_environment
generate_tpr
run_simulation

echo "=== Final Output Files ==="
ls -lh "${OUTPUT_PREFIX}".*
