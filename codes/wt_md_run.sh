#!/bin/bash

# Load GROMACS
source /usr/local/gromacs/bin/GMXRC

# Define directories and files
WT_DIR="ARHGAP5_WT_MD"
#MUT_DIR="ARHGAP5_MUT_MD"
WT_PDB="ARHGAP5_truncated_wt.pdb"
#MUT_PDB="ARHGAP5_truncated_mut.pdb"

# Function to create necessary MDP files
create_mdp_files() {
    local OUTPUT_DIR=$1
    cd $OUTPUT_DIR || { echo "Error: Cannot change to directory $OUTPUT_DIR"; exit 1; }

    # Energy minimization mdp
    cat << EOF > em.mdp
integrator              = steep
nsteps                  = 50000
emtol                   = 1000.0
emstep                  = 0.01
nstxout                 = 100
nstlog                  = 100
cutoff-scheme          = Verlet
nstlist                = 10
rlist                  = 1.0
coulombtype            = PME
rcoulomb               = 1.0
rvdw                   = 1.0
pbc                    = xyz
EOF

    # NVT equilibration mdp
    cat << EOF > nvt.mdp
title                   = NVT equilibration
define                  = -DPOSRES
integrator              = md
nsteps                  = 50000
dt                      = 0.002
nstxout                 = 500
nstvout                 = 500
nstfout                 = 500
nstlog                  = 500
nstenergy               = 500
nstxout-compressed      = 500
continuation            = no
constraint_algorithm    = lincs
constraints             = h-bonds
lincs_iter             = 1
lincs_order            = 4
cutoff-scheme          = Verlet
ns_type                 = grid
nstlist                 = 10
rlist                   = 1.0
coulombtype             = PME
rcoulomb                = 1.0
rvdw                    = 1.0
pbc                     = xyz
tcoupl                  = V-rescale
tc-grps                 = Protein Non-Protein
tau_t                   = 0.1 0.1
ref_t                   = 300 300
pcoupl                  = no
EOF

    # NPT equilibration mdp
    cat << EOF > npt.mdp
title                   = NPT equilibration
define                  = -DPOSRES
integrator              = md
nsteps                  = 50000
dt                      = 0.002
nstxout                 = 500
nstvout                 = 500
nstfout                 = 500
nstlog                  = 500
nstenergy               = 500
nstxout-compressed      = 500
continuation            = yes
constraint_algorithm    = lincs
constraints             = h-bonds
lincs_iter             = 1
lincs_order            = 4
cutoff-scheme          = Verlet
ns_type                 = grid
nstlist                 = 10
rlist                   = 1.0
coulombtype             = PME
rcoulomb                = 1.0
rvdw                    = 1.0
pbc                     = xyz
tcoupl                  = V-rescale
tc-grps                 = Protein Non-Protein
tau_t                   = 0.1 0.1
ref_t                   = 300 300
pcoupl                  = Parrinello-Rahman
pcoupltype             = isotropic
tau_p                   = 2.0
ref_p                   = 1.0
compressibility        = 4.5e-5
EOF

    # Production MD mdp
    cat << EOF > md.mdp
title                   = Production MD
integrator              = md
nsteps                  = 50000  ; 10 ps
dt                      = 0.002  ; 2 fs timestep
nstxout                 = 500    ; Save coordinates every 1 ps
nstvout                 = 500    ; Save velocities every 1 ps
nstfout                 = 500    ; Save forces every 1 ps
nstlog                  = 500    ; Save log every 1 ps
nstenergy               = 500    ; Save energies every 1 ps
nstxout-compressed      = 500    ; Save compressed output every 1 ps
continuation            = yes
constraint_algorithm    = lincs
constraints             = h-bonds
lincs_iter              = 1
lincs_order             = 4
cutoff-scheme           = Verlet
ns_type                 = grid
nstlist                 = 10
rlist                   = 1.0    ; Reduced to improve speed
coulombtype             = PME
rcoulomb                = 1.0    ; Reduced for faster PME
rvdw                    = 1.0    ; Reduced for efficiency
pbc                     = xyz
tcoupl                  = V-rescale
tc-grps                 = Protein Non-Protein
tau_t                   = 0.1 0.1
ref_t                   = 300 300
pcoupl                  = Parrinello-Rahman
pcoupltype              = isotropic
tau_p                   = 2.0
ref_p                   = 1.0
compressibility         = 4.5e-5
gen_vel                 = no
fourierspacing          = 0.16   ; Coarser grid for memory efficiency
fourier_nx              = 64
fourier_ny              = 64
fourier_nz              = 64
EOF

    cd ..
}

# Function to prepare and run all MD steps
setup_and_run_md() {
    local PDB_FILE=$1
    local OUTPUT_DIR=$2
 

    echo "Setting up MD simulation in $OUTPUT_DIR"
    mkdir -p $OUTPUT_DIR
    cd $OUTPUT_DIR || { echo "Error: Cannot change to directory $OUTPUT_DIR"; exit 1; }

    # Copy PDB file
    cp ../$PDB_FILE ./input.pdb || { echo "Error: Failed to copy PDB file $PDB_FILE"; exit 1; }

    # Verify input.pdb exists and has content
    if [ ! -s input.pdb ]; then
        echo "Error: input.pdb is missing or empty"
        exit 1
    fi

    # Generate topology (try different approaches if needed)
    echo "Generating topology..."
    if ! echo "6" | gmx pdb2gmx -f input.pdb -o processed.gro -water tip3p -ff amber99sb-ildn -ignh -ter; then
        echo "First topology generation attempt failed, trying alternative..."
        if ! echo "6" | gmx pdb2gmx -f input.pdb -o processed.gro -water tip3p -ff amber99sb-ildn -ignh; then
            echo "Error: Failed to generate topology"
            exit 1
        fi
    fi

    # Define box
    echo "Setting up simulation box..."
    gmx editconf -f processed.gro -o boxed.gro -c -d 1.2 -bt cubic || { echo "Error: editconf failed"; exit 1; }

    # Solvate
    echo "Adding solvent..."
    gmx solvate -cp boxed.gro -cs spc216.gro -o solvated.gro -p topol.top || { echo "Error: solvate failed"; exit 1; }

    # Add ions
    echo "Adding ions..."
    gmx grompp -f em.mdp -c solvated.gro -p topol.top -o ions.tpr -maxwarn 1 || { echo "Error: grompp (ions) failed"; exit 1; }
    echo "13" | gmx genion -s ions.tpr -o neutralized.gro -p topol.top -pname NA -nname CL -neutral || { echo "Error: genion failed"; exit 1; }

    # Energy minimization
    echo "Running energy minimization..."
    gmx grompp -f em.mdp -c neutralized.gro -p topol.top -o em.tpr -maxwarn 1 || { echo "Error: grompp (EM) failed"; exit 1; }
    gmx mdrun -v -deffnm em || { echo "Error: EM mdrun failed"; exit 1; }

    # NVT equilibration
    echo "Running NVT equilibration..."
    gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr -maxwarn 1 || { echo "Error: grompp (NVT) failed"; exit 1; }
    gmx mdrun -v -deffnm nvt || { echo "Error: NVT mdrun failed"; exit 1; }

    # NPT equilibration
    echo "Running NPT equilibration..."
    gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr -maxwarn 1 || { echo "Error: grompp (NPT) failed"; exit 1; }
    gmx mdrun -v -deffnm npt || { echo "Error: NPT mdrun failed"; exit 1; }

    # Production MD
    echo "Running production MD..."
    gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md.tpr -maxwarn 1 || { echo "Error: grompp (MD) failed"; exit 1; }
    gmx mdrun -v -deffnm md || { echo "Error: Production mdrun failed"; exit 1; }

    cd ..
}

# Create directories if they don't exist
mkdir -p $WT_DIR $MUT_DIR || { echo "Error: Cannot create directories"; exit 1; }

# Create MDP files
create_mdp_files $WT_DIR
create_mdp_files $MUT_DIR

# Run setup and MD for wild type
echo "Processing wild type structure..."
setup_and_run_md $WT_PDB $WT_DIR false

# Run setup and MD for mutant
echo "Processing mutant structure..."
setup_and_run_md $MUT_PDB $MUT_DIR true

echo "MD simulations completed for both wild type and mutant structures."
