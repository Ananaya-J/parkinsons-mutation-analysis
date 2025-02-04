#!/bin/bash

# Load GROMACS (adjust this based on your installation)
source /usr/local/gromacs/bin/GMXRC

# Define directories and files
#WT_DIR="ARHGAP5_WT_MD"
MUT_DIR="ARHGAP5_MUT_MD"
#WT_PDB="ARHGAP5_WT.pdb"
MUT_PDB="truncated.pdb"

# Create directories
mkdir -p $WT_DIR $MUT_DIR

# Function to create necessary MDP files
create_mdp_files() {
    local OUTPUT_DIR=$1
    cd $OUTPUT_DIR

    cat << EOF > em.mdp
integrator = steep
nsteps = 50000
emtol = 1000.0
emstep = 0.01
nstlist = 10
cutoff-scheme = Verlet
rlist = 1.0
coulombtype = PME
rcoulomb = 1.0
rvdw = 1.0
pbc = xyz
EOF

    cat << EOF > nvt.mdp
title                   = NVT equilibration
integrator              = md
nsteps                  = 25000  ; 50 ps
dt                      = 0.002  ; 2 fs timestep
nstxout                 = 500    ; Save coordinates every 1 ps
nstvout                 = 500    ; Save velocities every 1 ps
nstfout                 = 500    ; Save forces every 1 ps
nstlog                  = 500    ; Save log every 1 ps
nstenergy               = 500    ; Save energies every 1 ps
nstxout-compressed      = 500    ; Save compressed output every 1 ps
continuation            = no
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
EOF

    cat << EOF > npt.mdp
title                   = NPT equilibration
integrator              = md
nsteps                  = 25000  ; 50 ps
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
EOF

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

# Function to set up and run MD simulation
setup_and_run_md() {
    local PDB_FILE=$1
    local OUTPUT_DIR=$2

    echo "Setting up MD simulation in $OUTPUT_DIR"
    cd $OUTPUT_DIR || exit 1



    # Generate topology
    gmx pdb2gmx -f truncated.pdb -o processed.gro -water spce -ff amber99sb-ildn

    # Define box
    gmx editconf -f processed.gro -o boxed.gro -c -d 1.0 -bt cubic

    # Solvate system
    gmx solvate -cp boxed.gro -cs spc216.gro -o solvated.gro -p topol.top

    # Generate ions.tpr
    gmx grompp -f em.mdp -c solvated.gro -p topol.top -o ions.tpr -maxwarn 1

    # Add counterions (Select group 13 for SOL)
    echo "13" | gmx genion -s ions.tpr -o neutralized.gro -p topol.top -pname NA -nname CL -neutral

    # Prepare energy minimization
    gmx grompp -f em.mdp -c neutralized.gro -p topol.top -o em.tpr
    gmx mdrun -v -deffnm em

    # NVT Equilibration
    gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
    gmx mdrun -v -deffnm nvt

    # NPT Equilibration
    gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
    gmx mdrun -v -deffnm npt

    # Production MD
    gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md.tpr
    gmx mdrun -v -deffnm md

    echo "MD simulation completed for $OUTPUT_DIR"
    cd ..
}

# Create MDP files in both directories
#create_mdp_files $WT_DIR
create_mdp_files $MUT_DIR

# Run MD setup and simulation for both wild-type and mutated proteins
#setup_and_run_md $WT_PDB $WT_DIR
setup_and_run_md $MUT_PDB $MUT_DIR

echo "All MD simulations completed successfully."
