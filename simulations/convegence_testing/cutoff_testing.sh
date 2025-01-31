#!/bin/sh

# Check if the user provided the number of processors
if [ $# -ne 1 ]; then
    echo "Usage: $0 <number_of_processors>"
    exit 1
fi

NPROC=$1  # Number of processors provided by user

# Ensure necessary directories exist
mkdir -p input_ecutoff output_ecutoff

# Loop over ecutoff values
for ECUTOFF in 20 25 30 35 40 45 50 55 60 65 70 75 80
do
    cat > input_ecutoff/ecutoff_$ECUTOFF.in << EOF
&CONTROL
  calculation = 'scf'
  restart_mode = 'from_scratch'
  outdir = './results/'
  prefix = 'GeTe'
  pseudo_dir = '../pseudo/'
  max_seconds = 3600
  nstep = 10000
  disk_io = 'none'
/
&SYSTEM
  ibrav = 5 ! rhombohedral primitive cell
  celldm(1) = 8.267861436 ! Bohr
  celldm(4) = 0.532481505 ! cosine between the two vectors of the Bravais
  nat = 2
  ntyp = 2
  ecutwfc = $ECUTOFF
  ecutrho = $(($ECUTOFF * 8))
  vdw_corr = 'grimme-d2'
/
&ELECTRONS
  conv_thr = 1.0d-8
  electron_maxstep = 10000
/
ATOMIC_SPECIES
Ge     72.64 ge_pbe_v1.4.uspp.F.UPF
Te     127.6 Te_pbe_v1.uspp.F.UPF
ATOMIC_POSITIONS crystal
Ge           0.0049400000       0.0049400000       0.0049400000 
Te           0.4740600000       0.4740600000       0.4740600000 
K_POINTS automatic
8 8 8 0 0 0
EOF

    # Run pw.x with specified processors and wait for completion before next job
    mpirun -np $NPROC ./pw.x -i input_ecutoff/ecutoff_$ECUTOFF.in > output_ecutoff/ecutoff_$ECUTOFF.out
    wait
done

