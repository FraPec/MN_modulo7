#!/bin/bash

# Check if the user provided the number of processors
if [ $# -ne 1 ]; then
    echo "Usage: $0 <number_of_processors>"
    exit 1
fi

NPROC=$1  # Number of processors provided by user

if [ -d "input_Nk" ]; then
    echo "Folder exists."
else
    mkdir input_Nk
    echo "Folder created."
fi

if [ -d "output_Nk" ]; then
    echo "Folder exists."
else
    mkdir output_Nk
    echo "Folder created."
fi

for NK in 2 4 6 8 10 12 14 16
do
cat > input_Nk/Nk_$NK.in << EOF
&CONTROL
  calculation = 'scf'
  restart_mode = 'from_scratch'
  outdir = './e_bands/'
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
  ecutwfc = 40
  ecutrho = 320 ! approximately 8 x ecutwfc
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
Ge           0.00494000       0.00494000       0.00494000 
Te           0.47406000       0.47406000       0.47406000 
K_POINTS automatic
$NK $NK $NK 0 0 0

EOF

mpirun -np $NPROC ./pw.x -i input_Nk/Nk_$NK.in > output_Nk/Nk_$NK.out

done
