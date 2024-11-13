#!/bin/bash

# Prepare topology for all starting structures
gmx_mpi pdb2gmx -f trmd_af_monomer.pdb -p topol.top -i posre.itp -ignh -o processed.pdb<< EOF
1
1
EOF

gmx_mpi editconf -f processed.pdb -o centered.pdb -c -d 1.2 -bt dodecahedron
gmx_mpi solvate -cp centered.pdb -o solv.pdb -cs spc216.gro -p topol.top
gmx_mpi grompp -f ../ions.mdp -c solv.pdb -p topol.top -o ions.tpr
gmx_mpi genion -s ions.tpr -o solv_ions.pdb -p topol.top -pname MG -nname CLA -neutral -pq 2 -conc 0.012 << EOF
SOL
EOF
gmx_mpi grompp -f ../minim.mdp -c solv_ions.pdb -p topol.top -o em.tpr


