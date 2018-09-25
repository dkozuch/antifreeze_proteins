#!/bin/bash

protein=1hg7
gmx=gmx_164_gpu_pd241

#build box and pdb2gmx (4 for tip4p ice)
${gmx} editconf -f ${protein}.pdb -resnr 0 -o ${protein}.pdb
${gmx} editconf -f ${protein}.pdb -d 1.2 -bt dodecahedron -o ${protein}_box.pdb
${gmx} pdb2gmx -f ${protein}_box.pdb -o ${protein}_box.gro -ignh -p topol -ff amber03w <<EOF
7
EOF

#solvate protein
${gmx} solvate -cp ${protein}_box.gro -cs tip4p2005.gro -p topol.top -o ${protein}_sol.gro

#add ions (13 specifies ions should replace SOL molecules)
${gmx} grompp -f ions.mdp -c ${protein}_sol.gro -p topol.top -o ions.tpr
${gmx} genion -s ions.tpr -o ${protein}_ions.gro -p topol.top -pname NA -nname CL -neutral yes <<EOF
13
EOF

#make index file
${gmx} make_ndx -f ${protein}_ions.gro -o index.ndx <<EOF
q
EOF

#minimize
${gmx} grompp -f minim.mdp -c ${protein}_ions.gro -p topol.top -o ${protein}_min
${gmx} trjconv -f ${protein}_ions.gro -s ${protein}_min -ur compact -center -pbc mol -o ${protein}_compact.gro <<EOF
2
0
EOF
#mpirun -n 4 ${gmx} mdrun -v -deffnm ${protein}_min

sh remove_backups.sh
