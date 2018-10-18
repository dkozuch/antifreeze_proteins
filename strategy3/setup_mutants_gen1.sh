#!/bin/bash

gmx=gmx_164_gpu_pd241
scoring=antifreeze_score_v7.40_ev.py

generation=1
prev_gen=$(echo "$generation-1" | bc)
best_prev_mut=0
mkdir $generation

protein=1hg7
mutant=K62V

sim_num=17
mutant_number=$(echo $mutant | sed 's/[^0-9]*//g') #collected number from string
new_res=$(echo "${mutant: -1}") #collects last character of string
new_res=${new_res^^} #converts to uppercase
name=${protein}_${generation}.${sim_num}
folder=${generation}/${sim_num}
echo $folder

rm -r $folder
mkdir $folder
cp run_files/* $folder
cd $folder

#copy gro file and sequence from previous best (gen 0 since this file is only for gen 1)
cp ../../${prev_gen}/${best_prev_mut}/${protein}_${prev_gen}.${best_prev_mut}.fasta.txt ./${protein}.fasta.txt
cp ../../${prev_gen}/${best_prev_mut}/${protein}_sim.gro ./${protein}.gro
cp ../../${prev_gen}/${best_prev_mut}/${protein}_sim.tpr ./${protein}_old.tpr
cp ../../$scoring ./
mpirun -n 1 $gmx trjconv -f ${protein}.gro -s ${protein}_old.tpr -o ${protein}.pdb <<EOF
1
EOF
sed -i 's/OC1/OXT/g' ${protein}.pdb
sed -i 's/OC2/O  /g' ${protein}.pdb

#get sequence
sequence=$(head -n 1 ${protein}.fasta.txt) #read fasta
echo $sequence

#check that we are mutating the right residue
old_res=$(echo "${mutant:0:1}") #get first character
old_res=${old_res^^} #capitalize
mutant_num_bash=$(echo "$mutant_number-1 | bc") #bash starts at 0 but sed starts at 1
old_res_from_seq=$(echo "${sequence:$mutant_num_bash:1}") #old res from fasta
echo "Old res from mut string: "$old_res"; Old res from seq: "$old_res_from_seq
if [ $old_res = $old_res_from_seq ]
then
echo "Sequence match: YES"
else
echo "Error: Mutation location identity does not match mutant name. This should not happen."
fi

#mutate sequence
new_sequence=$(echo $sequence | sed s/./$new_res/$mutant_number)
echo $new_sequence
echo $new_sequence > ${name}.fasta.txt

#generate new pdb
Scwrl4 -i ${protein}.pdb -s ${name}.fasta.txt -o ${name}.pdb

sed -i "/protein=/c\protein=$name" prepare_protein.sh
sh prepare_protein.sh

mv slurm_${protein}_265K.cmd slurm_${name}.cmd
sed -i "/protein=/c\protein=$name" slurm_${name}.cmd

mv pbs_${protein}_265K.pbs pbs_${name}.pbs
sed -i "/protein=/c\protein=$name" pbs_${name}.pbs
sed -i "/#PBS -N/c\#PBS -N $name" pbs_${name}.pbs

#sbatch slurm_${name}.cmd
#qsub pbs_${name}.pbs

cd ../../

