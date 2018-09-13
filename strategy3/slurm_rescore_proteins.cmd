#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH -t 1:00:00
# sends mail when process begins, and
# when it ends. Make sure you define your email
# address.
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-type=fail
#SBATCH --mail-user=dkozuch@princeton.edu

module load anaconda
module load intel
module load openmpi

ngens=25
nmuts=17
protein=1ucs
python_file=antifreeze_score_v7.32_ev.py
kd=1.1
b=30
e=60

for i in $(seq 2 $ngens); do

	for j in $(seq 1 $nmuts); do
	
		echo "Scoring Gen."${i}" Mut. "${j}
		dir=${i}/${j}/
		proteinij=${protein}_${i}.${j}
		cp $python_file $dir
		cd $dir
		python $python_file $proteinij $kd $b $e > log_${python_file}_${proteinij}.txt
		
		cd ../../
		
	done
	
done
