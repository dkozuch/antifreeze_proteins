#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH --cpus-per-task=1
# --sockets-per-node=1
#SBATCH --gres=gpu:1
#SBATCH --gres-flags=enforce-binding
#SBATCH -t 24:00:00
# sends mail when process begins, and
# when it ends. Make sure you define your email
# address.
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-type=fail
#SBATCH --mail-user=dkozuch@princeton.edu


#Load necessary modules
module load openmpi
module load intel

nt=1
export OMP_NUM_THREADS=$nt

protein=1hg7
gmx=gmx_183
water_dynamics_file=getWaterDynamics.py
scoring_file=antifreeze_score_v7.40_ev

bns=40
ens=60
kdist=1.1

b=$(echo "$bns*1000" | bc)
e=$(echo "$ens*1000" | bc)

#get water dynamics
cutoff=8
hbl_folder=HBL_0t${cutoff}A_b${bns}e${ens}ns
sed -i "/^protein_name = /c\protein_name = "\""${protein}"\""" $water_dynamics_file
sed -i "/^b =/c\b = "\""${bns}"\""" $water_dynamics_file
sed -i "/^e =/c\e = "\""${ens}"\""" $water_dynamics_file
sed -i "/^cutoff =/c\cutoff = "${cutoff}"" $water_dynamics_file
echo "Running "$water_dynamics_file" ..."
python $water_dynamics_file > tmp_waterDynamics.log 2>&1
mkdir $hbl_folder
mv hbl_* $hbl_folder

sh remove_backups.sh
echo "Done"

#get score
python ${scoring_file}.py ${protein} $kdist $bns $ens > log_${scoring_file}_${protein}.txt


