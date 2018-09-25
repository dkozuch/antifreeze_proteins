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
gmx=gmx_164_gpu_pd241
water_dynamics_file=getWaterDynamics.py
scoring_file=antifreeze_score_v7.40_ev

bns=40
ens=60
kdist=1.1

b=$(echo "$bns*1000" | bc)
e=$(echo "$ens*1000" | bc)

#Remove old logs
rm mdrun_*
rm grompp_*

#minim
echo "Performing minimization"
${gmx} grompp -f minim -c ${protein}_ions.gro -p topol.top -n index.ndx -o ${protein}_min > grompp_min.txt 2>&1
srun ${gmx} mdrun -deffnm ${protein}_min -ntomp $nt > mdrun_min.txt 2>&1
        
sh remove_backups.sh

#eq
echo "Performing equilibration"
${gmx} grompp -f eq -c ${protein}_min.gro -p topol.top -n index.ndx -o ${protein}_eq -maxwarn 1 > grompp_eq.txt 2>&1
srun ${gmx} mdrun -deffnm ${protein}_eq -tunepme yes -ntomp $nt > mdrun_eq.txt 2>&1

sh remove_backups.sh

echo "Performing simulation"
#sim
${gmx} grompp -f sim -c ${protein}_eq.gro -t ${protein}_eq -p topol.top -n index.ndx -o ${protein}_sim > grompp_sim.txt 2>&1
srun ${gmx} mdrun -deffnm ${protein}_sim -tunepme yes -ntomp $nt > mdrun_sim.txt 2>&1

echo "Analyzing simulation"
#get fitted trajectory
$gmx trjconv -f ${protein}_sim -s ${protein}_eq -center -ur compact -pbc mol -b $b -e $e -o ${protein}_sim_b${bns}e${ens}ns <<EOF
2
0
EOF

#get SASA
$gmx sasa -f ${protein}_sim_b${bns}e${ens}ns -s ${protein}_eq -b $b -e $e -or sasa_r_b${bns}e${ens}ns -oa sasa_a_b${bns}e${ens}ns <<EOF
1
EOF

#get average coordinates
$gmx rmsf -f ${protein}_sim_b${bns}e${ens}ns -s ${protein}_eq -ox ${protein}_sim_b${bns}e${ens}ns_avg.pdb <<EOF
1
EOF

sh remove_backups.sh

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


