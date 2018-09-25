import numpy as np
import os
import sys

exec_flag = False
if len(sys.argv) > 1:
	if sys.argv[1] == "execute":
		exec_flag = True

log_file = "0/log_antifreeze_score_v7.40_ev_1hg7.txt"
fasta_file = "0/1hg7.fasta.txt"

aa_lib = ["A", "C", "D", "E", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V"] #no F,W,Y so we can run on 18 node therminator (aromatics too big anyway)
print aa_lib

def get_binding_resids(log_file):
	with open(log_file,"r") as f:
		for line in f:
			if " ".join(line.split()[0:4]) == "Sorted residues for mutation:":
				mut_resid_string = line.split(": ")[-1]
				mut_resids = [int(x) for x in mut_resid_string.split(", ")]
	return mut_resids

def get_fasta_seq(fasta_file):
	with open(fasta_file,"r") as f:
		for line in f:
			fasta = list(line)
	fasta = fasta[:-1] #remove "\n" at end of line
	return fasta
	
def generate_mutants(mut_resid,fasta,aa_lib):
	count = 1
	for i in aa_lib:
		if mut_resid != i:
			mutant = str(fasta[mut_resid]) + str(mut_resid+1) + str(i)  #add one since sed starts at 1
			print mutant
			os.system("sed -i '/mutant=/c\mutant=" + mutant + "' setup_mutants_gen1.sh")
			os.system("sed -i '/sim_num=/c\sim_num=" + str(count) + "' setup_mutants_gen1.sh")
			if exec_flag == True:
				os.system("sh setup_mutants_gen1.sh")
			count = count + 1


#######################################
mut_resids = get_binding_resids(log_file)
mut_resid = mut_resids[0]
print("Mutation location: " + str(mut_resid))

fasta = get_fasta_seq(fasta_file)
print fasta

generate_mutants(mut_resid,fasta,aa_lib)

