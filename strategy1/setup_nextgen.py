import numpy as np
import os
import sys

exec_flag = False
if len(sys.argv) > 2:
	if sys.argv[2] == "execute":
		exec_flag = True

next_gen = int(sys.argv[1])
prev_gen = next_gen - 1
protein="1ucs"

os.system("sed -i '/generation=/c\generation=" + str(next_gen) + "' setup_mutants.sh")

aa_lib = ["A","C", "D", "E", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V"] #no F,W,Y so we can run on 18 node therminator (aromatics too big anyway)
num_mut = len(aa_lib)
print aa_lib

def get_best_mutant():
	scores = []
	for i in range(1,num_mut+1):
		check = 0
		with open(str(prev_gen)+"/"+str(i)+"/"+"log_antifreeze_score_v7.20_"+protein+"_"+str(prev_gen)+"."+str(i)+".txt",'r') as f:
			for line in f:
				if " ".join(line.split()[0:2]) == "Fitted Score:":
					check = 1
					scores.append([i,float(line.split()[2])])
		if check == 0:	
			print "WARNING: No score found for Gen " + str(prev_gen) + " Mut " + str(i) + ". Using score of 0"
			scores.append([i,0])
	scores = np.array(scores)
	print "Scores:"
	print scores
	best_mut = int(scores[np.argmax(scores[:,1]),0])
	print("Best mutant: " + str(best_mut))
	print("Best score: " + str(np.max(scores[:,1])))
	return best_mut

def get_binding_resids(log_file):
	with open(log_file,"r") as f:
		for line in f:
			if " ".join(line.split()[0:4]) == "Residue selected for mutation:":
				mut_resid = int(line.split()[-1])
	return mut_resid

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
			mutant = str(fasta[mut_resid]) + str(mut_resid+1) + str(i)  #add one since bash starts at 1
			print mutant
			os.system("sed -i '/mutant=/c\mutant=" + mutant + "' setup_mutants.sh")
			os.system("sed -i '/sim_num=/c\sim_num=" + str(count) + "' setup_mutants.sh")
			if exec_flag == True:
				os.system("sh setup_mutants.sh")
			count = count + 1


#######################################
best_mut = get_best_mutant()
os.system("sed -i '/best_prev_mut=/c\'best_prev_mut'=" + str(best_mut) + "' setup_mutants.sh")

log_file = str(prev_gen)+"/"+str(best_mut)+"/"+"log_antifreeze_score_v7.20_"+protein+"_"+str(prev_gen)+"."+str(best_mut)+".txt"
mut_resid = get_binding_resids(log_file)
print("New mutation location: " + str(mut_resid+1)) #use index 1 because this is what bash will use

fasta_file = str(prev_gen)+"/"+str(best_mut)+"/"+ protein+"_"+str(prev_gen)+"."+str(best_mut)+".fasta.txt"
fasta = get_fasta_seq(fasta_file)
print fasta

generate_mutants(mut_resid,fasta,aa_lib)

