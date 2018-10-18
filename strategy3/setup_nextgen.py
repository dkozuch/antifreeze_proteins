import numpy as np
import os
import sys

protein = "1hg7"
python_file="log_antifreeze_score_v7.40_ev"

exec_flag = False
if len(sys.argv) > 2:
	if sys.argv[2] == "execute":
		exec_flag = True

next_gen = int(sys.argv[1])
num_mut = 17
os.system("sed -i '/generation=/c\generation=" + str(next_gen) + "' setup_mutants.sh")

aa_lib = ["A","C", "D", "E", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V"] #no F,W,Y so we can run on 18 node therminator (aromatics too big anyway)
num_mut = len(aa_lib)
#print aa_lib

def log_name(gen,mut=False,python_file=python_file,protein=protein):
	log_file = str(gen)+"/"+str(mut)+"/"+python_file +"_"+protein+"_"+str(gen)+"."+str(mut)+".txt"
	return log_file

def get_best_mutant(gen):
	print("Analyzing Gen: " + str(gen))
	scores = []
	if gen == 0:
		check = 0
		with open(log_name(0,0),'r') as f:
			for line in f:
				if " ".join(line.split()[0:2]) == "Fitted Score:":
					check = 1
					scores.append([1,float(line.split()[2])])
		if check == 0:
			print "WARNING: No score found for Gen " + str(gen) + ". Using score of 0"
			scores.append([i,0])
	else:
		for i in range(1,num_mut+1):
			check = 0
			with open(log_name(gen,i),'r') as f:
				for line in f:
					if " ".join(line.split()[0:2]) == "Fitted Score:":
						check = 1
						scores.append([i,float(line.split()[2])])
			if check == 0:
				print "WARNING: No score found for Gen " + str(gen) + " Mut " + str(i) + ". Using score of 0"
				scores.append([i,0])
	scores = np.array(scores)
	#print "Scores:"
	#print scores
	best_mut = int(scores[np.argmax(scores[:,1]),0])
	#print("Best mutant: " + str(best_mut))
	#print("Best score: " + str(np.max(scores[:,1])))
	return [best_mut,np.max(scores[:,1])]

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
			mutant = str(fasta[mut_resid]) + str(mut_resid+1) + str(i)  #add one since bash starts at 1
			print mutant
			os.system("sed -i '/mutant=/c\mutant=" + mutant + "' setup_mutants.sh")
			os.system("sed -i '/sim_num=/c\sim_num=" + str(count) + "' setup_mutants.sh")
			if exec_flag == True:
				os.system("sh setup_mutants.sh")
			count = count + 1


#######################################

#check all previous generations for best score
scores = []
for i in range(0,next_gen):
	scores.append(get_best_mutant(i))
scores = np.array(scores)
print scores
best_gen = int(np.argmax(scores[:,1]))
best_mut = int(scores[best_gen,0])
print("Best generation: " + str(best_gen) + ", Best mutant: " + str(best_mut))
os.system("sed -i '/prev_gen=/c\'prev_gen'=" + str(best_gen) + "' setup_mutants.sh")
os.system("sed -i '/best_prev_mut=/c\'best_prev_mut'=" + str(best_mut) + "' setup_mutants.sh")

#need sorted mutation residues
log_file = log_name(best_gen,mut=best_mut)
mut_resids = get_binding_resids(log_file)
print("Possible mutation locations (sorted): ")
print mut_resids

#get all tried sequences
seqs = []
for i in range(0,next_gen):
	if i == 0:
		fasta_file = str(i)+"/0/"+protein+"_0.0.fasta.txt"
		seq = get_fasta_seq(fasta_file)
		seqs.append(seq)
	else:
		for j in range(1,num_mut+1):
			fasta_file = str(i)+"/"+str(j)+"/"+ protein+"_"+str(i)+"."+str(j)+".fasta.txt"
			seq = get_fasta_seq(fasta_file)
			seqs.append(seq)

#get fasta for best
if best_gen != 0:
	best_fasta_file = str(best_gen)+"/"+str(best_mut)+"/"+ protein+"_"+str(best_gen)+"."+str(best_mut)+".fasta.txt"
else:
	best_fasta_file = str(best_gen)+"/0/"+ protein+"_0.0.fasta.txt"
best_seq = get_fasta_seq(best_fasta_file)
	
#test if mutation position already tried - THIS ASSUMES that ALL residues in aa_lib are tried EACH TIME
#not gaureneteed to work if this is not the case
for i in mut_resids:
	new_seq = list(best_seq)
	check = 0
	for j in aa_lib:
		if j != new_seq[i]:
			new_seq[i] = j
			check = 1
		if check == 1:
			break
	if new_seq in seqs:
		print "Mutation for original Gen " + str(best_gen) + " at residue " + str(i) + " already attempted. Proceeding to next position..."
	else:
		print "Mutation for original Gen " + str(best_gen) + " at residue " + str(i) + " is unique; using this position for next generation."
		mut_resid = i
		break

print("Mutants... (index +1 for bash)")
generate_mutants(mut_resid,best_seq,aa_lib)

