import numpy as np
import os
np.set_printoptions(suppress=True) #suppress scientific notation in numpy printing
#import matplotlib.pyplot as plt

protein="1ucs"
log_file = "log_antifreeze_score_v7.20_" + protein + "_"

#get fitted scores from log files
scores = []
for gen in range(0,9):
	if gen == 0:
		file = str(gen) + "/" + log_file + str(gen) + ".txt"
		with open(file,'r') as f:
			for line in f:
				if line.split()[0:2] == ["Fitted","Score:"]:
					score = float(line.split()[2])
					score_zero = score
	else:
		gen_scores = []
		for j in range(1,18):
			check = 0
			file = str(gen) + "/" + str(j) + "/" + log_file + str(gen) + "." + str(j) + ".txt"
			with open(file,'r') as f:
				for line in f:
					if line.split()[0:2] == ["Fitted","Score:"]:
						check = 1
						score = float(line.split()[2])
						gen_scores.append(score)
			if check == 0:
				print "Warning: no fitted score found for Gen: " + str(gen) + " Mut: " + str(j) + "; using score of 0"
				gen_scores.append(0)
		scores.append(gen_scores)

#find best scores
scores = np.array(scores)
np.savetxt("scores.txt",scores)
max_scores = np.max(scores,axis=1)
max_scores = np.insert(max_scores, 0, score_zero) #add gen 0 to list
print("Max scores:" + str(max_scores))
#plt.scatter(range(0,len(max_scores)),max_scores)
#plt.show()

#find best mutants
max_pos = np.argmax(scores, axis=1) + 1 #need to start at 1, not 0
print("Max mutants: " + str(max_pos)) 

#get data for best mutants
native = "0/" + log_file + "0.txt"
with open(native,'r') as f:
	for line in f:
		if line.split()[0:6] == ["Number", "of", "residues", "in", "binding", "face:"]:
			nr = float(line.split()[6])
		if line.split()[0:2] == ["Binding","HBL:"]:
			hblb = float(line.split()[2])
		if line.split()[0:2] == ["Nonbinding","HBL:"]:
			hbln = float(line.split()[2])
		if line.split()[0:2] == ["Projected","area:"]:
			area = float(line.split()[2])
		if line.split()[0:2] == ["Fitted","Score:"]:
			score = float(line.split()[2])
max_fitness = score
best_gen = 0
data = [[0,nr,hblb,hbln,area,score,max_fitness]]
for i in range(0,len(max_pos)):
	gen = i+1
	mutant = max_pos[i]
	file = str(gen) + "/" + str(mutant) + "/" + log_file + str(gen) + "." + str(mutant) + ".txt"
	#print file
	with open(file,'r') as f:
		for line in f:
			if line.split()[0:6] == ["Number", "of", "residues", "in", "binding", "face:"]:
				nr = float(line.split()[6])
			if line.split()[0:2] == ["Binding","HBL:"]:
				hblb = float(line.split()[2])
			if line.split()[0:2] == ["Nonbinding","HBL:"]:
				hbln = float(line.split()[2])
			if line.split()[0:2] == ["Projected","area:"]:
				area = float(line.split()[2])
			if line.split()[0:2] == ["Fitted","Score:"]:
				score = float(line.split()[2])
	if score >= max_fitness:
		max_fitness = score
		best_gen = gen
	data.append([gen,nr,hblb,hbln,area,score,max_fitness])
print("Max Fitness: " + str(max_fitness))
print("Best Gen/Mut: " + str(best_gen) + "/" + str(max_pos[best_gen-1]))

data = np.array(data)
print("Data Table: GEN, NR, HBLB, HBLN, A, S, MAX S")
print(data)
np.savetxt("log_data.txt",data)
#add label to data file
os.system('echo "# Generation, Number of Residues in Binding Face, HBL of Binding Face, HBL of Nonbinding Face, Area of Binding Face, Score, Max Score"|cat - log_data.txt > /tmp/out && mv /tmp/out log_data.txt')


#get sequences for native and best mutants
native_fasta = "0/" + str(protein) + ".fasta.txt"
with open(native_fasta,'r') as f:
	for line in f:
		native_seq = list(line)[:-1]
seqs = [native_seq]
for i in range(0,len(max_pos)):
	gen = i+1
	mutant = max_pos[i]
	fasta_file = str(gen) + "/" + str(mutant) + "/" + protein + "_" + str(gen) + "." + str(mutant) + ".fasta.txt"
	with open(fasta_file,'r') as f:
		for line in f:
			seq = list(line)[:-1]
	seqs.append(seq)
seqs = np.array(seqs)
#print seqs

#print mutations made
mut_list=[] 
for i in range(1,len(seqs[:,0])):
        for j in range(0,len(seqs[i,:])):
                if seqs[i,j] != seqs[i-1,j]:
                        mut_list.append([i,j,seqs[i-1,j],seqs[i,j]])
mut_list = np.array(mut_list)
print("Mutation list: " + str(mut_list))
np.savetxt("mutations.txt",mut_list,fmt="%s")

#get final mutations for best overall mutant
print("Mutations needed for best overall mutant: ")
for i in range(0,len(seqs[0])):
	if seqs[0,i] != seqs[best_gen,i]:
		print [i,seqs[0,i],seqs[best_gen,i]]







			
