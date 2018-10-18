print "Loading packages..."
import numpy as np
import os
from joblib import Parallel, delayed
import ast
import progressbar
from scipy.interpolate import griddata
import cPickle as pickle
np.set_printoptions(suppress=True) #suppress scientific notation in numpy printing
import matplotlib.pyplot as plt
print "Packages loaded..."
import sys


#define global
protein="1hg7"
log_file = "log_antifreeze_score_v7.40_ev_" + protein + "_"
ngens=int(sys.argv[1])
aa_lib = ["A","C", "D", "E", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V"] #no F,W,Y so we can run on 18 node therminator (aromatics too big anyway)
nmuts = len(aa_lib)
#import neural net data
#nn_data = np.load("nn_output.npy")

def save_object(obj, filename):
    with open(filename, 'wb') as output:  # Overwrites any existing file.
        pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)
		
def load_saved_object(filename):
	with open(filename, 'rb') as input:
		obj = pickle.load(input)
	return obj

class mutant:

	'''This class should hold all data and functions for a particular mutant'''

	def __init__(self,gen,mut):
		self.gen = gen
		self.mut = mut
		self.log_filename = str(gen)+"/"+str(mut)+"/"+log_file+str(gen)+"."+str(mut)+".txt"
		self.seq_filename = str(gen)+"/"+str(mut)+"/"+protein+"_"+str(gen)+"."+str(mut)+".fasta.txt"
			
	def get_log_data(self):
		check = 0
		with open(self.log_filename,'r') as f:
			for line in f:
				linesplit = line.split(": ")
				if len(linesplit) == 2:
					id, val = linesplit
					if id == "Residues (index) in binding face":
						self.binding_residues = ast.literal_eval(val) #ast converts string representation of list to actual list
					if id == "Binding HBL":
						self.lb = float(val)
					if id == "Nonbinding HBL":
						self.ln = float(val)
					if id == "Projected area":
						self.area = float(val)
					if id == "Fitted Score":
						self.fit_score = float(val)
						check = 1
		#if no binding face
		if check == 0:
			self.fit_score = 0
	
	#already score with nn					
	#def get_nn_score(self):
	#	point = [self.area,self.lb,self.ln]
	#	self.nn_score = np.max([0,griddata(nn_data[:,0:3],nn_data[:,3],point)]) #don't predict negative thermal hysteresis
	#	if np.isnan(self.nn_score): #convert nan to zero
	#		self.nn_score = 0
			
	def get_sequence(self):
		'''Get fasta seq from file
		sequence should be string in first line of file
		cut white space characters with rstrip'''
		with open(self.seq_filename,'r') as f:
			self.seq = list(f.readline().rstrip())
			
	def get_mutations(self,ref):
		mutations = []
		for i in range(0,len(self.seq)):
			if self.seq[i] != ref[i]:
				mutations.append([i,ref[i],self.seq[i]])
		self.mutations = mutations
				
def fill_mutant(gen,mut,get_mutations=False,ref=None):
	'''create single mutant object and load its data'''
	muti = mutant(gen,mut)
	muti.get_log_data()
	#muti.get_nn_score() #don't need this, already scored with nn
	muti.get_sequence()
	if get_mutations:
		muti.get_mutations(ref)
	return muti
					
def load_mutants(ngens,nmuts):
	mutants = []
	for gen in progressbar.progressbar(range(0,ngens+1)):
		#print "Loading mutants for gen: "+str(gen)
		if gen == 0:
			gen_mutants = [fill_mutant(gen,0)]
			ref = gen_mutants[0].seq
		else:
			inputs = range(1,nmuts+1)
			#parallelize loading mutants since takes forever with neural network evaluation
			#gen_mutants = Parallel(n_jobs=8)(delayed(fill_mutant)(gen,i,get_mutations=True,ref=ref) for i in inputs)
			gen_mutants = [ fill_mutant(gen,i,get_mutations=True,ref=ref) for i in inputs ]
		mutants.append(gen_mutants)
	return mutants
	
	
mutants = load_mutants(ngens,nmuts)
save_object(mutants,"analysis_obj_mutants.pkl")

#mutants = load_saved_object("analysis_obj_mutants.pkl")

#get max score per gen
max_scores = []
for i in range(0,len(mutants)):
	gen_fit_scores = [ x.fit_score for x in mutants[i] ]
	gen_max_score = np.max(gen_fit_scores)
	max_scores.append([i,gen_max_score])
max_scores = np.array(max_scores)

#track max score overall
max_scores_overall = [max_scores[0]]
for i in range(1,len(max_scores)):
	if max_scores[i,1] > max_scores_overall[i-1][1]:
		max_scores_overall.append([i,max_scores[i,1]])
	else:
		max_scores_overall.append([i,max_scores_overall[i-1][1]])
max_scores_overall = np.array(max_scores_overall)

plt.plot(max_scores[:,0],max_scores[:,1],'o',markersize=7)
plt.plot(max_scores_overall[:,0],max_scores_overall[:,1],'--o',markersize=1)
plt.xlabel("Generation")
plt.ylabel("Predicted $\Delta T_C$")
plt.legend(["Per generation","Overall"])
plt.show()


#look scores associated with mutating to each residue

		
	
			
		
		
			
