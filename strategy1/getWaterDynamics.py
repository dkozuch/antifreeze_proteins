import numpy as np
import time
import MDAnalysis
from MDAnalysis.analysis.waterdynamics import HydrogenBondLifetimes as HBL
import matplotlib.pyplot as plt
import os
from joblib import Parallel, delayed
import multiprocessing

protein_name = "1ucs_1.1"
b = "30"
e = "40"

def calcHBL(i):
	#water dynamics analysis
	selection1 = "resname SOL and around " + str(max_cutoff) + " (resnum "+ str(i) + ")"
	selection2 = selection1
	
	HBL_analysis = HBL(u, selection1, selection2, 0, 1000, dtmax)
	HBL_analysis.run()
	
	#save to txt file
	out = np.hstack((np.array(range(0,dtmax))[:, np.newaxis], HBL_analysis.timeseries))
	np.savetxt("hbl_" + str(i) + ".txt",out)
		
tpr = protein_name+"_sim.tpr"
trj = protein_name+"_sim_b"+str(b)+"e"+str(e)+"ns.xtc"

u = MDAnalysis.Universe(tpr,trj)
protein = u.select_atoms("protein")
nres = protein.n_residues
print("Number of residues: " + str(nres))

dtmax = 100
min_cutoff = 0
max_cutoff = 8

inputs = range(0,nres)
Parallel(n_jobs=16)(delayed(calcHBL)(i) for i in inputs)



