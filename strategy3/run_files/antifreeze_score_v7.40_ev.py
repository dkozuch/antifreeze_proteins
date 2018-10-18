'''Comments
Algorithm for scoring the antifreeze activity of a protein using geometry and hydrogen bond lifetime
Same as v20 except returns list of possible mutations ORDERED by hydrogen bond lifetime (smallest to largest)
'''

print("Importing necessary packages...")

import numpy as np
import math
import os
import sys
from joblib import Parallel, delayed
import time
from sklearn.decomposition import PCA
from scipy.spatial import ConvexHull
import subprocess
from scipy.interpolate import griddata
#from scipy import interpolate
#from scipy.optimize import curve_fit

print("All packages loaded")
print("Using analysis file: " + sys.argv[0])

processors=16

protein=sys.argv[1]
begin=sys.argv[3]
end=sys.argv[4]
print("Performing analysis for protein: " + protein)
pdb_file = protein+"_sim_b"+str(begin)+"e"+str(end)+"ns_avg.pdb"
sasa_r_file = "sasa_r_b"+str(begin)+"e"+str(end)+"ns.xvg"
sasa_a_file = "sasa_a_b"+str(begin)+"e"+str(end)+"ns.xvg"

k_neighbor = 20 #just for speeding up search
k_min_neighbors = 2 #for ensuring contiguous binding face
k_neighbor_near = 10 #for ensuring contiguous binding face
k_plane_dist = float(sys.argv[2]) #for finding res in binding face
k_split_dist = k_plane_dist #for finding res in binding face
k_in_plane = 8 #for ensuring binding face bigg enough
k_split = 0
k_sasa_r = 0.01
k_sasa_a = 0.001

#constants from fitting
c1 = 1.01275
c2 = 0.474931
c3 = 0.0431212
c4 = -0.0643819

hbl_folder="HBL_0t8A_b"+str(begin)+"e"+str(end)+"ns"
hbl_cutoff = 0.1
hbl_function = 1 #1 is continuous, 2 is intermitant

print("# Time start: " + begin)
print("# Time end: " + end)
print("# k_neighbor: " + str(k_neighbor))
print("# k_min_neighbors: " + str(k_min_neighbors))
print("# k_neighbor_near: " + str(k_neighbor_near))
print("# k_plane_dist: " + str(k_plane_dist))
print("# k_split_dist: " + str(k_split_dist))
print("# k_in_plane: " + str(k_in_plane))
print("# k_split: " + str(k_split))
print("# k_sasa_r: " + str(k_sasa_r))
print("# k_sasa_a: " + str(k_sasa_a))

print("# Score constants: " + str([c1,c2,c3,c4]))

print("# hbl_folder: " + str(hbl_folder))
print("# hbl_cutoff: " + str(hbl_cutoff))

def pad(l,n):
	'''l is a list of a string, n is the allowed length of the resulting string'''
	if len(l) == n:
		return l
	if len(l) < n:
		pads = " "*(n - len(l))
		return pads + l
	if len(l) > n:
		print "Error: string longer that allotted space"
		return "Error: string longer that allotted space"
	
def round_special(n,l):
	'''n is the number being rounded, l is the number of characters allowed in the number'''
	length = len(list(str(int(round(n,0)))))
	if n > -1 and n < 0: #account for fact that -0.42 rounds to 0 instead of -0 
		dec = l - length - 2
	else:
		dec = l - length - 1
	return "".join(list(str(round(n,dec)))[:l])
	
def distance(x0, x1, box):
	'''distance between two 3D coordinates taking into account periodic conditions'''
	delta = np.abs(x0 - x1)
	delta = np.where(delta > 0.5 * box, delta - box, delta)
	return np.sqrt((delta ** 2).sum(axis=-1))	

def angle(v1, v2):
	'''Returns the angle in radians between vectors 'v1' and 'v2' '''
	cosang = np.dot(v1, v2)
	sinang = np.linalg.norm(np.cross(v1, v2))
	return np.arctan2(sinang, cosang)
	
def point_plane_dist(point,plane):
	'''calculates minimimum distance from point to a plane; 
	return SIGNED number (useful for deciding which side of a plane a point is on'''
	dist = (np.sum(plane[0:3]*point)+plane[3])/float(np.linalg.norm(plane[0:3]))
	return dist
	
def raw_PCA(data, correlation = False, sort = True):
	'''Does PCA on data, currently used for plane fitting
	https://stackoverflow.com/questions/38754668/plane-fitting-in-a-3d-point-cloud'''
	
	mean = np.mean(data, axis=0)
	data_adjust = data - mean
	#: the data is transposed due to np.cov/corrcoef syntax
	if correlation:
		matrix = np.corrcoef(data_adjust.T)
	else:
		matrix = np.cov(data_adjust.T) 
	eigenvalues, eigenvectors = np.linalg.eig(matrix)
	if sort:
		#: sort eigenvalues and eigenvectors
		sort = eigenvalues.argsort()[::-1]
		eigenvalues = eigenvalues[sort]
		eigenvectors = eigenvectors[:,sort]
	return eigenvalues, eigenvectors

def best_fitting_plane(points, equation=False, calc_error=False):
	'''fit a plane to 3D points, requires PCA function defined above'''	
	w, v = raw_PCA(points)
	# the normal of the plane is the last eigenvector
	normal = v[:,2]
	# get a point from the plane
	point = np.mean(points, axis=0)
	if equation and calc_error:
		a, b, c = normal
		d = -(np.dot(normal, point))
		error = 0
		for i in points:
			dist = abs(point_plane_dist(i,np.array([a,b,c,d])))
			error = error + dist
		error = error/float(np.shape(points)[0])
		return np.array([a, b, c, d, error])
	elif equation:
		a, b, c = normal
		d = -(np.dot(normal, point))
		if np.any(np.iscomplex(np.array([a, b, c, d]))):
			print("Error: encountered imaginary number fitting plane to points: " + str(points))
			sys.exit()
		return np.array([a, b, c, d])
	else:
		return [point, normal]		
		
def get_surface_resids(sasa_file):
	'''returns the surface residues, calcualted from the gromacs SASA funciton'''

	#read sasa file (gromacs)
	comment_characters = ["#","@","@TYPE"]
	sasa_list = []
	with open(sasa_file,"r") as f:
			for line in f:
					lineSplit = line.split()
					if lineSplit[0] not in comment_characters:
							sasa_list.append([int(lineSplit[0]),float(lineSplit[1])])
	sasa_list = np.array(sasa_list)

	#list of residues at least equal to k_sasa
	#introduce functionality to accept residue numbers starting at 1
	surf_resids = []
	for i in range(0,np.shape(sasa_list)[0]):
			if sasa_list[i,1] >= k_sasa_r:
				surf_resids.append(sasa_list[i,0])
				
	return [surf_resids,sasa_list]
	
def get_surface_atoms(sasa_file):
	'''returns the surface atoms, calcualted from the gromacs SASA funciton'''
	#read sasa file (gromacs)
	comment_characters = ["#","@","@TYPE"]
	sasa_list = []
	with open(sasa_file,"r") as f:
			for line in f:
					lineSplit = line.split()
					if lineSplit[0] not in comment_characters:
							sasa_list.append([int(lineSplit[0]),float(lineSplit[1])])
	sasa_list = np.array(sasa_list)

	#list of atoms at least equal to k_sasa
	surf_atoms = []
	for i in range(0,np.shape(sasa_list)[0]):
		if sasa_list[i,1] >= k_sasa_a:
			surf_atoms.append(sasa_list[i,0])
			
	return surf_atoms
	
def edit_pdb_sasa(sasa_file,surf_atoms,pdb_file):
	'''edit pdb file with occupancy column being sasa'''
	fileout = protein + "_sasa_a_b"+str(begin)+"e"+str(end)+"ns.pdb"
	os.system("rm " + fileout)
	
	#read sasa file (gromacs)
	comment_characters = ["#","@","@TYPE"]
	sasa_list = []
	with open(sasa_file,"r") as f:
			for line in f:
					lineSplit = line.split()
					if lineSplit[0] not in comment_characters:
							sasa_list.append([int(lineSplit[0]),float(lineSplit[1])])
	sasa_list = np.array(sasa_list)

	with open(pdb_file, "r") as f:
		for line in f:
			if line.split()[0] == "ATOM":
				atom_num = int(line.split()[1])
				if atom_num in surf_atoms:
					tag = sasa_list[np.where(sasa_list[:,0] == atom_num)[0][0],1]
					tag = round_special(tag,5)
				else:
					tag = round_special(-1,5)
				linelist=list(line)
				linelist[56:59] = list(pad(tag,5))
				newline = "".join(linelist)
				with open(fileout,"a") as w:
					w.write(newline)
			else:
				with open(fileout,"a") as w:
					w.write(line)
					
def edit_pdb_sasa_r(sasa_file,pdb_file):
	'''edit pdb file with occupancy column being sasa of residues'''
	
	#read sasa file (gromacs)
	comment_characters = ["#","@","@TYPE"]
	sasa_list = []
	with open(sasa_file,"r") as f:
			for line in f:
					lineSplit = line.split()
					if lineSplit[0] not in comment_characters:
							sasa_list.append([int(lineSplit[0]),float(lineSplit[1])])
	sasa_list = np.array(sasa_list)
	
	fileout = protein + "_sasa_r_b"+str(begin)+"e"+str(end)+"ns.pdb"
	os.system("rm " + fileout)

	with open(pdb_file, "r") as f:
		for line in f:
			if line.split()[0] == "ATOM":
				resid = int(line.split()[4])
				tag = sasa_list[np.where(sasa_list[:,0] == resid)[0][0],1]
				tag = round_special(tag,5)
				linelist=list(line)
				linelist[56:59] = list(pad(tag,5))
				newline = "".join(linelist)
				with open(fileout,"a") as w:
					w.write(newline)
			else:
				with open(fileout,"a") as w:
					w.write(line)
		
def get_coordinates(pdb_file,surf_atoms,surf_resids):
	'''returns the coordinates of the residues, calculated as the geomtric center of atoms in the residue'''

	#get coordinates of surface residues
	res_coords = []
	with open(pdb_file,"r") as f:
		res_1 = -1
		for line in f:
			if line.split()[0] == "ATOM":
				res_2 = int(line.split()[4])
				atom_num = int(line.split()[1])
				if atom_num in surf_atoms and res_2 in surf_resids:
					coord = [float(line.split()[5]), float(line.split()[6]), float(line.split()[7])]
					#test if first residue
					if res_2 != res_1 and res_1 == -1:
						res_coordi = [coord]
						res_1 = res_2
					#test if changing residues
					elif res_2 != res_1 and res_1 != -1:
						res_coords.append(np.append(np.array([res_1]),np.mean(res_coordi,axis=0)))
						res_coordi = [coord]
						res_1 = res_2
					#test if continuing same residue
					elif res_2 == res_1:
						res_coordi.append(coord)
		#append last residue
		res_coords.append(np.append(np.array([res_1]),np.mean(res_coordi,axis=0)))
	res_coords = np.array(res_coords)
	
	return res_coords
	
def get_surf_atom_coordinates(pdb_file,surf_atoms):
	'''coordinates of atoms surf_atoms'''

	#allowed atoms
	#allowed_atoms = ["C","N","S","O"]

	#get coordinates of surface residues
	surf_atom_coords = []
	with open(pdb_file,"r") as f:
		for line in f:
			if line.split()[0] == "ATOM": #and line.split()[2] != "CA" and line.split()[2] != "N": #try to avoid backbone atoms
				res = int(line.split()[4])
				atom_num = int(line.split()[1])
				if atom_num in surf_atoms:
					coord = [float(line.split()[5]), float(line.split()[6]), float(line.split()[7])]
					surf_atom_coords.append(coord)
	
	return surf_atom_coords
	
def get_max_face_coordinates(pdb_file,surf_atoms,max_face_residues):
	'''coordinates of residues in ice binding face'''

	#get coordinates of surface residues
	max_face_coords = []
	with open(pdb_file,"r") as f:
		for line in f:
			if line.split()[0] == "ATOM":
				res = int(line.split()[4])
				atom_num = int(line.split()[1])
				if atom_num in surf_atoms and res in max_face_residues:
					coord = [float(line.split()[5]), float(line.split()[6]), float(line.split()[7])]
					max_face_coords.append(coord)
	
	return max_face_coords
	
def get_projected_area(max_face_residues):
	'''area of ice binding face '''

	#get solvent exposed atom coordinates of max plane
	max_face_coordinates = get_max_face_coordinates(pdb_file,surf_atoms,max_face_residues)
	max_face_coordinates = np.array(max_face_coordinates)

	pca = PCA(n_components=2)
	projected_max_face_coords = pca.fit_transform(max_face_coordinates)
	#plt.scatter(projected_max_face_coords[:,0],projected_max_face_coords[:,1])
	#plt.show()

	points = projected_max_face_coords
	points = points - np.min(points) #make sure all positive for area calculation
	hull = ConvexHull(points)
	#plt.plot(points[:,0], points[:,1], 'o')
	#for simplex in hull.simplices:
	#    plt.plot(points[simplex, 0], points[simplex, 1], 'k-')
	#plt.gca().set_aspect('equal', adjustable='box')
	#plt.show()

	#cannot use hull.area - gives weird answer, use numpy polygon method
	#method from https://stackoverflow.com/questions/24467972/calculate-area-of-polygon-given-x-y-coordinates
	#also called Shoelace formula 
	x = []
	y=[]
	for i in hull.vertices:
		x.append(points[i,0])
		y.append(points[i,1])
	x = np.array(x)
	y = np.array(y)
	area = (0.5*np.abs(np.dot(x,np.roll(y,1))-np.dot(y,np.roll(x,1))))/100 #divide by 100 to go to nm^2 like sasa measurement
	return area
	
	
	
def test_face(i):
	'''test face i for geomtric criteria'''

	points = np.array([res_coords[perm_list[i,0],1:4],res_coords[perm_list[i,1],1:4],res_coords[perm_list[i,2],1:4]])
	#test points are close to each other
	test_dists = [distance(points[0],points[1],100),distance(points[0],points[2],100),distance(points[1],points[2],100)]
	#print test_dists
	if np.max(test_dists) > k_neighbor:
		return
	#fit plane to set of points
	plane = best_fitting_plane(points, True)
	
	#find how many points close to the plane and which side they are on
	points_in_plane = 0
	neg_side = 0
	pos_side = 0
	face_clusteri = []
	distances = []
	for a in range(0,np.shape(res_coords)[0]):
		dista = point_plane_dist(res_coords[a,1:4],plane[0:4])
		distances.append(abs(dista))
		if abs(dista) < k_plane_dist:
			points_in_plane = points_in_plane + 1
			face_clusteri.append(res_coords[a,0])
		if abs(dista) > k_split_dist and dista < 0:
			neg_side = neg_side + 1
		if abs(dista) > k_split_dist and dista > 0:
			pos_side = pos_side + 1
		if np.min([neg_side,pos_side]) > k_split:
			break
	
	#make sure not selecting internal face and enough residues in face
	if np.min([neg_side,pos_side]) > k_split: 
		return
		
	#used to make cleaned face with outliers removed but allowed for weird planes that were not really planes
	#just delete plane if not all resiudes have right number of neighbors
	#face_clusteri_clean = []
	for a in face_clusteri:
		neighbor_dists = []
		for b in face_clusteri:
			if a != b:
				coord_a = res_coords[np.where(res_coords[:,0] == a)[0][0],1:4]
				coord_b = res_coords[np.where(res_coords[:,0] == b)[0][0],1:4]
				dist_ab = distance(coord_a,coord_b,100)
				neighbor_dists.append(dist_ab)
		
		#ensure k_min_neighbors within k_neighbor_near
		num_neighbors = np.sum( np.array(neighbor_dists) < k_neighbor_near )
		if num_neighbors < k_min_neighbors:
			#face_clusteri_clean.append(a)
			return

	#ensure enough residues still left
	if len(face_clusteri) >= k_in_plane:
		height = np.mean(np.sort(np.array(distances))[-5:-1])
		return [face_clusteri,height]
	else:
		return

def find_faces(res_coords):
	'''semi-random search for plane; kind of expensive since searches all combinations of 3 points in set
	fit plane to all possible combinations of 3 points 
	and test whether plane captures more points, with almost all other points on one side'''

	#find all possible permutations
	global perm_list  #avoid passing to parallel process by defining as global
	perm_list = []
	for i in range(0,np.shape(res_coords)[0]):
		for j in range(0,np.shape(res_coords)[0]):
			for k in range(0,np.shape(res_coords)[0]):
				perm_list.append((i,j,k))
	perm_list = np.unique(np.sort(np.array(perm_list),axis=1),axis=0) #delete duplicates irrespective of ordering ( [1,2,3] = [2,1,3] )
	perm_list = np.array([v for v in perm_list if len(set(v)) == len(v)]) #delete doubles and triples ( [ 2,2,1] is not allowed )
	print("Number of permutations to search for plane matching: " + str(np.shape(perm_list)[0]))

	print("Searching plane permutations...")
	
	#parallel execution
	face_search = Parallel(n_jobs=processors)(delayed(test_face)(i) for i in range(0,np.shape(perm_list)[0]))
	
	face_search = np.array(face_search)
	face_search = face_search[face_search != np.array(None)]
	face_clusters = [ x[0] for x in face_search ]
	face_clusters = np.array(face_clusters)
	heights = [ x[1] for x in face_search ]
	
	print("Number of viable faces detected geometrically: " + str(len(face_clusters)))
	
	return [face_clusters,heights]

def get_hbl_times(n_residues):
	'''get hydrogen bond lifetime decay time from HBL files calculated using mdanalysis'''
	folder=hbl_folder
	prefix="/hbl_"
	decayTimeList=[]
	for i in range(0,n_residues):
	
		filename=folder + prefix + str(i) + ".txt"
		hbl_data = np.loadtxt(filename)
		x = hbl_data[:,0]
		y = hbl_data[:,hbl_function]
		
		check=0
		for t in range(0,len(x)):
			if y[t] <= hbl_cutoff:
				decayTimei = x[t]*10
				check = 1
				break
		if check == 0:
			decayTimei = 0
			print("WARNING: HBL for residue " + str(i) + " never decayed below cutoff. Setting decay time to 0")
		decayTimeList.append([i,decayTimei])		

	decayTimeList = np.array(decayTimeList)
	np.savetxt("decayTimeList_b"+str(begin)+"e"+str(end)+"ns.txt",decayTimeList)
	return decayTimeList

def edit_pdb_hbl(decayTimeList,pdb_file):
	'''edit pdb file with occupancy column being the decay time'''
	
	fileout = protein + "_hbl_b"+str(begin)+"e"+str(end)+"ns.pdb"
	os.system("rm " + fileout)

	with open(pdb_file, "r") as f:
		for line in f:
			if line.split()[0] == "ATOM":
				resid = int(line.split()[4])
				tag = decayTimeList[np.where(decayTimeList[:,0] == resid)[0][0],1]
				tag = round_special(tag,5)
				linelist=list(line)
				linelist[56:59] = list(pad(tag,5))
				newline = "".join(linelist)
				with open(fileout,"a") as w:
					w.write(newline)
			else:
				with open(fileout,"a") as w:
					w.write(line)
					
def score_face(i):
	'''get scoring for face i'''
	
	#get binding hbl
	binding_hbl_listi = []
	for j in range(0,len(face_clusters[i])):
		binding_hbl_value = decayTimeList[np.where(decayTimeList[:,0] == face_clusters[i][j])[0][0],1]
		binding_hbl_listi.append(binding_hbl_value)		
	binding_hbl = np.mean(binding_hbl_listi)
		
	#get binding SASA
	binding_sasa_listi = []
	for j in range(0,len(face_clusters[i])):
		binding_sasa_value = sasa_resid_values[np.where(sasa_resid_values[:,0] == face_clusters[i][j])[0][0],1]
		binding_sasa_listi.append(binding_sasa_value)	
	binding_sasa = np.sum(binding_sasa_listi)

	weighted_binding_hbl = np.sum(np.multiply(binding_hbl_listi,binding_sasa_listi))/binding_sasa
	
	#get non-binding hbl
	nonbinding_hbl_listi = []
	for j in range(0,len(surf_resids)):
		if surf_resids[j] not in face_clusters[i]:
			hbl_value = decayTimeList[np.where(decayTimeList[:,0] == surf_resids[j])[0][0],1]
			nonbinding_hbl_listi.append(hbl_value)
	nonbinding_hbl = np.mean(nonbinding_hbl_listi)
	
	#get non-binding SASA
	nonbinding_sasa_listi = []
	for j in range(0,len(surf_resids)):
		if surf_resids[j] not in face_clusters[i]:
			sasa_value = sasa_resid_values[np.where(sasa_resid_values[:,0] == surf_resids[j])[0][0],1]
			nonbinding_sasa_listi.append(sasa_value)
	nonbinding_sasa = np.sum(nonbinding_sasa_listi)
	
	weighted_nonbinding_hbl = np.sum(np.multiply(nonbinding_hbl_listi,nonbinding_sasa_listi))/nonbinding_sasa
	
	#get projected area
	face_area = get_projected_area(face_clusters[i])
	
	#score = c1 + c2*binding_sasa + c3*(weighted_binding_hbl - weighted_nonbinding_hbl) + c4*heights[i]
	return [weighted_binding_hbl,weighted_nonbinding_hbl,binding_sasa,nonbinding_sasa,face_area]
	
	
def get_hbl_faces(decayTimeList,face_clusters,heights,surf_resids,sasa_resid_values):
	'''get hbl for each face and the non-binding hbl that excludes that face
	choose the plane with the best score as the max plane and calculate final score'''
	
	result_list = Parallel(n_jobs=processors)(delayed(score_face)(i) for i in range(0,np.shape(face_clusters)[0]))
	
	#normalize
	result_list = np.array(result_list)
	normalized_area = result_list[:,4]/np.abs(np.max(result_list[:,4]))
	hbl_delta = result_list[:,0] - result_list[:,1]
	normalized_hbl_delta = hbl_delta/np.abs(np.max(hbl_delta))
	scores = normalized_area + normalized_hbl_delta
	max_plane = np.argmax(scores)

	height = heights[max_plane]
	final_binding_hbl = result_list[max_plane,0]
	final_nonbinding_hbl = result_list[max_plane,1]
	final_binding_sasa = result_list[max_plane,2]
	final_nonbinding_sasa = result_list[max_plane,3]
	final_face_area = result_list[max_plane,4]
	max_score = scores[max_plane]
	
	return [max_plane,max_score,final_binding_hbl,final_nonbinding_hbl,height,final_binding_sasa,final_nonbinding_sasa,final_face_area]
	
	
def edit_pdb_hbl_faces(max_face_residues,pdb_file):
	'''edit pdb file with occupancy column being the average HBL of the face a residue was assigend to'''
	fileout = protein + "_hbl_faces_b"+str(begin)+"e"+str(end)+"ns.pdb"
	os.system("rm " + fileout)

	with open(pdb_file, "r") as f:
		for line in f:
			if line.split()[0] == "ATOM":
				resid = int("".join(line[23:26]))
				if resid in max_face_residues:
					tag = round_special(1,5)
				else:
					tag = round_special(-1,5)
				linelist=list(line)
				linelist[56:59] = list(pad(tag,5))
				newline = "".join(linelist)
				with open(fileout,"a") as w:
					w.write(newline)
			else:
				with open(fileout,"a") as w:
					w.write(line)

def get_res_mutation(max_face_residues,res_coords,decayTimeList):
	'''find residues in binding face and NEAR binding face for possible mutation, identify residue with lowest HBL
	returns list sorted by decayTime'''
	mut_res_list = list(max_face_residues[:]) #stupid python

#	Allow mutations near binding face
	for a in max_face_residues:
		for b in res_coords[:,0]:
			if a != b and b not in mut_res_list:
				coord_a = res_coords[np.where(res_coords[:,0] == a)[0][0],1:4]
				coord_b = res_coords[np.where(res_coords[:,0] == b)[0][0],1:4]
				dist_ab = distance(coord_a,coord_b,100)
				if dist_ab <= k_neighbor_near:
					mut_res_list.append(b)
	mut_hbl_values = []
	for r in mut_res_list:
		binding_hbl_value = decayTimeList[np.where(decayTimeList[:,0] == r)[0][0],1]
		mut_hbl_values.append([r,binding_hbl_value])
	mut_hbl_values = np.array(mut_hbl_values)
	sorted_mut_res = mut_hbl_values[mut_hbl_values[:,1].argsort()]
	return sorted_mut_res[:,0]


#######################################################################################################
#Run code

[surf_resids,sasa_resid_values] = get_surface_resids(sasa_r_file)
surf_atoms = get_surface_atoms(sasa_a_file)
edit_pdb_sasa(sasa_a_file,surf_atoms,pdb_file)
edit_pdb_sasa_r(sasa_r_file,pdb_file)
res_coords = get_coordinates(pdb_file,surf_atoms,surf_resids)
n_residues = int(np.max(res_coords[:,0]))+1
print("Number of residues found: " + str(n_residues))

decayTimeList = get_hbl_times(n_residues)
edit_pdb_hbl(decayTimeList,pdb_file)

start_time = time.time()
[face_clusters,heights] = find_faces(res_coords)
end_time = time.time()
print("Face search took: " + str(round(end_time - start_time,3)) + " seconds")

start_time = time.time()
[max_plane,max_score,final_binding_hbl,final_nonbinding_hbl,height,final_binding_sasa,final_nonbinding_sasa,final_face_area] = get_hbl_faces(decayTimeList,face_clusters,heights,surf_resids,sasa_resid_values)
max_face_residues = face_clusters[max_plane]
end_time = time.time()
print("Face scoring took: " + str(round(end_time - start_time,3)) + " seconds")
	
#fitted_score = c1 + c2*final_face_area + c3*final_binding_hbl + c4*final_nonbinding_hbl
#exec_string = "./evaluate_nn.wl "+str(final_face_area)+" "+str(final_binding_hbl)+" "+str(final_nonbinding_hbl)
#print exec_string
#fitted_score = float(subprocess.check_output(exec_string,cwd="../nn_files/",shell=True))

#interpolate fitted score from nerual network results - anything less than 0 or Nan is convertec to zero
nn_file = "nn_output_ext.npy"
print "Neural net grid from file: "+nn_file
nn_data = np.load(nn_file)
point = [final_face_area,final_binding_hbl,final_nonbinding_hbl]
fitted_score = np.max([0,griddata(nn_data[:,0:3],nn_data[:,3],point)])
if np.isnan(fitted_score):
	fitted_score = 0
	print("Interpolation returned NaN")

print("Residues (index) in binding face: " + ', '.join([str(int(x)) for x in max_face_residues]))
print("Number of residues in binding face: " + str(len(max_face_residues)))
print("Binding HBL: " + str(final_binding_hbl))
print("Nonbinding HBL: " + str(final_nonbinding_hbl))
print("Binding SASA: " + str(final_binding_sasa))
print("Nonbinding SASA: " + str(final_nonbinding_sasa))
print("Projected area: " + str(final_face_area))
print("Height: " + str(height))
print("Score: " + str(max_score))
print("Fitted Score: " + str(fitted_score))

edit_pdb_hbl_faces(max_face_residues,pdb_file)

sorted_mut_res = get_res_mutation(max_face_residues,res_coords,decayTimeList)
print("Sorted residues for mutation: " + ', '.join([str(int(x)) for x in sorted_mut_res]))















