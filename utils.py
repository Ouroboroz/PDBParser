import numpy as np
import matplotlib.pyplot as plt


def get_chain_bond_lengths(chain):
	'''
	Returns bond lengths of a chain as a list

	Parameters:
		chain : Chain :: Chain structure of residues
	'''
	residues = chain.get_residues()
	bond_lengths = []
	for residue_index in range(len(residues)-1):
		current_residue_coordinate = residues[residue_index].get_coordinate()
		next_residue_coordinate = residues[residue_index+1].get_coordinate()
		bond_length = np.linalg.norm(current_residue_coordinate - next_residue_coordinate)
		bond_lengths.append(bond_length)
	return bond_lengths

def get_chain_bond_angles(chain):
	'''
	Returns residue angles of a chain as a list

	Parameters:
		chain : Chain :: Chain structure of residues
	'''
	residues = chain.get_residues()
	bond_angles = []
	for residue_index in range(len(residues)-2):
		residue_1_coordinate = residues[residue_index].get_coordinate()
		residue_2_coordinate = residues[residue_index+1].get_coordinate()
		residue_3_coordinate = residues[residue_index+2].get_coordinate()
		bond_angle = calculate_vector_angle(residue_1_coordinate, residue_2_coordinate, 
												residue_3_coordinate)
		bond_angles.append(bond_angle)
	return bond_angles

def get_chain_torsion_angles(chain):
	'''
	Returns torsion angles of a chain as a list

	Parameters:
		chain : Chain :: Chain structure of residues
	'''
	residues = chain.get_residues()
	torsion_angles = []
	for residue_index in range(len(residues)-3):
		residue_1_coordinate = residues[residue_index].get_coordinate()
		residue_2_coordinate = residues[residue_index+1].get_coordinate()
		residue_3_coordinate = residues[residue_index+2].get_coordinate()
		residue_4_coordinate = residues[residue_index+3].get_coordinate()
		torsion_angle = calculate_torsion_angle(residue_1_coordinate, residue_2_coordinate,
												residue_3_coordinate, residue_4_coordinate)
		torsion_angles.append(torsion_angle)
	return torsion_angles

def plot_bond_angles_distribution(chain, hist_bins=None, hist_range=None, hist_density=False):
	'''
	Plots residue angles of a chain as a histogram

	Parameters:
		chain : Chain :: Chains structure of residues
		hist_bins : int or sequence or str :: Same as plt bins
		hist_range : tuple or None, default: None :: Same as plt range
		hist_density: bool, default: False :: Same as plt density
	'''
	bond_angles = get_chain_bond_angles(chain)
	plt.hist(bond_angles, bins=hist_bins, range=hist_range, density=hist_density)

def plot_bond_length_distribution(chain, hist_bins=None, hist_range=None, hist_density=False):
	'''
	Plots bond lengths of a chain as a histogram

	Parameters:
		chain : Chain :: Chains structure of residues
		hist_bins : int or sequence or str :: Same as plt bins
		hist_range : tuple or None, default: None :: Same as plt range
		hist_density: bool, default: False :: Same as plt density
	'''
	bond_lengths = get_chain_bond_lengths(chain)
	plt.hist(bond_lengths, bins=hist_bins, range=hist_range, density=hist_density)

def plot_torsion_angles_distribution(chain, hist_bins=None, hist_range=None, hist_density=False):
	'''
	Plots torsion angles of a chain as a histogram

	Parameters:
		chain : Chain :: Chains structure of residues
		hist_bins : int or sequence or str :: Same as plt bins
		hist_range : tuple or None, default: None :: Same as plt range
		hist_density: bool, default: False :: Same as plt density
	'''
	torsion_angles = get_chain_torsion_angles(chain)
	plt.hist(torsion_angles, bins=hist_bins, range=hist_range, density=hist_density)

def calculate_vector_angle(point_1, point_2, point_3):
	'''
	Calculates angle between vectors
	(point_1, point_2) and (point_2, point_3)
	Returns angle in radians

	Parameters:
		point_1,2,3 : numpy array :: Points 1,2,3 to define two vectors
	'''
	vector_1 = point_1 - point_2
	vector_2 = point_3 - point_2
	cosine_angle = np.dot(vector_1, vector_2) / (np.linalg.norm(vector_1) * np.linalg.norm(vector_2))
	angle = np.arccos(cosine_angle)

	return angle

def normalize(vector):
	'''
	Normalizes the vector
	Returns unit vector 

	Parameters:
		vector : np array :: vector to normalize
	'''
	return vector/np.linalg.norm(vector)

def calculate_torsion_angle(point_1, point_2, point_3, point_4):
	'''
	Calculates torsion angle between 4 points

	Parameters:
		point_1,2,3,4 : np array :: 4 points that form a torsion angle
	'''
	b1 = point_2 - point_1
	b2 = point_3 - point_2
	b3 = point_4 - point_3
	n1 = normalize(np.cross(b1,b2))
	n2 = normalize(np.cross(b2,b3))
	m1 = np.cross(n1,normalize(b2))
	x = np.dot(n1,n2)
	y = np.dot(m1,n2)
	return np.arctan2(y,x)

