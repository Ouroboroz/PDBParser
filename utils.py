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

def get_chain_residue_angles(chain):
	'''
	Returns residue angles of a chain as a list

	Parameters:
		chain : Chain :: Chain structure of residues
	'''
	residues = chain.get_residues()
	residue_angles = []
	for residue_index in range(len(residues)-2):
		residue_1_coordinate = residues[residue_index].get_coordinate()
		residue_2_coordinate = residues[residue_index+1].get_coordinate()
		residue_3_coordinate = residues[residue_index+2].get_coordinate()
		residue_angle = calculate_vector_angle(residue_1_coordinate, residue_2_coordinate, 
												residue_3_coordinate)
		residue_angles.append(residue_angle)
	return residue_angles

def plot_residue_angles_distribution(chain, hist_bins=None, hist_range=None, hist_density=False):
	'''
	Plots residue angles of a chain as a histogram

	Parameters:
		chain : Chain :: Chains structure of residues
		hist_bins : int or sequence or str :: Same as plt bins
		hist_range : tuple or None, default: None :: Same as plt range
		hist_density: bool, default: False :: Same as plt density
	'''
	residue_angles = get_chain_residue_angles(chain)
	plt.hist(residue_angles, bins=hist_bins, range=hist_range, density=hist_density)
