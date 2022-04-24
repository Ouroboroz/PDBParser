import pickle
from pathlib import Path
import numpy as np
from structures import Complex, Chain, Residue, Atom

class PDBParser:
	'''
	Parser for .pdb files into a Complex structure
	'''
	def __init__(self,verbose=False,local=True):
		'''
		Initializes PDBParser object
		Call parse_pdb to parse a certain file

		Parameters:
			verbose : bool, default - False :: Debug output toggle during parsing
			local : bool, default - True :: Calculate local system for residues

		'''
		self.verbose = verbose
		self.local = local

	def parse_pdb(self, path_to_pdb):
		'''
		Opens .pdb file and
		parses lines into Complex structure

		Parameters:
			path_to_pdb : string :: path to a .pdb file
		'''
		path = Path(path_to_pdb)
		assert path.suffix == '.pdb', 'PDB file needs to end in .pdb'
		
		complex_file_name = path.stem
		#complex file name should be in format (complex id)_(chainID)_(chainID)
		if len(complex_file_name.split('_')) != 3:
			raise FileNameError('File name must be in the format (complex id)_(chainID)_(chainID)')
		complex_name = complex_file_name.split('_')[0]
		complex_ligand_id = complex_file_name.split('_')[1]
		complex_protein_id = complex_file_name.split('_')[2]
		complex_structure = Complex(complex_name, proteinID=complex_protein_id,
									ligandID=complex_ligand_id)
		
		current_chainID = None
		current_resSeq = None
		current_residue_atom_list = []
		residue = Residue('init','residue')

		with path.open(encoding="utf-8") as file:
			for line in file:
				line = line.rstrip()
				parsed_pdb_row = parse_pdb_row(line)

				if current_chainID != parsed_pdb_row['chainID']:
					current_chainID = parsed_pdb_row['chainID']
					chain_structure = Chain(parsed_pdb_row['chainID'])
					complex_structure.add_chain(chain_structure)

				if current_resSeq != parsed_pdb_row['resSeq']:
					current_resSeq = parsed_pdb_row['resSeq']
					residue.add_atoms(current_residue_atom_list)
					current_residue_atom_list = []
					residue = Residue(parsed_pdb_row['name'],
										int(parsed_pdb_row['resSeq']),self.local)
					chain_structure.add_residue(residue)

				coordinate = np.array([
										float(parsed_pdb_row['x']),
										float(parsed_pdb_row['y']),
										float(parsed_pdb_row['z'])
									 ])
				
				atom = Atom(parsed_pdb_row['name'],parsed_pdb_row['element'], coordinate,
									parsed_pdb_row['altLoc'] == 'CA')
				current_residue_atom_list.append(atom)
		return complex_structure

	def parse_pdb_to_pickle(self,path_to_pdb, pickle_path=None):
		'''
		Parses .pdb into a pickle that stores Complex structure

		Parameters:
			path_to_pdb : string :: Path to .pdb file
			pickle_path : string, default - None :: Optional path to store pickle
		'''
		if pickle_path is None:
			path = Path(path_to_pdb)
			complex_file_name = path.stem
			pickle_path = './processedComplex/'+complex_file_name+'.pkl'
		complex_structure = self.parse_pdb(path_to_pdb)
		path = Path(pickle_path)
		with path.open() as f:
			pickle.dump(complex_structure, f)

def parse_pdb_from_pickle(pickle_path):
	'''
	Gets Complex structure from pickle file

	Parameters:
		pickle_path : string :: Path to pickle file
	'''
	path = Path(pickle_path)
	with path.open() as f:
		complex_structure = pickle.dump(f)
	return complex_structure


def parse_pdb_row(pdb_row):
	'''
	Parses a row of a .pdb file
	Specifically ATOM/HETATOM row

	Parameters:
		pdb_row : string :: Row of .pdb file
	'''
	indices = {'ATOM':(0,3),
			'HETATOM':(0,5),
			'serial':(6,10),
			'name':(12,15),
			'altLoc':16,
			'resName':(17,19),
			'chainID':21,
			'resSeq':(22,25),
			'x':(30,37),
			'y':(38,45),
			'z':(46,53),
			'element':(76,77)}
	row_length = len(pdb_row)
	if len(pdb_row) < 78:
		raise PDBError(f"Row has {row_length} characters but needs atleast 78")
	parsed_pdb_row = {}
	for item, indice_tup in indices.items():
		parsed_pdb_row[item] = pdb_row[indice_tup[0]:indice_tup[1]+1]
	return parsed_pdb_row

class FileNameError(Exception):
	''' 
	Exception thrown when filename incorrect
	'''

class PDBError(Exception):
	'''
	Exception thrown when there is a .pdb file related error
	'''
