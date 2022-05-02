import os
import pickle
from pathlib import Path
import argparse
import numpy as np
from structures import Complex, Chain, Residue, Atom

class PDBParser:
	'''
	Parser for .pdb files into a Complex structure
	'''
	def __init__(self,verbose=False, local=True):
		'''
		Initializes PDBParser object
		Call parse_pdb to parse a certain file

		Parameters:
			verbose : int, default - False :: Debug output toggle during parsing
						0 - No debug statements
						1 - Major debug statements
						2 - All debug statements
			local : bool, default - True :: Calculate local system for residue
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
		if path.suffix !='.pdb':
			raise FileNameError('PDB file needs to end in .pdb')
		
		complex_file_name = path.stem
		#complex file name should be in format (complex id)_(chainID)_(chainID)
		if len(complex_file_name.split('_')) != 3:
			raise FileNameError('File name must be in the format (complex id)_(chainID)_(chainID)')
		complex_name = complex_file_name.split('_')[0]
		complex_ligand_id = complex_file_name.split('_')[1]
		complex_protein_id = complex_file_name.split('_')[2]
		complex_structure = Complex(complex_name, proteinID=complex_protein_id,
									ligandID=complex_ligand_id)
		if self.verbose:
			print(f"Parsing {complex_name} with chains {complex_protein_id} and {complex_ligand_id}")

		current_chainID = None
		current_resSeq = None
		current_residue_atom_list = []
		residue = Residue('init','residue',local=False)

		with path.open(encoding="utf-8") as file:
			for line in file:
				line = line.rstrip()
				if self.verbose > 1:
					print('  Parsing', line)

				parsed_pdb_row = parse_pdb_row(line)

				if current_chainID != parsed_pdb_row['chainID']:
					if self.verbose > 2:
						print(f"    Starting Chain: {parsed_pdb_row['chainID']}")
					current_chainID = parsed_pdb_row['chainID']
					chain_structure = Chain(parsed_pdb_row['chainID'])
					complex_structure.add_chain(chain_structure)

				if current_resSeq != parsed_pdb_row['resSeq']:
					if self.verbose > 2:
						print(f"    Starting Residue: {parsed_pdb_row['resSeq']}")
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
				if self.verbose > 3:
					print(f"      Starting Atom: {parsed_pdb_row['name']} ({parsed_pdb_row['altLoc']}) \
								with element {parsed_pdb_row['element']} at coordinate {parsed_pdb_row}")
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
			os.makedirs('./processedComplex',exist_ok=True)
			pickle_path = './processedComplex/'+complex_file_name+'.pkl'
		if self.verbose:
			print(f'Parsing {path_to_pdb} to {pickle_path}')
		complex_structure = self.parse_pdb(path_to_pdb)
		path = Path(pickle_path)
		with path.open("wb") as f:
			pickle.dump(complex_structure, f)

def parse_pdb_from_pickle(pickle_path, verbose=False):
	'''
	Gets Complex structure from pickle file

	Parameters:
		pickle_path : string :: Path to pickle file
	'''
	path = Path(pickle_path)
	if verbose:
		print(f"Unpickling {pickle_path} into Complex")
	with path.open("r") as f:
		complex_structure = pickle.load(f)
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
			'altLoc':(16,16),
			'resName':(17,19),
			'chainID':(21,21),
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
		parsed_pdb_row[item] = pdb_row[indice_tup[0]:indice_tup[1]+1].strip()
	return parsed_pdb_row

class FileNameError(Exception):
	''' 
	Exception thrown when filename incorrect
	'''

class PDBError(Exception):
	'''
	Exception thrown when there is a .pdb file related error
	'''

if __name__ == "__main__":
	parser = argparse.ArgumentParser(prog='pdbparser',
                                    usage='%(prog)s [options] path',
                                    description='Process .pdb file to a pickle')
	parser.add_argument('Path',
                       metavar='path',
                       type=str,
                       help='the path to file')
	parser.add_argument('--verbose', '-v', action='count',
						help='set to have debug outputs')
	parser.add_argument('--local', '-l', action='store_true',
						help='set to calculate local system')
	parser.add_argument('--pickle', '-p',
						type=str,
						help='path to pickle file')

	args = parser.parse_args()
	pdb_file_path_arg = args.Path
	verbose_arg = args.verbose
	local_arg = args.local
	pickle_path_arg = args.pickle
	pdbparser = PDBParser(verbose=verbose_arg, local=local_arg)

	pdbparser.parse_pdb_to_pickle(path_to_pdb=pdb_file_path_arg, pickle_path=pickle_path_arg)




