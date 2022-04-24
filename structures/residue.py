import math
import numpy as np

class Residue:
	'''
	Structure that holds the representation of a Residue
	'''
	def __init__(self,name, resSeq, local=True):
		'''
		Initializes a Residue structure
		A residue is defined by the name, chain ID, sequence number

		Parameters:
			name : str :: name of the residue
			chainID : str :: chain identifier
			resSeq : int :: residue sequence number
			local : bool, default - True :: boolean if local system needs to be calculated
		'''
		self.name = name
		self.resSeq = resSeq
		self.alpha_carbon = None
		self.coordinate = None
		# carboxyl group (carbon, oxygen)
		self.carboxyl = [None,None]
		# amino group (nitrogen, hydrogen)
		self.amino = [None, None]
		# (u, t, v)
		self.u = None
		self.t = None
		self.v = None
		self.n = None
		self.atoms = []

		self.local = local

	def add_alpha_carbon(self, carbon=None):
		'''
		Sets the alpha carbon of the residue
		Alpha carbon acts as the coordinate of the residue

		Parameters:
			carbon : Atom, default - None :: Carbon atom that represents the alpha carbon
		'''
		if carbon.name != 'CA':
			raise ResidueError('Alpha Carbon must have name CA')
		self.alpha_carbon = carbon

	def add_carboxyl(self, carbon=None, oxygen=None):
		'''
		Sets the carboxyl group of the residue

		Parameters:
			carbon : Atom, default - None  :: Carbon atom of the carboxyl
			oxygen : Atom, default - None  :: Oxygen atom of the carboxyl
		'''
		if carbon.name != 'C':
			raise ResidueError('Carboxyl Carbon must have name C')
		if oxygen.name != 'O':
			raise ResidueError('Carboxyl Oxygen must have name O')
		if carbon is not None:
			self.carboxyl[0] = carbon
		if oxygen is not None:
			self.carboxyl[1] = oxygen

	def add_amino(self, nitrogen=None, hydrogen=None):
		'''
		Sets the amino group of the residue

		Parameters:
			nitrogen : Atom, default - None  :: Nitrogen atom of the amino
			oxygen : Atom, default - None  :: Oxygen atom of the amino
		'''
		if nitrogen is not None:
			self.amino[0] = nitrogen
		if hydrogen is not None:
			self.amino[1] = hydrogen

	def check_valid_residue(self):
		'''
		Asserts that the residue has a valid configuration
		'''
		if self.alpha_carbon is None:
			raise ResidueError('Missing Alpha Carbon')
		if self.carboxyl[0] is None:
			raise ResidueError('Missing Carboxyl Carbon')
		if self.carboxyl[1] is None:
			raise ResidueError('Missing Carboxyl Oxygen')
		if self.amino[0] is None:
			raise ResidueError('Missing Amino Nitrogen')
		if self.amino[1] is None:
			raise ResidueError('Missing Amino Oxygen')

	def calculate_system(self):
		'''
		Calculates the local structure system of the residue
		Need to call this after all atoms are added
		'''
		self.check_valid_residue()
		self.u = self.carboxyl[0].set_vec(self.alpha_carbon)
		self.t = self.amino[0].set_vec(self.alpha_carbon)
		u_cross_t = np.cross(self,self.u,self.t)
		self.n = u_cross_t/np.linalg.norm(u_cross_t)
		self.v = np.cross(self.n,self.u)

	def get_coordinate(self):
		'''
		Returns the coordinate of the residue
		Equivalent to the coordinates of the alpha carbon
		'''
		if self.alpha_carbon is None:
			raise ResidueError('Missing Alpha Carbon')
		return self.coordinate

	def get_system(self):
		'''
		Returns the local structure system (n,u,v)
		'''
		return self.n, self.u, self.v 

	def get_segid(self):
		'''
		Returns the sequence number the residue is
		'''
		return self.resSeq

	def get_group(self, name):
		'''
		Returns the group with the name parameter
		Will be a single atom if alpha_carbon, else a tuple
		returns None name is not valid

		Parameters:
			name : string :: Name of the group 
		'''
		if name == 'alpha':
			return self.alpha_carbon
		if name == 'carboxyl':
			return self.carboxyl
		if name == 'amino':
			return self.amino
		return None

	def get_atoms(self):
		'''
		Generator that returns all of the atoms of the residue
		In the order of alpha carbon, carboxyl atoms, amino atoms et al
		'''
		self.check_valid_residue()
		yield self.alpha_carbon
		for atom in self.carboxyl:
			yield atom
		for atom in self.amino:
			yield atom
		for atom in self.atoms:
			yield atom

	def add_atoms(self, atoms):
		'''
		Adds list of atoms to the residue
		Will automatically calculate local system if needed
		Returns number of atoms that are not C, O, N, H and CA
		'''
		unknown_atom_count = 0
		for atom in atoms:
			if atom.name == 'C':
				self.add_carboxyl(carbon=atom)
				continue
			if atom.name == 'O':
				self.add_carboxyl(oxygen=atom)
				continue
			if atom.name == 'N':
				self.add_amino(nitrogen=atom)
				continue
			if atom.name == 'H':
				self.add_amino(hydrogen=atom)
				continue
			if atom.name == 'CA':
				self.add_alpha_carbon(carbon=atom)
			self.atoms.append(atom)
			unknown_atom_count += 1
		if self.local:
			self.calculate_system()
		return unknown_atom_count == 0

	def get_relative_position_edge_features(self, coordinate):
		'''
		Calculates the relative position edge features of the residue

		Parameters:
			coordinate : numpy array :: coordinate of an residue
		'''
		system_tup = (self.n.T, self.u.T, self.v.t)
		return np.vstack(system_tup)@(coordinate - self.coordinate)

	def get_relative_orientation_edge_features(self,system):
		'''
		Calculates the relative orientation edge features of the residue

		Parameters:
			system : tuple :: tuple of system vectors
		'''
		n = system[0]
		u = system[1]
		v = system[2]
		system_tup = (n,u,v)
		return np.vstack(system_tup)@self.n, np.vstack(system_tup)@self.u, np.vstack(system_tup)@self.v

	def get_distance_edge_features(self, coordinate, sigmas):
		'''
		Calculates the distance edge features of the residue

		Parameters:
			coordinate : numpy array :: coordinate of an residue
			sigmas : list like :: list of sigmas for calculations
		'''
		f = []
		for sigma in sigmas:
			f.append(math.exp(np.linalg.norm(coordinate - self.coordinate)**2/(2*sigma**2)))
		return f

class ResidueError(Exception):
	'''
	Exception for all Residue errors
	'''
