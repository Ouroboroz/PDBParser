class Residue:
	def __init__(self,name, resSeq):
		'''
		Abstraction of a residue in a protein
		A residue is defined by the name, chain ID, sequence number

		Parameters:
			name : str :: name of the residue
			chainID : str :: chain identifier
			resSeq : int :: residue sequence number
		'''
		self.name = name
		self.resSeq = resSeq
		self.alpha_carbon = alpha_carbon
		self.coordinate = alpha_carbon.coordinate
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

	def add_alpha_carbon(self, carbon=None):
		self.alpha_carbon = carbon

	def add_carboxyl(self, carbon=None, oxygen=None):
		if carbon is not None:
			self.carboxyl[0] = carbon
		if oxygen is not None:
			self.carboxyl[1] = oxygen

	def add_amino(self, nitrogen=None, hydrogen=None):
		if nitrogen is not None:
			self.amino[0] = nitrogen
		if hydrogen is not None:
			self.amino[1] = hydrogen

	def calculate_local_system(self):
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
		self.u = carbon.set_vec(self.alpha_carbon)
		self.t = nitrogen.set_vec(self.alpha_carbon)
		u_cross_t = np.cross(self,u,self.t)
		self.n = u_cross_t/np.linalg.norm(u_cross_t)
		self.v = np.cross(self.n,self.u)

	def get_coordinate(self):
		return self.coordinate

	def get_system(self):
		return self.n, self.u, self.v 

	def get_id(self):
		return self.chainID

	def get_segid(self):
		return self.resSeq

	def get_group(self, name):
		if name == 'alpha':
			return self.alpha_carbon
		if name == 'carboxyl':
			return self.carboxyl
		elif name == 'amino':
			return self.amino
		return None

	def get_atoms(self):
		yield self.alpha_carbon
		for atom in self.carboxyl:
			yield atom
		for atom in self.amino:
			yield atom
		for atom in self.atoms:
			yield atom

	def add_atoms(self, atoms):
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
		return unknown_atom_count == 0

	def get_relative_position_edge_features(self, coordinate):
		system_tup = (self.n.T, self.u.T, self.v.t)
		return np.vstack(system_tup)@(coordinate - self.coordinate)

	def get_relative_orientation_edge_features(self,system):
		n = system[0]
		u = system[1]
		v = system[2]
		system_tup = (n,u,v)
		return np.vstack(system_tup)@n, np.vstack(system_tup)@u, np.vstack(system_tup)@v

	def get_distance_edge_features(self, coordinate, sigmas):
		f = []
		for sigma in sigmas:
			f.append(math.exp(np.linalg.norm(coordinate - self.coordinate)**2/(2*sigma**2)))
		return f

class ResidueError(Exception):
	pass