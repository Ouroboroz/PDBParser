class Atom:
	def __init__(self,name,element,coordinate,alpha = False):
		'''
		Abstraction of an atom in the protein's residue

		Parameters:
			name : str :: Name of the atom
			element : str :: Element of the atom
			coordinate : np array of dim 3 :: Coordinate of the atom in R3
			alpha : bool ::  boolean if atom is the alpha carbon

		'''
		self.name = name
		self.element = element
		assert(len(coordinate) == 3)
		self.coordinate = coordinate
		self.alpha = alpha
		self.vec = None

	def get_coordinate(self):
		return self.coordinate

	def get_vec(self):
		return self.vec

	def set_vec(self, alpha_carbon):
		vec = self.coordinate - alpha_carbon.coordinate
		self.vec = vec/np.linalg.norm(vec)
		return self.vec