import numpy as np

class Atom:
	'''
	Structures that holds the representation of an atom
	'''
	def __init__(self,name,element,coordinate,alpha = False):
		'''
		Initializes an Atom object

		Parameters:
			name : str :: Name of the atom
			element : str :: Element of the atom
			coordinate : np array of dim 3 :: Coordinate of the atom in R3
			alpha : bool, default - False ::  boolean if atom is the alpha carbon

		'''
		self.name = name
		self.element = element
		if len(coordinate) != 3:
			raise ValueError('Coordinate must be in 3D')
		self.coordinate = coordinate
		self.alpha = alpha
		self.vec = None

	def get_coordinate(self):
		'''
		Returns coordinate of the atom
		'''
		return self.coordinate

	def get_vec(self):
		'''
		Returns the vector of the atom to the alpha carbon
		'''
		return self.vec

	def set_vec(self, alpha_carbon):
		'''
		Sets the vector of the atom to the alpha carbon
		Returns the vector

		Parameters:
			alpha_carbon : Atom :: Alpha carbon atom
		'''
		vec = self.coordinate - alpha_carbon.coordinate
		self.vec = vec/np.linalg.norm(vec)
		return self.vec
