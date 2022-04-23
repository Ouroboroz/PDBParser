class Complex:
	'''
	Structure that holds the representation of a complex
	'''
	def __init__(self, name, proteinID, ligandID):
		'''
		Initializes a Complex structure
		A complex is defined by its name, protein and ligand id

		Parameters:
			name : str :: name of the residue
			proteinID : string :: ID of the protein chain
			ligandID : string :: ID of the ligand chain
		'''
		self.name = name
		self.proteinID = proteinID
		self.ligandID = ligandID
		self.protein = None
		self.ligand = None

	def add_chain(self, chain):
		'''
		Appends a chain to the complex
		Return True if chain matches ID, False else

		Parameters;
			chain : Chain :: Chain structure
		'''
		if chain.chainID == self.proteinID:
			self.protein = chain
			return True
		if chain.chainID == self.ligandID:
			self.ligand = chain
			return True
		return False
