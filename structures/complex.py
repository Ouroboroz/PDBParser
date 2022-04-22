class Complex:
	def __init__(self, name, proteinID, ligandID):
		self.name = name
		self.proteinID = proteinID
		self.ligandID = ligandID
		self.protein = None
		self.ligand = None

	def add_chain(self, chain):
		if chain.chainID == self.proteinID:
			self.protein = chain
			return True
		elif chain.chainID == self.ligandID:
			self.ligand = chain
			return True
		return False