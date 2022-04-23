class Chain:
	'''
	Structure that holds the representation of a chain in a complex
	'''
	def __init__(self, chainID):
		'''
		Initialzes a Chain structure

		Parameters:
			chainID : str :: ID of the chain
		'''
		self. chainID = chainID
		self.residueSequence = {}

	def add_residue(self, residue):
		'''
		Adds a Residue structure to the chain

		Parameters:
			residue: Residue :: Residue to add to the chain
		'''
		if residue.resSeq in self.residueSequence:
			raise ChainError('Residue already in the chain')
		self.residueSequence[residue.resSeq] = residue
		return True

	def get_residue_dict(self):
		'''
		Get dictionary representation of the residue
		resSeq : Residue
		'''
		return self.residueSequence
	
	def get_residue_list(self):
		'''
		Get list representation of the residue
		Order by resSeq
		'''
		return [self.residueSequence[i] for i in sorted(self.residueSequence.keys())]

	def get_residue_sequence(self):
		'''
		Get the sequence of residue numbers ordered
		'''
		return sorted(self.residueSequence.keys())
	
	def get_residues(self):
		'''
		Alternative name for get residue list
		'''
		return self.get_residue_list()

class ChainError(Exception):
	'''
	Exception for Chain errors
	'''
