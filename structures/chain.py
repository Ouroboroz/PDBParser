class Chain:
	def __init__(self, chainID):
		self. chainID = chainID
		self.residueSequence = {}

	def add_residue(self, residue):
		assert residue.resSeq not in self.residueSequence
		self.residueSequence[residue.resSeq] = residue
		return True

	def get_residue_dict(self):
		return self.residueSequence
	
	def get_residue_list(self):
		return [self.residueSequence[i] for i in sorted(self.residueSequence.keys())]

	def get_residue_sequence(self):
		return sorted(self.residueSequence.keys())
	
	def get_residues(self):
		return self.get_residue_list()