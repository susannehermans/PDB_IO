class PDBLine:
	def __init__(self, line):
		self.line = line

	def __str__(self) -> str:
		return self.line
    
    #@classmethod
    def from_line(cls, line):
        return cls(line)
