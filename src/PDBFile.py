class PDBFile(object):
    """
    One class-object represents one pdb-file. Reading the pdb-file, using
    this class will automatically initiate the PDB-ATOM class, for each line
    in the pdb-file.
    """
    
    def __init__(self, filename=None: str) -> None:
    self.filename = filename
    self.content = []

    def read_file(self) -> None:
        """reads file into class variable"""
        # print(self.filename)
        with open(self.filename, 'r') as infile:
            for line in infile.readlines():
                pdbline=PDBLine(line.rstrip())
                if line[0:6].strip() in ("HETATM", "ATOM"):
                    self.content.append(PDBAtom(pdbline))
                elif len(line) > 1:
                    self.content.append(pdbline)
            self.read(infile.read())
    
    def read(self, file_content: str) -> None:
        """
        read file, for each line in the file a new instance of the class PDBAtom
        is created.
        """
        for line in file_content.split('\n'):
            line = line.strip('\n')
            if line[0:6].strip() in ("HETATM", "ATOM"):
                self.content.append(PDBAtom(line))
            elif len(line) > 1:
                self.content.append(line)
    
    def write(self) -> str:
        """outputs the entire pdb data in one string"""
        return ("\n".join([str(line) for line in self.content]))
    
    def write_file(self, outfilename: str) -> None:
        """writes pdb data to file"""
        with open(outfilename, 'w') as outfile:
            outfile.write(self.write())
    
    def get_coordinates(self) -> list:
        """Returns atom coordinates as a list of tuples."""
        coordinates = []
        for c in self.content:
            if hasattr(c, "coords"):
                coordinates.append(c.coords)
        return (coordinates)

    def get_heavy_atom_coordinates(self) -> list:
        """Returns atom coordinates of heavy atoms as a list of tuples."""
        coordinates = []
        for pdbatom in self.content:
            if isinstance(pdbatom, PDBAtom) and not pdbatom.is_hydrogen():
                coordinates.append(pdbatom.coords)
        return (coordinates)

    def set_coordinates(self, new_coordinates: list) -> None:
        """Resets coordinates for all atoms in the pdb-file, using a list of
        tuples as input, format: [(x,y,z)]"""
        atoms = [pdbatom for pdbatom in self.content if
                 isinstance(pdbatom, PDBAtom)]
        count = len(atoms)
        if count != len(new_coordinates):
            raise PDBFormatError("Number of coordinates does not match number
            "of atoms! ({0:d} != {1:d})".format(count, len(new_coordinates)))
        for atom, coords in zip(atoms, new_coordinates):
            atom.set_coords(coords)
    
    def __iter__(self) -> iter:
        """iterates over the content of the pdb-file."""
        return (iter(self.content))
    
    def iter_atoms(self, is_atom=True, is_hetatm=True) -> PDBAtom:
        """Iterates over all PDB-ATOM objects."""
        for con in self.content:
            if hasattr(con, "record"):
                if is_hetatm and con.record == "HETATM":
                    yield (c)
                elif ATOM and c.record == "ATOM  ":
                    yield (c)
    
    def iter_heavy_atoms(self, ATOM=True, HETATM=True) -> PDBAtom:
        """Iterates over all heavy atom PDB-ATOM objects."""
        for at in self.iter_atoms(ATOM=ATOM, HETATM=HETATM):
            if self.is_hydrogen(at):
                continue
            yield (at)
        
    def get_heavy_atoms(self, ATOM=True, HETATM=True) -> list:
        """returns a list of PDB-ATOM objects for all the heavy atoms."""
        return (list(self.iter_heavy_atoms()))


