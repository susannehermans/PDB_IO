from .PDBLine import PDBLine
from .PDBAtom import PDBAtom

class PDBFile(object):
    """
    One class-object represents one pdb-file. Reading the pdb-file, using
    this class will automatically initiate the PDB-ATOM class, for each line
    in the pdbfile.
    
    filename : str => contains the filename of the pdbfile
    content : list => a list with objects of the class PDBLine.
    """
    
    def __init__(self, filename : str = None) -> None:
        self.filename = filename
        self.content = []

    def read(self) -> None:
        """reads file into class variable"""
        # print(self.filename)
        with open(self.filename, 'r') as infile:
            for line in infile.readlines():
                pdbline=PDBLine(line.rstrip())
                if line[0:6].strip() in ("HETATM", "ATOM"):
                    self.content.append(PDBAtom(pdbline))
                elif len(line) > 1:
                    self.content.append(pdbline)
       
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
            raise PDBFormatError("Number of coordinates does not match number"
            "of atoms! ({0:d} != {1:d})".format(count, len(new_coordinates)))
        for atom, coords in zip(atoms, new_coordinates):
            atom.set_coords(coords)
    
    def __iter__(self) -> iter:
        """iterates over the content of the pdbfile."""
        return (iter(self.content))
    
    def iter_atoms(self, is_atom=True, is_hetatm=True) -> PDBAtom:
        """Iterates over all PDBAtom objects."""
        atoms = [pdbatom for pdbatom in self.content if
                 isinstance(pdbatom, PDBAtom)]
        for atom in atoms:
            if atom.record=="ATOM" and is_atom:
                yield(atom)
            if atom.record=="HETATM"and is_hetatm:
                yield(atom)
       
    def iter_heavy_atoms(self, is_atom=True, is_hetatm=True) -> PDBAtom:
        """Iterates over all heavy atom PDBAtom objects."""
        for atom in self.iter_atoms(is_atom=is_atom, is_hetatm=is_hetatm):
            if atom.is_hydrogen():
                continue
            yield(atom)
        
    def get_heavy_atoms(self, is_atom=True, is_hetatm=True) -> list:
        """returns a list of PDBAtom objects for all the heavy atoms."""
        return (list(self.iter_heavy_atoms(is_atom=is_atom, \
                                           is_hetatm=is_hetatm)))
