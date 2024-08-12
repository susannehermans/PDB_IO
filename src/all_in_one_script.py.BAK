import sys

class PDBFormatError(ValueError):
	""" An Error called if any input violates PDB format conventions"""
	pass

class PDBLine:
    def __init__(self, line):
        self.line = line

    def __str__(self) -> str:
        return self.line
    
    #@classmethod
    def from_line(cls, line):
        return cls(line)

#@dataclass inherits from PDBLine
class PDBAtom(PDBLine):
    """
    Represents ATOM or HETATM entry as defined in PDB format description v3.30
    """
    
    def __init__(self, PDBLine):
        """
        Loads the data of the PDBLine into variables, one class object
        represents one Atom. NOTE: PDB-files shouldn't be splitted into columns,
        because columns might not have a separator in between.
        """
        
        self.record : str
        self.serial : int
        self.name :str
        self.altLoc : str
        self.resName : str
        self.chainID : str
        self.resSeq : str
        self.iCode : str
        self.coords : tuple
        self.occupancy : str
        self.tempFactor : str
        self.element : str
        self.charge : str
                
        self.from_line(PDBLine)

    #@classmethod
    def from_line(self, PDBLine):
        line = PDBLine.line.rstrip()
        self.record = line[0:6]  # one of ATOM or HETATM
        self.serial = int(line[6:11])
        self.name = line[12:16]
        self.altLoc = line[16:17]
        self.resName = line[17:20]
        self.chainID = line[21:22]
        self.resSeq = int(line[22:26])
        self.iCode = line[26:27]
        self.coords = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        self.occupancy = float(line[54:60])
        self.tempFactor = float(line[60:66])
        self.element = line[76:78]
        self.charge = line[78:80]
        return ()
        
    def __str__(self) -> str:
        """returns a string with atom data"""
        as_string = (
            "{record:6s}{serial:>5d} {name:4s}{altLoc:1s}{resName:3s} "
            "{chainID:1s}{resSeq:4d}{iCode:1s}   {x:>8.3f}{y:>8.3f}{z:>8.3f}"
            "{occupancy:>6.2f}{tempFactor:>6.2f}          {element:2s}{charge:2s}"
        ).format(
                record=self.record,
                serial=self.serial,
                name=self.name,
                altLoc=self.altLoc,
                resName=self.resName,
                chainID=self.chainID,
                resSeq=self.resSeq,
                iCode=self.iCode,
                x=self.coords[0],
                y=self.coords[1],
                z=self.coords[2],
                occupancy=self.occupancy,
                tempFactor=self.tempFactor,
                element=self.element,
                charge=self.charge)
        return (as_string)
        
    def get_coords(self) -> tuple:
        """returns the coordinates for this atom as a tuple."""
        return (self.coords)
    
    def set_coords(self, coords: tuple) -> None:
        """resets the coordinates of this atom."""
        if len(coords) != len(self.coords):
            raise PDBFormatError("Dimensionality of coordinates does not "
            "match! - ({0:d} != {1:d})".format(len(coords), len(self.coords)))
        self.coords = coords
    
    def is_hydrogen(self) -> bool:
        """Returns True if the atom is a hydrogen, otherwise False."""
        if self.element.strip() == "H":
            return (True)
        elif self.name.startswith("H") and not (
                self.name.startswith("He") or self.name.startswith("Hg")):
            return (True)
        else:
            return (False)


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


### MAIN SCRIPT
def __main__():
    """
    Detects the ligand "HETATM" of the pdb-file, and return ligand name
    and number. This is needed as input for CNA analyses.
    The output is a text file, in the work-directory, which contains:
    <pdbfile> <residue name> - <residue number>.
    If there are multiple ligands in the file, there will be multiple columns.
    The input is a text file with the name of a pdb-file in each line.
    """
    
    output_file = "residue_numbers.list"
    file_input = "/home/hermans/FLG_2024/prep3_structures/pdbfiles.list"
    
    # Read list of files which need to be processed
    with open(file_input, "r") as file_input_content:
        pdblist = []
        for line in file_input_content.readlines():
            pdblist.append(line.strip())
    
    # Get residue name and number of the ligands in the pdbfiles
    with open(output_file, "w") as output:
        for pdbfile_name in pdblist:
            # Read pdb file
            print("Processing '{}' ...".format(pdbfile_name))
            pdbfile = PDBFile(pdbfile_name)
            pdbfile.read()
            
            # iterate over HETATMs not ATOMs.
            resNames = set()
            for atom in pdbfile.iter_atoms(is_atom=False):
                # gets set of all residue names and numbers of all HETATMs
                resNames.add((atom.resName, atom.resSeq))
            resNames = list(resNames)
            # write output to file
            output.write(pdbfile_name)
            for tup in resNames:
                output.write(f"\t\t{tup[0]} - {tup[1]}")
            output.write("\n")
    
    print("\nFinished ... output was written to: residue_number.list\n")

if __name__ == "__main__":
    sys.exit(__main__())

