from typing import Tuple
from .PDBLine import PDBLine
from .PDBFormatError import PDBFormatError

class PDBAtom(PDBLine):
    """
    Represents ATOM or HETATM entry as defined in PDB format description v3.30
    """
    
    def __init__(self, line : str):
        """
        Loads the data of the PDBLine into variables, one class object
        represents one Atom. NOTE: PDB-files shouldn't be splitted into columns,
        because columns might not have a separator in between.
        """
        self.line = line.rstrip()
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
        
    def get_coords(self) -> Tuple[float]:
        """returns the coordinates for this atom as a tuple."""
        return (self.coords)
    
    def set_coords(self, coords: Tuple) -> None:
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
