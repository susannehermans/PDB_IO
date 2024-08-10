# PDB_IO

This repository is a custom PDB reader/writer with all the functionality that I need for 
building custom small-molecule ligands. This repository has been developed in Python3.7. 
by *Susanne Hermans*

## Installation

You can install PDB_IO directly from GitHub:

```sh
# Create a virtual environment (optional)
python -m venv pdbio
source pdbio/bin/activate

# Install directly from GitHub
pip install git+https://github.com/susannehermans/PDB_IO
```

## Usage

A simple import example is shown below. Other examples can be found in the `examples/` folder.

```python
# Import the library
from pdb_io import main
from pdb_io import PDBAtom

## Added features

### Version 1.0.1

* PDB_IO is a PDB file reader and writer which stores the atomdata of the PDB file into
class variables for easy usage. Similarly, I have developed a Mol2 file reader and writer
which reads the files into class variables which can be accessed in exactly the same way.
(Will be published soon!)

## Misc

The software is provided to you under the MIT license (see file `LICENSE.txt`).

## Reference

Please cite this github repository:

> S.M.A. Hermans (2024) PDB_IO repository https://github.com/susannehermans/PDB_IO
