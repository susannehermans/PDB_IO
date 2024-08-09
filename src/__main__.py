import sys


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


