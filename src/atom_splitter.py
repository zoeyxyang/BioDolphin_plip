'''
BioPython functions to get pdbs that only lipids and their surrounding protein residues(based on atom index)
'''


import os
from Bio import PDB
from Bio.PDB import PDBParser, PDBIO, Select, MMCIFParser
from Bio.PDB.mmcifio import MMCIFIO



class AtomSplitter:
    def __init__(self, out_dir='./data/assembly/pdbs_small_selected/'):
        """ Create parsing and writing objects, specify output directory. """
        self.out_dir = out_dir

    def make_cif(self, cif_path, protein_atoms):
        """ Create a new cif file containing only the specified residues and lipids.
        Returns the path to the created file.
        
        pdb_path: the path of the original pdb/cif file
        protein_chains, lipid_chains, lipid_resnames: the selections to retrieve
        """
                
        # Input/output definition
        (pdb_dir, pdb_fn) = os.path.split(cif_path)
        pdb_id = pdb_fn[:4]
        out_name = pdb_id + '.cif'
        out_path = os.path.join(self.out_dir, out_name)
        print(f'saving to {out_path}')
        
        # Obtain structures
        struct =  MMCIFParser().get_structure(pdb_id, cif_path)  #get the original structure (cif)
        io = MMCIFIO()
        io.set_structure(struct)
        io.save(out_path, select=AtomSelect(protein_atoms)) #save the structure with the selected chains

        return out_path
    





# To save lipids and the surrounding atoms
class AtomSelect(PDB.Select):
    def __init__(self, protein_atoms):
        self.protein_atoms = protein_atoms
  
    def accept_chain(self, chain):
        return True

    def accept_residue(self, residue):
        return True
        
    def accept_atom(self, atom):
        parentResidue = atom.get_parent()
        
        #print(f'atom serial:{atom.get_serial_number()}')
        
        if atom.get_serial_number() in self.protein_atoms: #protein residues that are in the close range of lipids
            return True
        elif "H" in parentResidue.get_id()[0]: #heteroatom (lipid)
            return True
        else:
            return False