'''
BioPython functions to get pdbs that only have one protein and one lipid specified
'''


import os
from Bio import PDB
from Bio.PDB import PDBParser, PDBIO, Select, MMCIFParser
from Bio.PDB.mmcifio import MMCIFIO

class ChainSplitter_one:
    def __init__(self, out_dir):
        """ Create parsing and writing objects, specify output directory. """
        #self.parser = PDB.PDBParser()
        #self.writer = PDB.PDBIO()
        self.out_dir = out_dir

    def make_pdb(self, pdbcif_path, bd_id, pdb_id, protein_chain, lipid_chain, lipid_resname, lipid_resnum, filetype='pdb'):
        """ Create a new PDB file containing only the specified chains.
        Returns the path to the created file.
        """
        # Get structure, write new file with only given protein and lipid
        if filetype == 'pdb':
            # output files
            out_name = bd_id + '.pdb'
            out_path = os.path.join(self.out_dir, out_name)
            
            #get pdb structure from the pdbcif_path (retrieve it if not existed)       
            struct = PDBParser().get_structure(pdb_id, pdbcif_path)
            io = PDB.PDBIO()
            io.set_structure(struct)
            io.save(out_path, select=ResidueSelect(protein_chain, lipid_chain, lipid_resname, lipid_resnum)) ## save the structure with the selected chains
        
        elif filetype == 'cif':
            # output files
            out_name = bd_id + '.cif'
            out_path = os.path.join(self.out_dir, out_name)

            #get cif structure from the pdbcif_path (retrieve it if not existed)       
            struct = MMCIFParser().get_structure(pdb_id, pdbcif_path)
            io = MMCIFIO()
            io.set_structure(struct)
            io.save(out_path, select=ResidueSelect(protein_chain, lipid_chain, lipid_resname, lipid_resnum)) ## save the structure with the selected chains

        return out_path





# To save a specific ligand in a specific chain (one to one pdb)
class ResidueSelect(PDB.Select):
    def __init__(self, protein_chain, ligand_chain, ligand_resname, ligand_resnum):
        self.protein_chain = protein_chain
        self.ligand_chain = ligand_chain
        self.chains = [protein_chain, ligand_chain]
        self.ligand_resname = ligand_resname
        #self.ligand_resnum = int(ligand_resnum) #TODO: check
        self.ligand_resnum = str(ligand_resnum)

    def accept_chain(self, chain):
        if chain.get_id() in self.chains:
            return True

    def accept_residue(self, residue):
                
        parentChain = residue.get_parent().get_id()

        # get protein:
        if (parentChain == self.protein_chain) and (residue.get_id()[0] == " "):
            
            return True

        # get lipid:
        if (parentChain == self.ligand_chain) and (str(residue.get_id()[1]) == self.ligand_resnum) and (residue.resname == self.ligand_resname): 
            #print('got lipid')
            return True
                
        # Reject all other atoms
        else:
            return False











                
