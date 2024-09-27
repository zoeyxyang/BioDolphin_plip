'''
This script is adapted from: https://stackoverflow.com/questions/11685716/how-to-extract-chains-from-a-pdb-file

----------------------------------------------------------------
To run the script: python3 chain_splitter.py pdb_test_list.txt ./

pdb_test_list.txt: file name of the pdbs and chains needed to be extracted
./: the directory of where the output pdbs are saved


'''

import os
from Bio import PDB

class ChainSplitter:
    def __init__(self, out_dir=None):
        """ Create parsing and writing objects, specify output directory. """
        self.parser = PDB.PDBParser()
        self.writer = PDB.PDBIO()
        if out_dir is None:
            out_dir = os.path.join(os.getcwd(), "chain_PDBs")
        self.out_dir = out_dir

    def make_pdb(self, pdb_path, protein_chains, lipid_chains, lipid_resnames, overwrite=False, struct=None):
        """ Create a new PDB file containing only the specified chains.

        Returns the path to the created file.

        :param pdb_path: full path to the crystal structure
        :param chain_letters: iterable of chain characters (case insensitive)
        :param overwrite: write over the output file if it exists
        """
        

        # Input/output files
        (pdb_dir, pdb_fn) = os.path.split(pdb_path)
        pdb_id = pdb_fn[3:7]
        #out_name = "%s_%s.pdb" % (pdb_id, "".join(protein_chains))
        out_name = pdb_id + '.pdb'
        out_path = os.path.join(self.out_dir, out_name)
        
        
        # Get structure, write new file with only given chains
        if struct is None:
            struct = self.parser.get_structure(pdb_id, pdb_path)
        self.writer.set_structure(struct)
        self.writer.save(out_path, select=SelectChains_ProtLipid(protein_chains, lipid_chains, lipid_resnames)) ## save the structure with the selected chains

        return out_path



class SelectChains_ProtLipid(PDB.Select):
    """ ## modifiied version for selecting  protein chains, plus certain lipid ligands  """
    def __init__(self, protein_chains, lipid_chains, lipid_resnames):
        '''
        protein_chains: a list of chain IDs of the proteins involved
        lipid_chains: a list of chain IDs of the lipids involved
        lipid_resnames: a list of residue names of the lipids involved
        '''
        self.protein_chains = protein_chains
        self.lipid_chains = lipid_chains
        self.lipid_resnames = lipid_resnames

    def accept_chain(self, chain): ## keep the protein chain and also the lipid chain
        #return (chain.get_id() in self.chain_letters) 
        if chain.get_id() in self.protein_chains:
            return True
        if chain.get_id() in self.lipid_chains:
            return True
        else:
            False

    def accept_residue(self, residue): 

        hetero_resnames = ["H_" + res for res in self.lipid_resnames]
        
        if residue.get_id()[0] in hetero_resnames: #accept if it is one of the lipids we want 
            return True
        if residue.get_id()[0] == " ": #blank is for standard amino acids (getting the proteins)
            return True
        else:
            False # delete any heteroatom that is not the one we want












if __name__ == "__main__":
    """ Parses PDB id's desired chains, and creates new PDB structures. """
    import sys
    if not len(sys.argv) == 3:
        print ("Usage: $ python %s 'pdb.txt' 'work_dir'" % __file__)
        sys.exit()

    pdb_textfn = sys.argv[1]
    outdir = sys.argv[2]

    pdbList = PDB.PDBList()
    splitter = ChainSplitter(outdir)

    #loop over pdb_proteinchain
    with open(pdb_textfn) as pdb_textfile:
        for line in pdb_textfile:
            strings_list = line.strip().split("\t")

            pdb_id = strings_list[0].lower() #pdb
            protein_chains = strings_list[1].split(",")
            lipid_chains = strings_list[2].split(",")
            lipid_resnames = strings_list[3].split(",")

            print(f'processing: pdb_id: {pdb_id}')

            pdb_fn = pdbList.retrieve_pdb_file(pdb_id, file_format='pdb') #fetch the pdb from the pdb_id
            
            try:
                splitter.make_pdb(pdb_fn, protein_chains, lipid_chains, lipid_resnames) #create new pdbs with certain chain
            except:
                print(f"biopython can't get the structure: {pdb_id}")
                with open("pdbs_nostruct.txt", "a") as f:
                    f.write(pdb_id+"\n")
                

                
