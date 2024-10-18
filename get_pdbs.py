'''
python script to get the interaction (pse files) by running BioPython and Plip
'''

import argparse
import pandas as pd
import warnings
import pathlib
from src.chain_splitter_assem import *
from src.chain_splitter_one import *
from pathlib import Path

warnings.filterwarnings("ignore")


# Functions for getting a list of elements to select in pdb:

'''
Get a text file of the pdbs and the chains and ccd selected
'''
def GetSelection(bdfile, mode='assembly'):
    # complex_PDB_ID, complex_Receptor_Chain, complex_Ligand_Chain, lipid_Ligand_ID_CCD

    if mode == 'assembly':
        filepath = "data/" + bdfile
        df = pd.read_csv(filepath, sep='\t')[["complex_PDB_ID", "complex_Receptor_Chain", "complex_Ligand_Chain", "lipid_Ligand_ID_CCD"]]
        df_merged = df.groupby(df.complex_PDB_ID.values).transform(lambda x: ','.join(x)).drop_duplicates()
        df_merged["complex_PDB_ID"] = df_merged["complex_PDB_ID"].apply(unique)
        df_merged["complex_Receptor_Chain"] = df_merged["complex_Receptor_Chain"].apply(unique)
        df_merged["complex_Ligand_Chain"] = df_merged["complex_Ligand_Chain"].apply(unique)
        df_merged["lipid_Ligand_ID_CCD"] = df_merged["lipid_Ligand_ID_CCD"].apply(unique)
        pathlib.Path('./data/assembly').mkdir(parents=True, exist_ok=True) 

        #TODO: eleminate the ones that exist already in assembly

        df_merged.to_csv("./data/assembly/selection_assembly.txt", sep='\t', index=False, header=False)
        return "./data/assembly/selection_assembly.txt"
    
    elif mode == 'one2one':
        filepath = "data/" + bdfile
        df = pd.read_csv(filepath, sep='\t')[["BioDolphinID", "complex_PDB_ID", "complex_Receptor_Chain", "complex_Ligand_Chain", "lipid_Ligand_ID_CCD", "complex_Residue_number_of_the_ligand"]]
        #each entry will need to generate one structure
        
        pathlib.Path('./data/one2one').mkdir(parents=True, exist_ok=True) 
        df.to_csv("./data/one2one/selection_one2one.txt", sep='\t', index=False, header=False)
        
        return "./data/one2one/selection_one2one.txt"


'''
Subfunction of GetSelection
'''
def unique(strings):
    unique_items = list(set(strings.split(",")))
    unique_strings = ",".join(unique_items)
    return unique_strings



# Funcitons to get the split pdbs with BioPython:

'''
Get the pdbs based on selection
'''
def GetSplitPDB(selectFile,  mode='assembly'):

    pdbList = PDB.PDBList()

    orgpdb_dir = "./data/pdbs_original" #save the original pdbs here

    if mode == 'assembly':
        outdir = "./data/assembly/pdbs_selected/" #save the selected pdbs here
        pathlib.Path(outdir).mkdir(parents=True, exist_ok=True) 
        splitter = ChainSplitter(outdir)
        with open(selectFile) as pdb_textfile:
            for line in pdb_textfile:
                strings_list = line.strip().split("\t")
                pdb_id = strings_list[0].lower() #pdb
                protein_chains = strings_list[1].split(",")
                lipid_chains = strings_list[2].split(",")
                lipid_resnames = strings_list[3].split(",")
                
                org_pdbfile = f"./data/assembly/pdbs_selected/{pdb_id}.pdb"
                if org_pdbpath.is_file():
                    print(f'split pdb for {pdb_id} exist, skipping this pdb')
                    continue
                else:
                    print(f'getting split pdb: {pdb_id}')

                pdb_fn = pdbList.retrieve_pdb_file(pdb_id, file_format='pdb', pdir=orgpdb_dir) #fetch the pdb from the pdb_id
                
                try:
                    splitter.make_pdb(pdb_fn, protein_chains, lipid_chains, lipid_resnames) #create new pdbs with certain chain
                    with open("./data/assembly/path_sele_all.txt", "a") as f:
                        f.write(outdir+pdb_id+".pdb\n")
                except:
                    print(f"biopython can't get the structure: {pdb_id}") #TODO: these big structures don't have pdb (only cif available)
                    with open("pdbs_nostruct.txt", "a") as f:
                        f.write(pdb_id+"\n")
                        
    elif mode == 'one2one':
        outdir = "./data/one2one/pdbs_selected/" #save the selected pdbs here
        pathlib.Path(outdir).mkdir(parents=True, exist_ok=True) 
        splitter = ChainSplitter_one(outdir)
        with open(selectFile) as BD_textfile:
            for line in BD_textfile:
                strings_list = line.strip().split("\t")
                bd_id = strings_list[0]
                pdb_id = strings_list[1]
                protein_chain = strings_list[2]
                lipid_chain = strings_list[3]
                lipid_resname = strings_list[4]
                if len(strings_list) == 6:
                    lipid_resnum = strings_list[5]
                    print(f'getting one2one structure for biodolphin ID: {bd_id}')
                else:
                    continue # if there is no resnum, skip that pdb (structure usually is too big and don't have pdb file installable)
                
               
                pdb_fn = pdbList.retrieve_pdb_file(pdb_id, file_format='pdb', pdir=orgpdb_dir) #fetch the pdb from the pdb_id
                
                try:
                    lipid_resnum = int(float(lipid_resnum))
                    splitter.make_pdb(pdb_fn, bd_id, pdb_id, protein_chain, lipid_chain, lipid_resname, lipid_resnum) #create new pdbs with certain protein and lipid
                    with open("./data/one2one/path_sele_all.txt", "a") as f:
                        f.write(outdir+bd_id+".pdb\n")
                except:
                    print(f"biopython can't get the structure of db ID: {bd_id}")
                    
                
            
            
                    

            




    


if __name__ == "__main__":
    # setup arguments
    parser = argparse.ArgumentParser(description='Add arguments')
    parser.add_argument('-d','--dataset', help='current dataset filename (.txt) in the data directory', default="BioDolphin_vr1.1.txt", type=str, required=False)
    parser.add_argument('--assembly', default=False, action='store_true')
    parser.add_argument('--one2one', default=False, action='store_true')

    args = parser.parse_args()

    BDFILE = args.dataset
    ASSEMBLY = args.assembly
    ONE2ONE = args.one2one


    # Get selection file and Run BioPython to get selected pdbs 
    
    if ASSEMBLY:
        # Get a file that specifies what to select:
        selectFilePath = GetSelection(BDFILE, mode='assembly') #Get a text file of the pdbs and the chains and ccd selected
        GetSplitPDB(selectFilePath, mode='assembly') #Get the pdbs based on selection


    elif ONE2ONE:
        selectFilePath = GetSelection(BDFILE, mode='one2one')  #Get a text file of the pdbs and the chains and ccd selected
        #selectFilePath = './data/one2one/selection_one2one.txt'
        GetSplitPDB(selectFilePath, mode='one2one')  #Get the pdbs based on selection

    else:
        raise Exception("no option chosen for pdb element selections")


    print('Finished running!')