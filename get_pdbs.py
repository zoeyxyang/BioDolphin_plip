'''
python script to get the interaction (pse files) by running BioPython and Plip
'''

import argparse
import pandas as pd
import warnings
import pathlib
from pathlib import Path
from src.chain_splitter_assem import *
from src.chain_splitter_one import *
from pathlib import Path
from Bio import PDB
from Bio.PDB import PDBParser, PDBIO, Select, MMCIFParser

warnings.filterwarnings("ignore")


# Functions for getting a list of elements to select in pdb:

'''
Get a text file of the pdbs and the chains and ccd that are present in the dataset 
(those are the list of all items with the elements that are supposed to be extracted for generating selected pdbs)
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

        df_merged.to_csv("./data/assembly/selection_assembly.txt", sep='\t', index=False, header=False)
        return "./data/assembly/selection_assembly.txt"
    
    elif mode == 'one2one':
        filepath = "result/" + bdfile # read from result file (since the entries are fixed at this point)
        df = pd.read_csv(filepath, sep='\t')[["BioDolphinID", "complex_PDB_ID", "complex_Receptor_Chain", "complex_Ligand_Chain", "lipid_Ligand_ID_CCD", "complex_Residue_number_of_the_ligand", "avail_struc"]]
        #each entry will need to generate one structure
        
        pathlib.Path('./data/one2one').mkdir(parents=True, exist_ok=True) 
        df.to_csv("./data/one2one/selection_one2one.txt", sep='\t', index=False, header=False)
        
        return "./data/one2one/selection_one2one.txt"


'''
Subfunction of GetSelection
'''
#Get the unique chains/ccdID within each pdb
def unique(strings):
    unique_items = list(set(strings.split(",")))
    unique_strings = ",".join(unique_items)
    return unique_strings



# Funcitons to get the split pdbs with BioPython:

'''
Get the pdbs based on selection
'''
def GetSplitPDB(selectFile,  mode='assembly', batch=None):

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
                
                sele_pdbfile = f"./data/assembly/pdbs_selected/{pdb_id}.pdb"
                sele_ciffile = f"./data/assembly/pdbs_selected/{pdb_id}.cif"

                if Path(sele_pdbfile).is_file():
                    print(f'split pdb for {pdb_id} exist, skipping this pdb')
                    continue
                
                elif Path(sele_ciffile).is_file():
                    print(f'split cif for {pdb_id} exist, skipping this cif')
                    continue
                
                else: #the split assembly pdb/cif doens't exist
                    print(f'getting split structure for: {pdb_id}')
                    try: #try getting pdb
                        pdb_fn = pdbList.retrieve_pdb_file(pdb_id, file_format='pdb', pdir=orgpdb_dir) #fetch the pdb from the pdb_id
                        splitter.make_pdb(pdb_fn, protein_chains, lipid_chains, lipid_resnames, filetype='pdb') #create new pdbs with certain chain
                        with open("./data/assembly/path_sele_all.txt", "a") as f: #path_sele_all will store the path of those selected pdbs for later plip runs
                            f.write(outdir+pdb_id+".pdb\n")
                    except:
                        try: #if no pdb, try getting cif instead
                            cif_fn = pdbList.retrieve_pdb_file(pdb_id, file_format='mmCif', pdir=orgpdb_dir) #fetch the cif from the pdb_id
                            splitter.make_pdb(cif_fn, protein_chains, lipid_chains, lipid_resnames, filetype='cif') 
                            with open("./data/assembly/path_sele_all.txt", "a") as f: #path_sele_all will store the path of those selected pdbs for later plip runs
                                f.write(outdir+pdb_id+".cif\n")
                            
                        except:
                            print(f"######biopython can't get the structure: {pdb_id}#########")
                            with open("pdbs_nostruct.txt", "a") as f:
                                f.write(pdb_id+"\n")
          
          
                       
                        
    elif mode == 'one2one':
        outdir = "./data/one2one/pdbs_selected/" #save the selected pdbs here
        pathlib.Path(outdir).mkdir(parents=True, exist_ok=True) 
        splitter = ChainSplitter_one(outdir)
        with open(selectFile) as BD_textfile:
            for line in BD_textfile:
                strings_list = line.strip().split("\t")
                assert len(strings_list) == 7
                bd_id = strings_list[0]
                pdb_id = strings_list[1]
                protein_chain = strings_list[2]
                lipid_chain = strings_list[3]
                lipid_resname = strings_list[4]
                lipid_resnum = strings_list[5]
                avail_struc = strings_list[6]
                
                if batch is not None:
                    if pdb_id[:2] != batch:
                        continue

                #if len(strings_list) == 6:
                #    lipid_resnum = strings_list[5]
                #    print(f'getting one2one structure for biodolphin ID: {bd_id}')
                #else:
                #    continue # if there is no resnum, skip that pdb (structure usually is too big and don't have pdb file installable)
                sele_pdbfile = f"./data/one2one/pdbs_selected/{bd_id}.pdb"
                sele_ciffile = f"./data/one2one/pdbs_selected/{bd_id}.cif"

                if Path(sele_pdbfile).is_file():
                    #print(f'split pdb for {bd_id} exist, skipping this pdb')
                    continue
                
                elif Path(sele_ciffile).is_file():
                    #print(f'split cif for {bd_id} exist, skipping this cif')
                    continue

                else: #the split one-2-one pdb/cif doens't exist, generate the structure
                    print(f'<<< getting split structure for: {bd_id}')
                    if avail_struc == "pdb":
                        pdb_fn = pdbList.retrieve_pdb_file(pdb_id, file_format='pdb', pdir=orgpdb_dir) #fetch the pdb from the pdb_id
                        try:
                            #lipid_resnum = int(float(lipid_resnum))
                            splitter.make_pdb(pdb_fn, bd_id, pdb_id, protein_chain, lipid_chain, lipid_resname, lipid_resnum, filetype="pdb") #create new pdbs with certain protein and lipid
                            with open("./data/one2one/path_sele_all.txt", "a") as f:
                                f.write(outdir+bd_id+".pdb\n")
                            #print(f'got split structure for: {bd_id}.pdb>>>')
                        except Exception as e:
                            print(f"biopython can't get the one-2-one pdb structure of BDID: {bd_id} due to {e}")

                    elif avail_struc == "cif":
                        cif_fn = pdbList.retrieve_pdb_file(pdb_id, file_format='mmCif', pdir=orgpdb_dir) #fetch the pdb from the pdb_id
                        try:
                            splitter.make_pdb(cif_fn, bd_id, pdb_id, protein_chain, lipid_chain, lipid_resname, lipid_resnum, filetype="cif") #create new pdbs with certain protein and lipid
                            with open("./data/one2one/path_sele_all.txt", "a") as f:
                                f.write(outdir+bd_id+".cif\n")
                            #print(f'got split structure for: {bd_id}.cif>>>')
                        except Exception as e:
                            print(f"biopython can't get the one-2-one cif structure of BDID: {bd_id} due to {e}")
                    
                    
                    
                
            
            
            
            




    


if __name__ == "__main__":
    # setup arguments
    parser = argparse.ArgumentParser(description='Add arguments')
    parser.add_argument('-d','--dataset', help='current dataset filename (.txt) in the data directory (result for one-2-one)', default="BioDolphin_vr1.1.txt", type=str, required=False)
    parser.add_argument('--assembly', default=False, action='store_true')
    parser.add_argument('--one2one', default=False, action='store_true')
    parser.add_argument('-b','--batch', help='pdb ID number to run on', default=None, type=str, required=False)

    args = parser.parse_args()

    BDFILE = args.dataset
    ASSEMBLY = args.assembly
    ONE2ONE = args.one2one
    BATCH = args.batch
    print(f'batch serial is:{BATCH}')


    # Get selection file and Run BioPython to get selected pdbs 
    
    if ASSEMBLY:
        # Get a file that specifies what to select:
        selectFilePath = GetSelection(BDFILE, mode='assembly') #Get a text file of the pdbs and the chains and ccd selected
        GetSplitPDB(selectFilePath, mode='assembly') #Get the pdbs based on selection


    elif ONE2ONE:
        selectFilePath = GetSelection(BDFILE, mode='one2one')  #Get a text file of the pdbs and the chains and ccd selected
        GetSplitPDB(selectFilePath, mode='one2one', batch=BATCH)  #Get the pdbs based on selection

    else:
        raise Exception("no option chosen for pdb element selections")


    print('Finished running!')