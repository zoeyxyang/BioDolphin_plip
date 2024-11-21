'''
python>
script to get a list of paths (of the selected assembly structure) that needs plip to run it on 

(1) Remove the obsolute structures
(2) Label pdbIDs as cif/pdb based on the sturcture available
(3) Read the dataset in data, and detect which pdbs don't have the plip results yet
(4) Save the list of paths as path_plip_all.txt

'''
import argparse
import pandas as pd
import os

def DetectStruc(pdbid):
    filetype = pdb2file[pdbid]
    return filetype
    
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Add arguments')
    parser.add_argument('-d','--dataset', help='current dataset filename (.txt) in the data directory', default="BioDolphin_vr1.1.txt", type=str, required=False)
    args = parser.parse_args()
    BDFILE = args.dataset

    # read biodolphin dataset
    print(f'Reading the dataset {BDFILE} from ./data')
    df = pd.read_csv(f'./data/{BDFILE}', sep='\t')
    
    # get a list of pdbs that don't have sturctures (usually they are obsolete pdbs that will have the new pdbs in the dataset already)
    with open("pdbs_nostruct.txt", "r") as f:
        data = f.read()
        obsolete_pdbs = list(set(data.split("\n")))
        
    # delete the entries with these obsolete pdbs
    df_size_before = df.shape[0]
    df = df[~df['complex_PDB_ID'].isin(obsolete_pdbs)]
    df_size_after = df.shape[0]
    
    print(f'finished removing {df_size_before - df_size_after} entries that have obsolete pdbs')
    
    # remove ghost lipid chains (potentially from incorrect lipid chains from pdb api)
    df_size_before = df.shape[0]
    df = df.dropna(subset=['complex_Residue_number_of_the_ligand'])
    df_size_after = df.shape[0]
    
    print(f'finished removing {df_size_before - df_size_after} entries that have ghost lipid chains')
    
    
    # Add new column to indicate if the sturcture has pdb file
    pdb_selected_path = './data/assembly/pdbs_selected/'
    
    pdb2file = {}
    for path in os.scandir(pdb_selected_path):
        if path.is_file():
            pdbid = path.name[:4]
            filetype = path.name[5:]
            pdb2file[pdbid] = filetype
            
    df['avail_struc'] = df['complex_PDB_ID'].apply(DetectStruc)
    print('finished labeling pdbIDs for structure files')
    
    
    # Detect assmebly pdb files that don't have the plip results yet
    plip_result_path = './data/assembly/plip_result/'
    plip_pdbs = os.listdir(plip_result_path)
    
    all_pdbs = list(set(df['complex_PDB_ID'].tolist()))
    
    with open('./data/assembly/path_plip_all.txt', 'w') as f:
        for pdbid in all_pdbs:
            if pdbid not in plip_pdbs:                
                if pdb2file[pdbid] == 'pdb':
                    print(f'no plip available for {pdbid} that has pdb structure file')
                    f.write(f'./data/assembly/pdbs_selected/{pdbid}.pdb\n')
                
    print(f'finished generating path_plip_all.txt file')
    
    # save the dataset that has the structure file tagged
    prefix = BDFILE.replace('.txt', '')
    df.to_csv(f'./data/{prefix}_strtag.txt',  sep ='\t', index=False)
    df.to_csv(f'./data/{prefix}_strtag.csv', index=False)
    
    
    
        
    
        
        

        
        
    

