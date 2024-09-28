'''
python > tag the biodolphin file so it has a column to show if the pdb has assembly interactions based on plip results
'''
import argparse
import pandas as pd
import pathlib
import warnings
warnings.filterwarnings("ignore")






if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Add arguments')
    parser.add_argument('-d','--dataset', help='current dataset filename (.txt) in the data directory', default="BioDolphin_vr1.1.txt", type=str, required=False)
    args = parser.parse_args()
    BDFILE = args.dataset

    
    #df = pd.read_csv(f'./data/{BDFILE}', sep='\t')
    df = pd.read_csv(f'./data/{BDFILE}', sep='\t')
    colnames = ['complex_PDB_ID', 'pdb_has_assembly']
    assem_tags = pd.read_csv('./data/assembly/assembly.txt', names=colnames, header=None)
    
    #df['pdb_has_assembly'] = False
    new_df = df.merge(assem_tags, on='complex_PDB_ID')
    
    prefix = BDFILE.replace('.txt', '')
    outdir = "./result/" #save the selected pdbs here
    pathlib.Path(outdir).mkdir(parents=True, exist_ok=True) 
    new_df.to_csv(f'./result/{BDFILE}',  sep ='\t', index=False)
    new_df.to_csv(f'./result/{prefix}.csv', index=False)
    
    print('finished!')
    
    
    
   
    
    
    
    