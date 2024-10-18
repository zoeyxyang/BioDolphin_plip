'''
python 

> tag the biodolphin file so it has a column to show if the pdb has assembly interactions based on plip results
> do a final formatting to make sure there are no redundant entries
> print out the statistics of the final dataset
'''
import argparse
import pandas as pd
import pathlib
import warnings
import ast

warnings.filterwarnings("ignore")



def FinalFormat(df):


    return df





def GetStat(df):
    #print(df.columns)
    # total interaction entries
    print(f'Number of interactions: {df.shape[0]}')
    
    # number of unique PDBs
    unique_pdbs = set(df['complex_PDB_ID'].to_list())
    print(f'number of unique PDBs: {len(unique_pdbs)}')

    # number of unique lipid
    df_unique_lipid = df.drop_duplicates(subset='lipid_Ligand_ID_CCD')
    print(f'number of unique lipids: {df_unique_lipid.shape[0]}')

    # number of unique proteins
    unique_protien_uniprot = set(df['protein_UniProt_ID'].to_list())
    print(f'number of unique proteins: {len(unique_protien_uniprot)}')

    # unique lipids :
    print('####################### lipid statistics (unique): ###########################')
    print('lipid category counts (unique):')
    print(df_unique_lipid["lipid_Lipidmaps_categories"].value_counts())

    print('lipid molecular weights distrubution (unique):')
    df_unique_lipid['bins'] = pd.cut(x=df_unique_lipid['lipid_Molecular_weight'], bins=[0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2100],
                    labels=['0-100', 
                            '100-200', 
                            '200-300',
                            '300-400', 
                            '400-500',
                            '500-600',
                            '600-700',
                            '700-800',
                            '800-900',
                            '900-1000',
                            '1000-1100',
                            '1100-1200',
                            '1200-1300',
                            '1300-1400',
                            '1400-1500',
                            '1500-1600',
                            '1600-1700',
                            '1700-1800',
                            '1800-1900',
                            '1900-2000',
                            '2000-2100'                   
                            ])

    print(df_unique_lipid['bins'].value_counts(sort=False, normalize=True))



    # all lipids:
    print('####################### lipid statistics (all): ###########################')
    print('lipid category counts (all):')
    print(df["lipid_Lipidmaps_categories"].value_counts())

    print('lipid molecular weights distribution (all):')
    df['bins'] = pd.cut(x=df['lipid_Molecular_weight'], bins=[0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2100],
                    labels=['0-100', 
                            '100-200', 
                            '200-300',
                            '300-400', 
                            '400-500',
                            '500-600',
                            '600-700',
                            '700-800',
                            '800-900',
                            '900-1000',
                            '1000-1100',
                            '1100-1200',
                            '1200-1300',
                            '1300-1400',
                            '1400-1500',
                            '1500-1600',
                            '1600-1700',
                            '1700-1800',
                            '1800-1900',
                            '1900-2000',
                            '2000-2100'                   
                            ])
    print(df['bins'].value_counts(sort=False, normalize=True))



    #  proteins:
    print('####################### protein statistics (unique): #######################')
    df['protein_MembraneType'] = df.apply(lambda x: MergeMem(x.protein_MembraneType_UniProt, x.protein_MembraneType_DeepLoc), axis=1)


    
    ## -> organism percentages
    df_unique_protein = df.drop_duplicates(subset='protein_UniProt_ID')
    df_unique_protein["protein_Organism"].value_counts().to_csv('stat/unique_org.csv')
    print(f'saved organism stats of unique proteins to stat/uniquepro_org.csv')

    ## -> membrane type
    print( df_unique_protein["protein_MembraneType"].value_counts())


    # all proteins:
    print('####################### protein statistics (all): #######################')

    ## -> organism percentages
    df["protein_Organism"].value_counts().to_csv('stat/allpro_org.csv')
    print(f'saved organism stats of unique proteins to stat/allpro_org.csv')

    ## -> membrane type
    print( df["protein_MembraneType"].value_counts())

    #TODO: print affinity ranges
    print('####################### affinity ranges: #######################')
    print(f'Kd(nM) range: {df['complex_avgAffinity_Kd(nM)'].min()} ~ {df['complex_avgAffinity_Kd(nM)'].max()}')
    print(f'Ki(nM) range: {df['complex_avgAffinity_Ki(nM)'].min()} ~ {df['complex_avgAffinity_Ki(nM)'].max()}')
    print(f'IC50(nM) range: {df['complex_avgAffinity_IC50(nM)'].min()} ~ {df['complex_avgAffinity_IC50(nM)'].max()}')
    print(f'EC50(nM) range: {df['complex_avgAffinity_EC50(nM)'].min()} ~ {df['complex_avgAffinity_EC50(nM)'].max()}')
    print(f'-logKd/Ki range: {df['complex_avgAffinity_-logKd/Ki'].min()} ~ {df['complex_avgAffinity_-logKd/Ki'].max()}')
    print(f'Ka(M^-1) range: {df['complex_avgAffinity_Ka(M^-1)'].min()} ~ {df['complex_avgAffinity_Ka(M^-1)'].max()}')

    


# subfunction of GetStat
def MergeMem(uni, deep):

    if pd.isna(uni):
        out = deep
    else:
        out = uni
    
    out_format = []
    for term in ast.literal_eval(out):
        if term == 'Lipid anchor':
            term = 'Lipid-anchor'
        out_format.append(term)
    return out_format










if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Add arguments')
    parser.add_argument('-d','--dataset', help='current dataset filename (.txt) in the data directory', default="BioDolphin_vr1.1.txt", type=str, required=False)
    args = parser.parse_args()
    BDFILE = args.dataset

    
    #df = pd.read_csv(f'./data/{BDFILE}', sep='\t')
    print(f'Reading the dataset {BDFILE} from ./data')
    df = pd.read_csv(f'./data/{BDFILE}', sep='\t')

    print(f'tagging the entries with assembly in pdbs')
    colnames = ['complex_PDB_ID', 'pdb_has_assembly']
    assem_tags = pd.read_csv('./data/assembly/assembly.txt', names=colnames, header=None)
    new_df = df.merge(assem_tags, on='complex_PDB_ID', how='left') #Note: default is inner merge so large struc without pdb may be deleted

    print(f'formatting the dataset')
    new_df = FinalFormat(new_df)

    print(f'getting the statistics of the dataset')
    GetStat(new_df)

    #TODO: resume the below to save results
    '''
    print(f'Saving the final result as {BDFILE} in ./result')
    prefix = BDFILE.replace('.txt', '')
    outdir = "./result/" #save the selected pdbs here
    pathlib.Path(outdir).mkdir(parents=True, exist_ok=True) 
    new_df.to_csv(f'./result/{BDFILE}',  sep ='\t', index=False)
    new_df.to_csv(f'./result/{prefix}.csv', index=False)
    '''
    
    print('finished!')
    
    
    
   
    
    
    
    