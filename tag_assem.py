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
    #order entries based on BioDolphinIDs
    df = df.sort_values(by=['BioDolphinID'])

    return df





def GetStat(df):
    #print(df.columns)
    f = open("./stat/statistics.txt", "w")

    # total interaction entries
    f.write(f'Number of interactions: {df.shape[0]}\n')

    
    # number of unique PDBs
    unique_pdbs = set(df['complex_PDB_ID'].to_list())
    f.write(f'number of unique PDBs: {len(unique_pdbs)}\n')

    # number of unique lipid
    df_unique_lipid = df.drop_duplicates(subset='lipid_Ligand_ID_CCD')
    f.write(f'number of unique lipids: {df_unique_lipid.shape[0]}\n')

    # number of unique proteins
    unique_protien_uniprot = set(df['protein_UniProt_ID'].to_list())
    f.write(f'number of unique proteins: {len(unique_protien_uniprot)}\n')

    # unique lipids :
    f.write('####################### lipid statistics (unique): ###########################\n')
    f.write('lipid category counts (unique):\n')
    f.write(df_unique_lipid["lipid_Lipidmaps_categories"].value_counts().to_string())
    f.write('\n')

    f.write('lipid molecular weights distrubution (unique):\n')
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

    f.write(df_unique_lipid['bins'].value_counts(sort=False, normalize=True).to_string())
    f.write('\n')
    
    avg = df_unique_lipid['lipid_Molecular_weight'].mean()
    f.write(f'average molecular weight of unique lipids = {avg}')
    f.write('\n')




    # all lipids:
    f.write('####################### lipid statistics (all): ###########################\n')
    f.write('lipid category counts (all):\n')
    f.write(df["lipid_Lipidmaps_categories"].value_counts().to_string())
    f.write('\n')

    f.write('lipid molecular weights distribution (all):\n')
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
    f.write(df['bins'].value_counts(sort=False, normalize=True).to_string())
    f.write('\n')
    
    avg = df['lipid_Molecular_weight'].mean()
    f.write(f'average molecular weight of all lipids = {avg}')
    f.write('\n')





    #  proteins:
    f.write('####################### protein statistics (unique): #######################\n')
    df['protein_MembraneType'] = df.apply(lambda x: MergeMem(x.protein_MembraneType_UniProt, x.protein_MembraneType_DeepLoc), axis=1)


    
    ## -> organism percentages
    df_unique_protein = df.drop_duplicates(subset='protein_UniProt_ID')
    df_unique_protein["protein_Organism"].value_counts().to_csv('stat/unique_org.csv')
    f.write(f'saved organism stats of unique proteins to stat/uniquepro_org.csv')

    ## -> membrane type
    f.write( df_unique_protein["protein_MembraneType"].value_counts().to_string())
    f.write('\n')

    ## -> protein family
    #f.write(df_unique_protein["protein_InterPro"].value_counts().to_string())
    #f.write('\n')


    # all proteins:
    f.write('####################### protein statistics (all): #######################\n')

    ## -> organism percentages
    df["protein_Organism"].value_counts().to_csv('stat/allpro_org.csv')
    f.write(f'saved organism stats of unique proteins to stat/allpro_org.csv\n')

    ## -> membrane type
    f.write( df["protein_MembraneType"].value_counts().to_string())
    f.write('\n')

    #print affinity ranges
    f.write('####################### affinity ranges: #######################\n')
    kd_min, kd_max = df['complex_avgAffinity_Kd(nM)'].min(), df['complex_avgAffinity_Kd(nM)'].max()
    f.write(f'Kd(nM) range: {kd_min} ~ {kd_max}\n')
    ki_min, kd_max = df['complex_avgAffinity_Ki(nM)'].min(), df['complex_avgAffinity_Ki(nM)'].max()
    f.write(f'Ki(nM) range: {ki_min} ~ {kd_max}\n')
    ic50_min, ic50_max = df['complex_avgAffinity_IC50(nM)'].min(), df['complex_avgAffinity_IC50(nM)'].max()
    f.write(f'IC50(nM) range: {ic50_min} ~ {ic50_max}\n')
    ec50_min, ec50_max = df['complex_avgAffinity_EC50(nM)'].min(), df['complex_avgAffinity_EC50(nM)'].max()
    f.write(f'EC50(nM) range: {ec50_min} ~ {ec50_max}\n')
    log_min, log_max = df['complex_avgAffinity_-logKd/Ki'].min(), df['complex_avgAffinity_-logKd/Ki'].max()
    f.write(f'-logKd/Ki range: {log_min} ~ {log_max}\n')
    ka_min, ka_max = df['complex_avgAffinity_Ka(M^-1)'].min(), df['complex_avgAffinity_Ka(M^-1)'].max()
    f.write(f'Ka(M^-1) range: {ka_min} ~ {ka_max}\n')

    f.close()

    # remove extra columns
    df.drop(columns=['bins', 'protein_MembraneType'], inplace=True)

    

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
    parser.add_argument('-d','--dataset', help='current dataset filename (.txt) in the data directory', default="BioDolphin_vr1.1_expand.txt", type=str, required=False)
    args = parser.parse_args()
    BDFILE = args.dataset
    prefix = BDFILE.replace('_expand.txt', '')

    
    #df = pd.read_csv(f'./data/{BDFILE}', sep='\t')
    print(f'Reading the dataset {BDFILE} from ./data')
    df = pd.read_csv(f'./data/{BDFILE}', sep='\t')

    print(f'tagging the entries with assembly in pdbs')
    colnames = ['complex_PDB_ID', 'pdb_has_assembly']
    assem_tags = pd.read_csv('./data/assembly/assembly.txt', names=colnames, header=None) #mapping pdbid to True/False based on if it has assembly interactions
    new_df = df.merge(assem_tags, on='complex_PDB_ID', how='left') 
    assert df.shape[0] == new_df.shape[0] # adding the assembly pdbs should not change the dataset size

    print(f'formatting the dataset')
    new_df = FinalFormat(new_df)

    print(f'getting the statistics of the dataset')
    GetStat(new_df)

    # save results
    
    print(f'Saving the final result as {prefix}.txt and {prefix}.csv in ./result')
    
    outdir = "./result/" #save the selected pdbs here
    pathlib.Path(outdir).mkdir(parents=True, exist_ok=True) 
    new_df.to_csv(f'./result/{prefix}.txt',  sep ='\t', index=False)
    new_df.to_csv(f'./result/{prefix}.csv', index=False)
    
    
    print('finished!')
    
    
