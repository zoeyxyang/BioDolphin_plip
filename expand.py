'''
python>
script to expand the dataset for cross chain interaction detected by plip

'''
import argparse
import pandas as pd
import os
import glob
import json



'''
Get a dictionary to map pdbID: {lipid_code: [protein_chains]}
'''
def GetPlip_Dict(plip_result_path):
    
    plip_entry_dict_path = './data/assembly/plip_entry_dict.json'
    
    if os.path.isfile(plip_entry_dict_path):
        print(f'reading plip_entry_dict.json from ./data/assembly/')
        with open(plip_entry_dict_path, 'r') as f:
            plip_entry_dict = json.load(f)
    
    else:
        print(f'start generating plip_entry_dict.json')
        plip_entry_dict = {}
        pdb_names = os.listdir(plip_result_path)
        for pdbID in pdb_names:
            print(f'process: {pdbID}')
            plip_entry_dict[pdbID] = Lipid2ProChains(pdbID, plip_result_path)
            
        with open('./data/assembly/plip_entry_dict.json', 'w') as f:
            json.dump(plip_entry_dict, f, ensure_ascii=False, indent=4)
            
    return plip_entry_dict


'''
Subfunction of GetPlip_Dict:
to get a directory to map all lipids to the protein chains they interacts within a pdbID
'''
def Lipid2ProChains(pdbID, plip_result_path):
    plip_result_folder = f'{plip_result_path}/{pdbID}'
    lipid2chain_dict = {}
    
    psepaths = glob.glob(f'{plip_result_folder}/*pse')
    psefiles = [os.path.splitext(os.path.basename(x))[0] for x in psepaths]
    lipidnames =  ["_".join(y) for y in ([x.split("_")[-3:] for x in psefiles])]
    for lipidname in lipidnames:
        lipid2chain_dict[lipidname] = []
        lipid_csv_list = glob.glob(f'{plip_result_folder}/*{lipidname}.csv') # a list of csv that is related to that lipid
        for lipid_csv in lipid_csv_list:
            protein_chains = pd.read_csv(lipid_csv)['RESCHAIN'].to_list()
            protein_chains = list(set(protein_chains))
            lipid2chain_dict[lipidname].extend(protein_chains)
            
        lipid2chain_dict[lipidname] = list(set(lipid2chain_dict[lipidname]))
        
    return lipid2chain_dict



'''
Get a dictionary to map pdbID: {lipid_code: [protein_chains]} from biodolphin dataset
'''
def GetBD_Dict(df):
    
    bd_entry_dict_path = './data/assembly/bd_entry_dict.json'
    
    if os.path.isfile(bd_entry_dict_path):
        print(f'reading bd_entry_dict.json from ./data/assembly/')
        with open(bd_entry_dict_path, 'r') as f:
            bd_entry_dict = json.load(f)

    else:
        df_sub = df[['complex_PDB_ID', 'complex_Receptor_Chain', 'complex_Ligand_Chain', 'complex_Residue_number_of_the_ligand', 'lipid_Ligand_ID_CCD']]
        bd_entry_dict = {}
        for _, row in df_sub.iterrows():
            pdbID, protein_chain, lipid_chain, lipid_ccd = row['complex_PDB_ID'], row['complex_Receptor_Chain'], row['complex_Ligand_Chain'], row['lipid_Ligand_ID_CCD']
            try:
                lipid_resnum = int(float(row['complex_Residue_number_of_the_ligand']))
            except:
                lipid_resnum = row['complex_Residue_number_of_the_ligand']
                
            
            lipid_code = f'{lipid_ccd}_{lipid_chain}_{lipid_resnum}'
            if pdbID not in bd_entry_dict:
                bd_entry_dict[pdbID] = {}
            if lipid_code not in bd_entry_dict[pdbID]:
                bd_entry_dict[pdbID][lipid_code] = []
            bd_entry_dict[pdbID][lipid_code].append(protein_chain)
            
        with open(bd_entry_dict_path, 'w') as f:
            json.dump(bd_entry_dict, f, ensure_ascii=False, indent=4)
            
    return bd_entry_dict


'''
Get a dictionary to map pdbID: {lipid_code: [protein_chains]} that map the entries that BioDolphin is missing from plip detected interactions
'''
def GetDiff_Dict(plip_entry_dict, bd_entry_dict):
    
    add_entry_dict_path = './data/assembly/add_entry_dict.json'
    
    if os.path.isfile(add_entry_dict_path):
        print(f'reading add_entry_dict.json from ./data/assembly/')
        with open(add_entry_dict_path, 'r') as f:
            add_entry_dict_clean = json.load(f)
    
    else:
        print(f'start generating add_entry_dict.json')
        add_entry_dict = plip_entry_dict.copy()
        # loop over items in plip_entry_dict and delete the ones that already exist in bd_entry_dict
        #count = 0
        for pdbID in add_entry_dict:
            for lipidcode in add_entry_dict[pdbID]:
                try:
                    protein_chains = bd_entry_dict[pdbID][lipidcode]
                    add_entry_dict[pdbID][lipidcode] = [chain for chain in add_entry_dict[pdbID][lipidcode] if chain not in protein_chains]
                except Exception as e:
                    #print(e)
                    if 'UNK' not in lipidcode:
                        pass
                        #print(f'processing pdb {pdbID} lipid {lipidcode} not found in biodolphin')
                        #count+=1
        #print(f'number of count = {count}')
        
        #clean empty items:
        
        add_entry_dict_clean = {}
        for pdbID in add_entry_dict:
            for lipidcode in add_entry_dict[pdbID]:
                pro_chains = add_entry_dict[pdbID][lipidcode]
                if len(pro_chains) != 0: #only adding protein chains if it is not empty
                    #print(f'pro_chains= {pro_chains}')
                    if pdbID not in add_entry_dict_clean:
                        add_entry_dict_clean[pdbID] = {}
                    if lipidcode not in add_entry_dict_clean[pdbID]:
                        add_entry_dict_clean[pdbID][lipidcode] = pro_chains 
                        
                
        with open(add_entry_dict_path, 'w') as f:
                json.dump(add_entry_dict_clean, f, ensure_ascii=False, indent=4)
           
            
    return add_entry_dict_clean 



'''
Get a dataframe of the entries detected by plip that needs to be added
'''
def Dataset_toAdd(add_entry_dict_clean, df):
    
    df_toAdd_path = './data/plip_newentries.csv'
    
    if os.path.isfile(df_toAdd_path):
        print(f'reading plip_newentries.csv from ./data/')
        df_toAdd = pd.read_csv(df_toAdd_path)
    
    else:
        print(f'generating add_entry_dict.json')
        entry_dfs = []
        # setting a list of column for the lipid information
        lipid_columns_list = [col for col in df.columns if col.startswith('lipid_')]
        lipid_columns_list.extend(["complex_Residue_number_of_the_ligand", "complex_Ligand_Chain", "complex_Ligand_Serial_Number"]) 
        # setting a list of column for the protein information
        protein_columns_list = [col for col in df.columns if col.startswith('protein_')]
        protein_columns_list.extend(["complex_PubMed_ID", "complex_Resolution", "complex_PDB_ID","complex_Receptor_Chain", "complex_Receptor_asymChain", "avail_struc"]) 
        for pdbID in add_entry_dict_clean:
            for lipidcode in add_entry_dict_clean[pdbID]:
                lipidcode_list = lipidcode.split('_')
                lipidCCD, lipidChain, lipidResnum = lipidcode_list[0], lipidcode_list[1], lipidcode_list[2]
                # get the information of this lipid
                lipid_data_df = df[(df['complex_PDB_ID'] == pdbID) & (df['complex_Ligand_Chain'] == lipidChain)  & (df['lipid_Ligand_ID_CCD']== lipidCCD) & (df['complex_Residue_number_of_the_ligand']== lipidResnum)]
                #lipid_data_df= lipid_data_df[(lipid_data_df['complex_Residue_number_of_the_ligand']== lipidResnum) | (lipid_data_df['complex_Residue_number_of_the_ligand']== str(float(lipidResnum)))]            
                #((df['complex_Residue_number_of_the_ligand']== int(lipidResnum)) | (df['complex_Residue_number_of_the_ligand']== float(lipidResnum)) | (df['complex_Residue_number_of_the_ligand']== lipidResnum))
                lipid_data_df = lipid_data_df[lipid_columns_list].head(1)
              
                
                # get the information of the protein chains
                for proChain in add_entry_dict_clean[pdbID][lipidcode]:
                    proChain = str(proChain)
                    print(f'getting new entry: pdb {pdbID}, lipid {lipidcode}, protein chain {proChain}')
                    protein_data_df = df[(df['complex_PDB_ID'] == pdbID) & (df['complex_Receptor_Chain'] == proChain)]
                    protein_data_df = protein_data_df[protein_columns_list].head(1)
                    # combine the dataframe into one row 
                    combined_data = pd.concat([lipid_data_df.reset_index(drop=True), protein_data_df.reset_index(drop=True)], axis=1)
                    
                    if lipid_data_df.shape[0] != 1:
                        print(f'####lipid data has {lipid_data_df.shape[0]} rows')
                        continue
                        #raise Exception("incorrect lipid data rows ")
                    if protein_data_df.shape[0] != 1:
                        print(f'####protein data has {protein_data_df.shape[0]} rows')
                        continue
                        #raise Exception("incorrect protein data rows ")
                    
                    entry_dfs.append(combined_data)
  

        df_toAdd = pd.concat(entry_dfs, axis=0, ignore_index=True)
        df_toAdd['entry_source'] = 'plip'
        df_toAdd.reset_index(drop=True, inplace=True)
        df_toAdd.to_csv("./data/plip_newentries.csv") 
       
    
    return df_toAdd


'''
Function to assign BioDolphin ID
'''
def getID(pdb, chain_rec, chain_lig, CCD_lig, serial_lig):
    ID = f'BD{pdb}-{chain_rec}-{chain_lig}-{CCD_lig}{serial_lig}'
    return ID




if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Add arguments')
    parser.add_argument('-d','--dataset', help='taggged dataset filename BioDolpin_vr1.1_strtag.txt in the data directory', default="BioDolphin_vr1.1_strtag.txt", type=str, required=False)
    args = parser.parse_args()
    BDFILE_STRTAG = args.dataset
    prefix=BDFILE_STRTAG.replace("_strtag.txt", "")
    
    print(f'reading the file: {BDFILE_STRTAG}')
    df = pd.read_csv(f'./data/{BDFILE_STRTAG}', sep='\t')
    
    print(f'prefix of the data is: {prefix}')
    
    
    # change data type for resnum
    def ChangeType(num):
        try:
            out = str(int(num))
        except:
            out = num   
        return out
    
    df['complex_Residue_number_of_the_ligand'] = df['complex_Residue_number_of_the_ligand'].apply(ChangeType)
    
    # Change other columns in df to be strings so the matching can be correct later
    df['complex_Ligand_Chain'] = df['complex_Ligand_Chain'].astype(str)
    df['lipid_Ligand_ID_CCD'] = df['lipid_Ligand_ID_CCD'].astype(str)
    df['complex_Receptor_Chain'] = df['complex_Receptor_Chain'].astype(str)

    # Get the expanded dataset
    plip_result_path = './data/assembly/plip_result'
    plip_entry_dict = GetPlip_Dict(plip_result_path) #-> get dictionary: plip_entry_dict.json
    bd_entry_dict = GetBD_Dict(df) #-> get dictionary: bd_entry_dict.json
    add_entry_dict_clean = GetDiff_Dict(plip_entry_dict, bd_entry_dict) #-> get dictionary: add_entry_dict.json
    add_df = Dataset_toAdd(add_entry_dict_clean, df) # new dataframe of the entries to add
    
    # Assign BDid
    add_df["BioDolphinID"] = add_df.apply(lambda x: getID(x.complex_PDB_ID, x.complex_Receptor_Chain, x.complex_Ligand_Chain, x.lipid_Ligand_ID_CCD, x.complex_Ligand_Serial_Number), axis=1)

    # Merge the dataframes and drop duplicates of BioDolphinID
    merged_df = pd.concat([df, add_df], axis=0) 
    merged_df = merged_df.drop_duplicates(subset='BioDolphinID')

    # Save the dataset
    merged_df.to_csv(f'./data/{prefix}_expand.txt',  sep ='\t', index=False)
    merged_df.to_csv(f'./data/{prefix}_expand.csv', index=False)
 
    
    
    