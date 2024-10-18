import csv
import json
import os
import glob
import pandas as pd
from collections import defaultdict





'''
Go through all plip pdb directories and write the paths 
'''
def AddJson_allpdbs(resultpath='../plip_result/', outfile='../psenames.json', BioDolphin_path='./result/BioDolphin_vr1.1.txt'):
    pdb_names = os.listdir(resultpath)
    #pdb_names= ['7x2c', '2lbv'] #for small test
    FullList = []
    BD2plip_dict = BD2plip_Mapping(BioDolphin_path=BioDolphin_path)


    #i = 0
    for pdb in pdb_names:
        pdbid_lower, pdbid_upper, lipidnames, BDmapping = AddJson_pdb(pdb, resultpath, BD2plip_dict) 
        if pdbid_lower is not None:
            FullList.append(GenerateDict(pdbid_lower, pdbid_upper, lipidnames, BDmapping))
        #i +=1
        #if i >10:
        #    break

    with open(outfile, "w") as f:
        json.dump(FullList, f, indent=4)  # indent for pretty formatting

    return FullList
    

    #with open(outfile, 'w') as f:
    #    for pdb in pdb_names:
    #        pdbid_lower, pdbid_upper, lipidnames, BDmapping = AddJson_pdb(pdb, resultpath) 
    #        if pdbid_lower is not None:
    #            FullList.append(GenerateDict(pdbid_lower, pdbid_upper, lipidnames, BDmapping))
                #f.write(f'{pdbid_lower},{pdbid_upper},{lipidnames}\n')

# A subfunction of AddJson_allpdbs
def GenerateDict(pdbid_lower, pdbid_upper, lipidnames, BDmapping):
    pdb_plip_dict = {}
    pdb_plip_dict['pdbid'] = pdbid_lower
    pdb_plip_dict['pdbid_upper'] = pdbid_upper
    pdb_plip_dict['interactions'] = lipidnames
    pdb_plip_dict['BDmapping'] = BDmapping
    return pdb_plip_dict


# Function to process interaction files in one pdb directory
'''
Given a pdbID and the path for the plip results, return the pdbID, pdbID upper case, lipid codes in the pdb directory, mapping of lipid codes to the related BioDolphin ID
'''
def AddJson_pdb(pdbid, dirpath, BD2plip_dict):
    print(f'processing pdb: {pdbid}')
    pdb_path = dirpath + pdbid + '/'
    psepaths = glob.glob(f'{pdb_path}/*pse')
    psefiles = [os.path.splitext(os.path.basename(x))[0] for x in psepaths]
    lipidnames =  ["_".join(y) for y in ([x.split("_")[-3:] for x in psefiles])]
    #print(lipidnames)
    if len(lipidnames) == 0:
        return None, None, None, None

    for lipid in lipidnames:
        csv_files = glob.glob(f'{pdb_path}/*{lipid}.csv')
        lipid_json_path = f'{pdb_path}/{lipid}.json'
        csv_files_to_json(csv_files, lipid_json_path)

    BDmapping = GetBDmapping(pdbid,lipidnames, BD2plip_dict)
    return pdbid, pdbid.upper(), lipidnames, BDmapping
#AddJson_pdb('1a05', dirpath='../plip_result/')



# A subfunction of AddJson_pdb
'''
Convert multiple CSV files to a JSON file
'''
def csv_files_to_json(csv_files, json_file_path):
    combined_data = {}

    for csv_file_path in csv_files:
        file_name = os.path.splitext(os.path.basename(csv_file_path))[0]
        namelist = file_name.split("_")[:-3]
        file_name =  "_".join(namelist)  
        
        with open(csv_file_path, mode='r', encoding='utf-8') as csv_file:
            # Read the CSV file
            csv_reader = csv.DictReader(csv_file)
            file_data = {}
            
            # Initialize the dictionary for this file
            for row in csv_reader:
                for key, value in row.items():
                    if key not in file_data:
                        file_data[key] = []
                    file_data[key].append(value)

        # Store the data from this CSV file under its filename
        combined_data[file_name] = file_data

    # Write to JSON file
    with open(json_file_path, mode='w', encoding='utf-8') as json_file:
        json.dump(combined_data, json_file, indent=4)


# A subfunction of AddJson_pdb #TODO
def GetBDmapping(pdbid,lipidnames, BD2plip_dict):
    #print(f'pdbid: {pdbid}')
    #print(f'lipidnames: {lipidnames}')
    
    BD2plip_dict_pdb = BD2plip_dict[pdbid]
    #print(f'BD2plip_dict_pdb: {BD2plip_dict_pdb}')
    BDmapping = {key: value for key, value in BD2plip_dict_pdb.items() if key in lipidnames}
    #print(f'BDmapping: {BDmapping}')
    return BDmapping




    


'''
Given the BioDolphin dataset, generate the directory mapping PDBid and lipid codes (plip format) to the corresponding BioDolphin IDs
'''
def BD2plip_Mapping(BioDolphin_path='./result/BioDolphin_vr1.1.txt'):
    df = pd.read_csv(BioDolphin_path, sep='\t')
    df = df[['BioDolphinID', 'complex_PDB_ID', 'lipid_Ligand_ID_CCD', 'complex_Residue_number_of_the_ligand', 'complex_Ligand_Chain']]
    df['lipid_code'] = df.apply(lambda x: getplipcode(x.complex_PDB_ID, x.lipid_Ligand_ID_CCD, x.complex_Ligand_Chain, x.complex_Residue_number_of_the_ligand), axis=1)
        
    BD2plip_dict = {}

    for i in range(len(df)):

        bdid, pdbid, lipidcode = df.iloc[i]['BioDolphinID'],df.iloc[i]['complex_PDB_ID'], df.iloc[i]['lipid_code']

        if lipidcode is None:
            continue

        if pdbid not in BD2plip_dict: #new pdb
            BD2plip_dict[pdbid] = {}
            
        if lipidcode not in BD2plip_dict[pdbid]:
            BD2plip_dict[pdbid][lipidcode] = []
        
        BD2plip_dict[pdbid][lipidcode].append(bdid)
    return BD2plip_dict




# subfunction of BD2plip_Mapping
def getplipcode(pdbid, lipidccd, lipidchain, lipidresnum):
    try:
        lipidresnum = int(float(lipidresnum))
    except Exception as e:
        if isinstance(lipidresnum, str):
            lipidresnum = lipidresnum
        else:
            print(pdbid)
            print(e)
            return None
    return f'{lipidccd}_{lipidchain}_{lipidresnum}'


    






if __name__ == '__main__':

    print(AddJson_allpdbs(resultpath='./data/assembly/plip_result/', outfile='./data/assembly/psenames.json', BioDolphin_path='./result/BioDolphin_vr1.1.txt'))
    #TODO: change this to json format
    #mapping = {'lipidA': [], 'lipidB':[]}
    #newdict = GenerateDict("1a05", "1A05", ["lipidA, lipidB"], mapping)
    #print(newdict)

    #print(BD2plip_Mapping())