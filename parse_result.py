'''
Python >
Parse the result txt file into csv table format
'''


import argparse
import pandas as pd
import os



def ParseResult(resultfile, output_dir, findAssem=False):
    csv_list = []
    COLLECT = 0
    ligandtype = None
    PDB_hasASSEM = False

    with open(resultfile) as file:
        for i, line in enumerate(file):
            line = line.rstrip()
            if i == 0:
                pdbID = line.split(" ")[-1][:4].lower()
            if "SMALLMOLECULE" in line:
                CCD_chain_resnum, ligandtype = ParseLigand(line)

            if ligandtype == "SMALLMOLECULE":

                if "Interacting chain(s)" in line:
                    interacting_chain = ParseInterChain(line)
                    if "," in interacting_chain and findAssem:
                        PDB_hasASSEM = True
                    
                elif "**" in line:
                    interacting_type = ParseInterType(line)
                    csv_list = []
                    NAME_CSV = f'{interacting_type}_{CCD_chain_resnum}.csv'
   
                elif "|" in line:
                    COLLECT -= 1
                    entry = ParseEntry(line)
                    csv_list.append(entry)

                elif "+" in line:
                    COLLECT += 1
                
                if line == "" and COLLECT == 1:
                    df = pd.DataFrame(csv_list[1:], columns=csv_list[0])
                    df.to_csv(output_dir + NAME_CSV, index=False)
                    COLLECT = 0
    if findAssem:
        return PDB_hasASSEM
    



'''
Subfunction for ParseResult >>>
'''

def ParseLigand(line):
    ligandinfo = line.split(" ")[0]
    ligandinfo = ligandinfo.replace(":", "_")
    ligandtype = line.split("-")[-1].strip()
    return ligandinfo, ligandtype


def ParseInterChain(line):
    protchains = line.split(" ")[-1]
    return protchains


def ParseInterType(line):
    interactiontype = line.replace("**", "")
    interactiontype = interactiontype.replace(" ","_")
    return interactiontype


def ParseEntry(line):
    entrylist = [x.strip() for x in line.split('|') if x.strip() != ""]
    return entrylist


'''
<<< Subfunction for ParseResult
'''



def WriteAssem(pdb, hasAssem, outfile):
    with open(outfile, 'a') as f:
        f.write(pdb + ',' + str(hasAssem) + "\n")
        
        


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Add arguments')
    parser.add_argument('--assembly', default=False, action='store_true')
    parser.add_argument('--one2one', default=False, action='store_true')

    args = parser.parse_args()
    
    if args.assembly:
        plip_result_path = "./data/assembly/plip_result/"
        dirlist = [item for item in os.listdir(plip_result_path)]
        for d_name in dirlist:
            #print(f'parsing the result file for directory: {d_name}')
            dir_path = plip_result_path + d_name + "/"
            result_file = dir_path + "report.txt"
            try:
                pdb_hasASSEM = ParseResult(result_file, output_dir=dir_path, findAssem=True)
                WriteAssem(d_name, pdb_hasASSEM, outfile="./data/assembly/assembly.txt")
            except:
                print(f'pdb: {d_name} has no result.txt file')
            #print(f'pdb: {d_name} has assembly: {pdb_hasASSEM}')
            
            #break


    if args.one2one:
        plip_result_path = "./data/one2one/plip_result/"
        pass #TODO: not implemeted yet


