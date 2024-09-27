
import pandas as pd
import argparse
import os






def GetPDBids(file, sep="\t"):
    df = pd.read_csv(file, sep=sep)
    pdbList = df.iloc[:, 0].tolist()

    return pdbList






if __name__ == "__main__":

    # setup arguments
    parser = argparse.ArgumentParser(description='Add arguments')
    parser.add_argument('--assembly', default=False, action='store_true')
    parser.add_argument('--one2one', default=False, action='store_true')

    args = parser.parse_args()

    ASSEMBLY = args.assembly
    ONE2ONE = args.one2one

    
    if ASSEMBLY:
        requirePDBList = set(GetPDBids('./data/assembly/selection_assembly.txt', sep="\t"))

        plip_result_path = './data/assembly/plip_result/'
        resultPDBList = set([item for item in os.listdir(plip_result_path)])

        missing = requirePDBList - resultPDBList

        print(missing)
        print(f'number of pdbs missing: {len(missing)}')
       
    