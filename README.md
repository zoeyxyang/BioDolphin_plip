# BioDolphin_plip
Run plip to get the interaction results of all structures in BioDolphin

## Set up the conda environment
`cd BioDolphin_plip` \
`conda env create -f environment.yaml` \
`conda activate plip_env` 

## Download plip software
download plip.simg in the BioDolphin_plip

## place data
Place the updated BioDolphin file in `data\`. 



## Runing assembly interactions:

### Run BioPython to get the pdbs with only proteins/lipids of interest
To get assembly interactions (lipids can interact with multiple proteins) within each pdb: \
`python get_pdbs.py --assembly` \
If running on a cluster, use the slurm script instead, example: `sbatch run_get_pdbs.slurm` \

--> This will produce pdbs with selected components in the directory `data/assembly/pdbs_selected/`\
--> `path_sele_all.txt` will store a list of paths pointing to all these pdbs with selected components



### Run plip to get the interaction files
To run plip on assembly interactions (lipids can interact with multiple proteins) within each pdb: \
`source run_plip_assem.sh` (note that this can only run if the num of structures are less than 1000)\
If running on a cluster, use the slurm script instead: `source run_plip_assem_parallel.sh` \
This script submit slurm scripts for the selected structures (1000 structures in one slurm script) \

--> This will produce plip results under `data/assembly/plip_result/`


### Format the result files into csv tables for webserver 
To get csv tables for all the result.txt files in plip_result/(pdbID): \
`python parse_result.py --assembly` \

--> This will produce assembly.txt file in ./data/assembly that gives True/False on whether that pdb has assembly interactions.



### Run scripts for checking missing plip results

`python check_result.py --assembly`






## Runing one to one interactions:

`python get_pdbs.py --one2one` \
(If running on a cluster, use the slurm script instead: `sbatch run_get_pdbs_one2on2.slurm`) \
`source run_plip_structures_one2one.sh` \
(If running on a cluster, use the slurm script instead: `sbatch run_plip_structures_one2one.slurm`) \

