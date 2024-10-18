# BioDolphin_plip
Run plip to get the interaction results of all structures in BioDolphin

## Set up the conda environment
`cd BioDolphin_plip` \
`conda env create -f environment.yaml` \
`conda activate plip_env` 

## Download plip software
download plip.simg in the BioDolphin_plip

## Place data
Place the updated BioDolphin file in `data\`. 



## Getting assembly structures and interactions:

**1. Getting structures with BioPython** \
To get assembly chains (lipids can interact with multiple proteins) within each pdb: \
`python get_pdbs.py --assembly` (or `sbatch run_get_pdbs.slurm` on cluster)\ 
--> This will produce pdbs with selected components in the directory `data/assembly/pdbs_selected/`\
--> `path_sele_all.txt` will store a list of paths pointing to all these pdbs with selected components 



<br/>



**2. Run plip to get the interaction files** \
To run plip on assembly interactions (lipids can interact with multiple proteins) within each pdb: \
`source run_plip_assem.sh` (note that this can only run if the num of structures are less than 1000)\
If running on a cluster, use the slurm script instead: `source run_plip_assem_parallel.sh` \
This script submit slurm scripts for the selected structures (1000 structures in one slurm script) \
--> This will produce plip results under `data/assembly/plip_result/`

<br/>

**3. Format the result files into for webserver** \
To get csv tables for all the result.txt files in plip_result/(pdbID): \
`python parse_result.py --assembly` \
--> This will produce assembly.txt file in ./data/assembly that gives True/False on whether that pdb has assembly interactions. (if the pdb is not avaialble, it will not be tagged)\
To get json files and produce files for finding paths:\
`python convert.py`


<br/>




**4.Run scripts for checking missing plip results** 

`python check_result.py --assembly`


<br/>


## Getting one-to-one Structures and Interactions:

**1. Getting structures with BioPython** \
To get one to one chains for each biodolphin entry within each pdb: \
`python get_pdbs.py --one2one` (or `sbatch run_get_pdbs_one.slurm` on cluster)\

#TODO: plip interactions for one2one pdbs





## Final format and statistics

`python tag_assem.py`


