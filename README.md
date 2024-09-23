# BioDolphin_plip
Run plip to get the interaction results of all structures in BioDolphin

### Set up the conda environment
`cd BioDolphin_plip` \
`conda env create -f environment.yaml` \
`conda activate plip_env` 


### place data
Place the updated BioDolphin file in `data\`. \
Make sure to put only one BioDolphin file so the program won't get confused.


### Run BioPython to get the pdbs with only proteins/lipids of interest
To get assembly interactions (lipids can interact with multiple proteins) within each pdb: \
`python get_pdbs.py --assembly` \
(If running on a cluster, use the slurm script instead: `sbatch run_get_pdbs_assembly.slurm`) \

To get one-to-one interactions: \
`python get_pdbs.py --one2one` \
(If running on a cluster, use the slurm script instead: `sbatch run_get_pdbs_one2on2.slurm`)


--> This will produce pdbs with specific components defined by the above scripts in the directory `data/pdbs_selected/assembly` or `data/pdbs_selected/one2one` \


### Run plip to get the interaction files

To run plip on assembly interactions (lipids can interact with multiple proteins) within each pdb: \
`source run_plip_structures_assembly.sh` \
(If running on a cluster, use the slurm script instead: `sbatch run_plip_structures_assembly.slurm`) \

To run plip on one-to-one interactions: \
`source run_plip_structures_one2one.sh` \
(If running on a cluster, use the slurm script instead: `sbatch run_plip_structures_one2one.slurm`) \