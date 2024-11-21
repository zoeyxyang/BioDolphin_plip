# BioDolphin_plip
Run plip to get the interaction results of all structures in BioDolphin

## Set up the conda environment
`cd BioDolphin_plip` \
`conda env create -f environment.yaml` \
`conda activate plip_env` 

## Download plip software
download plip.simg in the BioDolphin_plip

## Place data
`mkdir data`
Place the updated BioDolphin file in `data\`. 



## Getting assembly structures and interactions:
Assembly structures: all lipid proteins elements in each pdb.

**1. Getting structures with BioPython** 

`python get_pdbs.py --assembly` (or `sbatch run_get_pdbs.slurm` on cluster)

--> This will produce pdbs with selected components in the directory `data/assembly/pdbs_selected/` 

--> `path_sele_all.txt` will store a list of paths pointing to all these selected pdbs/cifs that don't exist previously before and will be created this time



<br/>

**2. (Skip) Convert cif to pdb** 

`python cif2pdb.py`

--> This currently will not work because pdb can't have multiple charector letters.


<br/>



**3. Run plip to get the interaction files** 

(1) Remove the obsolute structures and entries with ghost lipid chains\
(2) Label pdbIDs as cif/pdb based on the sturcture available\
(3) Read the dataset in data, and detect which pdbs don't have the plip results yet(cif files will not be run)\
(4) Save the list of paths as path_plip_all.txt

`python prepare_plip.py -d BioDolphin_vr1.1.txt` 

--> This will produce the file `path_plip_all.txt` under `data\assembly\`to tell plip which pdbs to run on

--> This will also produce `BioDolphin_vr*.strtag.txt or .csv` with an extra column that indicates if the pdbid has pdb or cif structure

(skip the following in this section if no new plip to run)\
`source run_plip_assem.sh` (note that this can only run if the num of structures are less than 1000)\
If the number of structures are more than 1000,  use the slurm script instead: `source run_plip_assem_parallel.sh` to run on a cluster \
This script submit slurm scripts for the selected structures (1000 structures in one slurm script)\
Those scripts reads `path_plip_all.txt` and run plip on the structures indicated there

--> This will produce plip results under `data/assembly/plip_result/`

To get csv tables for all the result.txt files in plip_result/(pdbID): \
`python parse_result.py --assembly` \
--> This will produce `assembly.txt` file in ./data/assembly that gives True/False on whether that pdb has assembly interactions. (if the pdb is not avaialble, it will not be tagged)\


<br/>


**4.Expand the csv file to include extra entries that are detected by plip**

`python expand.py`\
--> This will produce `BioDolphin_vr1.1_expand.txt` and `BioDolphin_vr1.1_expand.csv` in `/data`.




**5. Format the result files for webserver** 

`python convert.py`\
--> This will generate json files in plip directories for reading as tables\
--> This will also generate `psenames.json` for pointing paths and assigning BioDolphin entries for interaction pages for webserver


<br/>


**6.(optional) Run scripts for checking missing plip results** 

`python check_result.py --assembly` (only structures with cif should be missing)

<br/>


**7. Final format and statistics**

`python tag_assem.py`

--> tag the assembly pdbs and give a final statistics
--> save the final results in `./result`


<br/>

## Getting one-to-one Structures (for molstar visualization):

**1. Getting structures with BioPython** \
To get one to one chains for each biodolphin entry within each pdb: \
`python get_pdbs.py --one2one` (or `sbatch run_get_pdbs_one.slurm` on cluster)\








