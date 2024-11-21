#!/bin/bash

#SBATCH -A bio230204  # Allocation name
#SBATCH --nodes=1        # Total # of nodes
#SBATCH --time=72:00:00   # Total run time limit (hh:mm:ss)
#SBATCH -J biopython     # Job name
#SBATCH -p shared      # Queue (partition) name
#SBATCH --mail-user=useremailaddress
#SBATCH --mail-type=all   # Send email to above address at begin and end of job
#SBATCH --mem-per-cpu=10G

#go in the work dir
cd "/anvil/projects/x-bio230204/Database/BioDolphin_plip"

#activate environment
module load anaconda/2021.05-py38
conda activate plip

#execute the python script
batch=$1 #the first two letter of pdbs to run on in this job
python get_pdbs.py -d BioDolphin_vr1.1.txt -b ${batch} --one2one > "slurm_get_pdbs_one${batch}.out"
