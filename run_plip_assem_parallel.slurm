#!/bin/bash
# FILENAME:  plip_runs

#SBATCH -A bio230204  # Allocation name
#SBATCH --nodes=1        # Total # of nodes
#SBATCH --time=72:00:00   # Total run time limit (hh:mm:ss)
#SBATCH -J plip_run     # Job name
#SBATCH -p shared      # Queue (partition) name
#SBATCH --mail-user=useremailaddress
#SBATCH --mail-type=all   # Send email to above address at begin and end of job
#SBATCH --mem-per-cpu=8G

#go in the work dir
cd "/anvil/projects/x-bio230204/Database/BioDolphin_plip"

#put the text file into an array
pathfile=$1
unset my_array
mapfile -t my_array < ${pathfile}

echo "running structures on the path:"
echo ${pathfile}

#run plip on the pdbs listed in $my_array
singularity run plip_v2.3.0.simg -f ${my_array[@]} -yvtxp --model 1 -o "./data/assembly/plip_result/"
