#!/bin/bash

pathfile='./data/assembly/path_sele_all.txt'
fileprefix='./data/assembly/path_sele'

split -l 1000 $pathfile $fileprefix -x --additional-suffix=".txt"


# loop over all path files to submit slurm script
unset paths_array
paths_array=($(ls ./data/assembly/path_sele*))
delete='./data/assembly/path_sele_all.txt'
paths_array=( "${paths_array[@]/$delete}" )

echo 'finished splitting the path files'



for path_curr in "${paths_array[@]}"
do
    sbatch run_plip_assem_parallel.slurm "$path_curr"
done

echo 'finished submitting the slurm scripts'