#!/bin/bash


#put the text file into an array
unset my_array
mapfile -t my_array < "./data/assembly/path_sele.txt" #reading all files at once

#run plip on the pdbs listed in $my_array
singularity run plip_v2.3.0.simg -f ${my_array[@]} -yvtxp --model 1 -o "./data/assembly/plip_result/"