#!/bin/bash


# Create an array with letters from a to z
numbers=({8..8})
letters=({a..z})

# Loop over the array
for num in "${numbers[@]}"; do
    for letter in "${letters[@]}"; do
        sbatch run_get_pdbs_one.slurm "$num$letter"
        echo "submit: ${num}${letter}"
    done
done