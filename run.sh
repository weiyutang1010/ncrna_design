#!/bin/bash

# echo "((((...))))." | ./main --mode ncrna_design --objective pyx_sampling --init uniform --verbose --lr 0.01 -k 1000 --step 2500 > results/temp/2.txt &
# echo "((((...))))." | ./main --mode ncrna_design --objective pyx_sampling --init targeted --verbose --lr 0.005 -k 1000 --step 2500 > results/temp/3.txt &
# echo "(((...)))" | ./main --mode ncrna_design --objective pyx_sampling --init targeted --verbose -k 1000 --step 5 > result.txt &
# echo "(((((((........)))))))" | ./main --mode ncrna_design --objective pyx_sampling --init uniform --verbose -k 10000 --step 5 > result.txt &

# Check if file argument is provided
if [ $# -eq 0 ]; then
    echo "Usage: $0 <datafile>"
    exit 1
fi

# Check if file exists
if [ ! -f "$1" ]; then
    echo "File $1 not found!"
    exit 1
fi

# Read file line by line
while IFS= read -r line; do
    # Split line by space
    puzzles=($line)

    # echo "${puzzles[1]}" | ./main --mode ncrna_design --objective pyx_sampling --init uniform --verbose --lr 0.005 -k 1000 --step 2500 > results/sampling_uniform/${puzzles[0]}.txt &
    # echo "${puzzles[1]}" | ./main --mode ncrna_design --objective pyx_sampling --init targeted --verbose --lr 0.005 -k 1000 --step 2500 > results/sampling_targeted/${puzzles[0]}.txt &
    # echo "${puzzles[1]}" | ./main --mode ncrna_design --objective pyx_sampling --init uniform --verbose --lr 0.005 -k 2500 --step 2500 > results/sampling_uniform_k2500/${puzzles[0]}.txt &
    # echo "${puzzles[1]}" | ./main --mode ncrna_design --objective pyx_sampling --init targeted --verbose --lr 0.005 -k 2500 --step 2500 > results/sampling_targeted_k2500/${puzzles[0]}.txt &
    echo "${puzzles[1]}" | ./main --mode ncrna_design --objective deltaG --init uniform --verbose --lr 0.005 --step 2500 > results/deltaG_uniform/${puzzles[0]}.txt &
    echo "${puzzles[1]}" | ./main --mode ncrna_design --objective deltaG --init targeted --verbose --lr 0.005 --step 2500 > results/deltaG_targeted/${puzzles[0]}.txt &
done < "$1"