#!/bin/bash

# Check if file argument is provided
if [ $# -ne 2 ]; then
    echo "Usage: $0 <data file> <result folder>"
    exit 1
fi

# Check if file exists
if [ ! -f "data/eterna/$1" ]; then
    echo "File data/eterna/$1 not found!"
    exit 1
fi

# Check if directory exists
if [ ! -d "results/$2" ]; then
    echo "Folder $2 not found!"
    exit 1
fi

# Check if directory exists
if [ ! -d "analysis/$2" ]; then
    mkdir analysis/$2
fi


# Read file line by line
while IFS= read -r line; do
    # Split line by space
    puzzles=($line)

    echo "${puzzles[0]} ${puzzles[1]}"
    python analysis.py --mode graph --folder $2 --file ${puzzles[0]} > ./analysis/$2/${puzzles[0]}.txt

    # echo "${puzzles[1]}" | ./main --mode ncrna_design --objective pyx_sampling --init uniform --verbose --lr 0.005 -k 1000 --step 2500 > results/sampling_uniform/${puzzles[0]}.txt &
    # echo "${puzzles[1]}" | ./main --mode ncrna_design --objective pyx_sampling --init targeted --verbose --lr 0.005 -k 1000 --step 2500 > results/sampling_targeted/${puzzles[0]}.txt &
done < "data/eterna/$1"