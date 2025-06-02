#!/bin/bash

# weiyu: should be run from the root directory
folder_1=eterna100_uniform
folder_2=eterna100_targeted

# Check if file argument is provided
if [ $# -lt 1 ]; then
    echo "Usage: $0 <data file>"
    exit 1
fi

# Check if data file exists
if [ ! -f "$1" ]; then
    echo "Data file $1 not found!"
    exit 1
fi

# create results directory
if [ ! -d "./results" ]; then
    mkdir ./results
    echo "Created directory \"results\""
fi

# create analysis directory
if [ ! -d "./analysis" ]; then
    mkdir ./analysis/
    echo "Directory ./analysis created."
fi

# create graph directory
if [ ! -d "./graphs" ]; then
    mkdir ./graphs/
    echo "Directory ./graphs created."
fi

if [ ! -d "./results/$folder_1" ]; then
    mkdir results/$folder_1
    echo "Directory ./results/$folder_1 created."
fi

if [ ! -d "./analysis/$folder_1" ]; then
    mkdir analysis/$folder_1
    echo "Directory ./analysis/$folder_1 created."
fi

if [ ! -d "./graphs/$folder_1" ]; then
    mkdir graphs/$folder_1
    echo "Directory ./graphs/$folder_1 created."
fi

if [ ! -d "./results/$folder_2" ]; then
    mkdir results/$folder_2
    echo "Directory ./results/$folder_2 created."
fi

if [ ! -d "./analysis/$folder_2" ]; then
    mkdir analysis/$folder_2
    echo "Directory ./analysis/$folder_2 created."
fi

if [ ! -d "./graphs/$folder_2" ]; then
    mkdir graphs/$folder_2
    echo "Directory ./graphs/$folder_2 created."
fi

# reset the built-in timer
SECONDS=0

# Read file line by line
while IFS= read -r line; do
    # Split line by space
    puzzles=($line)

    # uniform intialization
    echo "${puzzles[1]}" | ./main --init uniform --boxplot --seed 10 > results/$folder_1/${puzzles[0]}.txt
    python analysis.py --folder "$folder_1" --file ${puzzles[0]}.txt > ./analysis/$folder_1/${puzzles[0]}.txt

    # epsilon-targeted (eps = 0.75) intialization
    echo "${puzzles[1]}" | ./main --init targeted --eps 0.75 --boxplot --seed 10 > results/$folder_2/${puzzles[0]}.txt
    python analysis.py --folder "$folder_2" --file ${puzzles[0]}.txt > ./analysis/$folder_2/${puzzles[0]}.txt
done < "$1"

# write the total runtime to a file
echo "Total runtime: ${SECONDS} seconds" > runtime.log