#!/bin/bash

# weiyu: should be run from the root directory

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

if [ ! -d "./results/eterna100_uniform" ]; then
    mkdir results/eterna100_uniform
    echo "Directory ./results/eterna100_uniform created."
fi

if [ ! -d "./analysis/eterna100_uniform" ]; then
    mkdir analysis/eterna100_uniform
    echo "Directory ./analysis/eterna100_uniform created."
fi

if [ ! -d "./graphs/eterna100_uniform" ]; then
    mkdir graphs/eterna100_uniform
    echo "Directory ./graphs/eterna100_uniform created."
fi

if [ ! -d "./results/eterna100_targeted" ]; then
    mkdir results/eterna100_targeted
    echo "Directory ./results/eterna100_targeted created."
fi

if [ ! -d "./analysis/eterna100_targeted" ]; then
    mkdir analysis/eterna100_targeted
    echo "Directory ./analysis/eterna100_targeted created."
fi

if [ ! -d "./graphs/eterna100_targeted" ]; then
    mkdir graphs/eterna100_targeted
    echo "Directory ./graphs/eterna100_targeted created."
fi

# reset the built-in timer
SECONDS=0

# Read file line by line
while IFS= read -r line; do
    # Split line by space
    puzzles=($line)

    # uniform intialization
    echo "${puzzles[1]}" | ./main --init uniform --boxplot --seed 10 > results/eterna100_uniform/${puzzles[0]}.txt
    python analysis.py --folder "eterna100_uniform" --file ${puzzles[0]}.txt > ./analysis/eterna100_uniform/${puzzles[0]}.txt

    # epsilon-targeted (eps = 0.75) intialization
    echo "${puzzles[1]}" | ./main --init targeted --eps 0.75 --boxplot --seed 10 > results/eterna100_targeted/${puzzles[0]}.txt
    python analysis.py --folder "eterna100_targeted" --file ${puzzles[0]}.txt > ./analysis/eterna100_targeted/${puzzles[0]}.txt
done < "$1"

# write the total runtime to a file
echo "Total runtime: ${SECONDS} seconds" > runtime.log