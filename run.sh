#!/bin/bash

# weiyu: should be run from the root directory

# Check if file argument is provided
if [ $# -lt 2 ]; then
    echo "Usage: $0 <result folder> <data file>"
    exit 1
fi

# Check if data file exists
if [ ! -f "$2" ]; then
    echo "Data file $2 not found!"
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

if [ ! -d "./results/$1" ]; then
    echo "Directory ./results/$1 created."
    mkdir results/$1
fi

if [ ! -d "./analysis/$1" ]; then
    echo "Directory ./analysis/$1 created."
    mkdir analysis/$1
fi

if [ ! -d "./graphs/$1" ]; then
    echo "Directory ./graphs/$1 created."
    mkdir graph/$1
fi

# Read file line by line
while IFS= read -r line; do
    # Split line by space
    puzzles=($line)

    # example 1: softmax mode optimizing for p(y* | x) with boxplot (~30 seconds)
    echo "${puzzles[1]}" | ./main --step 200 --boxplot --num_thread 4 > results/$1/${puzzles[0]}.txt

    # example 2: softmax mode optimizing for NED with boxplot 
    # echo "${puzzles[1]}" | ./main --step 200 --obj ned --is_lazy --boxplot --num_thread 4 > results/$1/${puzzles[0]}.txt
    
    python analysis.py --folder "$1" --file ${puzzles[0]}.txt > ./analysis/$1/${puzzles[0]}.txt

done < "$2"

