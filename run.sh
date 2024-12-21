#!/bin/bash

# weiyu: should be run from the root directory

# Check if file argument is provided
if [ $# -lt 1 ]; then
    echo "Usage: $0 <data file>"
    exit 1
fi

# Check if file exists
if [ ! -f "data/$1.txt" ]; then
    echo "File data/$1 not found!"
    exit 1
fi

# create results directory
if [ ! -d "results" ]; then
    mkdir results
    echo "Created directory \"results\""
fi

# Read file line by line
while IFS= read -r line; do
    # Split line by space
    puzzles=($line)

    # example 1: softmax mode optimizing for p(y* | x) with boxplot (max 2000 steps)
    # echo "${puzzles[1]}" | ./main --step 2000 --boxplot > results/example/${puzzles[0]}.txt

    # example 2: softmax mode optimizing for NED with boxplot (max 2000 steps)
    # echo "${puzzles[1]}" | ./main --step 2000 --obj ned --boxplot > results/example_2/${puzzles[0]}.txt

    # example 3: projection mode optimizing for p(y* | x) with lr decay (max 2000 steps)
    echo "${puzzles[1]}" | ./main --step 2000 --projection --nesterov --lr_decay --adaptive_lr --boxplot > results/example_3/${puzzles[0]}.txt

done < "data/$1.txt"

