#!/bin/bash

# weiyu: only works if running from the root directory

# Check if file argument is provided
if [ $# -eq 0 ]; then
    echo "Usage: $0 <data file> <result folder>"
    exit 1
fi

# Check if file exists
if [ ! -f "data/$1" ]; then
    echo "File data/$1 not found!"
    exit 1
fi

# create results directory
if [ ! -d "results" ]; then
    mkdir results
    echo "Created directory \"results\""
fi

# create results directory
if [ ! -d "results/$2" ]; then
    mkdir results/$2
    echo "Created directory \"results/$2\""
fi

# Read file line by line
while IFS= read -r line; do
    # Split line by space
    puzzles=($line)

    # format: replace [] below
    # echo "${puzzles[1]}" | ./main [args] > results/[folder]/${puzzles[0]}.txt

    # example: softmax mode
    echo "${puzzles[1]}" | ./main --step 500 --boxplot > results/$2/${puzzles[0]}.txt

    # example_2: projection mode with adaptive lr decay
    # echo "${puzzles[1]}" | ./main --step 500 --projection --init targeted --eps 1.0 --lr 0.01 --lr_decay --adaptive_lr --lr_decay_rate 0.5 > results/$2/${puzzles[0]}.txt
done < "data/$1"

