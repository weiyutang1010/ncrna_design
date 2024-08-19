#!/bin/bash

# Check if result folder is provided
if [ $# -eq 0 ]; then
    echo "Usage: $0 <result folder> <data file (optional)>"
    exit 1
fi

# Check if directory exists
if [ ! -d "results/$1" ]; then
    echo "Folder $1 not found!"
    exit 1
fi

if [ $# -eq 1 ]; then
    # Check if directory exists
    if [ ! -d "./analysis/$1" ]; then
        echo "Directory ./analysis/$1 created."
        mkdir analysis/$1
    fi

    for file in "./results/$1"/*; do
        # Print the file name
        filename=$(basename "$file")
        filename_without_extension="${filename%.*}"

        if [ -f "./results/$1/$filename_without_extension.txt" ]; then
            # Print the file name
            echo "Processing ./results/$1/$filename_without_extension.txt"
            python analysis.py --folder "$1" --file $filename_without_extension.txt > ./analysis/$1/$filename_without_extension.txt &
        fi
    done
fi

if [ $# -eq 2 ]; then
    if [ ! -f "./data/eterna/$2.txt" ]; then
        echo "data file $2.txt not found!"
        exit 1
    fi

    # Check if directory exists
    if [ ! -d "./analysis/$1" ]; then
        echo "Directory ./analysis/$1 created."
        mkdir analysis/$1
    fi

    while IFS= read -r line; do
        puzzles=($line)
        file="${puzzles[0]}"

        if [ -f "./results/$1/$file.txt" ]; then
            # Print the file name
            echo "Processing ./results/$1/$file.txt"
            python analysis.py --folder "$1" --file $file.txt > ./analysis/$1/$file.txt &
        fi
    done < "data/eterna/$2.txt"
fi
