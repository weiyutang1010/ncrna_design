#!/bin/bash

# Check if file argument is provided
if [ $# -ne 2 ]; then
    echo "Usage: $0 <result folder> <data file>"
    exit 1
fi

if [ ! -f "./data/eterna/$2.txt" ]; then
    echo "data file $2.txt not found!"
    exit 1
fi

# Check if directory exists
if [ ! -d "results/$1" ]; then
    echo "Folder $1 not found!"
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
        echo "./results/$1/$file.txt"
        python analysis.py --folder "$1" --file $file > ./analysis/$1/$file.txt &
    fi
done < "data/eterna/$2.txt"

#     # varying sample size
#     if [ ! -d "./results/$1_p$puzzle" ]; then
#         echo "Folder $1 not found!"
#         exit 1
#     fi

#     if [ ! -d "./analysis/$1_p$puzzle" ]; then
#         mkdir "./analysis/$1_p$puzzle"
#     fi

#     for file in "./results/$1_p$puzzle"/*; do
#         if [ -f "$file" ]; then
#             # Print the file name
#             filename=$(basename "$file")
#             filename_without_extension="${filename%.*}"
#             # Print the filename
#             echo "$filename_without_extension"

#             python analysis.py --folder "$1_p$puzzle" --file $filename_without_extension > ./analysis/$1_p$puzzle/$filename_without_extension.txt &
#         fi
#     done
# done

# start=0
# end=0

# for (( seed=start; seed<=end; seed++ )); do
#     if [ ! -d "./results/$1_$seed" ]; then
#         echo "./results/$1_$seed is missing"
#         exit 1
#     fi

#     if [ ! -d "./analysis/$1_$seed" ]; then
#         mkdir "./analysis/$1_$seed"
#     fi

#     while IFS= read -r line; do
#         puzzles=($line)
#         file="${puzzles[0]}"

#         if [ -f "./results/$1_$seed/$file.txt" ]; then
#             # Print the file name
#             echo "./results/$1_$seed/$file.txt"
#             python analysis.py --folder "$1_$seed" --file $file > ./analysis/$1_$seed/$file.txt &
#         fi
#     done < "data/eterna/$2.txt"
# done

# Resample Iter
# if [ ! -d "./results/$1" ]; then
#     echo "./results/$1 is missing"
#     exit 1
# fi

# if [ ! -d "./analysis/$1" ]; then
#     mkdir "./analysis/$1"
# fi

# for folder_path in "./results/$1"/*; do
#     folder=$(basename "$folder_path")
#     if [ ! -d "./analysis/$1/$folder" ]; then
#         mkdir "./analysis/$1/$folder"
#     fi

#     echo $folder
#     for file in "./results/$1/$folder"/*; do
#         filename=$(basename "$file")
#         filename_without_extension="${filename%.*}"
        
#         # Print the filename
#         echo "$filename_without_extension"

#         python analysis.py --folder "$1/$folder" --file $filename_without_extension > ./analysis/$1/$folder/$filename_without_extension.txt &
#     done
# done