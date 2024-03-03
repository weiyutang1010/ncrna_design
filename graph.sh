#!/bin/bash

# Check if file argument is provided
if [ $# -ne 1 ]; then
    echo "Usage: $0 <result folder>"
    exit 1
fi

# # Check if file exists
# if [ ! -f "data/eterna/$1.txt" ]; then
#     echo "File data/eterna/$1.txt not found!"
#     exit 1
# fi

# Check if directory exists
# if [ ! -d "results/$1" ]; then
#     echo "Folder $1 not found!"
#     exit 1
# fi

# Check if directory exists
# if [ ! -d "./analysis/$1" ]; then
#     echo "Directory ./analysis/$1 created."
#     mkdir analysis/$1
# fi


# Read file line by line
# while IFS= read -r line; do
#     # Split line by space
#     puzzles=($line)

#     echo "${puzzles[0]} ${puzzles[1]}"
#     python analysis.py --mode graph --folder $2 --file ${puzzles[0]} > ./analysis/$2/${puzzles[0]}.txt &

# done < "data/eterna/$1.txt"

# puzzles=(1 3 8 10 11 15 20 23 26 30 33 40 41 47 57 65 66 88)

# for puzzle in "${puzzles[@]}"; do
#     # for file in "./results/$1"/*; do
#     #     if [ -f "$file" ]; then
#     #         # Print the file name
#     #         filename=$(basename "$file")
#     #         filename_without_extension="${filename%.*}"
#     #         # Print the filename
#     #         echo "$filename_without_extension"

#     #         python analysis.py --mode graph --folder $1 --file $filename_without_extension > ./analysis/$1/$filename_without_extension.txt &
#     #     fi
#     # done

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

#             python analysis.py --mode graph --folder "$1_p$puzzle" --file $filename_without_extension > ./analysis/$1_p$puzzle/$filename_without_extension.txt &
#         fi
#     done
# done

start=0
end=9

for (( seed=start; seed<=end; seed++ )); do
    if [ ! -d "./results/$1_$seed" ]; then
        echo "./results/$1_$seed is missing"
        exit 1
    fi

    if [ ! -d "./analysis/$1_$seed" ]; then
        mkdir "./analysis/$1_$seed"
    fi

    for file in "./results/$1_$seed"/*; do
        if [ -f "$file" ]; then
            # Print the file name
            filename=$(basename "$file")
            filename_without_extension="${filename%.*}"
            # Print the filename
            echo "$filename_without_extension"

            python analysis.py --mode graph --folder "$1_$seed" --file $filename_without_extension > ./analysis/$1_$seed/$filename_without_extension.txt &
        fi
    done
done
