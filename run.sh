#!/bin/bash


# echo "((((...))))." | ./main --mode ncrna_design --obj pyx_exact --step 1500 --lr 0.005 --verbose --init targeted > results/exact_targeted/8.txt &
# echo "(((((.....)))))" | ./main --mode ncrna_design --obj pyx_exact --step 1500 --lr 0.005 --verbose --init targeted > results/exact_targeted/1.txt &

# Check if file argument is provided
if [ $# -eq 0 ]; then
    echo "Usage: $0 <datafile>"
    exit 1
fi

# Check if file exists
if [ ! -f "data/eterna/$1.txt" ]; then
    echo "File $1 not found!"
    exit 1
fi

# Read file line by line
# while IFS= read -r line; do
#     # Split line by space
#     puzzles=($line)

#     # echo "${puzzles[1]}" | ./main --mode ncrna_design --objective pyx_sampling --init uniform --verbose --lr 0.005 -k 1000 --step 2500 > results/sampling_uniform/${puzzles[0]}.txt &
#     # echo "${puzzles[1]}" | ./main --mode ncrna_design --objective pyx_sampling --init targeted --verbose --lr 0.005 -k 1000 --step 2500 > results/sampling_targeted/${puzzles[0]}.txt &
#     # echo "${puzzles[1]}" | ./main --mode ncrna_design --objective pyx_sampling --init uniform --verbose --lr 0.005 -k 2500 --step 2500 > results/sampling_uniform_k2500/${puzzles[0]}.txt &
#     # echo "${puzzles[1]}" | ./main --mode ncrna_design --objective pyx_sampling --init uniform --verbose --lr 0.01 -k 2500 --step 3500 > results/sampling_uniform_ex/${puzzles[0]}.txt &
#     # echo "${puzzles[1]}" | ./main --mode ncrna_design --objective pyx_sampling --init uniform --verbose --lr 0.005 -k 2500 --step 2500 > results/sampling_uniform_k2500_fixed/${puzzles[0]}.txt &
#     # echo "${puzzles[1]}" | ./main --mode ncrna_design --objective pyx_sampling --init targeted --verbose --lr 0.005 -k 2500 --step 2500 > results/sampling_targeted_k2500/${puzzles[0]}.txt &
#     # echo "${puzzles[1]}" | ./main --mode ncrna_design --objective pyx_sampling --init targeted --verbose --lr 0.01 -k 2500 --step 3500 > results/sampling_targeted_ex/${puzzles[0]}.txt &
#     # echo "${puzzles[1]}" | ./main --mode ncrna_design --objective pyx_sampling --init targeted --verbose --lr 0.005 -k 2500 --step 2500 > results/sampling_targeted_k2500_fixed/${puzzles[0]}.txt &
#     # echo "${puzzles[1]}" | ./main --mode ncrna_design --objective deltaG --init uniform --verbose --lr 0.005 --step 5000 > results/deltaG_uniform/${puzzles[0]}.txt &
#     # echo "${puzzles[1]}" | ./main --mode ncrna_design --objective deltaG --init targeted --verbose --lr 0.005 --step 5000 > results/deltaG_targeted/${puzzles[0]}.txt &

#     # echo "${puzzles[1]}" | ./main --mode test_gradient --objective deltaG --init targeted
#     echo "${puzzles[1]}" | ./main --mode ncrna_design --objective pyx_sampling --init sequence --lr 0.01 -k 1600 --step 2500 --eps 0.00 --seq ${puzzles[2]} --verbose > results/sampling_sequence_00/${puzzles[0]}.txt &
#     echo "${puzzles[1]}" | ./main --mode ncrna_design --objective pyx_sampling --init sequence --lr 0.01 -k 1600 --step 2500 --eps 0.05 --seq ${puzzles[2]} > results/sampling_sequence_05/${puzzles[0]}.txt &
#     # echo "${puzzles[1]}" | ./main --mode ncrna_design --objective pyx_sampling --init sequence --lr 0.01 -k 1600 --step 2500 --eps 0.10 --seq ${puzzles[2]} > results/sampling_sequence_10/${puzzles[0]}.txt &
#     # echo "${puzzles[1]}" | ./main --mode ncrna_design --objective pyx_sampling --init sequence --lr 0.01 -k 1600 --step 2500 --eps 0.15 --seq ${puzzles[2]} > results/sampling_sequence_15/${puzzles[0]}.txt &
# done < "data/eterna/$1.txt"

# varying k
# sample_sizes=(200 400 800 1600 3200 6400)

# while IFS= read -r line; do
#     # Split line by space
#     puzzles=($line)

#     if [ ! -d "./results/sample_size_p${puzzles[0]}" ]; then
#         mkdir "./results/sample_size_p${puzzles[0]}"
#     fi

#     for k in "${sample_sizes[@]}"; do
#         echo "${puzzles[1]}" | ./main --mode ncrna_design --obj pyx_sampling --step 2500 --lr 0.005 --k $k --init uniform > results/sample_size_p${puzzles[0]}/k$k.txt &
#     done
# done < "data/eterna/$1.txt"

# varying resampling iteration
# resample_iter=(1 2 4 8 16 32)

# while IFS= read -r line; do
#     # Split line by space
#     puzzles=($line)

#     if [ ! -d "./results/resample_iter_p${puzzles[0]}" ]; then
#         mkdir "./results/resample_iter_p${puzzles[0]}"
#     fi

#     for it in "${resample_iter[@]}"; do
#         echo "${puzzles[1]}" | ./main --mode ncrna_design --obj pyx_sampling --step 2500 --lr 0.005 --k 2500 --iter $it --init uniform > results/resample_iter_p${puzzles[0]}/iter$it.txt &
#     done
# done < "data/eterna/$1.txt"

# Random Initializations


# Run with seed random initializations (using xargs)
max_concurrent=16

execute_program() {
    read id struct seed <<< "${1}"

    if [ ! -d "./results/sampling_random_$seed" ]; then
        mkdir "./results/sampling_random_$seed"
    fi

    # if [ ! -d "./results/sampling_epsilon_$seed" ]; then
    #     mkdir "./results/sampling_epsilon_$seed"
    # fi

    # if [ ! -f "./results/sampling_random_$seed/${id}.txt" ]; then
        echo "${struct}" | ./main --mode ncrna_design --obj pyx_sampling --init random --step 2500 --lr 0.01 --k 1600 --seed "${seed}" > "./results/sampling_random_${seed}/${id}.txt"
    # fi

    # if [ ! -f "./results/sampling_epsilon_$seed/${id}.txt" ]; then
    #     echo "${struct}" | ./main --mode ncrna_design --obj pyx_sampling --init epsilon --step 2500 --lr 0.01 --k 1600 --seed "${seed}" > "./results/sampling_epsilon_${seed}/${id}.txt"
    # fi
}

export -f execute_program

declare -a lines
while IFS= read -r line; do
    lines+=("$line")
done < "data/eterna/$1_seed.txt"

printf "%s\0" "${lines[@]}" | xargs -0 -n 1 -P "$max_concurrent" bash -c 'execute_program "$@"' _



