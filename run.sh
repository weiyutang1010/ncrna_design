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
while IFS= read -r line; do
    # Split line by space
    puzzles=($line)

    # echo "${puzzles[1]}" | ./main --mode ncrna_design --objective pyx_sampling --init uniform --verbose --lr 0.005 -k 1000 --step 2500 > results/sampling_uniform/${puzzles[0]}.txt &
    # echo "${puzzles[1]}" | ./main --mode ncrna_design --objective pyx_sampling --init targeted --verbose --lr 0.005 -k 1000 --step 2500 > results/sampling_targeted/${puzzles[0]}.txt &
    # echo "${puzzles[1]}" | ./main --mode ncrna_design --objective pyx_sampling --init uniform --verbose --lr 0.005 -k 2500 --step 2500 > results/sampling_uniform_k2500/${puzzles[0]}.txt &
    # echo "${puzzles[1]}" | ./main --mode ncrna_design --objective pyx_sampling --init uniform --verbose --lr 0.01 -k 2500 --step 3500 > results/sampling_uniform_ex/${puzzles[0]}.txt &
    # echo "${puzzles[1]}" | ./main --mode ncrna_design --objective pyx_sampling --init uniform --verbose --lr 0.005 -k 2500 --step 2500 > results/sampling_uniform_k2500_fixed/${puzzles[0]}.txt &
    # echo "${puzzles[1]}" | ./main --mode ncrna_design --objective pyx_sampling --init targeted --verbose --lr 0.005 -k 2500 --step 2500 > results/sampling_targeted_k2500/${puzzles[0]}.txt &
    # echo "${puzzles[1]}" | ./main --mode ncrna_design --objective pyx_sampling --init targeted --verbose --lr 0.01 -k 2500 --step 3500 > results/sampling_targeted_ex/${puzzles[0]}.txt &
    # echo "${puzzles[1]}" | ./main --mode ncrna_design --objective pyx_sampling --init targeted --verbose --lr 0.005 -k 2500 --step 2500 > results/sampling_targeted_k2500_fixed/${puzzles[0]}.txt &
    # echo "${puzzles[1]}" | ./main --mode ncrna_design --objective deltaG --init uniform --verbose --lr 0.005 --step 5000 > results/deltaG_uniform/${puzzles[0]}.txt &
    # echo "${puzzles[1]}" | ./main --mode ncrna_design --objective deltaG --init targeted --verbose --lr 0.005 --step 5000 > results/deltaG_targeted/${puzzles[0]}.txt &

    # echo "${puzzles[1]}" | ./main --mode test_gradient --objective deltaG --init targeted
    # echo "${puzzles[1]}" | ./main --mode ncrna_design --objective pyx_sampling --init sequence --lr 0.01 -k 1600 --step 2500 --eps 0.00 --seq ${puzzles[2]} --verbose > results/sampling_sequence_00/${puzzles[0]}.txt &
    # echo "${puzzles[1]}" | ./main --mode ncrna_design --objective pyx_sampling --init sequence --lr 0.01 -k 1600 --step 2500 --eps 0.05 --seq ${puzzles[2]} > results/sampling_sequence_05/${puzzles[0]}.txt &
    # echo "${puzzles[1]}" | ./main --mode ncrna_design --objective pyx_sampling --init sequence --lr 0.01 -k 1600 --step 2500 --eps 0.10 --seq ${puzzles[2]} > results/sampling_sequence_10/${puzzles[0]}.txt &
    # echo "${puzzles[1]}" | ./main --mode ncrna_design --objective pyx_sampling --init sequence --lr 0.01 -k 1600 --step 2500 --eps 0.15 --seq ${puzzles[2]} > results/sampling_sequence_15/${puzzles[0]}.txt &
    # echo "${puzzles[1]}" | ./main --mode ncrna_design --objective pyx_sampling --init random --lr 0.01 -k 1600 --step 2500 > results/sampling_time/${puzzles[0]}.txt
    
    # echo "${puzzles[1]}" | ./main --mode ncrna_design --objective pyx_sampling --init uniform --test -k 1250 --step 1500 > /dev/null 2> ./variance/p${puzzles[0]}.txt &
    
    # echo "${puzzles[1]}" | ./main --objective pyx_direct --init uniform -k 2500 --lr 0.001 > results/sampling_pyx_uniform/${puzzles[0]}.txt &
    # echo "${puzzles[1]}" | ./main --objective pyx_direct --init targeted -k 2500 --lr 0.001 > results/sampling_pyx_targeted/${puzzles[0]}.txt &

    # echo "${puzzles[1]}" | ./main --objective pyx_direct --init uniform --best_k 2500 -k 2500 --step 2000 --lr 0.01 > results/boxplot_uniform/${puzzles[0]}.txt &
    # echo "${puzzles[1]}" | ./main --objective pyx_direct --init targeted --best_k 2500 -k 2500 --step 2000 --lr 0.01 > results/boxplot_targeted/${puzzles[0]}.txt &

    # 
    echo "${puzzles[1]}" | ./main --init uniform_sm -k 2500 --step 2000 --lr 0.001 > results/sampling_pyx_uniform_sm/${puzzles[0]}.txt &
    echo "${puzzles[1]}" | ./main --init uniform_sm -k 2500 --step 2000 --lr 0.1 --softmax > results/sampling_pyx_uniform_sm_softmax/${puzzles[0]}.txt &
    # echo "${puzzles[1]}" | ./main --init targeted_sm -k 2500 --lr 0.001 > results/sampling_pyx_targeted_sm/${puzzles[0]}.txt &
done < "data/eterna/$1.txt"

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
# resample_iter=(16)

# while IFS= read -r line; do
#     # Split line by space
#     puzzles=($line)

#     if [ ! -d "./results/resample_iter/resample_iter_p${puzzles[0]}_uniform" ]; then
#         mkdir "./results/resample_iter/resample_iter_p${puzzles[0]}_uniform"
#     fi

#     # if [ ! -d "./results/resample_iter/resample_iter_p${puzzles[0]}_targeted" ]; then
#     #     mkdir "./results/resample_iter/resample_iter_p${puzzles[0]}_targeted"
#     # fi

#     for it in "${resample_iter[@]}"; do
#         echo "${puzzles[1]}" | ./main --test --lr 0.001 -k 2500 --iter $it --init uniform > ./results/resample_iter/resample_iter_p${puzzles[0]}_uniform/iter$it.txt &
#         # echo "${puzzles[1]}" | ./main --test --lr 0.001 -k 2500 --iter $it --init targeted > ./results/resample_iter/resample_iter_p${puzzles[0]}_targeted/iter$it.txt &
#     done
# done < "data/eterna/$1.txt"

# Random Initializations

# Run with seed random initializations (using xargs)
# max_concurrent=24

# execute_program() {
#     read id struct seed <<< "${1}"

#     if [ ! -d "./results/sampling_pyx_random_sm_$seed" ]; then
#         mkdir "./results/sampling_pyx_random_sm_$seed"
#     fi

#     if [ ! -f "./results/sampling_pyx_random_sm_$seed/${id}.txt" ]; then
#         # n50
#         # echo "${struct}" | ./main --init random_sm --lr 0.001 --k 1500 --seed "${seed}" > "./results/sampling_pyx_random_sm_${seed}/${id}.txt"
        
#         # n100
#         echo "${struct}" | ./main --init random_sm --lr 0.001 --k 2500 --seed "${seed}" > "./results/sampling_pyx_random_sm_${seed}/${id}.txt"

#         # n100
#         # echo "${struct}" | ./main --obj pyx_direct --init random --lr 0.001 --k 2500 --seed "${seed}" > "./results/sampling_pyx_random_${seed}/${id}.txt"
#     fi

#     # read id struct eps idx <<< "${1}"

#     # if [ ! -d "./results/sampling_epsilon_$idx" ]; then
#     #     mkdir "./results/sampling_epsilon_$idx"
#     # fi

#     # if [ ! -f "./results/sampling_epsilon_$idx/${id}.txt" ]; then
#     #     echo "${struct}" | ./main --mode ncrna_design --obj pyx_sampling --init epsilon --step 2500 --lr 0.01 --k 1600 --eps "${eps}" > "./results/sampling_epsilon_${idx}/${id}.txt"
#     # fi
# }

# export -f execute_program

# declare -a lines
# while IFS= read -r line; do
#     lines+=("$line")
# done < "data/eterna/$1_seed.txt"
# # done < "data/eterna/$1_eps.txt"

# printf "%s\0" "${lines[@]}" | xargs -0 -n 1 -P "$max_concurrent" bash -c 'execute_program "$@"' _


# lr
# p61
# echo ".....((.((.((...)).((...))...)))).((.((.((...)).((...))...))))....." | ./main --obj pyx_direct --init targeted --step 2000 -k 2500 --lr 0.01 > results/lr/p61_targeted_lr_010.txt &
# echo ".....((.((.((...)).((...))...)))).((.((.((...)).((...))...))))....." | ./main --obj pyx_direct --init targeted --step 2000 -k 2500 --lr 0.005 > results/lr/p61_targeted_lr_005.txt &
# echo ".....((.((.((...)).((...))...)))).((.((.((...)).((...))...))))....." | ./main --obj pyx_direct --init targeted --step 2000 -k 2500 --lr 0.001 > results/lr/p61_targeted_lr_001.txt &

# echo ".....((.((.((...)).((...))...)))).((.((.((...)).((...))...))))....." | ./main --obj pyx_direct --init uniform --step 2000 -k 2500 --lr 0.01 > results/lr/p61_uniform_lr_010.txt &
# echo ".....((.((.((...)).((...))...)))).((.((.((...)).((...))...))))....." | ./main --obj pyx_direct --init uniform --step 2000 -k 2500 --lr 0.005 > results/lr/p61_uniform_lr_005.txt &
# echo ".....((.((.((...)).((...))...)))).((.((.((...)).((...))...))))....." | ./main --obj pyx_direct --init uniform --step 2000 -k 2500 --lr 0.001 > results/lr/p61_uniform_lr_001.txt &

# p71
# echo "....((....((....((....((..((..((...))..))....))......)).......)).........))............." | ./main --obj pyx_direct --init targeted --step 2000 -k 2500 --lr 0.01 > results/lr/p71_targeted_lr_010.txt &
# echo "....((....((....((....((..((..((...))..))....))......)).......)).........))............." | ./main --obj pyx_direct --init targeted --step 2000 -k 2500 --lr 0.005 > results/lr/p71_targeted_lr_005.txt &
# echo "....((....((....((....((..((..((...))..))....))......)).......)).........))............." | ./main --obj pyx_direct --init targeted --step 2000 -k 2500 --lr 0.001 > results/lr/p71_targeted_lr_001.txt &

# echo "....((....((....((....((..((..((...))..))....))......)).......)).........))............." | ./main --obj pyx_direct --init uniform --step 2000 -k 2500 --lr 0.01 > results/lr/p71_uniform_lr_010.txt &
# echo "....((....((....((....((..((..((...))..))....))......)).......)).........))............." | ./main --obj pyx_direct --init uniform --step 2000 -k 2500 --lr 0.005 > results/lr/p71_uniform_lr_005.txt &
# echo "....((....((....((....((..((..((...))..))....))......)).......)).........))............." | ./main --obj pyx_direct --init uniform --step 2000 -k 2500 --lr 0.001 > results/lr/p71_uniform_lr_001.txt &

# p14
# echo ".....((...((....))...((....))...((....))...((....))...((....))...((....))...))...................." | ./main --obj pyx_direct --init targeted --step 2000 -k 2500 --lr 0.01 > results/lr/p14_targeted_lr_010.txt &
# echo ".....((...((....))...((....))...((....))...((....))...((....))...((....))...))...................." | ./main --obj pyx_direct --init targeted --step 2000 -k 2500 --lr 0.005 > results/lr/p14_targeted_lr_005.txt &
# echo ".....((...((....))...((....))...((....))...((....))...((....))...((....))...))...................." | ./main --obj pyx_direct --init targeted --step 2000 -k 2500 --lr 0.001 > results/lr/p14_targeted_lr_001.txt &

# echo ".....((...((....))...((....))...((....))...((....))...((....))...((....))...))...................." | ./main --obj pyx_direct --init uniform --step 2000 -k 2500 --lr 0.01 > results/lr/p14_uniform_lr_010.txt &
# echo ".....((...((....))...((....))...((....))...((....))...((....))...((....))...))...................." | ./main --obj pyx_direct --init uniform --step 2000 -k 2500 --lr 0.005 > results/lr/p14_uniform_lr_005.txt &
# echo ".....((...((....))...((....))...((....))...((....))...((....))...((....))...))...................." | ./main --obj pyx_direct --init uniform --step 2000 -k 2500 --lr 0.001 > results/lr/p14_uniform_lr_001.txt &