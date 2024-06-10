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
    # echo "${puzzles[1]}" | ./main --init uniform_sm -k 2500 --step 2000 --lr 0.001 > results/sampling_pyx_uniform_sm/${puzzles[0]}.txt &
    # echo "${puzzles[1]}" | ./main --init uniform_sm --lr 0.01 --softmax > results/sampling_pyx_uniform_sm_softmax/${puzzles[0]}.txt &
    # echo "${puzzles[1]}" | ./main --init uniform_sm --lr 0.0001 --softmax --adam > results/sampling_pyx_uniform_sm_softmax_adam_0001/${puzzles[0]}.txt &
    # echo "${puzzles[1]}" | ./main --init uniform_sm --lr 0.001 --softmax --adam > results/sampling_pyx_uniform_sm_softmax_adam/${puzzles[0]}.txt &
    # echo "${puzzles[1]}" | ./main --init uniform_sm --lr 0.1 --softmax > results/sampling_pyx_uniform_sm_softmax_1/${puzzles[0]}.txt &
    
    # echo "${puzzles[1]}" | ./main --init targeted_sm --eps 0.50 --lr 0.01 --softmax > results/sampling_pyx_targeted_sm_softmax_eps_050/${puzzles[0]}.txt &
    # echo "${puzzles[1]}" | ./main --init targeted_sm --eps 0.75 --lr 0.01 --softmax > results/sampling_pyx_targeted_sm_softmax_eps_075/${puzzles[0]}.txt &
    # echo "${puzzles[1]}" | ./main --init targeted_sm --eps 0.90 --lr 0.01 --softmax > results/sampling_pyx_targeted_sm_softmax_eps_090/${puzzles[0]}.txt &
    
    # echo "${puzzles[1]}" | ./main --init targeted_sm --eps 0.75 --lr 0.01 --softmax --step 2000 > results/sampling_pyx_targeted_sm_softmax_eps_075_no_adam/${puzzles[0]}.txt &
    # echo "${puzzles[1]}" | ./main --init targeted_sm --eps 0.75 --lr 0.01 --softmax --adam > results/sampling_pyx_targeted_sm_softmax_eps_075_adam/${puzzles[0]}.txt &
    
    # echo "${puzzles[1]}" | ./main --init targeted_sm --eps 0.75 --lr 0.1 --softmax --adam --step 2000 > results/sampling_pyx_targeted_sm_softmax_eps_075_adam_lr_1/${puzzles[0]}.txt &
    # echo "${puzzles[1]}" | ./main --init targeted_sm --eps 0.75 --lr 0.001 --softmax --adam --step 2000 > results/sampling_pyx_targeted_sm_softmax_eps_075_adam_lr_001/${puzzles[0]}.txt &
    
    # echo "${puzzles[1]}" | ./main --test --init targeted_sm --eps 0.75 --lr 0.01 --softmax --adam --step 2000 > results/sampling_pyx_targeted_sm_softmax_eps_075_adam_lr_decay/49_lr_01_no_decay.txt &
    # echo "${puzzles[1]}" | ./main --test --init targeted_sm --eps 0.75 --lr 0.01 --softmax --adam --step 2000 --lr_decay --lr_decay_rate 0.96 --staircase 200 > results/sampling_pyx_targeted_sm_softmax_eps_075_adam_lr_decay/49_lr_01_decay_96_stair_200.txt &
    # echo "${puzzles[1]}" | ./main --test --init targeted_sm --eps 0.75 --lr 0.01 --softmax --adam --step 2000 --lr_decay --lr_decay_rate 0.92 --staircase 200 > results/sampling_pyx_targeted_sm_softmax_eps_075_adam_lr_decay/49_lr_01_decay_92_stair_200.txt &

    # echo "${puzzles[1]}" | ./main --test --init targeted_sm --eps 0.75 --lr 0.1 --softmax --adam --step 2000 > results/sampling_pyx_targeted_sm_softmax_eps_075_adam_lr_decay/49_lr_1_no_decay.txt &
    # echo "${puzzles[1]}" | ./main --test --init targeted_sm --eps 0.75 --lr 0.1 --softmax --adam --step 2000 --lr_decay --lr_decay_rate 0.96 --staircase 200 > results/sampling_pyx_targeted_sm_softmax_eps_075_adam_lr_decay/49_lr_1_decay_96_stair_200.txt &
    # echo "${puzzles[1]}" | ./main --test --init targeted_sm --eps 0.75 --lr 0.1 --softmax --adam --step 2000 --lr_decay --lr_decay_rate 0.92 --staircase 200 > results/sampling_pyx_targeted_sm_softmax_eps_075_adam_lr_decay/49_lr_1_decay_92_stair_200.txt &

    # echo "${puzzles[1]}" | ./main --init targeted_sm --eps 0.05 --lr 0.1 --softmax > results/sampling_pyx_uniform_sm_softmax_eps_05/${puzzles[0]}.txt &
    # echo "${puzzles[1]}" | ./main --init targeted_sm --eps 0.10 --lr 0.1 --softmax > results/sampling_pyx_uniform_sm_softmax_eps_10/${puzzles[0]}.txt &
    # echo "${puzzles[1]}" | ./main --init targeted_sm --eps 0.20 --lr 0.1 --softmax > results/sampling_pyx_uniform_sm_softmax_eps_20/${puzzles[0]}.txt &
    # echo "${puzzles[1]}" | ./main --init targeted_sm > results/sampling_pyx_targeted_sm/${puzzles[0]}.txt &

    # echo "${puzzles[1]}" | ./main --init targeted_sm --eps 0.75 --lr 0.01 --softmax --adam --lr_decay --lr_decay_rate 0.96 --staircase --lr_decay_step 100 > results/sampling_pyx_targeted_sm_softmax_eps_075_adam_lr_decay/${puzzles[0]}.txt &
    # echo "${puzzles[1]}" | ./main --init uniform_sm --lr 0.01 --softmax --adam --lr_decay --lr_decay_rate 0.96 --staircase --lr_decay_step 50 > results/sampling_pyx_uniform_softmax_adam_lr_decay/${puzzles[0]}.txt &
    # echo "${puzzles[1]}" | ./main --init uniform_sm --lr 0.01 -k 5000 --softmax --adam > results/sampling_pyx_uniform_sm_softmax_adam/${puzzles[0]}.txt &
    
    # Time test
    # echo "${puzzles[1]}" | ./main --test --init uniform_sm --lr 0.01 -k 2500 --softmax --adam --step 5 > results/sampling_pyx_uniform_sm_softmax_adam_time/${puzzles[0]}.txt
    # echo "${puzzles[1]}" | ./main --test --init targeted_sm --lr 0.01 -k 2500 --softmax --adam --step 5 > results/sampling_pyx_targeted_sm_softmax_adam_time/${puzzles[0]}.txt
    # echo "${puzzles[1]}" | ./main --init uniform_sm --lr 0.01 -k 2500 --softmax --adam --step 5 > results/sampling_pyx_uniform_sm_softmax_adam_time_parallel/${puzzles[0]}.txt
    # echo "${puzzles[1]}" | ./main --init targeted_sm --eps 0.75 --lr 0.01 -k 2500 --softmax --adam --step 5 > results/sampling_pyx_targeted_sm_softmax_adam_time_parallel/${puzzles[0]}.txt
    # echo "${puzzles[1]}" | ./main --init targeted_sm --eps 0.75 --lr 0.01 -k 2500 --softmax --adam --step 2000 > results/sampling_pyx_targeted_sm_softmax_adam_time_parallel_cache/${puzzles[0]}.txt
    # echo "${puzzles[1]}" | ./main --init uniform_sm --lr 0.01 -k 2500 --softmax --adam --step 2000 > results/sampling_pyx_uniform_sm_softmax_adam_time_parallel_cache/${puzzles[0]}.txt
    
    # Nesterov
    # echo "${puzzles[1]}" | ./main --init targeted_sm --nesterov > results/sampling_pyx_targeted_sm_nesterov/${puzzles[0]}.txt &
    # echo "${puzzles[1]}" | ./main --init uniform_sm --nesterov > results/sampling_pyx_uniform_sm_nesterov/${puzzles[0]}.txt &
    # echo "${puzzles[1]}" | ./main --init targeted_sm --nesterov --step 2000 > results/sampling_pyx_targeted_sm_nesterov/${puzzles[0]}.txt &
    # echo "${puzzles[1]}" | ./main --init targeted_sm --step 2000 > results/sampling_pyx_targeted_sm_no_nesterov/${puzzles[0]}.txt &
    
    # Nesterov + lr_decay
    # echo "${puzzles[1]}" | ./main --init targeted_sm --nesterov --lr_decay --lr_decay_rate 0.96 --staircase --lr_decay_step 50 > results/sampling_pyx_targeted_sm_nesterov_lr_decay/${puzzles[0]}.txt &
    # echo "${puzzles[1]}" | ./main --init targeted_sm --nesterov --lr_decay --lr_decay_rate 0.92 --staircase --lr_decay_step 50 > results/sampling_pyx_targeted_sm_nesterov_lr_decay/${puzzles[0]}.txt &

    # Nesterov + adaptive lr_decay
    # echo "${puzzles[1]}" | ./main --init targeted_sm --nesterov --lr_decay --lr_decay_rate 0.92 > results/sampling_pyx_targeted_sm_nesterov_adaptive_lr_decay_1/${puzzles[0]}.txt &
    # echo "${puzzles[1]}" | ./main --init targeted_sm --nesterov --lr_decay --lr_decay_rate 0.50 > results/sampling_pyx_targeted_sm_nesterov_adaptive_lr_decay_50/${puzzles[0]}.txt &
    # Puzzle 74
    # echo "${puzzles[1]}" | ./main --init targeted_sm --initial_lr 0.0001 --nesterov --lr_decay --lr_decay_rate 0.50 > results/sampling_pyx_targeted_sm_nesterov_adaptive_lr_decay_50/${puzzles[0]}.txt &

    # Kmers Analysis
    # echo "${puzzles[1]}" | ./main --test --init targeted_sm --eps 0.75 --lr 0.01 --softmax --adam --step 2000 > results/sampling_pyx_targeted_sm_softmax_eps_075_adam_kmers/${puzzles[0]}.txt

    # NED
    # echo "${puzzles[1]}" | ./main --obj ned --init targeted_sm --eps 0.75 --lr 0.01 --softmax --adam > results/sampling_ned_targeted_sm_softmax_eps_075_adam/${puzzles[0]}.txt &
    # echo "${puzzles[1]}" | ./main --obj log_ned --init targeted_sm --eps 0.75 --lr 0.01 --softmax --adam > results/sampling_log_ned_targeted_sm_softmax_eps_075_adam/${puzzles[0]}.txt &
    # echo "${puzzles[1]}" | ./main --obj log_ned --init uniform_sm --lr 0.01 --softmax --adam > results/sampling_log_ned_uniform_sm_softmax_adam/${puzzles[0]}.txt &
    # echo "${puzzles[1]}" | ./main --obj ned --init targeted_sm --eps 0.75 --lr 0.005 --softmax --adam > results/sampling_ned_targeted_sm_softmax_eps_075_adam_lr_005/${puzzles[0]}.txt &
    # echo "${puzzles[1]}" | ./main --obj log_ned --init targeted_sm --eps 0.75 --lr 0.005 --softmax --adam > results/sampling_log_ned_targeted_sm_softmax_eps_075_adam_lr_005/${puzzles[0]}.txt &

    # main log NED script
    # echo "${puzzles[1]}" | ./main --obj log_ned --init targeted_sm --eps 0.75 > results/sampling_log_ned_targeted_sm_softmax_eps_075_adam/${puzzles[0]}.txt &
    # echo "${puzzles[1]}" | ./main --obj log_ned --init targeted_sm --eps 0.75 --is_lazy -b 150 > results/sampling_log_ned_targeted_sm_softmax_eps_075_adam/${puzzles[0]}.txt &
    # echo "${puzzles[1]}" | ./main --obj log_ned --init uniform_sm > results/sampling_log_ned_uniform_sm_softmax_adam/${puzzles[0]}.txt &
    # use sampling_log_ned_targeted_sm_softmax_eps_075_adam_copy for shorter puzzles

    # structural distance
    # echo "${puzzles[1]}" | ./main --obj dist --init uniform_sm > results/sampling_dist_uniform_sm_softmax_adam/${puzzles[0]}.txt &
    # echo "${puzzles[1]}" | ./main --obj dist --init uniform_sm -k 1000 > results/sampling_dist_uniform_sm_softmax_adam/${puzzles[0]}.txt & 
    
    # nemo objective
    # echo "${puzzles[1]}" | ./main --obj dist --init uniform_sm > results/sampling_dist_nemo_uniform_sm_softmax_adam/${puzzles[0]}.txt &
    # echo "${puzzles[1]}" | ./main --obj dist --init targeted_sm --eps 0.75 > results/sampling_dist_nemo_targeted_sm_eps_075_softmax_adam/${puzzles[0]}.txt &
    
    # free energy diff
    # echo "${puzzles[1]}" | ./main --obj ediff --init targeted_sm --eps 0.75 > results/sampling_ediff_nemo_targeted_sm_eps_075_softmax_adam/${puzzles[0]}.txt &
    # echo "${puzzles[1]}" | ./main --obj ediff --init targeted_sm --eps 0.75 -b 100 > results/sampling_ediff_nemo_targeted_sm_eps_075_softmax_adam/${puzzles[0]}.txt &

    # small beamsize
    # echo "${puzzles[1]}" | ./main --init targeted_sm --eps 0.75 --lr 0.01 --softmax --adam -b 10 --best_k 10 --step 50 > results/small_beam_size/${puzzles[0]}.txt &
    
    # lazy vs non lazy mode
    # echo "${puzzles[1]}" | ./main --obj log_ned --init uniform_sm --step 5 > results/sampling_log_ned_uniform_sm_softmax_adam_time/${puzzles[0]}.txt
    # echo "${puzzles[1]}" | ./main --obj log_ned --init uniform_sm --step 5 --is_lazy  > results/sampling_log_ned_uniform_sm_softmax_adam_lazy_time/${puzzles[0]}.txt

    # ablation study
    # echo "${puzzles[1]}" | ./main --init targeted --eps 0.75 --lr 0.01 --softmax --adam > results/sampling_pyx_targeted_softmax_eps_075_adam/${puzzles[0]}.txt &
    # echo "${puzzles[1]}" | ./main --init targeted_sm --eps 0.75 --lr 0.01 --no_trimismatch > results/sampling_pyx_targeted_softmax_eps_075_adam_no_trimismatch/${puzzles[0]}.txt &

    # composite objective function
    echo "${puzzles[1]}" | ./main --test --obj comp > results/sampling_comp_targeted_sm_eps_075_softmax_adam/${puzzles[0]}.txt &
done < "data/eterna/$1.txt"

