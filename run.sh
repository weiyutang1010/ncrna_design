#!/bin/bash

echo "((((...))))." | ./main --mode ncrna_design --objective pyx_sampling --init uniform --verbose --lr 0.01 -k 1000 --step 2500 > results/temp/2.txt &
echo "((((...))))." | ./main --mode ncrna_design --objective pyx_sampling --init targeted --verbose --lr 0.005 -k 1000 --step 2500 > results/temp/3.txt &
# echo "(((...)))" | ./main --mode ncrna_design --objective pyx_sampling --init targeted --verbose -k 1000 --step 5 > result.txt &
# echo "(((((((........)))))))" | ./main --mode ncrna_design --objective pyx_sampling --init uniform --verbose -k 10000 --step 5 > result.txt &