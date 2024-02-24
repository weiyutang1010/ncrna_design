#!/bin/bash

echo "(((...)))" | ./main --mode ncrna_design --objective pyx_sampling --init targeted --verbose -k 10 --step 5 > result.txt &