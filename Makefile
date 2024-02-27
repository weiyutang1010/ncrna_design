################################
# Makefile
#
# author: He Zhang
# edited by: 03/2019
################################

CC=g++
DEPS=src/Outside.cpp src/Inside.cpp src/main.h src/Utils/energy_parameter.h src/Utils/feature_weight.h src/Utils/intl11.h src/Utils/intl21.h src/Utils/intl22.h src/Utils/utility_v.h src/Utils/utility.h src/Sampling.cpp src/LinearPartition.cpp src/LinearPartition.h
CFLAGS=-std=c++17 -O3 -fopenmp
.PHONY : clean main
objects=bin/main

main: src/main.cpp $(DEPS) 
		chmod +x main
		mkdir -p bin
		$(CC) src/main.cpp $(CFLAGS) -Dlpv -DSPECIAL_HP -o bin/main

clean:
	-rm $(objects)