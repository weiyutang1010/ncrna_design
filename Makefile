################################
# Makefile
#
# author: He Zhang
# edited by: 03/2019
################################

CC=g++
DEPS=src/Outside.cpp src/ExpectedPartition.h src/Utils/energy_parameter.h src/Utils/feature_weight.h src/Utils/intl11.h src/Utils/intl21.h src/Utils/intl22.h src/Utils/utility_v.h src/Utils/utility.h
CFLAGS=-std=c++11 -O3
.PHONY : clean expectedpartition
objects=bin/expectedpartition

expectedpartition: src/ExpectedPartition.cpp $(DEPS) 
		chmod +x expectedpartition
		mkdir -p bin
		$(CC) src/ExpectedPartition.cpp $(CFLAGS) -Dlpv -o bin/expectedpartition

clean:
	-rm $(objects)