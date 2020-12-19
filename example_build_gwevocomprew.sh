#!/bin/bash

g++ -o gwevo26r.run gwintrew.cpp /phenod/data/pniksa/mypref/lib/libgsl.a /phenod/data/pniksa/mypref/lib/libgslcblas.a \
-I/phenod/data/pniksa/mypref/include -lm -I/phenod/data/pniksa/cubature-1.0.2 -I/phenod/data/pniksa/boost_1_61_0/ \
-L/phenod/data/pniksa/boost_1_61_0/libs -std=c++11 -O1 -fopenmp \
-lpthread -g -"ulimit -c unlimited"
set OMP_NUM_THREADS=64
./gwevo26r.run
