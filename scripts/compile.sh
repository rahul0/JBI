#!/bin/bash

rm *.o
rm core*
rm 1a
g++ -O3  -I$gsl_home/include -I$m_inc  -c 1a.cpp -fopenmp
g++ -fopenmp -O3 -o 1a 1a.o  -lgsl -lgslcblas  -lm -L$gsl_home/lib -L$m_lib $m_lib/libMinuit2.a
