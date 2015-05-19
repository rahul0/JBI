#!/bin/bash

rm *.o
rm core*
rm gibbs_linked_underlining 
g++ -O3    -I$gsl_home/include -I$m_inc  -c gibbs_linked_underlining.cpp -fopenmp
g++ -fopenmp -O3   -o gibbs_linked_underlining gibbs_linked_underlining.o  -lgsl -lgslcblas  -lm -L$gsl_home/lib -L$m_lib $m_lib/libMinuit2.a
