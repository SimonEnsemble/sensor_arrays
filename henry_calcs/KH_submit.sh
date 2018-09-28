#!/bin/bash

# use current working directory for input and output
# default is to use the users home directory
#$ -cwd

# name this job
#$ -N WONKOB_iqeqcharged.cif_C2H6_250

#$ -pe thread 4 # use 4 threads/cores

# send stdout and stderror to this file
#$ -o jobz/WONKOB_iqeqcharged.cif_C2H6.o
#$ -e jobz/WONKOB_iqeqcharged.cif_C2H6.e

# select queue - if needed; mime5 is SimonEnsemble priority queue but is restrictive.
##$ -q mime5

# print date and time
date
julia -p 4 run_henry.jl C2H6 WONKOB_iqeqcharged.cif 298.000000 UFF.csv 250
