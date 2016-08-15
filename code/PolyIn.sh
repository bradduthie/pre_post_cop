#!/bin/bash
# SGE SUBMISSION SCRIPT 
# This is an embarrassingly parallel simulation
# Submit to 1 core -- should print to one output file
# I estimate that this program will take less than 24 hrs to run
#
#$ -V -cwd -l h_rt=24:00:00 -o out.$JOB_ID
./PolyIn
