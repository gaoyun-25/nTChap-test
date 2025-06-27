#!/bin/bash
# set -x

eval "$(conda shell.bash hook)"
conda activate Haplosim_py2
haplogenerator.py -f $1 -o $2 --model poisson -s $3 -p $4 -v -m $5 --dosage $6
# echo "-f $1 -o $2 --model poisson  -s $3 -p $4 -v -m $5 --dosage $6"
# haplogenerator.py -h 
