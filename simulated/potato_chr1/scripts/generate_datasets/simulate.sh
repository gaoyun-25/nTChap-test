#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate repeat
benchmarkingParametersFile="/home/gaoyun/poly/code/phasing/final_nTChap/process_data/simulated/potato_chr1/scripts/benchmarks/benchmarkingParameters_potato_chr1.tsv"
my_code="/home/gaoyun/poly/code/phasing/final_nTChap/process_data/simulated/potato_chr1/scripts/generate_datasets/"

python $my_code"simulate.py" $benchmarkingParametersFile $my_code

