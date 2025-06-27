#!/bin/bash

# # The pipeline of nPhase to generate Saccharomyces cerevisiae datasets
# echo "download data"
# python3.8 downloadData.py benchmarkingParameters.tsv && echo "success!"

# echo "generate groundTruth"
# # every long reads map;every short reads map, then call_variants
# python3.8 processReads.py benchmarkingParameters.tsv groundTruth && echo "success!"

# echo "generates virtual hybrids"
# # given coverage downsample every fastq and then merge(short; long reads)
# python3.8 hybridGenerator.py benchmarkingParameters.tsv && echo "success!"

# echo "Maps short reads and long reads, variants calls and subsets short read VCFs, with and without indels"
# python3.8 processReads.py benchmarkingParameters.tsv virtualHybrids && echo "success!"
