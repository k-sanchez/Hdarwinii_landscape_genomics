#!/bin/bash

# Absolute path to the stacks 2 protocol directory
work=~/OneDrive/PRODUCCION/ongoing/Hdar_posdoc/data/RADseq/stacks

# Loop over the M values
# iterates over every value between 1 and 6
for M in {1..6}; do
    # create a new output directory per Stacks run
    out=$work/param_opt_all/denovo_M${M}
    mkdir -p $out
    # move into the new directory
    cd $out
    # stacks command
    denovo_map.pl \
        --samples $work/process_radtags/ \
        --popmap $work/info/paramopt_subset_all.tsv \
        -o $out \
        -M $M \
        -n $M \
        -T 6 \
        --min-samples-per-pop 0.8 \
        -X "ustacks:--force-diff-len"
done
