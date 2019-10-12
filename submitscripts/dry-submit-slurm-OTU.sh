#!/bin/bash

snakemake -np \
        --jobs 100 -s Snakefile-otu \
        --use-conda --cluster-config submitscripts/cluster.yaml --cluster "sbatch --parsable --qos=unlim --partition={cluster.queue} --job-name=tagseq.{rule}.{wildcards} --mem={cluster.mem}gb --time={cluster.time} --ntasks={cluster.threads} --nodes={cluster.nodes}"

