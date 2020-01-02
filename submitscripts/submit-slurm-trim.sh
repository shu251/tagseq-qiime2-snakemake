#!/bin/bash

snakemake   \
        --jobs 100 -s Snakefile-otu \
        --use-conda --until multiqc --cluster-config submitscripts/cluster.yaml --cluster "sbatch --parsable --partition={cluster.queue} --job-name=otu.{rule}.{wildcards} --mem={cluster.mem}gb --time={cluster.time} --ntasks={cluster.threads} --nodes={cluster.nodes}"

