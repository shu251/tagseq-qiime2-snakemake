# Snakemake file to generate ASVs from trimmed reads
configfile: "config.yaml"

import io 
import os
import pandas as pd
import pathlib
from snakemake.exceptions import print_exception, WorkflowError

#----SET VARIABLES----#

PROJ = config["proj_name"]
INPUTDIR = config["raw_data"]
SCRATCH = config["scratch"]

SAMPLE_TABLE = pd.read_table(config["input_seq"])
SAMPLE_LIST = list(SAMPLE_TABLE.Run)
SAMPLE_NAME = list(SAMPLE_TABLE.SampleName)

SAMPLE_SUFFIX = config["sample_suffix"]
OUTPUTDIR = config["outputDIR"]

# ASV specific:
MANIFEST = config["manifest"]
PRIMER_F = config["primerF"]
PRIMER_R = config["primerR"]


#----DEFINE RULES----#

rule all:
  input:
    import_qza = expand("{base}/asv/{project}-PE-demux.qza", base = SCRATCH, project = PROJ),
    rm_primers = expand("{base}/asv/{project}-PE-demux-noprimer.qza", base = SCRATCH, project = PROJ),
    table = expand("{base}/asv/{project}-asv-table.qza", base = SCRATCH, project = PROJ),
    rep = expand("{base}/asv/{project}-rep-seqs.qza", base = SCRATCH, project = PROJ),
    stats = expand("{base}/asv/{project}-stats-dada2.qza", base = SCRATCH, project = PROJ),

rule import_qiime:
  input:
    manifest = MANIFEST
  output:
    SCRATCH + "/asv/{project}-PE-demux.qza"
  log:
    SCRATCH + "/logs/qiime2/{project}_q2.log"
  shell:
   "qiime tools import \
       --type 'SampleData[PairedEndSequencesWithQuality]' \
       --input-path {input.manifest} \
       --output-path {output} \
       --input-format PairedEndFastqManifestPhred33"
#
# As of 6/30/2019 only using raw input fastq - need to change to trimmed reads as input for manifest 
#

rule rm_primers:
  input:
    raw = SCRATCH + "/asv/{project}-PE-demux.qza"
  output:
    SCRATCH + "/asv/{project}-PE-demux-noprimer.qza"
  log:
    SCRATCH + "/logs/qiime2/{project}_primer_q2.log"
  shell:
    "qiime cutadapt trim-paired \
       --i-demultiplexed-sequences {input.raw} \
       --p-front-f {config[primerF]} \
       --p-front-r {config[primerR]} \
       --p-error-rate {config[primer_err]} \
       --p-overlap {config[primer_overlap]} \
       --o-trimmed-sequences {output}"

#
# I'd like to add a move to main directory step here, rather than working from scratch? advisable?
#

rule dada2:
  input:
    SCRATCH + "/asv/{project}-PE-demux-noprimer.qza"
  output:
    table = SCRATCH + "/asv/{project}-asv-table.qza",
    rep = SCRATCH + "/asv/{project}-rep-seqs.qza",
    stats = SCRATCH + "/asv/{project}-stats-dada2.qza"
  log:
    SCRATCH + "/logs/qiime2/{project}_dada2_q2.log"
  shell:
    "qiime dada2 denoise-paired \
  	--i-demultiplexed-seqs {input} \
	--p-trunc-q {config[truncation_err]} \
	--p-trunc-len-f {config[truncation_len-f]} \
	--p-trunc-len-r {config[truncation_len-r]} \
	--p-max-ee {config[quality_err]} \
	--p-n-reads-learn {config[training]} \
	--p-chimera-method {config[chimera]} \
	--o-table {output.table} \
	--o-representative-sequences {output.rep} \
	--o-denoising-stats {output.stats}"
