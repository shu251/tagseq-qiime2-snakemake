# Snakemake file to QC tag-sequence reads
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

#----DEFINE RULES----#

rule all:
  input:
    html = expand("{scratch}/qc/fastqc/{sample}_{suf}_{num}_fastqc.html", scratch = SCRATCH, sample = SAMPLE_LIST, suf = SAMPLE_SUFFIX, num = [1,2]),

rule fastqc:
  input:    
    expand("{rawfastq}/{sample}_{suf}_{num}.fastq.gz", rawfastq = INPUTDIR, scratch = SCRATCH, sample = SAMPLE_LIST, suf = SAMPLE_SUFFIX, num = [1,2])
  output:
    html = expand("{scratch}/qc/fastqc/{sample}_{suf}_{num}_fastqc.html", scratch = SCRATCH, sample = SAMPLE_LIST, suf = SAMPLE_SUFFIX, num = [1,2]),
    zip = expand("{scratch}/qc/fastqc/{sample}_{suf}_{num}_fastqc.zip", scratch = SCRATCH, sample = SAMPLE_LIST, suf = SAMPLE_SUFFIX, num = [1,2])
  params: ""
  log:
    expand("{scratch}/qc/fastqc/logs/{sample}_{suf}_{num}_fastqc.log", scratch = SCRATCH, sample = SAMPLE_LIST, suf = SAMPLE_SUFFIX, num = [1,2])
  wrapper:
    "0.35.1/bio/fastqc"

