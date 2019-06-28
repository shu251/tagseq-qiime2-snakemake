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
    # Trimmed read output
    r1 = expand("{scratch}/trimmed/{sample}_1.fastq.gz", scratch = SCRATCH, sample = SAMPLE_LIST),
    r2 = expand("{scratch}/trimmed/{sample}_2.fastq.gz", scratch = SCRATCH, sample = SAMPLE_LIST),
    # reads where trimming entirely removed the mate
    r1_unpaired = expand("{scratch}/trimmed/{sample}_1.unpaired.fastq.gz", scratch = SCRATCH, sample = SAMPLE_LIST),
    r2_unpaired = expand("{scratch}/trimmed/{sample}_2.unpaired.fastq.gz", scratch = SCRATCH, sample = SAMPLE_LIST),

rule fastqc:
  input:    
    expand("{rawfastq}/{sample}_{suf}_{num}.fastq.gz", rawfastq = INPUTDIR, sample = SAMPLE_LIST, suf = SAMPLE_SUFFIX, num = [1,2])
  output:
    html = expand("{scratch}/qc/fastqc/{sample}_{suf}_{num}_fastqc.html", scratch = SCRATCH, sample = SAMPLE_LIST, suf = SAMPLE_SUFFIX, num = [1,2]),
    zip = expand("{scratch}/qc/fastqc/{sample}_{suf}_{num}_fastqc.zip", scratch = SCRATCH, sample = SAMPLE_LIST, suf = SAMPLE_SUFFIX, num = [1,2])
  params: ""
  log:
    expand("{scratch}/qc/fastqc/logs/{sample}_{suf}_{num}_fastqc.log", scratch = SCRATCH, sample = SAMPLE_LIST, suf = SAMPLE_SUFFIX, num = [1,2])
  wrapper:
    "0.35.1/bio/fastqc"
# Error in how snakemake finds this env?

rule trimmomatic_pe:
  input:
    r1 = expand("{rawfastq}/{sample}_{suf}_1.fastq.gz", rawfastq = INPUTDIR, sample = SAMPLE_LIST, suf = SAMPLE_SUFFIX),
    r2 = expand("{rawfastq}/{sample}_{suf}_2.fastq.gz", rawfastq = INPUTDIR, sample = SAMPLE_LIST, suf = SAMPLE_SUFFIX)
  output:
    r1 = "{scratch}/trimmed/{sample}_1.fastq.gz",
    r2 = "{scratch}/trimmed/{sample}_2.fastq.gz",
    # reads where trimming entirely removed the mate
    r1_unpaired = "{scratch}/trimmed/{sample}_1.unpaired.fastq.gz",
    r2_unpaired = "{scratch}/trimmed/{sample}_2.unpaired.fastq.gz"
  log:
   "{scratch}/trimmed/logs/trimmomatic/{sample}.log"
  params:
    trimmer = ["LEADING:2", "TRAILING:2", "SLIDINGWINDOW:4:2", "MINLEN:25"],
    extra = ""
  wrapper:
   "0.35.1/bio/trimmomatic/pe"
