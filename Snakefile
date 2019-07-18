# Snakemake file - input raw fastq reads to generate Amplicon Sequence Variants
# SHu 07-2019
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
OUTPUTDIR = config["outputDIR"]

# Use glob statement to find all samples in 'raw_data' directory
SAMPLE_LIST,NUMS = glob_wildcards("INPUTDIR/{sample}_{num}.fastq.gz")
# Unique the output variables from glob_wildcards
SAMPLE_SET = set(SAMPLE_LIST)
SET_NUMS = set(NUMS)

# Create final manifest file for qiime2
MANIFEST = pd.read_csv(config["manifest"]) #Import manifest
MANIFEST['filename'] = MANIFEST['absolute-filepath'].str.split('/').str[-1] #add new column with only file name
MANIFEST.rename(columns = {'absolute-filepath':'rawreads-filepath'}, inplace = True)
PATH_TRIMMED = "trimmed" # name of directory with trimmed reads
NEWPATH = os.path.join(SCRATCH, PATH_TRIMMED)
MANIFEST['filename'] = MANIFEST['filename'].str.replace(".fastq.gz", "_trim.fastq.gz")
MANIFEST['absolute-filepath'] = NEWPATH+ "/" + MANIFEST['filename']    
MANIFEST[['sample-id','absolute-filepath','direction']].set_index('sample-id').to_csv('manifest-trimmed.txt')
MANIFEST_FINAL = config["manifest-trimmed"]

# Database information to assign taxonomy
DB_classifier = config["database"]

# All qiime2 artifact files
ARTIFACT = glob_wildcards("SCRATCH/qiime2/asv/PROJ-{arti}.qza")

#----DEFINE RULES----#

rule all:
  input:
    # fastqc output before trimming
    html = expand("{scratch}/fastqc/{sample}_{num}_fastqc.html", scratch = SCRATCH, sample=SAMPLE_SET, num=SET_NUMS),
    zip = expand("{scratch}/fastqc/{sample}_{num}_fastqc.zip", scratch = SCRATCH, sample=SAMPLE_SET, num=SET_NUMS),
    orig_html = SCRATCH + "/fastqc/raw_multiqc.html",
    orig_stats = SCRATCH + "/fastqc/raw_multiqc_general_stats.txt",
    # Trimmed data output
    trimmedData = expand("{scratch}/trimmed/{sample}_{num}_trim.fastq.gz", scratch= SCRATCH, sample=SAMPLE_SET, num=SET_NUMS), 
    html_trim = expand("{scratch}/fastqc/{sample}_{num}_trimmed_fastqc.html", scratch= SCRATCH, sample=SAMPLE_SET, num=SET_NUMS),
    zip_trim = expand("{scratch}/fastqc/{sample}_{num}_trimmed_fastqc.zip", scratch= SCRATCH, sample=SAMPLE_SET, num=SET_NUMS),
    trim_html = SCRATCH + "/fastqc/trimmed_multiqc.html", #next change to include proj name
    trim_stats = SCRATCH + "/fastqc/trimmed_multiqc_general_stats.txt",
    # QIIME2 outputs
    q2_import = SCRATCH + "/qiime2/asv/" + PROJ + "-PE-demux.qza",
    q2_primerRM = SCRATCH + "/qiime2/asv/" + PROJ + "-PE-demux-noprimer.qza",
    table = SCRATCH + "/qiime2/asv/" + PROJ + "-asv-table.qza",
    rep = SCRATCH + "/qiime2/asv/" + PROJ + "-rep-seqs.qza",
    stats = SCRATCH + "/qiime2/asv/" + PROJ + "-stats-dada2.qza",
    sklearn = SCRATCH + "/qiime2/asv/" + PROJ +	"-tax_sklearn.qza",
    biom = SCRATCH + "/qiime2/asv/table/feature-table.biom",
    table_tsv = SCRATCH + "/qiime2/asv/" + PROJ + "-asv-table.tsv",
    table_tax = SCRATCH + "/qiime2/asv/tax_dir/taxonomy.tsv",
    # q2 visualization outputs
    # qzv = expand("{scratch}/qiime2/asv/viz/{proj}-{arti}.qzv", scratch = SCRATCH, proj = PROJ, arti = ARTIFACT)

rule fastqc:
  input:    
    INPUTDIR + "/{sample}_{num}.fastq.gz"
  output:
    html = SCRATCH + "/fastqc/{sample}_{num}_fastqc.html",
    zip = SCRATCH + "/fastqc/{sample}_{num}_fastqc.zip"
  params: ""
  log:
    SCRATCH + "/logs/fastqc/{sample}_{num}.log"
  wrapper:
    "0.35.2/bio/fastqc"

rule trimmomatic_pe:
  input:
    r1 = INPUTDIR + "/{sample}_1.fastq.gz",
    r2 = INPUTDIR + "/{sample}_2.fastq.gz"
  output:
    r1 = SCRATCH + "/trimmed/{sample}_1_trim.fastq.gz",
    r2 = SCRATCH + "/trimmed/{sample}_2_trim.fastq.gz",
    # reads where trimming entirely removed the mate
    r1_unpaired = SCRATCH + "/trimmed/{sample}_1.unpaired.fastq.gz",
    r2_unpaired = SCRATCH + "/trimmed/{sample}_2.unpaired.fastq.gz"
  log:
    SCRATCH + "/trimmed/logs/trimmomatic/{sample}.log"
  params:
    trimmer = ["LEADING:2", "TRAILING:2", "SLIDINGWINDOW:4:2", "MINLEN:25"],
    extra = ""
  wrapper:
    "0.35.2/bio/trimmomatic/pe"

rule fastqc_trim:
  input:
    SCRATCH + "/trimmed/{sample}_{num}_trim.fastq.gz"
  output:
    html = SCRATCH + "/fastqc/{sample}_{num}_trimmed_fastqc.html",
    zip = SCRATCH + "/fastqc/{sample}_{num}_trimmed_fastqc.zip"
  params: ""
  log:
    SCRATCH + "/logs/fastqc/{sample}_{num}_trimmed.log"
  wrapper:
    "0.35.2/bio/fastqc"

rule multiqc:
  input:
    orig = expand("{scratch}/fastqc/{sample}_{num}_fastqc.zip", scratch= SCRATCH, sample=SAMPLE_SET, num=SET_NUMS),
    trimmed = expand("{scratch}/fastqc/{sample}_{num}_trimmed_fastqc.zip", scratch= SCRATCH, sample=SAMPLE_SET, num=SET_NUMS)
  output:
    orig_html = SCRATCH + "/fastqc/raw_multiqc.html", 
    orig_stats = SCRATCH + "/fastqc/raw_multiqc_general_stats.txt",
    trim_html = SCRATCH + "/fastqc/trimmed_multiqc.html", 
    trim_stats = SCRATCH + "/fastqc/trimmed_multiqc_general_stats.txt"
  conda:
   "envs/multiqc-env.yaml"
  shell: 
    """
    multiqc -n multiqc.html {input.orig} #run multiqc
    mv multiqc.html {output.orig_html} #rename html
    mv multiqc_data/multiqc_general_stats.txt {output.orig_stats} #move and rename stats
    rm -rf multiqc_data #clean-up
    #repeat for trimmed data
    multiqc -n multiqc.html {input.trimmed} #run multiqc
    mv multiqc.html {output.trim_html} #rename html
    mv multiqc_data/multiqc_general_stats.txt {output.trim_stats} #move and rename stats
    rm -rf multiqc_data	#clean-up
    """ 

rule import_qiime:
  input:
    MANIFEST_FINAL
  output:
    q2_import = SCRATCH + "/qiime2/asv/" + PROJ + "-PE-demux.qza"
  log:
    SCRATCH + "/qiime2/logs/" + PROJ + "_q2.log"
  conda:
    "envs/qiime2-2019.4.yaml"
  shell:
   "qiime tools import \
       --type 'SampleData[PairedEndSequencesWithQuality]' \
       --input-path {input} \
       --output-path {output.q2_import} \
       --input-format PairedEndFastqManifestPhred33"

rule rm_primers:
  input:
    q2_import = SCRATCH + "/qiime2/asv/" + PROJ + "-PE-demux.qza"
  output:
    q2_primerRM = SCRATCH + "/qiime2/asv/" + PROJ + "-PE-demux-noprimer.qza"
  log:
    SCRATCH + "/qiime2/logs/" + PROJ + "_primer_q2.log"
  conda:
    "envs/qiime2-2019.4.yaml"
  shell:
    "qiime cutadapt trim-paired \
       --i-demultiplexed-sequences {input.q2_import} \
       --p-front-f {config[primerF]} \
       --p-front-r {config[primerR]} \
       --p-error-rate {config[primer_err]} \
       --p-overlap {config[primer_overlap]} \
       --o-trimmed-sequences {output.q2_primerRM}"

rule dada2:
  input:
    q2_primerRM = SCRATCH + "/qiime2/asv/" + PROJ + "-PE-demux-noprimer.qza"
  output:
    table = SCRATCH + "/qiime2/asv/" + PROJ + "-asv-table.qza",
    rep = SCRATCH + "/qiime2/asv/" + PROJ + "-rep-seqs.qza",
    stats = SCRATCH + "/qiime2/asv/" + PROJ + "-stats-dada2.qza"
  log:
    SCRATCH + "/qiime2/logs/" + PROJ + "_dada2_q2.log"
  conda:
    "envs/qiime2-2019.4.yaml"
  shell:
    "qiime dada2 denoise-paired \
        --i-demultiplexed-seqs {input.q2_primerRM} \
        --p-trunc-q {config[truncation_err]} \
        --p-trunc-len-f {config[truncation_len-f]} \
        --p-trunc-len-r {config[truncation_len-r]} \
        --p-max-ee {config[quality_err]} \
        --p-n-reads-learn {config[training]} \
        --p-chimera-method {config[chimera]} \
        --o-table {output.table} \
        --o-representative-sequences {output.rep} \
        --o-denoising-stats {output.stats}"

rule assign_tax:
  input:
    rep = SCRATCH + "/qiime2/asv/" + PROJ + "-rep-seqs.qza",
    db_classified = DB_classifier
  output:
    sklearn = SCRATCH + "/qiime2/asv/" + PROJ + "-tax_sklearn.qza"
  log:
    SCRATCH + "/qiime2/logs/" + PROJ + "_sklearn_q2.log"
  conda:
    "envs/qiime2-2019.4.yaml"
  shell:
    "qiime feature-classifier classify-sklearn \
	--i-classifier {input.db_classified} \
	--i-reads {input.rep} \
	--o-classification {output.sklearn}"

rule gen_table:
  input:
    table = SCRATCH + "/qiime2/asv/" + PROJ + "-asv-table.qza"
  output:
    table_tsv = SCRATCH + "/qiime2/asv/table/feature-table.biom"
  log:
    SCRATCH + "/qiime2/logs/" + PROJ + "_exportBIOM_q2.log"
  conda:
    "envs/qiime2-2019.4.yaml"
  params:
    directory(SCRATCH + "/qiime2/asv/table")
  shell:
    "qiime tools export --input-path {input.table} --output-path {params}"

rule convert:
  input:
    SCRATCH + "/qiime2/asv/table/feature-table.biom"
  output:
    SCRATCH + "/qiime2/asv/" + PROJ + "-asv-table.tsv"
  log:
    SCRATCH + "/qiime2/logs/" + PROJ + "_exportTSV_q2.log"
  conda:
    "envs/qiime2-2019.4.yaml"
  shell:
    "biom convert -i {input} -o {output} --to-tsv"

rule gen_tax:
  input:
    sklearn = SCRATCH + "/qiime2/asv/" + PROJ + "-tax_sklearn.qza"
  output:
    table_tax = SCRATCH + "/qiime2/asv/tax_dir/taxonomy.tsv"
  log:
    SCRATCH + "/qiime2/logs/" + PROJ + "_exportTAXTSV_q2.log"
  conda:
    "envs/qiime2-2019.4.yaml"
  params:
    directory(SCRATCH + "/qiime2/asv/tax_dir")
  shell:
    "qiime tools export --input-path {input.sklearn} --output-path {params}"

#rule viz:
#  input:
##    qza = SCRATCH + "/qiime2/asv/{proj}-{arti}.qza"
#    qza = expand("{scratch}/qiime2/asv/{proj}-{arti}.qza", scratch = SCRATCH, proj = PROJ, arti = ARTIFACT)
#  output:
##    qzv = SCRATCH + "/qiime2/asv/viz/{proj}-{arti}.qzv"
#    qzv = expand("{scratch}/qiime2/asv/viz/{proj}-{arti}.qzv", scratch = SCRATCH, proj = PROJ, arti = ARTIFACT)
#  log:
#    SCRATCH + "/qiime2/logs/" + PROJ + "_viz_q2.log"
#  conda:
#    "envs/qiime2-2019.4.yaml"
#  shell:
#    "qiime demux summarize --i-data {input.qza} --o-visualization {output.qzv}"
