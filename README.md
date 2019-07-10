# tagseq-qiime2-snakemake
Pipeline to run qiime2 with snakemake



# Set up

* make working directory of raw_data (sra explorer tool), intro dataset as test
* Bioproject: PRJNA393172 - test data from Hu et al. 2018
* Go to [SRA explorer](https://ewels.github.io/sra-explorer/) to search with BioProject, select bash script for downloading FastQ files.
* Make sure database is set up  - link to other repo
* Generate manifest file

## Set up working directory & conda environment
```
git clone https://github.com/shu251/tagseq-qiime2-snakemake.git
cd tagseq-qiime2-snakemake

# Create conda environment
conda env create --name snake-tagseq --file envs/snake.yaml 

# Enter environment
source activate snake-tagseq

# Check versions and correct environment set up
snakemake --version
qiime info
```

## Download tag-sequencing data from SRA

```
cd tagseq-qiime2-snakemake/raw_data/
# Copy output from SRA explorer into a new text file 'download.sh'
nano download.sh

# Download all sequences
bash download.sh
```

## Create required files for input into QIIME2 and snakemake
The ```write-manifest.R``` script will input the list of fastq files it finds in ```raw_data/``` to generate a QIIME2 specific manifest file (manifest.txt) and a text file with the SRR ID list for input into snakemake and a corresponding list of sample names.
```
# Run R script
Rscript write-manifest.R

# outputs:
## manifest.txt
## SampleList.txt

#
# Note to update this later to input trimmed reads
#


```


## Run quality control and trimming
Modify config.yaml file to include path to raw sequences and where output trimming and QC files will be written to.

Snakemake pipeline performs fastqc on all raw reads. Then uses trimmomatic to remove barcodes and repeats the fastqc post-trimming. Finally, uses data from fastqc runs to generate a multiqc report.


```
snakemake --use-conda
# Error message related to how snakemake is reading in conda (finding?) conda environment

#
# Issue with conda - snakemake versions not compatible
#


```


## Generate ASVs
Run write-manifest.R script to input trimmed reads instead of raw fastq reads.  

Update ```config.yaml``` with preferences for how DADA2 will trim, denoise, chimera check, and perform ASV determination.

For more information ```qiime dada2 denoise-paired --help```

```
# Modify the config.yml file to add your preferences for each flag:
#Truncate sequences at the 3' end when sequence quality may decrease
--p-trunc-len-f #forward read truncation
--p-trunc-len-r #reverse read truncation

#Trim forward and reverse reads at 5' end based on quality
--p-trim-left-f #forward read trim
--p-trim-left-r #reverse read trim

#Quality threshold for above trimming
--p-trunc-q

#Any reads with expected error higher than this value will be removed
--p-max-ee

#Choice of chimera removal - Choices('pooled', 'none', 'consensus')
--p-chimera-method

#Number of reads to consider in training set for error model
--p-n-reads-learn # defaul is 1 million
```

Output from dada2 pipeline is 3 QIIME2 artifacts:
* A feature table - or 'OTU/ASV table' ```--o-table```
* A table of the representative sequences for the feature table ```--o-representative-sequences```
* Statistics related to denoising and ASV determination ```--o-denoising-stats```


## Visualization QIIME2 run so far
Next add the summarization features to get .qzv files.. see previous documentation on these visualizations

## Assign taxonomy
Assign taxonomy - link to db-build repo and introduce 2 ways to perform this
  Add to snakefile, provide the option for the user to specify sklearn or the consensus-vsearch option?


## Generate output tables
show commands - add to Snakefile

## Compile taxonomy and count tables
R script for this



### TROUBLESHOOTING
* conda and snakemake version control? only applies to wrapper usage, but anticipate it will become a larger issue
* mv to separate directory? wasn't there an issue before
