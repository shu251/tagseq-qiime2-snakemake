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
```
_Output files_  
* ```manifest.txt```: this is the manifest file you are looking for. However, the file path leads to the raw fastq reads before trimming. We will fix this later in the Snakefile. Alternatively, please manually alter OR create this file so you fill in your desired sample names in the first column, the full path to the fastq reads you want to make ASVs from, and the last column states the direction of the paired read
* ```SampleList.txt```: simply a list of the SRR IDs and the sample names you actually get from the SraRunInfo.csv file you can download from SRA/NCBI.


## Modify your config file
Open and take a look at ```config.yaml```. Below is a breakdown of the config file as you will need to modify this to fit your needs before running Snakemake with your data. 
This config.yaml file is set up to run with the test sequences downloaded from above
```
proj_name: Diel_18S # Replace this with your project name. This will be what the final ASV qiime2 artifact and ASV table is named.
scratch:  /vortexfs1/scratch/sarahhu # Change this to where you want all your outputs to be written. Since I am working on a HPC, I write all of this to scratch and move it later.
outputDIR: /vortexfs1/omics/huber/shu/tagseq-qiime2-snakemake
primerF: CCAGCASCYGCGGTAATTCC # Forward primer sequence
primerR: ACTTTCGTTCTTGATYRA # Reverse primer sequence
manifest: manifest.txt #Input of manifest file that you need to provide. Or use the R script written above to generate this file
sample_names: SampleList.txt #Sample list output from the R script used above
manifest-trimmed: manifest-trimmed.txt # Final manifest file, the Snakemake pipeline will create this file, as the first few steps generate the trimmed fastq files which we ACTUALLY want to use in the QIIME2 part, so the Snakemake file is written to input your manifest file and modify the filepath

# QIIME2-specific flags
## Primer removal
primer_err: 0.4
primer_overlap: 3


## DADA2 - ASV flags
truncation_err: 2
truncation_len-f: 200
truncation_len-r: 200
quality_err: 2
training: 1000 #should be set higher for a non-test dataset
chimera: pooled
```

## Execute Snakemake to run everything



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
