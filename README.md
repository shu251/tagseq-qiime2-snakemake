# Pipeline to run qiime2 with snakemake
_July 2019_
_Sarah Hu_

## Set up before starting
* Run below with dataset downloaded from SRA, or on your own paired end fastq files
* Generate a .csv file to link SRR IDs (or raw fastq sequence names) to actual sample names (if applicable)
* Build database for assigning taxonomy, [directions here](https://github.com/shu251/db-build-microeuks)
* Familiarize yourself with conda environments if you are not already [Explanation using R](https://alexanderlabwhoi.github.io/post/anaconda-r-sarah/)

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
## Output should be: 5.5.2
```

### Place all raw fastq files in one place

_Two options to use test data_:
* Download test data from Hu et al. 2018 using [SRA explorer](https://ewels.github.io/sra-explorer/) and search for Bioproject: PRJNA393172. Copy the output of bash script to a new file in ```tagseq-qiime2-snakemake/raw_data/``` and execute to download.
* Use provided bash scripts to download full or subset of sequences (this script was generated using the SRA explorer)

```
## To use provided bash script:
# migrate to raw_data directory
cd raw_data

# Download subset of sequences (n = 10)
bash download-subset.sh

# or download the full set:
# bash download.sh

# Exit raw_data directory
cd ..
```

**If using your own fastq files** place in ```raw_data``` directory.

## Create required files for input into Snakefile
The ```write-manifest.R``` script will input the list of fastq files it finds in ```raw_data/``` to generate a QIIME2 specific manifest file (manifest.txt) and a text file with the SRR ID list for input into snakemake and a corresponding list of sample names. Ensure you are in an environment to run and R script. [Follow these directions for more information](https://alexanderlabwhoi.github.io/post/anaconda-r-sarah/)
```
# Run R script (make sure you are in an R-enabled environment)
Rscript write-manifest.R
```
_Output files_  
* ```manifest.txt```: this is the manifest file you are looking for. However, the file path leads to the raw fastq reads before trimming. We will fix this later in the Snakefile. Alternatively, please manually alter OR create this file so you fill in your desired sample names in the first column, the full path to the fastq reads you want to make ASVs from, and the last column states the direction of the paired read
* ```SampleList.txt```: simply a list of the SRR IDs and the sample names you actually get from the SraRunInfo.csv file you can download from SRA/NCBI.

_Why do I need this SampleList.txt file?_ The SRR fastq files I downloaded directly from SRA do not have actual sample names. So after [going to this page](https://www.ncbi.nlm.nih.gov/sra/?term=PRJNA393172), click on the *'Send to:'* drop down menu on the upper right of the webpage and click *File*. Then select *'RunInfo'* format and create the file. This will download a .csv file that has sample information and all the other metadata the sequence submitter provided to SRA. The manifest.txt file input into qiime2 requires this information and it will make your life easier to change your file names to usable sample names for your study, rather than a series of random letters and numbers.


## Modify your config file
Open and take a look at ```config.yaml```. Below is a breakdown of the config file as you will need to modify this to fit your needs before running Snakemake with your data. 
This config.yaml file is set up to run with the test sequences downloaded from above:

```proj_name: Diel_18S``` Replace this with your project name. This will be what the final ASV qiime2 artifact and ASV table is named.
```scratch:  /vortexfs1/scratch/sarahhu``` Change this to where you want all your outputs to be written. Since I am working on a HPC, I write all of this to scratch and move it later.
```outputDIR: /vortexfs1/omics/huber/shu/tagseq-qiime2-snakemake``` Change this to your eventual output directory, file will get moved later from scratch to this output
```primerF: CCAGCASCYGCGGTAATTCC``` Forward primer sequence, in this case I'm using the V4 hypervariable region
```primerR: ACTTTCGTTCTTGATYRA``` Reverse primer sequence
```manifest: manifest.txt``` Input of manifest file that you need to provide. Or use the R script written above to generate this file
```sample_names: SampleList.txt``` Sample list output from the R script used above
```manifest-trimmed: manifest-trimmed.txt``` Final manifest file, the Snakemake pipeline will create this file, as the first few steps generate the trimmed fastq files which we ACTUALLY want to use in the QIIME2 part, so the Snakemake file is written to input your manifest file and modify the filepath

_The config file_ also includes parameters you set when generating ASVs.
Modify primer removal step with error and required overlap
```
primer_err: 0.4
primer_overlap: 3
```

For the DADA2 step to determine ASVs, you also need to specify the forward and reverse read truncation
```
truncation_err: 2
truncation_len-f: 200
truncation_len-r: 200

#Truncate sequences at the 3' end when sequence quality may decrease
--p-trunc-len-f #forward read truncation
--p-trunc-len-r #reverse read truncation

#Trim forward and reverse reads at 5' end based on quality
--p-trim-left-f #forward read trim
--p-trim-left-r #reverse read trim


quality_err: 2
training: 1000 #should be set higher for a non-test dataset
chimera: pooled
#Quality threshold for above trimming
--p-trunc-q

#Any reads with expected error higher than this value will be removed
--p-max-ee

#Choice of chimera removal - Choices('pooled', 'none', 'consensus')
--p-chimera-method

#Number of reads to consider in training set for error model
--p-n-reads-learn # defaul is 1 million
```

## Visualization QIIME2 run so far
Next add the summarization features to get .qzv files.. see previous documentation on these visualizations

## Assign taxonomy
Assign taxonomy - link to db-build repo and introduce 2 ways to perform this
  Add to snakefile, provide the option for the user to specify sklearn or the consensus-vsearch option?


## Generate output tables
show commands - add to Snakefile

## Compile taxonomy and count tables
R script for this



### To do:
* *Error with Child-Parent directory-file issue on the last few rules to export qiime2 information*
* move things to output dir
* enter provided R environment and run R script to compile output tax and asv information, link to other tutorial
* check to see how directory stuff is set up.. ok?
* add visualization rule to generate all the qzv files...
* How to select or skip steps? i.e. don't repeat dada2? dont do triming and fastqc step?
* How to add an option? Alternate snakefiles in other directories for (a) vsearch tax assignment and (b) OTU clustering?
