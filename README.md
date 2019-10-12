# Pipeline to run qiime2 with snakemake
_Sept 2019_
_Sarah Hu_

Run [QIIME2](https://www.nature.com/articles/s41587-019-0209-9) using [snakemake](https://academic.oup.com/bioinformatics/article/28/19/2520/290322). Requires input of raw .fastq reads. Below describes steps to set up environment and run a test dataset.
This Snakemake pipeline can be easily scaled up to larger datasets and includes scripts to submit Snakemake jobs to slurm.  

## _Pipeline_
1. Generate manifest.txt file from a directory containing all raw .fastq files (R)
2. Read in raw .fastq files and perform *fastqc*, *Trimmomatic*, and repeat *fastqc* on newly trimmed reads (snakemake)
3. Modify manifest.txt file so input files are the trimmed .fastq reads
4. Import all trimmed reads as QIIME2 artifact
5. Remove primers with cutadapt

_Amplicon Sequence Variants_  

6. Run DADA2, which performs additional quality trimming, filtering, and denoising steps. Also removes chimeric sequences. Finally, determines *Amplicon Sequence Variants*
7. Assigns taxonomy using QIIME2 feature classifer
8. Generates ASV count and taxonomy tables
9. Compiles ASV count + taxonomy table (R)

_Operational Taxonomic Units_   

10. Merged paired end reads and filter reads by qscore
11. Dereplicate sequences
12. Cluster into OTUs (open reference & closed reference currently written)
13. Determine chimeras using uchime with the reference database
14. Remove chimeras from representative sequences and cluster table
15. Assign taxonomy using QIIME2 feature classifier
16. Generate OTU count and taxonomy tables
17. Compile OTU count + taxonomy table
***
## Before starting
* If you're new to snakemake and/or qiime2, run the below using the provided test data. If you want to learn more about qiime2, [see tutorials on their website](https://docs.qiime2.org/2019.7/). And if you're new to snakemake, [learn more here](https://snakemake.readthedocs.io/en/stable/) or follow a recommended [tutorial](https://github.com/ctb/2019-snakemake-ucdavis).
* Before running qiime2, build your preferred database for assigning taxonomy, [directions here for a microeuk db](https://github.com/shu251/db-build-microeuks)
* Familiarize yourself with conda environments if you are not already [Explanation using R](https://alexanderlabwhoi.github.io/post/anaconda-r-sarah/)
* Be familiar with the qiime2 options to determine Amplicon Sequence Variants

## 1. Set up working directory & conda environment

Create and launch a conda environment to run this pipeline.

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

### _1.1 Place all raw fastq files in one place_

_Two options to use test data_:
* Download test data from Hu et al. 2018 using [SRA explorer](https://ewels.github.io/sra-explorer/) and search for *Bioproject: PRJNA393172*. Copy the output of bash script to a new file in ```tagseq-qiime2-snakemake/raw_data/``` and execute to download.
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

### **If using your own fastq files**  
* Place them in their own ```raw_data``` directory (see step to update config.yaml to tell snakemake where these are located).
* Make sure they are labeled so the last numbers/letters in the file names designate the read pair. Input fastq files should be labeled either: Sample01_treatment1_R1.fastq.gz or Sample01_treatment1_1.fastq.gz

## 2. Create required files for input into Snakefile

The ```write-manifest.R``` script will input the list of fastq files it finds in ```raw_data/``` to generate a QIIME2 specific manifest file (manifest.txt) and a text file with the SRR ID list for input into snakemake and a corresponding list of sample names. Ensure you are in an environment to run an R script. [Follow these directions for more information](https://alexanderlabwhoi.github.io/post/anaconda-r-sarah/). To execute this script, enter an R environment or if already enabled, use ```Rscript``` (see below).

*Run the R code*
```
# Run R script (make sure you are in an R-enabled environment)
Rscript write-manifest.R
```
_Output files_
* ```manifest.txt```: this is the manifest file you are looking for. However, the file path leads to the raw fastq reads before trimming. We will fix this later in the Snakefile. Alternatively$
* ```SampleList.txt```: simply a list of the SRR IDs and the sample names you actually get from the SraRunInfo.csv file you can download from SRA/NCBI.

* Run below with dataset downloaded from SRA, or on your own paired end fastq files
* Generate a .csv file to link SRR IDs (or raw fastq sequence names) to actual sample names (if applicable)
 
_Why do I need this manifest.txt & SampleList.txt file?_   
The SRR fastq files I downloaded directly from SRA do not have actual sample names. So after [searching the NCBI/SRA database here](https://www.ncbi.nlm.nih.gov/sra/?term=PRJNA393172), I need to download a file that will link these uninformative SRRXXX IDs to actual sample names. One way to do this is to click on the *'Send to:'* drop down menu on the upper right of the webpage and click *File*. Then select *'RunInfo'* format and create the file. This will download a .csv file that has sample information and all the other metadata the sequence submitter provided to SRA. This information can be used for the ```sample-id``` column in the ```manifest.txt``` file which is required to run QIIME2. _To run this on your own data_, either make your own ```manifest.txt``` (see below) or download this .csv file in the same way. Alternatively, you can modify the ```write-manifest.R``` code to import a table that links your SRR IDs to sample names that you specify.

*Make your own manifest.txt*
Format your own txt file (using a text editor or a CSV file) exactly this way. Every line must be comma separated, with sample-id, followed by the path of fastq file, followed by either "forward" or "reverse". 
```
sample-id,absolute-filepath,direction
sample1,$PWD/raw_seqs_dir/Test01_full_L001_R1.fastq.gz,forward
sample1,$PWD/raw_seqs_dir/Test01_full_L001_R2.fastq.gz,reverse
sample2,$PWD/raw_seqs_dir/Test02_full_L001_R1.fastq.gz,forward
sample2,$PWD/raw_seqs_dir/Test02_full_L001_R2.fastq.gz,reverse
```
* Replace $PWD with your path
* The fastq files can be gziped. 
* List all of your fastq files. 
* Save the file and name it "manifest.txt".


## 3. Modify your config file to run snakemake
Now, take a look at ```config.yaml```. Below is a breakdown of the parameters you need to revise in your config file. Edit this file to fit your computer set-up. 
This config.yaml file is set up to run with the test sequences downloaded from above:

```proj_name: Diel_18S``` Replace this with your project name. This will be what the final ASV qiime2 artifact and ASV table is named.  
```raw_data: /vortexfs1/omics/huber/shu/tagseq-qiime2-snakemake/raw_data``` Point config file to location of your raw fastq reads, in this case a directory called 'raw_data'   
```scratch:  /vortexfs1/scratch/sarahhu``` Change this to where you want all your outputs to be written. Since I am working on a HPC, I write all of this to scratch and move it later.   
```manifest: manifest.txt``` Input of manifest file that you need to provide, generated using R script or created by you.   
```sample_names: SampleList.txt``` Sample list output from the R script used above   
```manifest-trimmed: manifest-trimmed.txt``` Final manifest file, the Snakemake pipeline will create this file, as the first few steps generate the trimmed fastq files which we ACTUALLY want to use in the QIIME2 part, so the Snakemake file is written to input your manifest file and modify the filepath   

_The config file_ also includes parameters you set when generating ASVs.
```primerF: CCAGCASCYGCGGTAATTCC``` Forward primer sequence, in this case I'm using the V4 hypervariable region  
```primerR: ACTTTCGTTCTTGATYRA``` Reverse primer sequence   

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

For clustering OTUs and the steps preceeding OTU clustering, set parameters in the config.yaml file:
```
# Merge paired end reads
minoverlap: 10
maxdiff: 4
minlength: 270

# Quality filtering of merged reads
minphred: 4
qualwindow: 3 

# Open reference OTU clustering
perc_id: 0.90
otu-thread: 1

# Chimera checking
chimera-thread: 1
```
If you're unsure about parameters, consult the QIIME2 reference manual.


Ahead of running this pipeline, prepare a database so the reference sequences can be assigned a taxonomy.
You can use [this pipeline](https://github.com/shu251/db-build-microeuks) to do this.
*In the config.yaml file*, change this line to direct the Snakefile to where your database is stored.
```
## Path to your database for assigning taxonomy
database: /vortexfs1/omics/huber/shu/db/pr2-db/V4-pr2_4.11.1-classifier.qza
```

## 4. Test ```snakemake``` dry run & execute full pipeline

Use *-s* to select between the Snakemake pipelines
```
snakemake -np -s Snakefile-asv
snakemake -np -s Snakefile-otu
# output should be all green and display no errors
``` 

To run the full pipeline make sure you enable the ```--use-conda``` flag. This is because snakemake uses the conda environments stored in ```envs/``` to execute rules.
```
# For ASVs
snakemake --use-conda -s Snakefile-asv

# For OTUs
snakemake --use-conda -s Snakefile-otu
```

## 5. Run on HPC with SLURM

[Read about executing snakemake on a cluster](https://snakemake.readthedocs.io/en/stable/executable.html) and another explanation on how to execute with a submit script can be found [here](https://hpc-carpentry.github.io/hpc-python/17-cluster/).    
Review the submit scripts available in ```submitscripts```. Files in this directory include another config file called ```cluster.yaml```, and two scripts to submit your snakemake pipeline to the cluster with and without the dry run option.   
First, open and modify the ```cluster.yaml``` to fit your machine. Then test run using the provided submit scripts.
```
# Make sure you are still in the snake-tagseq conda environment

## For ASVs:
bash submitscripts/dry-submit-slurm-ASV.sh 

## For OTUs:
bash submitscripts/dry-submit-slurm-OTU.sh
```
Outputs should all be green and will report how many jobs and rules snakemake plans to run. This will allow you to evaluate any error and syntax messages.  

Once ready, use the submit-slurm.sh script to submit the jobs to slurm. Run with screen, tmux, or nohup.
```
# Full run for ASVs:
bash submitscripts/submit-slurm-ASV.sh

# Full run for OTUs:
bash submitscripts/submit-slurm-OTU.sh
```
Run the above in *tmux* or *screen*, as it will print out the jobs it submits to SLURM. You will be able to monitor the progress this way.

## 6. Output from pipeline

* **ASV table:** ```[PROJ]-asv-table.tsv``` includes samples by columns and ASVs by row. Values represent number of sequences per ASV (row) in a given sample (column).  
* **Assigned taxonomy:** ```tax_dir/taxonomy.tsv``` represents the full taxonomic name for each ASV. The same ASV identifer (string of letters and numbers under 'Feature ID') as the asv-table.  

To combine these table, you can run the R script from the ```../../qiime2/asv/``` directory.
```
# Migrate to directory with .qza and .tsv outputs from snakemake pipeline
cd ../../../qiime2/asv/

# Ensure R is enabled
# conda activate r_3.5.1 # to activate the R conda environment

# Path to R script
Rscript /vortexfs1/omics/huber/shu/tagseq-qiime2-snakemake/make-asv-table.R
```

* **CountTable-wtax-DATE.txt** - ASV table with counts per sample and the taxonomic identities
* **Output_stats.txt** - quick stats on results, including how many sequences per sample, and how many ASVs with only 1 or 2 sequences are in final results

Find an introduction to R for processing ASV or OTU tables [here](https://github.com/shu251/PreliminaryFigures_V4_tagseq).

## 7. qiime2 visualization option

QIIME2 offers away to visualize the data types (artifact files) using [an interative viewer](https://docs.qiime2.org/2019.4/concepts/#data-files-visualizations). An output directory of these ```.qzv``` files will be created at the end of the Snakefile run. These files can be brought locally and [drag and dropped to here](https://view.qiime2.org). In the above QC steps, whenever a .qza file was worked on, you had the option to run this:
```
# ../../qiime2/asv/ #In this directory

# Activate a qiime2 environment (see instructions on QIIME2 website)
conda activate qiime2-2019.4

# Insert any of the .qza artifact files generated
qiime demux summarize --i-data PROJECT-STEP.qza --o-visualization PROJECT-STEP.qzv
```
***

## References
* Köster, Johannes and Rahmann, Sven. “Snakemake - A scalable bioinformatics workflow engine”. Bioinformatics 2012. https://doi.org/10.1093/bioinformatics/bts480.
* Bolyen E. _et al._ 2019. Reproducible, interactive, scalable and extensible microbiome data science using QIIME 2. Nature Biotechnology 37: 852–857. https://doi.org/10.1038/s41587-019-0209-9
* Wingett SW, Andrews S. FastQ Screen: A tool for multi-genome mapping and quality control. F1000Res. 2018 Aug 24 [revised 2018 Jan 1];7:1338. doi: 10.12688/f1000research.15931.2. eCollection - https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
* Ewels P., _et al._ 2016. MultiQC: summarize analysis results for multiple tools and samples in a single report. https://doi.org/10.1093/bioinformatics/btw354.
* Bolger, A. M., Lohse, M., & Usadel, B. (2014). Trimmomatic: A flexible trimmer for Illumina Sequence Data. Bioinformatics, btu170.
* Martin, M., 2011. Cutadapt Removes Adapter Sequences From High-Throughput Sequencing Reads. https://doi.org/10.14806/ej.17.1.200.
* Callahan BJ, McMurdie PJ, Rosen MJ, Han AW, Johnson AJ, Holmes SP. DADA2: High-resolution sample inference from Illumina amplicon data. Nat Methods. 2016;13(7):581–583. doi:10.1038/nmeth.3869
* R Core Team (2017). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/
***


#### To do:
* Update so R script is integrated into snakemake and creates tables with taxonomy for both ASV and OTU pipelines
* Add in final step to generate a fasta file with fasta headers which correspond to the reference sequences

_last updated 09-24-2019_

