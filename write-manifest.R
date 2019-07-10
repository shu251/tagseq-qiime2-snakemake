# Read in SRA sample IDs and list of raw fastq files to generate a manifest file ready for qiime2 import

# Required libraries
library(reshape2);library(tidyverse)

# Import SRA run, extract run id and sample names so we can link meaningful sample names with fastq IDs
samplename <- read.csv("raw_data/SraRunInfo.csv")
paths <- as.data.frame(list.files(path = "raw_data", pattern = ".fastq.gz", full.names = TRUE))
colnames(paths)[1]<-c("absolute-filepath")

# Extract Run information
paths_run <- paths %>% mutate(Run = str_replace(`absolute-filepath`, "raw_data/SRR(\\d*?)_(\\w+).fastq.gz","SRR\\1"))

# Join with sample name from Sra table
paths_runinfo <- inner_join(paths_run, samplename, by="Run")
paths_runinfo$pre_path <- paste(getwd(), "/", sep="") # add prefix so final path will be absolut
paths_runinfo$PATH <- paste(paths_runinfo$pre_path, paths_runinfo$`absolute-filepath`, sep="")

manifest <- data.frame(paths_runinfo$SampleName, paths_runinfo$PATH)

colnames(manifest)[1:2]<-c("sample-id","absolute-filepath")

# if _1.fastq.gz fill new column with forward, else reverse
manifest$direction<-ifelse(grepl("_1.fastq.gz", manifest$`absolute-filepath`),"forward", "reverse")

# Write output as a manifest file
write.table(manifest, "manifest.txt",quote=FALSE,col.names=TRUE,row.names=FALSE,sep=",")


# Generate sample list:
sample_list1 <- data.frame(paths_runinfo$SampleName, paths_runinfo$Run)
colnames(sample_list1)[1:2]<-c("SampleName", "Run")
sample_list2 <- unique(sample_list1)
write.table(sample_list2, "SampleList.txt",quote=FALSE,col.names=TRUE,row.names=FALSE,sep="\t")

# First line must include the column headers
# Example:
# sample-id,absolute-filepath,direction
# sample2,/path/to/file/sample2_1.fastq.gz,forward
# sampl3,/path/to/file/sample3_2.fastq.gz,reverse

