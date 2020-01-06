# Read in SRA sample IDs and list of raw fastq files to generate a manifest file ready for qiime2 import

# Required libraries
library(reshape2);library(tidyverse)

# Import all files in currect directory
paths <- as.data.frame(list.files(pattern = ".fastq.gz", full.names = TRUE))
colnames(paths)[1]<-"FASTQ"

# Get local path & add to dataframe
path_to_files <- getwd()
paths$PATH <- path_to_files

# Extract sample ID
paths_run <- paths %>% mutate(SAMPLEID = str_replace(FASTQ, "(\\w*?)_(\\w+).fastq.gz","\\1"))
## ^See wildcard options on this line to modify how R script pulls out your sample IDs from fastq files

paths_run$SAMPLEID <- gsub("./","", paths_run$SAMPLEID)
paths_run$FASTQ <- gsub("./","", paths_run$FASTQ)

# Write full path
paths_run$FULL_PATH <- paste(paths_run$PATH, paths_run$FASTQ, sep="/")

# Write manifest:
manifest_orig <- data.frame(paths_run$SAMPLEID, paths_run$FULL_PATH)
colnames(manifest_orig)[1:2]<-c("sample-id","absolute-filepath")

# If R1_001.fastq.gz or _1.fastq.gz ending, fill new column with forward, else reverse
manifest_orig$direction <- ifelse(grepl("R1_001.fastq.gz|_1.fastq.gz|R1.fastq.gz", manifest_orig$`absolute-filepath`), "forward", "reverse")
## ^Works with these options as file name suffix (R1_001.fastq.gz|_1.fastq.gz|R1.fastq.gz), modify if you have something else

# Write output as a manifest file
write.table(manifest_orig, "manifest-orig.txt",quote=FALSE,col.names=TRUE,row.names=FALSE,sep=",")


