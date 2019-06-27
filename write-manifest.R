# Read in SRA sample IDs and list of raw fastq files to generate a manifest file ready for qiime2 import
# Required libraries
library(reshape2)

# Extract raw fastq.gz file list, first the absolute paths and then only the file names
## Uses the raw_data directory to find all fastq.gz files as input
files<-as.data.frame(list.files(path="raw_data", pattern= ".fastq.gz"))
paths<-as.data.frame(list.files(path="raw_data", pattern= ".fastq.gz", full.names=TRUE))
join<-cbind(paths,files)
colnames(join)[1:2]<-c("filepath", "tmp")

# Parse sample-ID and file path
out<-colsplit(join$tmp, "_", c("Run","proj"))
out2<-colsplit(out$proj, "\\.", c("sample-id", "ReadPair"))
df<-data.frame(out2$`sample-id`, join$filepath)
colnames(df)[1:2]<-c("sample-id","absolute-filepath")

# Generate direction column
df$direction<-ifelse(grepl("_1.fastq.gz", df$`absolute-filepath`),"forward", "reverse")
write.table(df, "manifest.txt",quote=FALSE,col.names=TRUE,row.names=FALSE,sep=",")


# First line must include the column headers
# Example:
# sample-id,absolute-filepath,direction
# sample2,/path/to/file/sample2_1.fastq.gz,forward
# sampl3,/path/to/file/sample3_2.fastq.gz,reverse

