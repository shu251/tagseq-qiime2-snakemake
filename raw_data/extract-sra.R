sra_info <- read.csv("SraRunInfo.csv")

# Generate column to download forward reads
sra_info$TMP_1 <- paste("curl -L ", sra_info$download_path, "/", sra_info$Run, "_1.fastq.gz -o ", sra_info$Run, "_", sra_info$SampleName, "_1.fastq.gz", sep="")

# Repeat for reverse reads
sra_info$TMP_2 <- paste("curl -L ", sra_info$download_path, "/", sra_info$Run, "_2.fastq.gz -o ", sra_info$Run, "_", sra_info$SampleName, "_2.fastq.gz", sep="")

# Write each of these outputs to a table
write.table(sra_info$TMP_1, file = "download_read1.sh", quote = FALSE, sep = "", row.names = FALSE, col.names = FALSE)

write.table(sra_info$TMP_2, file = "download_read2.sh", quote = FALSE, sep = "", row.names = FALSE, col.names = FALSE)
