### Import qiime2 outputs to generate ASV table & print preliminary ASV results
library(dplyr)

# Find output tsv files with wildcards:
count_tsv <- list.files(pattern = "asv-table.tsv")
tax_tsv <- list.files(pattern = "taxonomy.tsv", path = "tax_dir/", full.names = TRUE)
# Import tsv files
count <- read.delim(count_tsv, header=FALSE, row.names=NULL, stringsAsFactors=FALSE, na.strings = "n/a")
tax <- read.delim(tax_tsv, header=TRUE, row.names=NULL, stringsAsFactors=FALSE, na.strings = "n/a")
# head(count)
# head(tax)


# Re-classify and re-format
colnames(count)<-as.character(unlist(count[2,]))
count1<-count[-c(1,2),]
colnames(count1)[1]<-"Feature.ID"
# head(count1[1:2,]) # check dataframe structure
x<-dim(count1)[2] # number of columns
# Convert class in count1 to numerics
count2<-count1
x<-dim(count2)[2] # number of columns
count2[2:x] <- lapply(count2[2:x], function(x) as.numeric(as.character(x)))
# Get base stats and generate a table called "dataset_info"
TOTAL_SEQUENCES<-apply(count2[2:x],2,sum) #number of sequences per sample
TOTAL_ASV<-colSums(count2[2:x]>0) # ASVs/OTUs per sample
SINGLE_ASV<-colSums(count2[2:x]==1) # ASVs/OTUs with only 1 seq
DOUBLE_ASV<-colSums(count2[2:x]==2) # ASVs/OTUs that have only 2 seqs
dataset_info<-data.frame(TOTAL_SEQUENCES, TOTAL_ASV,
                         SINGLE_ASV, DOUBLE_ASV)
# Write stats as output table:
write.table(dataset_info, file="Output_stats.txt", quote=FALSE, row.names=TRUE, col.names=TRUE)


# Join count table with taxonomy information
counts_wtax <- left_join(count2, tax, by = "Feature.ID")
# Write to new table - date included
write.table(counts_wtax, file = (paste("CountTable-wtax-", Sys.Date(),".txt",sep="")), quote=FALSE, row.names=FALSE, col.names=TRUE)
