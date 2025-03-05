library(data.table)
library(ggplot2)

output_folder="/home/cnorton5/scr4_hlee308/cnorton5/old_nanopore/overlap_hg19"
gene_mapping="/home/cnorton5/data_hlee308/cnorton5/ref/hg19_other/tx2gene_hg19.67.csv"

table1 <- fread(paste0(output_folder, "/overlap_matrix.tsv"), sep="\t")
mapping <- fread(gene_mapping, sep=",")

table2 <- merge(table1, mapping, by.x="Gene", by.y="transcript_id")

#Let's put the last column (gene_name) in the first column
table2 <- table2[, c(7, 1:6)]

#Sort by the first column
table2 <- table2[order(table2$gene_name),]

#Rename the columns
colnames(table2) <- c("Gene", "ENSEMBL", "UP", "EXONS", "INTRONS", "DOWN", "ANY")

#Drop the ENSEMBL column
table2 <- table2[, -2]

#Now, we need to collapse the rows that have the same gene name, and sum the values in columns 3-7
table2 <- table2[, lapply(.SD, sum), by=Gene]

#Then, replace any non-zero values with 1
table2[, 2:6] <- lapply(table2[, 2:6], function(x) ifelse(x > 0, 1, 0))

write.table(table2, file=paste0(output_folder, "/overlap_matrix_gene.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

#Let's get the sum of each column 3-7
sums <- colSums(table2[, 2:6])







