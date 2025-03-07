if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
#BiocManager::install(c("bsseq"))

library(data.table)
library(bsseq)
library(BiocParallel)

#This will need to be started within a conda environment with a newer version of R
# than what the cluster has available. In my case, I used "conda activate basicR"

#Load the data
base_dir<- "/home/cnorton5/scr4_hlee308/cnorton5/old_nanopore"
output_folder=paste0(base_dir, "/DMR_Results_hg19")
input_folder1=paste0(base_dir,"/ONT_HG01_hg19/methylbed")
input_folder2=paste0(base_dir,"/ONT_HG02_hg19/methylbed")


files <- c(
  paste0(input_folder1, "/hg19_barcode01.sorted.bed.gz"),
  paste0(input_folder1, "/hg19_barcode04.sorted.bed.gz"),
  paste0(input_folder2, "/EXP-NBD104_barcode01.bed.gz"),
  paste0(input_folder2, "/EXP-NBD104_barcode04.bed.gz")
)


bsseq_nano <- bsseq::read.modkit(files, rmZeroCov = TRUE, strandCollapse=TRUE)

bsseq_nano <- sort(bsseq_nano)
#Let's save the object
saveRDS(bsseq_nano, file=paste0(output_folder, "/bsseq_nano.RData"))

#Read in the object
bsseq_nano <- readRDS(file=paste0(output_folder, "/bsseq_nano.RData"))

#Smooth the object
## TODO WHY WON"T THIS SMOOTH
## maybe try: less workers (4), smaller gap, smaller h
bsseq_nano_smooth <- BSmooth(BSseq = bsseq_nano, ns=70, h=1000, maxGap = 10^8, BPPARAM = MulticoreParam(workers = 8, progressbar = TRUE))

#Now that the object is smoothed, let's save it
saveRDS(bsseq_nano_smooth, file=paste0(output_folder, "/bsseq_nano_smooth.RData"))

#Read the smoothed object
bsseq_nano_smooth <- readRDS(file=paste0(output_folder, "/bsseq_nano_smooth.RData"))

#Remove the original object
rm(bsseq_nano)

#Now let's inspect the samples within the object
samples <- sampleNames(bsseq_nano_smooth)

#Let's assign labels
control <- grep("barcode01", samples, value = TRUE)
test <- grep("barcode04", samples, value = TRUE)

#DMR finder
tstat <- BSmooth.tstat(bsseq_nano_smooth, test, control, estimate.var="same", verbose=FALSE, local.correct=FALSE, mc.cores=8)
regions <- dmrFinder(tstat, stat="tstat", cutoff=c(-3,3), maxGap=300, verbose=TRUE)
regions <- regions[regions$n > 4]

#Let's save the regions
write.table(regions, file=paste0(output_folder, "/soft_DMR_regions.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

#Let's get value counts of regions chr column
table(regions$chr)
