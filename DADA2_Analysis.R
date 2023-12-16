# Based on the tutorial: https://benjjneb.github.io/dada2/tutorial.html
# Adapted by Felipe Vaz Peres (felipe.vzps@gmail.com)

# Cool references:
# https://www.nature.com/articles/nrmicro3451.pdf
# https://www.nature.com/articles/nmicrobiol201649.pdf
# https://www.nature.com/articles/ismej2017119
# https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0227434
# https://www.mdpi.com/2076-2607/10/10/1961

# Dataset:
# Six soil samples from Amazonian forest and pasture.
# Sequencing was performed on the Illumina MiSeq platform in PE 250bp mode.
# Primers used were 515F (Parada et al. 2016) and 806R (Apprill et al. 2015) -> 19 and 20 bp, respectively.
# https://earthmicrobiome.org/protocols-and-standards/16s/

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2")

library(dada2)
package.version("dada2")

library(ShortRead)
package.version("ShortRead")

getwd()
setwd("/home/felipe/Documents/sample-inference-from-amplicon/rawData")
path <- "/home/felipe/Documents/sample-inference-from-amplicon/rawData"
list.files(path)

fnFs <- sort(list.files(path, pattern="_R1.fastq.gz"))
fnRs <- sort(list.files(path, pattern="_R2.fastq.gz"))

# Extracting sample names, assuming file names follow the format: SAMPLENAME_XXX.fastq
sampleNames <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# Specify the full path for fnFs and fnRs
fnFs <- file.path(path, fnFs)
fnRs <- file.path(path, fnRs)
fnFs[1:6]
fnRs[1:6]

# Most Illumina sequencing data shows a trend of decreasing average quality towards the end of sequencing reads
# Check the quality (Quality Score - Phred)
plotQualityProfile(fnFs[1:6])
plotQualityProfile(fnRs[1:6])

# Filter forward and reverse reads based on parameters determined from the Phred Score
filt_path <- file.path(path, "filteredData") # Place filtered files in a subdirectory
filtFs <- file.path(filt_path, paste0(sampleNames, 1:6, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sampleNames, 1:6, "_R_filt.fastq.gz"))
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(220,210),
                     maxN=0, maxEE=c(2,2), trimLeft=c(19,20), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)
head(out)
# truncLen -> Truncate reads after truncLen bases. Shorter sequences are discarded
# rm.phix -> Remove bases that do not match A, T, C, or G
# compress -> Filtered files will be compressed into gzip
# multithread -> Allows parallel tasks

# Reviewing quality, now with filtered data
plotQualityProfile(filtFs[1:6])
plotQualityProfile(filtRs[1:6])

# View sequencing error rates
errF <- learnErrors(filtFs, multithread=FALSE)
# 62565270 total bases in 311270 reads from 6 samples will be used for learning the error rates
errR <- learnErrors(filtRs, multithread=FALSE)
# 59141300 total bases in 311270 reads from 6 samples will be used for learning the error rates.

# Visualize sequencing error rates by bases
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

# Dereplication (It is the identification of unique sequences, so that only one copy of each sequence is considered)
# Allows a significant reduction in computer memory requirements
derepFs <- derepFastq(filtFs)
derepRs <- derepFastq(filtRs)

# Identify dereplicated samples
names(derepFs) <- sampleNames
names(derepRs) <- sampleNames

# Inference on samples - how many unique sequences exist in each sample?
dadaFs <- dada(derepFs, err=errF, multithread=FALSE)
#Sample 1 - 55786 reads in 25387 unique sequences.
#Sample 2 - 47722 reads in 21020 unique sequences.
#Sample 3 - 60818 reads in 24824 unique sequences.
#Sample 4 - 68871 reads in 37305 unique sequences.
#Sample 5 - 31178 reads in 19172 unique sequences.
#Sample 6 - 46895 reads in 27997 unique sequences.

dadaRs <- dada(derepRs, err=errR, multithread=FALSE)
#Sample 1 - 55786 reads in 28997 unique sequences.
#Sample 2 - 47722 reads in 25155 unique sequences.
#Sample 3 - 60818 reads in 31351 unique sequences.
#Sample 4 - 68871 reads in 42366 unique sequences.
#Sample 5 - 31178 reads in 19944 unique sequences.
#Sample 6 - 46895 reads in 28863 unique sequences.

# Merge paired reads - combine forward and reverse sequences
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
#45782 paired-reads (in 557 unique pairings) successfully merged out of 52871 (in 2356 pairings) input.
#38333 paired-reads (in 484 unique pairings) successfully merged out of 45017 (in 2000 pairings) input.
#50321 paired-reads (in 616 unique pairings) successfully merged out of 57908 (in 2545 pairings) input.
#48824 paired-reads (in 1050 unique pairings) successfully merged out of 62683 (in 4408 pairings) input.
#21239 paired-reads (in 498 unique pairings) successfully merged out of 27512 (in 1937 pairings) input.
#33108 paired-reads (in 739 unique pairings) successfully merged out of 42009 (in 2924 pairings) input.

# Build the sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
# Samples ASVs
# 6        2179

# Inspect the distribution of sequence lengths
table(nchar(getSequences(seqtab)))
# Most sequences have 253bp (1954 sequences)
# 238  240  248  252  253  254  255  256  257  263  264  271 
# 1    1    1    38   1954 164  14   1    1    1    1    2

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=FALSE, verbose=TRUE)
# Identified 54 bimeras out of 2179 input sequences.

dim(seqtab.nochim)
# 6 2125

# method="consensus" -> Samples are checked independently for chimeras
sum(seqtab.nochim)/sum(seqtab)
# 0.981663 (less than 2% chimeras - its ok!)

# Now inspect the distribution of sequence lengths with chimeras removed
table(nchar(getSequences(seqtab.nochim)))
# 238  240  248  252  253  254  255  256  257  263  264  271 
# 1    1    1    38   1905 160  13  1    1    1    1    2

# Track sequence quantities throughout the pipeline - important to see the final number of sequences
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sampleNames
head(track)
#     input filtered denoisedF denoisedR merged nonchim
# A1 108005    55786     54220     53897  45782   44818
# A2  90233    47722     46110     46064  38333   37368
# A3 117013    60818     59185     59044  50321   49122
# A7 133600    68871     65140     65104  48824   48070
# A8  60867    31178     28789     28872  21239   21029
# A9  93554    46895     43851     43844  33108   32843

write.table(track, "results/track.txt", sep="\t", quote=F)

# Access taxonomy
# Download can be done here: https://zenodo.org/record/4587955
taxa <- assignTaxonomy(seqtab.nochim, "C:/Users/PC/Desktop/bioinfo_microbiana/silva_nr99_v138.1_train_set.fa", multithread=FALSE)
write.table(taxa, "results/taxa.txt", sep="\t", quote=F)

# Extract sequences for easy visualization
taxa.print <- taxa
rownames(taxa.print) <- NULL
head(taxa.print)
#       Kingdom    Phylum           Class                 Order               Family               Genus
# [1,] "Bacteria" "Proteobacteria" "Alphaproteobacteria" "Rhizobiales"       "Xanthobacteraceae"  NA   
# [2,] "Bacteria" "Proteobacteria" "Alphaproteobacteria" "Rhizobiales"       "Xanthobacteraceae"  NA   
# [3,] "Archaea"  "Crenarchaeota"  "Nitrososphaeria"     "Nitrososphaerales" "Nitrososphaeraceae" NA   
# [4,] "Bacteria" "Proteobacteria" "Alphaproteobacteria" "Rhizobiales"       "Xanthobacteraceae"  NA   
# [5,] "Archaea"  "Crenarchaeota"  "Nitrososphaeria"     "Nitrososphaerales" "Nitrososphaeraceae" NA   
# [6,] "Archaea"  "Crenarchaeota"  "Nitrososphaeria"     "Nitrososphaerales" "Nitrososphaeraceae" NA

# Simplify headers with names (ASV_1, ASV_2...)
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")
for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

# Create two tables with identifications and quantifications of ASVs
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "ASVs.fa")

# count table
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "results/ASVs_counts.txt", sep="\t", quote=F)

#tax table
asv_tax <- taxa
row.names(asv_tax) <- sub(">", "", asv_headers)
write.table(asv_tax, "results/ASVs_taxonomy.txt", sep="\t", quote=F)

# Organize tables for analysis
counts <- read.table("results/ASVs_counts.txt", header = TRUE, sep = "\t", dec = ".")
taxonomy <- read.table("results/ASVs_taxonomy.txt", header = TRUE, sep = "\t", dec = ".")

# Merge tables
counts_taxonomy <- as.data.frame(c(taxonomy, counts)) #concat
write.table(counts_taxonomy, "results/counts_taxonomy.txt", sep="\t", quote=F)

# Bonus - rarefaction curve
install.packages("vegan")
library(vegan); packageVersion("vegan")

rarefaction <- t(counts)
S <- specnumber(rarefaction)
raremax <- min(rowSums(rarefaction))
rarecurve(rarefaction, step = 1, sample = raremax, col = "blue", cex = 0.6)