
## For the webpage

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Genomics pipeline tutorial
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Author: Laura Twomey
# Date: January 2023
# Version: 1.0

# Setting up environment ===================================================

# Clean environment
rm(list = ls(all.names = TRUE)) # will clear all objects including hidden objects
gc() # free up memory and report the memory usage
options(max.print = .Machine$integer.max, scipen = 999, stringsAsFactors = F, dplyr.summarise.inform = F) # avoid truncated output in R console and scientific notation

# Loading relevant libraries 
library(tidyverse) # includes ggplot2, for data visualisation. dplyr, for data manipulation.
library(Rqc) # for QC 
library(pander) # to open HTML files
library(QuasR) # for filtering and trimming reads
library(BiocIO) # for import function

# Set input path
path <- "~/Biostatsquid/Scripts/Genomics/"
setwd(path)

list.files(path)
set.seed(42)

#sample data folder
folder <- (system.file(package = "QuasR", "extdata"))
dir(folder) # will show the contents of the folder

# Import data ===================================================
# Copy the chip files to our working directory
filenames <- list.files(folder)[grepl("chip|hg19sub.fa", list.files(folder))]
file.copy(paste(folder, filenames, sep = '/'), ".", recursive = TRUE)
chip_filenames <- filenames[grepl('^chip.*bz2$', filenames)]
list.files(path)

qcRes <- rqc(path = path, # input path
             pattern = "^chip.*bz2$", # files that start with chip and end with bz2
             openBrowser = FALSE, 
             outdir = path) # output path

# QC ============================================================

## 1. Sequence quality per base/cycle ===========================
rqcCycleQualityBoxPlot(qcRes)

## 2. Sequence content per base/cycle ===========================
rqcCycleBaseCallsLinePlot(qcRes)

## 3. Read frequency plot ===========================
rqcReadFrequencyPlot(qcRes)

## 4. Other QC metrics ===========================================
rqc_report <- rqcReport(qcRes, outdir = path, file = "rqc_report")
openFileInOS(rqc_report)

# Filtering and trimming reads ============================================================
fastqFiles <- paste(folder, chip_filenames, sep = '/')

# defined processed fastq file names
outfiles <- paste(tempfile(pattern = c("processed_1_", "processed_2_"), 
                           tmpdir = gsub('\\/$', '', path)), ".fastq", sep = "")

# process fastq files
preprocessReads(fastqFiles, outfiles, 
                nBases = 1, # remove reads that have more than 1 N
                truncateEndBases = 3, # trim 3 bases from the end of the reads 
                minLength = 10) # remove reads shorter than 10 base-pairs 

# We can also use the ShortRead package to filter reads in ways that are 
# not possible using the QuasR::preprocessReads() function. 
for(file in 1:length(fastqFiles)){
  
  #for debug
  #file <- 2
  fastqFile <- fastqFiles[file]
  
  # set up streaming with block size 1000
  f <- FastqStreamer(fastqFile, readerBlockSize = 1000)
  
  # we set up a while loop to call yield() function to go through the file
  # every time we call the yield() function 1000 read portion of the file will be read successively. 
  while(length(fq <- yield(f))) {
    
    # remove reads where all quality scores are < 20
    # get quality scores per base as a matrix
    qPerBase <- as(quality(fq), "matrix")
    
    # get number of bases per read that have Q score < 20
    qcount <- rowSums(qPerBase <= 20)
    
    # Number of reads where all Phred scores >= 20
    print(paste('# reads whith all Phred scores >= 20 for file', as.character(gsub('.*chip', 'chip', fastqFile))))
    print(fq[qcount == 0])
    
    writeFastq(fq[qcount == 0],
               paste(gsub('\\.fq.bz2', '', fastqFile), "Qfiltered", ".fq.bz2", sep = "_"),
               mode = "a")  # write fastq file with mode="a", so every new block is written out to the same file
  }
}


# Mapping/aligning reads to the genome ============================================================

# genome file in fasta format
genomeFile <- "hg19sub.fa"

# text file containing sample names and fastq file paths
#sampleFile <- "extdata/samples_chip_single.txt"
filenames <- list.files(path)[grepl('chip.*fq.bz2|processed', list.files(path))]

align_txt <- cbind(filenames, paste0('Sample', rep(1:length(filenames))))
colnames(align_txt) <- c('FileName', 'SampleName')
write.table(align_txt, 'align_txt.txt', col.names = T, quote = F, row.names = F, sep = "\t")


sampleFile <- "align_txt.txt"
proj <- qAlign(sampleFile, genomeFile)

# For each input sequence file, there will be one bam file with alignments 
# against the reference genome, and one for each auxiliary target sequence 
# with alignments of reads without genome hits.
list.files(pattern = ".bam$")
# Each alignment file is accompanied by two additional files with suffixes .bai and .txt:
list.files(pattern = "^chip_1_1_")[1:3]

# Plot alignment statistics
qQCReport(proj, pdfFilename = "qc_report.pdf")

# Alignment statistics
alignmentStats(proj)

# Export genome wig file from alignments
# For visualization in a genome browser, alignment coverage along the genome 
# can be exported to a (compressed) wig file using the qExportWig function. 
qExportWig(proj, binsize = 100L, scaling = TRUE, collapseBySample = TRUE)













# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Genomics pipeline tutorial
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Author: Laura Twomey
# Date: January 2023
# Version: 1.0

# Setting up environment ===================================================

# Clean environment
rm(list = ls(all.names = TRUE)) # will clear all objects including hidden objects
gc() # free up memory and report the memory usage
options(max.print = .Machine$integer.max, scipen = 999, stringsAsFactors = F, dplyr.summarise.inform = F) # avoid truncated output in R console and scientific notation

# Loading relevant libraries 
library(tidyverse) # includes ggplot2, for data visualisation. dplyr, for data manipulation.
library(Rqc) # for QC 
#library(fastqcr) # for other QC metrics - Java
library(pander) # to open HTML files
library(QuasR) # for filtering and trimming reads
#library(Gviz) # for gene alignment viz
library(BiocIO) # for import function

# Set input path
path <- "~/Biostatsquid/Scripts/Genomics/"
setwd(path)

list.files(path)
set.seed(42)

#sample data folder
folder <- (system.file(package = "QuasR", "extdata"))
dir(folder) # will show the contents of the folder

# Import data ===================================================
# feeds fastq.qz files found in path to quality check function
# qcRes <- rqc(path = system.file(package = "ShortRead", "extdata/E-MTAB-1147"), 
#           pattern = ".fastq.gz", 
#           openBrowser = FALSE)

# Copy the chip files to our working directory
filenames <- list.files(folder)[grepl("chip|hg19sub.fa", list.files(folder))]
file.copy(paste(folder, filenames, sep = '/'), ".", recursive = TRUE)
chip_filenames <- filenames[grepl('^chip.*bz2$', filenames)]
list.files(path)

qcRes <- rqc(path = path,
          pattern = "^chip.*bz2$",
          openBrowser = FALSE, 
          outdir = path)

# QC ============================================================

## 1. Sequence quality per base/cycle ===========================
rqcCycleQualityBoxPlot(qcRes)

## 2. Sequence content per base/cycle ===========================
rqcCycleBaseCallsLinePlot(qcRes)

## 3. Read frequency plot ===========================
rqcReadFrequencyPlot(qcRes)

## 4. Other QC metrics ===========================================

# # Only works if MacOS or Linux
# # install the FASTQC java tool
# fastqc_install()
# # call FASTQC and record the resulting statistics
# # in fastqc_results folder
# fastqc(fq.dir = folder,qc.dir = "fastqc_results")

rqc_report <- rqcReport(qcRes, outdir = path, file = "rqc_report")
openFileInOS(rqc_report)

# Filtering and trimming reads ============================================================
# obtain a list of fastq file paths
# fastqFiles <- system.file(package = "ShortRead",
#                           "extdata/E-MTAB-1147",
#                           c("ERR127302_1_subset.fastq.gz",
#                             "ERR127302_2_subset.fastq.gz")
# )

# fastqFiles <- system.file(package = "QuasR",
#                           "extdata/",
#                           dir(folder)[grepl('^chip', dir(folder))]
#                           ) 

fastqFiles <- paste(folder, chip_filenames, sep = '/')
  
# defined processed fastq file names
outfiles <- paste(tempfile(pattern = c("processed_1_", "processed_2_"), 
                           tmpdir = gsub('\\/$', '', path)), ".fastq", sep = "")

# process fastq files
preprocessReads(fastqFiles, outfiles, 
                nBases = 1, # remove reads that have more than 1 N
                truncateEndBases = 3, # trim 3 bases from the end of the reads 
                Lpattern = "ACCCGGGA", # Remove ACCCGGGA pattern if it occurs at the start
                minLength = 10) # remove reads shorter than 10 base-pairs 

# We can also use the ShortRead package to filter reads in ways that are 
# not possible using the QuasR::preprocessReads() function. 
# obtain a list of fastq file paths
# fastqFile <- system.file(package="ShortRead",
#                          "extdata/E-MTAB-1147",
#                          "ERR127302_1_subset.fastq.gz")
# 
# # read fastq file
# fq = readFastq(fastqFile)
# 
# # get quality scores per base as a matrix
# qPerBase = as(quality(fq), "matrix")
# 
# # get number of bases per read that have quality score below 20
# # we use this
# qcount = rowSums( qPerBase <= 20) 
# 
# # Number of reads where all Phred scores >= 20
# fq[qcount == 0]
# 
# # write out fastq file with only reads where all 
# # quality scores per base are above 20
# writeFastq(fq[qcount == 0], 
#            paste(fastqFile, "Qfiltered", sep="_")) 

# set up streaming with block size 1000
# every time we call the yield() function 1000 read portion
# of the file will be read successively. 

for(file in 1:length(fastqFiles)){

  #for debug
  #file <- 2
  fastqFile <- fastqFiles[file]

  f <- FastqStreamer(fastqFile, readerBlockSize = 1000)
  
  # we set up a while loop to call yield() function to go through the file
  while(length(fq <- yield(f))) {

    # remove reads where all quality scores are < 20
    # get quality scores per base as a matrix
    qPerBase <- as(quality(fq), "matrix")

    # get number of bases per read that have Q score < 20
    qcount <- rowSums(qPerBase <= 20)
    
    # Number of reads where all Phred scores >= 20
    print(paste('# reads whith all Phred scores >= 20 for file', as.character(gsub('.*chip', 'chip', fastqFile))))
    print(fq[qcount == 0])

    writeFastq(fq[qcount == 0],
               paste(gsub('\\.fq.bz2', '', fastqFile), "Qfiltered.fq.bz2", sep = "_"),
               mode = "a")  # write fastq file with mode="a", so every new block is written out to the same file
  }
}


# Mapping/aligning reads to the genome ============================================================

# genome file in fasta format
genomeFile <- "hg19sub.fa"

# text file containing sample names and fastq file paths
#sampleFile <- "extdata/samples_chip_single.txt"
filenames <- list.files(path)[grepl('chip.*fq.bz2|processed', list.files(path))]

align_txt <- cbind(filenames, paste0('Sample', rep(1:length(filenames))))
colnames(align_txt) <- c('FileName', 'SampleName')
write.table(align_txt, 'align_txt.txt', col.names = T, quote = F, row.names = F, sep = "\t")


sampleFile <- "align_txt.txt"
proj <- qAlign(sampleFile, genomeFile)

# For each input sequence file, there will be one bam file with alignments 
# against the reference genome, and one for each auxiliary target sequence 
# with alignments of reads without genome hits.
list.files(pattern = ".bam$")
# Each alignment file is accompanied by two additional files with suffixes .bai and .txt:
list.files(pattern = "^chip_1_1_")[1:3]

# Plot alignment statistics
qQCReport(proj, pdfFilename = "qc_report.pdf")

# Alignment statistics
alignmentStats(proj)

# Export genome wig file from alignments
# For visualization in a genome browser, alignment coverage along the genome 
# can be exported to a (compressed) wig file using the qExportWig function. 
qExportWig(proj, binsize = 100L, scaling = TRUE, collapseBySample = TRUE)








# Count alignements using gCount
annotFile <- "extdata/hg19sub_annotation.gtf"
chrLen <- scanFaIndex(genomeFile)
chrominfo <- data.frame(chrom = as.character(seqnames(chrLen)),
                        length = width(chrLen),
                        is_circular = rep(FALSE, length(chrLen)))
txdb <- makeTxDbFromGFF(file = annotFile, format = "gtf",
                        chrominfo = chrominfo,
                        dataSource = "Ensembl",
                        organism = "Homo sapiens")

gr1 <- import("Sample1.wig.gz")
gr2 <- import("Sample2.wig.gz")
axisTrack <- GenomeAxisTrack()
dTrack1 <- DataTrack(range = gr1, name = "Sample 1", type = "h")
dTrack2 <- DataTrack(range = gr2, name = "Sample 2", type = "h")
txTrack <- GeneRegionTrack(txdb, name = "Transcripts", showId = TRUE)
plotTracks(list(axisTrack, dTrack1, dTrack2, txTrack),
           chromosome = "chr3", extend.left = 1000)
