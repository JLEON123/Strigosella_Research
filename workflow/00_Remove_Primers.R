# -----------------------------------------------------------------------------#
# Meta-amplicon analysis recipe
# Processing raw amplicon reads
# Author: Geoffrey Zahn
# Edited: Josh Leon
# Software versions:  R v 4.1.1
#                     tidyverse v 1.3.2
#                     dada2 v 1.20.0
#                     decontam v 1.12.0
#                     phyloseq v 1.36.0
#                     purrr v 0.3.5
#                     Biostrings v 2.60.2
#                     patchwork v 1.1.2
#                     ShortRead v 1.50.0
# -----------------------------------------------------------------------------#

#############################################################
#### This script calls cutadapt to remove any primers    ####
#### You must have cutadapt installed on your system     ####
#### and present in your PATH. See cutadapt installation ####
#### documents for instructions.                         ####
#############################################################

# PACKAGES, SCRIPTS, AND SETUP ####
library(tidyverse); packageVersion("tidyverse")
library(dada2); packageVersion("dada2")
library(purrr); packageVersion("purrr")
library(Biostrings); packageVersion("Biostrings")
library(ShortRead); packageVersion("ShortRead")

# PARSE FILE PATHS ####

# File parsing
path <- "../Sequence_data/Strigosella_raw_reads" # CHANGE to the subdirectory containing your demultiplexed fastq files
filtpath <- "./Pre_filtered_reads" # CHANGE to the subdirectory where you want your cutadapted reads to live


# your filenames might have a different pattern for determining FWD and REV reads
# Change "pattern" to accomodate your file names
fnFs <- sort(list.files(path, pattern = "R1", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "R2", full.names = TRUE))


# CHECK FOR AND REMOVE PRIMER SITES WITH CUTADAPT ####

############# You will need to change these two values to match your data #############

# Here, you should supply the primer sequences used during PCR
# The example data in this workflow was generated using the 515f - 806r primer pair
# to amplify the V4 region of bacterial 16S rDNA
FWD <- "CCTACGGGNGGCWGCAG" # Sequence of FWD primer
REV <- "GACTACNVGGGTMTCTAATCC"  # Sequence of REV primer
######################################################################################################

allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients; REV.orients

# Prefilter to remove reads with ambiguous (N) bases ####
fnFs.filtN <- file.path(filtpath, basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- file.path(filtpath, basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE) # on Windows, set multithread = FALSE

# Discover primer matches, regardless of orientation ####
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))

# Run cutadapt ####
# If the following command returns an error, you do not have cutadapt installed correctly
system2("cutadapt", args = "--version")

path.cut <- file.path("./cut_adapt_reads")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC) 
# Run Cutadapt
for(i in seq_along(fnFs)) {
  system2("cutadapt", args = c(R1.flags, R2.flags, "-n", 2, "--minimum-length 100", # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i])) # input files
}

# sanity check
# This should show no occurences in any of the orientations now
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))
