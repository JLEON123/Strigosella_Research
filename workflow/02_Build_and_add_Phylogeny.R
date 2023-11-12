# -----------------------------------------------------------------------------#
# Meta-amplicon analysis recipe
# Building and adding a phylogeny to the cleaned phyloseq object
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
#                     DECIPHER v 2.20.0
#                     vegan v 2.6.4
#                     phangorn v 2.10.0
#                     msa v 1.24.0
#                     ape v 5.6.2
#                     seqinr v 4.2.16
# -----------------------------------------------------------------------------#

#################################################################################
#                               Main workflow                                   #
# Perform multiple sequence alignment of all ASVs, build distance matrix,       # 
# construct and refine a phylogenetic tree, add the tree to the phyloseq object #
#           With larger data sets, this can be a long process...                #
# Further, proper phylogenetics is beyond the scope of this tutorial.           #
#################################################################################

# Packages and functions ####
library(tidyverse); packageVersion("tidyverse")
library(phyloseq); packageVersion("phyloseq")
library(vegan); packageVersion("vegan")
library(phangorn); packageVersion("phangorn")
library(msa); packageVersion("msa")
library(ape); packageVersion("ape")
library(seqinr); packageVersion("seqinr")
library(DECIPHER); packageVersion("DECIPHER")
# Read in phyloseq object from first script output ####
ps <- readRDS("./Final_phyloseq_object.RDS")
sample_data(ps)

# tax_table <- as.data.frame(ps@tax_table)
# simplify ASV names
seqs <- rownames(tax_table(ps))
names(seqs) <- paste0("ASV_",1:length(seqs)) # This propagates to the tip labels of the tree

# Multiple sequence alignment  ####
# With msa, takes extremely long (my data took ~43 hours)
#alignment <- msa(seqs,method = "Muscle", type = "dna",verbose = TRUE,order = "input",maxiters = 10)
# save progress 
#saveRDS(alignment,"./cache/16S_dna_alignment_muscle.RDS")
msa_align <- readRDS("./cache/16S_dna_alignment_muscle.RDS")

# With Decipher
# Create function to convert seqs to DNAStrings
To_DNA <- function(seq) { 
  return (DNAString(seq))
}
seqs_DEC <- lapply(seqs, To_DNA)
# Convert the list to a DNAstringSet
seqs_DEC <- DNAStringSet(seqs_DEC)
# nucleotide sequences need to be in the same orientation
# if they are not, then they can be reoriented (optional)
seqs_DEC <- OrientNucleotides(seqs_DEC)
# Run the alignment
#Took < 1 hour with my data
Alignment_DEC <- AlignSeqs(myXStringSet = seqs_DEC,iterations = 10, verbose = TRUE)
# view the alignment in a browser (optional)
#BrowseSeqs(Alignment_DEC, highlight=0)

# write the alignment to a new FASTA file
writeXStringSet(Alignment_DEC,
                file="./cache/16S_DECIPHER_Alignment.fasta")

# Convert to DNAbin
Alignment_DEC <- as.DNAbin(Alignment_DEC)
# Convert to phangorn format
phang.align = as.phyDat(Alignment_DEC)


# Model testing
mt <- modelTest(phang.align)
mt$Model

# distance - maximum likelihood ####
dm <- dist.ml(phang.align,model = "JC69")


# save progress
saveRDS(dm,"./cache/16S_ML_Distance.RDS")

# Initial neighbor-joining and upgma trees ####
treeNJ <- NJ(dm) # Note, tip order != sequence order
treeUPGMA <- upgma(dm)

# find parsimony scores for each tree
parsimony(treeNJ, phang.align)
parsimony(treeUPGMA, phang.align)

# search through treespace with NNI and SPR algorithm
optim <- optim.parsimony(treeNJ,phang.align)


# save progress
saveRDS(treeNJ, "./cache/16S_treeNJ.RDS")
saveRDS(optim, "./cache/16S_tree_optim.RDS")



# Maximum likelihood ####
fit = pml(treeNJ, data=phang.align)
fitJC <- optim.pml(fit, model = "JC", rearrangement = "stochastic")


# edit tip labels
fit$tree$tip.label <- seqs


# save trees
saveRDS(fit,"./cache/16S_fit_treeNJ.RDS")



# add tree to phyloseq object ####
ps2 <- phyloseq(tax_table(tax_table(ps)),
                otu_table(otu_table(ps)),
                sample_data(sample_data(ps)),
                phy_tree(fit$tree))


# Save updated phyloseq object with tree
saveRDS(ps2, "./cache/ps_not-cleaned_w_tree.RDS")



