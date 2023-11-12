# -----------------------------------------------------------------------------#
# Meta-amplicon analysis recipe
# Cleaning up the phyloseq object
# Author: Geoffrey Zahn
# Software versions:  R v 3.6.3
#                     tidyverse v 1.3.0
#                     phyloseq v 1.30.0
#                     ShortRead v 1.44.3
#                     Biostrings v 2.54.0
# -----------------------------------------------------------------------------#

#################################################################################
#                               Main workflow                                   #
#        Remove non-bacteria, chloroplast and mitochondrial sequences           #
#                                                                               #
#################################################################################

# PACKAGES, SCRIPTS, AND SETUP ####
library(tidyverse); packageVersion("tidyverse")
library(phyloseq); packageVersion("phyloseq")
library(ShortRead); packageVersion("ShortRead")
library(Biostrings); packageVersion("Biostrings")

# load phyloseq object with phylogenetic tree ####
ps <- readRDS("./cache/ps_not-cleaned.RDS") # change to non-phylogeny stuff

# Find and remove non-bacteria ####
ps_nonbact <- subset_taxa(ps, Kingdom != "Bacteria")

# quick plot to look at kingdom-level taxonomy
ps %>% transform_sample_counts(function(x){x/sum(x)}) %>%
  plot_bar(fill="Kingdom")
ggsave("./graphs/16S_Kingdom-Level_Taxonomic_Proportions.png",dpi=500) # save figure for later use

# same plot, but non-bacteria, for sanity check
ps_nonbact %>% 
  transform_sample_counts(function(x){x/sum(x)}) %>%
  plot_bar(fill="Kingdom")


# REMOVE NON-BACTERIA, CHLOROPLASTS, MITOCHONDRIA, and empty samples/taxa ####
# 8903 ASVs to start
ps <- subset_taxa(ps, Kingdom == "Bacteria")
tax <- tax_table(ps)

# 8871 ASVs after removing non-bacterial sequences (32 Archaea ASVs)

ps <- subset_taxa(ps,Class != "Chloroplast")

# 7640 ASVs after removing Chloroplast sequences (1231 Chloroplast ASVs)

ps <- subset_taxa(ps, taxa_sums(ps) > 0) #7058
ps <- subset_samples(ps, sample_sums(ps) > 0)

# 7058 ASVs after removing ASVs with 0 reads (582); No empty samples


# Save DNA sequences apart from rownames (from subsetted ps object)
seqs <- taxa_names(ps)
seqs <- DNAStringSet(seqs)
saveRDS(seqs,"./cache/16S_ASV_reference_sequences.RDS")

# Save RDS object for cleaned up Phyloseq object
saveRDS(ps, file = "./cache/clean_phyloseq_object.RDS")

# Changing the sam data because the original is dumb
current_sam_data <- sample_data(ps)
Good_sam_data <- readRDS("./data/True_Meta.RDS")
Good_sam_data <- sample_data(Good_sam_data)
# Need to change the otu table to match new sam data
otu <- otu_table(ps)
row.names(otu) <- row.names(Good_sam_data)
identical(row.names(Good_sam_data),row.names(otu))

tax <- tax_table(ps)

ps <- phyloseq(otu, tax, Good_sam_data)

saveRDS(ps, file = "./Final_phyloseq_object.RDS")


