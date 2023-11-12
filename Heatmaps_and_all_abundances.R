library(phyloseq)
library(tidyverse)

theme_set(theme_bw())

ps <- readRDS("./Final_phyloseq_object.RDS")
ps_actino <- subset_taxa(ps, Phylum == "Actinobacteria")

ps_act_merged <- merge_samples(ps_actino, "Structue")
ps_actino_merged.prop <- transform_sample_counts(ps_actino_merged, function(x) x / sum(x) )
