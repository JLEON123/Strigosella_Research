# -----------------------------------------------------------------------------#
# Meta-amplicon analysis recipe
# Exploring cleaned data using phyloseq and corncob packages
# Author: Geoffrey Zahn
# Edited: Josh Leon
# Software versions:  R v 4.1.1
#                     tidyverse v 1.3.2
#                     phyloseq v 1.36.0
#                     purrr v 0.3.5
#                     Biostrings v 2.60.2
#                     corncob v 0.3.0
#                     vegan v 2.6.4
#                     patchwork v 1.1.2
# -----------------------------------------------------------------------------#

# PACKAGES, SCRIPTS, AND SETUP ####
library(tidyverse); packageVersion("tidyverse")
library(phyloseq); packageVersion("phyloseq")
library(purrr); packageVersion("purrr")
library(Biostrings); packageVersion("Biostrings")
library(corncob); packageVersion("corncob")
library(vegan); packageVersion("vegan")
library(patchwork); packageVersion("patchwork")
library(janitor) ; packageVersion("janitor")

source("./scripts/bbdml_helper.R")
source("./scripts/theme2.R")
#################################################################################
#                               Main workflow                                   #
#  Explore alpha and beta diversity, visualize data set, test hypotheses,       #
#  search for differentially abundant taxa                                      #
#                                                                               #
#################################################################################


# Load cleaned phyloseq object ####
ps <- readRDS("./Final_phyloseq_object.RDS")

# merge taxa at genus level
ps_genus <- ps %>% 
  tax_glom(taxrank = "Genus")

# Clean up ASV names to show taxonomy
ASV_names <- otu_table(ps) %>% colnames()
ASV_taxa <- otu_to_taxonomy(ASV_names,ps,level = c("Phylum","Class","Order","Family","Genus"))

genus_names <- otu_table(ps_genus) %>% colnames()
genus_taxa <- otu_to_taxonomy(genus_names,ps_genus,level = c("Phylum","Class","Order","Family","Genus"))


# Getting to know your phyloseq data ####

# number of taxa
ntaxa(ps)

# number of samples
nsamples(ps)

# sample names
sample_names(ps)
rank_names(ps)
# taxa names
taxa_names(ps)

# overview of taxonomy
tax_table(ps)
tax_table(ps)[,1] %>% table()
tax_table(ps)[,2] %>% table()
tax_table(ps)[,3] %>% table()
tax_table(ps)[,4] %>% table()
tax_table(ps)[,5] %>% table()
tax_table(ps)[,6] %>% table()

# sample metadata
sample_data(ps) %>% as_tibble()

# ASV table
otu_table(ps)

# how many sequences observed in each sample?
otu_table(ps) %>% rowSums()

# how many times was each taxon observed?
otu_table(ps) %>% colSums()

# how many different samples was each taxon found in?
asv <- otu_table(ps) %>% as("matrix") %>% as.data.frame() # convert to matrix before you can convert to data frame
asv[asv>0] <- 1 # convert to presence/absence
colSums(asv) # sum of presence (present = 1; absent = 0)

# what was the most widespread taxon (not abundance)
widespread_seq <- names(asv)[which(colSums(asv) == max(colSums(asv)))] # this gives long sequence
tax_table(ps)[widespread_seq,] # this pull that row from ASV table

# what was most abundant (raw counts) taxon?
abund_seq <- which(otu_table(ps) %>% colSums() == max(otu_table(ps) %>% colSums()))
tax_table(ps)[abund_seq,]
otu_table(ps)[,abund_seq]

# access the phylogenetic tree
phy_tree(ps)
plot_tree(ps,color="Phylum")
# this tree needs some refining, but it's will work for this tutorial



# Alpha diversity metrics ####
alpha <- estimate_richness(ps) %>% 
  select(Observed,Shannon, Simpson)

# since we have no singletons (consequence of DADA2), we cannot use Chao1 estimate!
# plot alpha diversity for every sample
plot_richness(ps, 
              measures = c("Observed","Shannon","Simpson"), 
              color = "Structure", 
              sortby = "Observed") +
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(face = 'bold', size = 20, color = "#000000"),
        axis.title.y = element_text(face = 'bold', size = 25, color = "#000000"),
        strip.background = element_rect(fill = "#000000"), 
        strip.text = element_text(color = "#FFFFFF", size = 20),
        legend.text = element_text(face = "bold", size = 20),
        legend.title = element_text(face = "bold", size = 25),
        legend.key = element_rect(size = 25),
        legend.title.align = 0.5) +
  scale_discrete_manual(values = pal,
                        aesthetics = 'color') 

ggsave("./graphs/alpha_diversity_point_plot_structure.png",dpi=300, width = 12, units = "in")
# plot, grouped by colony color with added boxplot
plot_richness(ps, 
              x = "Structure",
              measures = c("Observed","Shannon","Simpson"),
              color = "Structure") + 
  geom_boxplot(aes(color = Structure, fill = Structure), 
               alpha = 0.5, show.legend = FALSE) +
  geom_point(size = 3, show.legend = FALSE) +
  scale_discrete_manual(values = pal,
                        aesthetics = c('color', "fill")) +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        strip.background = element_rect(fill = "#000000"), 
        strip.text = element_text(color = "#FFFFFF", size = 20), 
        axis.text.x = element_text(angle = .5, hjust = .5),
        legend.background = element_blank(),
        legend.text = element_blank(),
        legend.title = element_blank(),
        legend.key = element_blank())

  
ggsave("./graphs/alpha_diversity_boxplot_structure.png",dpi=300,width = 12, unit = "in")

# Add data to alpha date frame
alpha$Sample <- row.names(alpha)
alpha$Structure <- ps@sam_data$Structure
alpha$Location <- ps@sam_data$Location

# model alpha
mod_observed <- glm(data=alpha,
                    formula = Observed ~ Location * Structure)

sink("./cache/alpha_mod_observed_summary.txt")
summary(mod_observed)
sink(NULL)

sink("./cache/alpha_mod_observed_summary.txt")
broom::tidy(mod_observed) %>% 
  filter(p.value<0.05)
sink(NULL)


# glm shannon
mod_shannon <- glm(data=alpha,
                   formula = Shannon ~ Location * Structure)
broom::tidy(mod_shannon) %>% 
  saveRDS("./cache/shannon_div_df.RDS")

sink("./cache/alpha_mod_shannon_summary.txt")
broom::tidy(mod_shannon)
sink(NULL)
install.packages("report")
report::report(mod_observed)

# transform raw counts to relative abundance ####
ps_ra <- transform_sample_counts(ps, fun = function(x){x/sum(x)})

# Beta-diversity ####

# Ordination
ordcap <-  ordinate(ps_ra, "CAP", "bray", ~Structure)


(
  ord1 <- plot_ordination(ps_ra,ordcap,color = "Structure", shape = "Location") +
    geom_point(size=4)  +
    theme(legend.position = "top",
          legend.key = element_rect(size = 25), 
          axis.line = element_line(color = "#000000")) +
    #stat_ellipse(show.legend = FALSE) +
    scale_discrete_manual(values = pal,
                          aesthetics = 'color')
)

ggsave("./graphs/CAP_Ordination_with_location.png",dpi=500,width = 12, unit = "in")
# try another ordination method
nmds <- ordinate(ps_ra,method = "NMDS")

(
ord2 <- plot_ordination(ps_ra,nmds,color = "Structure") +
    geom_point(size=4)  + 
    theme(legend.position = "none") +
    stat_ellipse() +
    labs(title = "NMDS - Bray") +
    scale_discrete_manual(values = pal,
                          aesthetics = 'color')
)


# also try with unifrac distance, which takes phylogeny into account
unifrac.dist <- UniFrac(ps_ra)
unifrac <- ordinate(ps_ra,method = "NMDS",distance = unifrac.dist)

(
ord3 <- plot_ordination(ps_ra,unifrac,color = "Structure") +
    geom_point(size=4) +
    theme(legend.position = "none") +
    stat_ellipse() +
    labs(title = "NMDS - Weighted Unifrac") +
    scale_discrete_manual(values = pal,
                          aesthetics = 'color')
)

# combine all plots into one figure for comparison
ord1 / ord2 / ord3
ggsave("./graphs/ordinations_by_structure.png",dpi=500,width = 12, unit = "in")

# permanova ####
# pull out components
asv <- otu_table(ps_ra) %>% as("matrix") %>% as.data.frame()
meta <- sample_data(ps_ra) %>% as.data.frame()



# run permanova model with Structure and Location as predictors (with interaction term included)
permanova.bray <- 
  vegan::adonis2(asv ~ meta$Structure * meta$Location ,method = "bray")%>% 
  broom::tidy()
permanova.bray


# save output to a file
sink("./cache/permANOVA_Bray_structure_Table.txt")
permanova.bray
sink(NULL)

# try with jaccard distance as well
permanova.jaccard <- 
  vegan::adonis2(asv ~ meta$Structure * meta$Location, method = "jaccard")%>% 
  broom::tidy()
permanova.jaccard

sink("./cache/permANOVA_structure_Table_jaccard.txt")
permanova.jaccard
sink(NULL)

# Differential abundance/dispersion tests ####

# Here, we use the corncob package to test for taxa that have differential abundance/dispersion between groups


# use non-transformed data!
set.seed(123)
da_analysis_Location <- differentialTest(formula = ~ Location, #abundance
                                         phi.formula = ~ 1, #dispersion
                                         formula_null = ~ 1, #mean
                                         phi.formula_null = ~ 1,
                                         test = "Wald", 
                                         boot = FALSE,
                                         data = ps_genus,
                                         fdr_cutoff = 0.05,
                                         full_output = TRUE)

plot(da_analysis_Location) +
  theme(legend.position = 'none')

# Get model info
mods <- da_analysis_Location$significant_models

capture_mods <- function(x){
  y <- x$coefficients %>% 
    as.data.frame() %>%
    janitor::clean_names()
  names(y) <- c("estimate","std_error","t_value","p_value")
  y <- y %>% 
    mutate(p_value = p_value %>% round(6))
  y <- y %>% 
    filter(row.names(y) %>% grepl(pattern="Location"))
  return(y)
}
mods[[1]]$coefficients
bbdml_mods <- map(mods,capture_mods)
# find the significant taxa
da_analysis_Location$significant_taxa
sig_taxa <- da_analysis_Location$significant_taxa %>% otu_to_taxonomy(data=ps_genus)

names(bbdml_mods) <- sig_taxa
joined_mods <- bbdml_mods %>% 
  purrr::reduce(full_join)
joined_mods$taxon <- names(bbdml_mods)
joined_mods <- joined_mods %>% 
  select(taxon,estimate,std_error,t_value,p_value)

joined_mods %>%
  saveRDS("./data/bbdml_significant_mod_tables.RDS")
# This is a helper function for plotting corncob results. It's found in "scripts/bbdml_helper.R" 
bbdml_obj <- multi_bbdml(da_analysis_Location,
                         ps_object = ps_genus,
                         mu_predictor = "Location",
                         phi_predictor = "Location",
                         taxlevels = 6)
multi
# # filter significant taxa to only those with absolute effect sizes > 1.5
# is that unit log odds???

find_mu <- function(x,mu=1.5){
  y <- x$b.mu[2]
  z <- abs(y) > 1.5
  print(z)
}

new_bbdml_obj <- bbdml_obj[map(bbdml_obj, find_mu) %>% unlist]
names(new_bbdml_obj)


# INDICSPECIES ####
comm <- ps_genus %>% 
  transform_sample_counts(function(x){x/sum(x)}) %>% 
  otu_table() %>% as("matrix")
base::colnames(comm) <- genus_taxa

groups <- ps_genus@sam_data$Location

# multipatt analysis
BiocManager::install("indicspecies")
indval <- indicspecies::multipatt(comm,
                                  groups,
                                  control = how(nperm=999))

summary(indval)

# pull genus names of taxa that indicate EITHER CARE or MC
x <- row.names(indval$A)
indicspecies_taxa <- x[which(indval$sign$p.value < 0.05)] %>% 
  str_split("_") %>% 
  map_chr(5)
# RANDOM FOREST ####
# use relative abundance transformations for Ranger
ps_genus <- ps_genus %>% 
  transform_sample_counts(function(x){x/sum(x)})

# get data ready for RF modeling
meta <- microbiome::meta(ps_genus) %>% 
  select(Location) %>% 
  mutate(Cap_Reef = case_when(Location == "Capitol Reef" ~ TRUE, # east is logical respons
                              Location== "Manning Canyon" ~ FALSE)) %>% 
  select(Cap_Reef)


asv <- otu_table(ps_genus) %>% 
  as("matrix") %>% 
  as.data.frame()
df <- 
  meta %>% 
  bind_cols(asv)

# train RF model on full data set
BiocManager::install("vip")
ranger_model <- ranger::ranger(Cap_Reef~., 
                               data = df, 
                               classification = TRUE, 
                               probability = TRUE,
                               importance = 'permutation')
# find top important taxa
top <- vip::vip(ranger_model,num_features=20) # find most important factors for success (survival)
top +
  theme(axis.text.y = element_blank())
vip_taxa <- corncob::otu_to_taxonomy(top$data$Variable,data = ps_genus) %>% unname()
vip_taxa

pred <- predict(ranger_model,df) # how does it do predicting itself (no cross-validation)
preds <- pred$predictions %>% 
  as.data.frame() %>% 
  bind_cols(df$Cap_Reef)
names(preds) <- c("prob_T","prob_F","Capitol_Reef")

preds %>% 
  ggplot(aes(x=prob_T,fill=Capitol_Reef)) +
  geom_density()
vip_taxa %>% 
  saveRDS("./data/ranger_vip_taxa.RDS")


# COMMON DIFFABUND TAXA ####
# find taxa that were detected by all methods
corncob_taxa <- sig_taxa %>% 
  str_split("_") %>% 
  map_chr(6)

# subset bbdml tests to just those found by multipatt as well
new_bbdml_obj <- new_bbdml_obj[which(names(new_bbdml_obj) %in% indicspecies_taxa)]

# also subset bbdml tests to those found by Ranger as well
vip_taxa <- vip_taxa %>% 
  str_split("_") %>% 
  map_chr(6)
new_bbdml_obj <- new_bbdml_obj[which(names(new_bbdml_obj) %in% vip_taxa)]

new_bbdml_obj %>% names() %>% 
  saveRDS("./data/final_significant_taxa.RDS")
new_bbdml_obj %>% 
  saveRDS("./data/final_significant_bbdml_list.RDS")

# another helper function found in the same file
plot_multi_bbdml(new_bbdml_obj,
                 color="Location", 
                 pointsize = 3)


# save plots as RDS objects
plots <- ls(pattern = "bbdml_plot_")
for(i in plots){
  x <- base::get(i)
  saveRDS(x,file=file.path("./cache",paste0(i,".RDS")))
}


# PLOT ####
# stick all plots together
plots <- ls(pattern = "^bbdml_plot_")

p1 <- bbdml_plot_1 +
  labs(color="Location") +
  theme(legend.position = 'none',
        axis.title.y = element_blank())
p2 <- bbdml_plot_2 +
  labs(color="Location") +
  theme(axis.title.y = element_blank(),
        legend.position = 'none')
p3 <- bbdml_plot_3 +
  labs(color="Location",
       y = "Relative Abundance") +
  theme(legend.position = 'none',
        axis.title.y = element_text(size = 15, face = "bold"))
p4 <- bbdml_plot_4 +
  labs(color="Location") +
  theme(axis.title.y = element_blank(),
        legend.position = 'none')
p5 <- bbdml_plot_5 +
  labs(color="Location") +
  theme(legend.position = 'none',
        axis.title.y = element_blank())
p6 <- bbdml_plot_6 +
  labs(color="Location") +
  theme(axis.title.y = element_blank(),
        legend.position = 'none') 





(p1 + p2) / (p3 + p4) / (p5 + p6) +
  plot_layout(guides = "collect") & 
  theme(legend.position = 'bottom', 
        legend.text = element_text(face = "bold", size = 15),
        axis.text.y = element_text(size = 15),
        legend.key.size = unit(1, 'in'), 
        legend.title = element_blank(), 
        plot.title = element_text(face = "bold", size = 20))

ggsave("./graphs/differential_abundance_sig_taxa.png", height = 8, width = 10,dpi=300)




