
### Codivergence analysis adapted following the method of Youngblut et al. 2019. Nat Comms ###

# Paul A. O'Brien
# E. paul.obrien@my.jcu.edu.au


# 01 Format data


## Set directory and load packages ----
setwd("~/Documents/R/Ramaciotti_16S_dataset/Codivergence_and_Neutral/Codivergence")

library(dplyr)
library(ape)
library(tidyr)
library(phyloseq)



## Import and format data ----
# Note: 16S reads pre-processed using QIIME2 (refer to phylosymbiosis scripts for details)
# QIIME2 artifacts include ASV table, taxonomy table and ASV phylogenetic tree

q_data <- "~/Documents/R/Ramaciotti_16S_dataset/Codivergence_and_Neutral/Data/Qiime2_artifacts/" # set folder

# Load asv table
asv <- read.table(paste0(q_data, "table-dada2-filtered.txt"), sep = '\t', row.names = 1, header = T, strip.white = T) %>% as.matrix()
colnames(asv) <- gsub("\\.", "-", colnames(asv)) # edit column names for downstream processing

# Load taxonomy
tax <- read.table(paste0(q_data, "taxonomy_silva.tsv"), sep = '\t', row.names = 1, header = T, strip.white = T)

# Load metadata
meta <- read.table(paste0(q_data, "metadata.tsv"), sep = '\t', row.names = 1, header = T, strip.white = T)
meta$SampleID <- row.names(meta) # Create sampleID column

# Combine soft coral and gorgonian to make Octocoral group
meta$Taxonomy <- gsub("Soft Coral", "Octocoral", meta$Taxonomy) 
meta$Taxonomy <- gsub("Gorgonian", "Octocoral", meta$Taxonomy)

# Load ASV phylogeny 
asv_tree <- read.tree(paste0(q_data, "qiime_insertion_tree_trimmed_all_samples.nwk"))

# Format taxonomy table
tax$Confidence <- NULL
tax <- separate(tax, Taxon, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
                    sep = ";", remove = TRUE, convert = FALSE, extra = "warn", fill = "warn")

tax <- as.data.frame(apply(tax, 2, function(x) gsub(pattern = "D_[0-6]__", replacement = "", x))) 

# Combine to make phyloseq object 
physeq <- phyloseq(
  otu_table(asv %>% as.matrix(), taxa_are_rows=TRUE),
  sample_data(meta),
  tax_table(tax %>% as.matrix()),
  phy_tree(asv_tree)
)
physeq


## Subset and normalise data ----

# Remove samples with no corresponding host phylogeny
physeq <- prune_samples(!sample_names(physeq) %in% c("Gor1", "Gor2", "Gor3","Pmas1", "Pmas2", "Pmas3", "Hist1", "Hist2", "Hist3", "Lobo1"), physeq)

# Subset to samples of interest
all <- subset_samples(physeq, Taxonomy=="Coral" | Taxonomy == "Octocoral" | Taxonomy == "Sponge")
all_sw <- subset_samples(physeq, Taxonomy=="Coral" | Taxonomy == "Octocoral" | Taxonomy == "Sponge" | Taxonomy == "Seawater") # with seawater included

# Save
saveRDS(all, "asv_non-rarefied_phyloseq")
saveRDS(all_sw, "asv_sw_non-rarefied_phyloseq")

# Create host groups
coral <- subset_samples(physeq, Taxonomy=="Coral")
octocoral <- subset_samples(physeq, Taxonomy=="Octocoral")
sponge <- subset_samples(physeq, Taxonomy=="Sponge")


# Normalise data 
rarefy <- function(phy) {
  phy = rarefy_even_depth(phy, 
                    sample.size = min(sample_sums(phy)),
                    rngseed = 5812, replace = FALSE, 
                    trimOTUs = TRUE, verbose = TRUE)
  return(phy)
}

all_rare <- rarefy(all)
all_sw_rare <- rarefy(all_sw)
coral_rare <- rarefy(coral)
octocoral_rare <- rarefy(octocoral)
sponge_rare <- rarefy(sponge)

# Join together
asv_rare_all <- list(coral_rare, octocoral_rare, sponge_rare, all_rare)

# Save
saveRDS(asv_rare_all, "asv_rarefied_phyloseq_all")
saveRDS(all_sw_rare, "asv_sw_rarefied_phyloseq_all")

asv_rare_all <- readRDS(paste0(filepath, "asv_rarefied_phyloseq_all"))
all_sw_rare <- readRDS(paste0(filepath, "asv_sw_rarefied_phyloseq_all"))

# Create metadata list
metadata <- list()
for (i in seq_along(asv_rare_all)) {
  metadata[[i]] <- as.data.frame(as.matrix(sample_data(asv_rare_all[[i]])))
}
saveRDS(metadata, "metadata_files")

# Metadata with seawater
metadata_sw <- as.data.frame(as.matrix(sample_data(all_sw_rare)))
saveRDS(metadata_sw, "metadata_sw")


## Import host tree and format ----

# Import tree
host_all <- read.tree("~/Documents/R/Ramaciotti_16S_dataset/Phylogenies/Final_host/all_species_Bayes_1M_HKY_topcon_newik.txt")
plot(host_all)

host_all <- drop.tip(host_all, tip = "Pol1_Ascidian") # remove as not in analysis

# Edit tip labels
tip_labs <- host_all$tip.label

tip_labs2 <- tip_labs %>% gsub("_C_singaporense", "", .) %>%  gsub("_S_pistillata", "", .) %>% gsub("_S_hysterix", "", .) %>%
  gsub("_P_damnicornis", "", .) %>% gsub("_P_verrucosa", "", .) %>% gsub("_D_heliopora", "", .) %>%
  gsub("_E_mammiformis", "", .) %>% gsub("_P_cylindrica", "", .) %>% gsub("_P_massive", "", .) %>%
  gsub("_Pachyseris_sp", "", .) %>% gsub("_A_formosa", "", .) %>% gsub("_A_hyacinthus", "", .) %>%
  gsub("_P_cactus", "", .) %>% gsub("_Xenia_sp", "", .) %>% gsub("_Briareum_sp2", "", .) %>%
  gsub("_Sinularia_sp", "", .) %>% gsub("_Lobophytum_sp", "", .) %>% gsub("_Briareum_sp", "", .) %>% 
  gsub("_Sarcophyton_sp", "", .) %>% gsub("_I_hippurus", "", .) %>% gsub("_Pinnigorgia_sp", "", .) %>% 
  gsub("_Cladiella_sp", "", .) %>% gsub("_Clavularia_sp", "", .) %>% gsub("_I_blue", "", .) %>% 
  gsub("_Coscinoderma_sp", "", .) %>% gsub("_I_ramosa", "", .) %>% gsub("_C_foliascens", "", .) %>% 
  gsub("_", "-", .) %>% gsub("Irc2", "Icr2", .) %>% gsub("Ahy3", "Ahy4", .) %>% gsub("Ahy2", "Ahy3", .) %>% gsub("Ahy1", "Ahy2", .) 
 
tips <- data.frame(tip_labs, tip_labs2)
tips

host_all$tip.label <- as.character(tips[,2])

# Check host data and ASV data match 
setdiff(host_all$tip.label, unique(sample_data(asv_rare_all[[4]])$SampleID))
setdiff(unique(sample_data(asv_rare_all[[4]])$SampleID), host_all$tip.label)

# Extract clades for host groups
c_tips <- host_all$tip.label[1:58]
o_tips <- host_all$tip.label[59:113]
s_tips <- host_all$tip.label[114:141]

coral <- keep.tip(host_all, c_tips)
octocoral <- keep.tip(host_all, o_tips)
sponge <- keep.tip(host_all, s_tips)

plot(coral)
plot(octocoral)
plot(sponge)

# Creat subsampled host tree for downstream analysis
metadata <- readRDS("metadata_files") 
metadata <- metadata[[4]]

# Save trees
write.tree(host_all, "all_samples_tree_final")
write.tree(coral, "coral_tree_final")
write.tree(octocoral, "octocoral_tree_final")
write.tree(sponge, "sponge_tree_final")



## Create core ASV dataset ----

# Filter OTUs to core microbiome (in this case, those which occur in > 50% of host species replicates) 
physeq <- readRDS("asv_rarefied_phyloseq_all") # using rarefied data (contains a list of coral, octocoral, sponge and all_samples datasets)
summary(physeq)

# First temporarily remove ASV tree
tree <- list()
tax <- list()
otu <- list()
sam <- list()
physeq_noTree <- list()
for(i in seq_along(physeq)) {
  tree[[i]] <- phy_tree(physeq[[i]])
  tax[[i]] <- tax_table(physeq[[i]])
  otu[[i]] <- otu_table(physeq[[i]])
  sam[[i]] <- sample_data(physeq[[i]])
  physeq_noTree[[i]] <- phyloseq(tax[[i]], otu[[i]], sam[[i]])
}

# Set core limit
core_cutoff = 0.5

# Get list of host species to subset core ASVs 
species_list <- list()
for(i in seq_along(physeq_noTree)) {
  species_list[[i]] <- as.character(unique(sample_data(physeq_noTree[[i]])$Species))
  print(species_list[i])
}


# Filter core microbes
#Coral
core_coral <- list()
for(i in seq_along(species_list[[1]])) {
  print(species_list[[1]][i])
  core_coral[[i]] <- subset_samples(physeq_noTree[[1]], Species == species_list[[1]][i])
  core_coral[[i]] <- filter_taxa(core_coral[[i]], function(x) sum(x > 0) / length(x) > core_cutoff, TRUE)
  print(core_coral[i])
}
names(core_coral) <- species_list[[1]]

# Octocoral
core_octocoral <- list()
for(i in seq_along(species_list[[2]])) {
  print(species_list[[2]][i])
  core_octocoral[[i]] <- subset_samples(physeq_noTree[[2]], Species == species_list[[2]][i])
  core_octocoral[[i]] <- filter_taxa(core_octocoral[[i]], function(x) sum(x > 0) / length(x) > core_cutoff, TRUE)
  print(core_octocoral[i])
}
names(core_octocoral) <- species_list[[2]]

# Sponge
core_sponge <- list()
for(i in seq_along(species_list[[3]])) {
  print(species_list[[3]][i])
  core_sponge[[i]] <- subset_samples(physeq_noTree[[3]], Species == species_list[[3]][i])
  core_sponge[[i]] <- filter_taxa(core_sponge[[i]], function(x) sum(x > 0) / length(x) > core_cutoff, TRUE)
  print(core_sponge[i])
}
names(core_sponge) <- species_list[[3]]

# All samples
core_all <- list()
for(i in seq_along(species_list[[4]])) {
  print(species_list[[4]][i])
  core_all[[i]] <- subset_samples(physeq_noTree[[4]], Species == species_list[[4]][i])
  core_all[[i]] <- filter_taxa(core_all[[i]], function(x) sum(x > 0) / length(x) > core_cutoff, TRUE)
  print(core_all[i])
}
names(core_all) <- species_list[[4]]

# Create a separate phyloseq object for each host species with core OTUs filtered
lapply(names(core_coral), function(x) assign(x, core_coral[[x]], envir = .GlobalEnv)) # assigns a value to each name in the list, the value being the index of that list
lapply(names(core_octocoral), function(x) assign(x, core_octocoral[[x]], envir = .GlobalEnv))
lapply(names(core_sponge), function(x) assign(x, core_sponge[[x]], envir = .GlobalEnv))
lapply(names(core_all), function(x) assign(x, core_all[[x]], envir = .GlobalEnv))


# Merge into one phyloseq file
coral <- merge_phyloseq(`A. hyacinthus`, `P. speciosa`, `P. cylindrica_CS`, `S. hystrix_CS`,
                        `D. heliopora`, `E. mammiformis`, `A. formosa`, `P. cactus`, `P. damicornis`,
                        `P. massive`, `P. cylindrica`, `S. hystrix`, `S. pistillata`, `P. verrucosa`, tree[[1]])

octocoral <- merge_phyloseq(`Briareum sp`, `Briareum sp2`, `Clavularia sp`, `Sinularia sp_CS`, `Sarcophyton sp_CS`,
                            `Heteroxenia sp`, `I. hippuris`, `Sinularia sp2`, `Pinnigorgia sp`, `Sarcophyton sp`,
                            `Sinularia sp`, `Cladiella sp`, tree[[2]])

sponge <- merge_phyloseq(`C. foliascens`, `Coscinoderma sp`, `I. ramosa_CS`, `I. ramosa`, `Ircinia sp`, `C. singaporensis`, tree[[3]])

all <- merge_phyloseq(`A. hyacinthus`, `P. speciosa`, `P. cylindrica_CS`, `S. hystrix_CS`,
                      `D. heliopora`, `E. mammiformis`, `A. formosa`, `P. cactus`, `P. damicornis`,
                      `P. massive`, `P. cylindrica`, `S. hystrix`, `S. pistillata`, `P. verrucosa`, 
                      `Briareum sp`, `Briareum sp2`, `Clavularia sp`, `Sinularia sp_CS`, `Sarcophyton sp_CS`,
                      `Heteroxenia sp`, `I. hippuris`, `Sinularia sp2`, `Pinnigorgia sp`, `Sarcophyton sp`,
                      `Sinularia sp`, `Cladiella sp`, `C. foliascens`, `Coscinoderma sp`, `I. ramosa_CS`, 
                      `I. ramosa`, `Ircinia sp`, `C. singaporensis`, tree[[4]])


# Join back together
physeq <- list(coral, octocoral, sponge, all)
summary(physeq)

# Get metadata
metadata <- list()
for (i in seq_along(physeq)) {
  metadata[[i]] <- as.data.frame(as.matrix(sample_data(physeq[[i]])))
}

# Get microbe taxonomy
taxonomy <- list()
for (i in seq_along(physeq)) {
  taxonomy[[i]] <- as.data.frame(as.matrix(tax_table(physeq[[i]]))) 
}

saveRDS(taxonomy, "asv_taxonomy_files")
saveRDS(metadata, "metadata_files")
saveRDS(physeq, "physeq_core_all")





