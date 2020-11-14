### Codivergence analysis adapted following the method of Youngblut et al. 2019. Nat Comms ###

# Paul A. O'Brien
# E. paul.obrien@my.jcu.edu.au


# 03 Run PACo model on each microbe genus and family to test for cophylogeny

## Set directory and load packages ----
setwd("~/Documents/R/Ramaciotti_16S_dataset/Codivergence_and_Neutral/Codivergence/PACo/")

library(dplyr)
library(tidyr)
library(ggplot2)
library(phyloseq)
library(ape)
library(paco)
library(doParallel)
library(future)
library(future.apply)
library(multcomp)


# Function (Youngblut et al. 2019)
rescale_dist_mtx = function(m){
  m = m %>% as.matrix
  labs = m %>% colnames
  n_row = m %>% nrow
  n_col = m %>% ncol
  x = m %>% as.vector 
  x = scales::rescale(x) 
  m = matrix(x, nrow=n_row, ncol=n_col)
  colnames(m) = labs
  rownames(m) = labs
  m = m %>% as.dist
  return(m)
}


## Load and filter data ----

# Phyloseq file - rarefied tables
physeq <- readRDS("asv_rarefied_phyloseq_all")
physeq

# Filter out low abundance ASVs
for(i in seq_along(physeq)) {
  physeq[[i]] <- filter_taxa(physeq[[i]], function(x) sum(x) > 10, TRUE)
}
physeq

# Host phylogeny 
tree_files <- list.files(pattern = "tree_final")[c(2,3,4,1)]
host_trees <- list()
for (f in seq_along(tree_files)) {
  host_trees[[f]] <- read.tree(tree_files[[f]])
}
summary(host_trees)

# Metadata 
metadata <- readRDS("metadata_files")

# Taxonomy
tax <- list()
for(i in seq_along(physeq)) {
  tax[[i]] <- as.data.frame(as.matrix(tax_table(physeq[[i]])))
}


## Subset by microbe genus and family ----

# Make a list of OTU tables for each genus and family subset 
# Create genera list
gen_list <- list()
for(i in seq_along(tax)) {
  gen_list[[i]] <- as.vector(unique(tax[[i]]$Genus))
  gen_list[[i]] <- gen_list[[i]][!is.na(gen_list[[i]])]
}

fam_list <- list()
for(i in seq_along(tax)) {
  fam_list[[i]] <- as.vector(unique(tax[[i]]$Family))
  fam_list[[i]] <- fam_list[[i]][!is.na(fam_list[[i]])]
}

# Pull out OTU tables for those genera/family
# Genera
c_gen_list <- list() # coral
for(i in seq_along(gen_list[[1]])) {
  print(gen_list[[1]][i])
  c_gen_list[[i]] <- subset_taxa(physeq[[1]], Genus == gen_list[[1]][i])
}

o_gen_list <- list() # octocoral
for(i in seq_along(gen_list[[2]])) {
  print(gen_list[[2]][i])
  o_gen_list[[i]] <- subset_taxa(physeq[[2]], Genus == gen_list[[2]][i])
}

s_gen_list <- list() # sponge
for(i in seq_along(gen_list[[3]])) {
  print(gen_list[[3]][i])
  s_gen_list[[i]] <- subset_taxa(physeq[[3]], Genus == gen_list[[3]][i])
}

a_gen_list <- list() # all
for(i in seq_along(gen_list[[4]])) {
  print(gen_list[[4]][i])
  a_gen_list[[i]] <- subset_taxa(physeq[[4]], Genus == gen_list[[4]][i])
}

names(c_gen_list) <- as.vector(gen_list[[1]])
names(o_gen_list) <- as.vector(gen_list[[2]])
names(s_gen_list) <- as.vector(gen_list[[3]])
names(a_gen_list) <- as.vector(gen_list[[4]])

# Family
c_fam_list <- list() # coral
for(i in seq_along(fam_list[[1]])) {
  print(fam_list[[1]][i])
  c_fam_list[[i]] <- subset_taxa(physeq[[1]], Family == fam_list[[1]][i])
}

o_fam_list <- list() # octocoral
for(i in seq_along(fam_list[[2]])) {
  print(fam_list[[2]][i])
  o_fam_list[[i]] <- subset_taxa(physeq[[2]], Family == fam_list[[2]][i])
}

s_fam_list <- list() # sponge
for(i in seq_along(fam_list[[3]])) {
  print(fam_list[[3]][i])
  s_fam_list[[i]] <- subset_taxa(physeq[[3]], Family == fam_list[[3]][i])
}

a_fam_list <- list() # all
for(i in seq_along(fam_list[[4]])) {
  print(fam_list[[4]][i])
  a_fam_list[[i]] <- subset_taxa(physeq[[4]], Family == fam_list[[4]][i])
}


names(c_fam_list) <- as.vector(fam_list[[1]])
names(o_fam_list) <- as.vector(fam_list[[2]])
names(s_fam_list) <- as.vector(fam_list[[3]])
names(a_fam_list) <- as.vector(fam_list[[4]])


# Remove genus / family low number of representatives to give more meaningful phylogeny 

# Genus 
# Coral
c_gen_list_s <- list()
c_gen_name <- vector()
for(i in seq_along(c_gen_list)) {
  if (nrow(tax_table(c_gen_list[[i]])) > 3) {
    c_gen_list_s[i] <- c_gen_list[i]
    c_gen_name[i] <- names(c_gen_list[i])
  } else {
    NULL
  }
}
c_gen_list <- c_gen_list_s[-which(sapply(c_gen_list_s, is.null))]
c_gen_name <- c_gen_name[!is.na(c_gen_name)]
names(c_gen_list) <- c_gen_name
saveRDS(c_gen_list, "coral_genera_physeq")


# Octocoral
o_gen_list_s <- list()
o_gen_name <- vector()
for(i in seq_along(o_gen_list)) {
  if (nrow(tax_table(o_gen_list[[i]])) > 3) {
    o_gen_list_s[i] <- o_gen_list[i]
    o_gen_name[i] <- names(o_gen_list[i])
  } else {
    NULL
  }
}
o_gen_list <- o_gen_list_s[-which(sapply(o_gen_list_s, is.null))]
o_gen_name <- o_gen_name[!is.na(o_gen_name)]
names(o_gen_list) <- o_gen_name
saveRDS(o_gen_list, "octocoral_genera_physeq")

# Sponge
s_gen_list_s <- list()
s_gen_name <- vector()
for(i in seq_along(s_gen_list)) {
  if (nrow(tax_table(s_gen_list[[i]])) > 3) {
    s_gen_list_s[i] <- s_gen_list[i]
    s_gen_name[i] <- names(s_gen_list[i])
  } else {
    NULL
  }
}
s_gen_list <- s_gen_list_s[-which(sapply(s_gen_list_s, is.null))]
s_gen_name <- s_gen_name[!is.na(s_gen_name)]
names(s_gen_list) <- s_gen_name
saveRDS(s_gen_list, "sponge_genera_physeq")


# All
a_gen_list_s <- list()
a_gen_name <- vector()
for(i in seq_along(a_gen_list)) {
  if (nrow(tax_table(a_gen_list[[i]])) > 9) {
    a_gen_list_s[i] <- a_gen_list[i]
    a_gen_name[i] <- names(a_gen_list[i])
  } else {
    NULL
  }
}
a_gen_list <- a_gen_list_s[-which(sapply(a_gen_list_s, is.null))]
a_gen_name <- a_gen_name[!is.na(a_gen_name)]
names(a_gen_list) <- a_gen_name
saveRDS(a_gen_list, "all_genera_physeq")


# Family
# Coral
c_fam_list_s <- list()
c_fam_name <- vector()
for(i in seq_along(c_fam_list)) {
  if (nrow(tax_table(c_fam_list[[i]])) > 3) {
    c_fam_list_s[i] <- c_fam_list[i]
    c_fam_name[i] <- names(c_fam_list[i])
  } else {
    NULL
  }
}
c_fam_list <- c_fam_list_s[-which(sapply(c_fam_list_s, is.null))]
c_fam_name <- c_fam_name[!is.na(c_fam_name)]
names(c_fam_list) <- c_fam_name
saveRDS(c_fam_list, "coral_family_physeq")


# Octocoral
o_fam_list_s <- list()
o_fam_name <- vector()
for(i in seq_along(o_fam_list)) {
  if (nrow(tax_table(o_fam_list[[i]])) > 3) {
    o_fam_list_s[i] <- o_fam_list[i]
    o_fam_name[i] <- names(o_fam_list[i])
  } else {
    NULL
  }
}
o_fam_list <- o_fam_list_s[-which(sapply(o_fam_list_s, is.null))]
o_fam_name <- o_fam_name[!is.na(o_fam_name)]
names(o_fam_list) <- o_fam_name
saveRDS(o_fam_list, "octocoral_family_physeq")


# Sponge
s_fam_list_s <- list()
s_fam_name <- vector()
for(i in seq_along(s_fam_list)) {
  if (nrow(tax_table(s_fam_list[[i]])) > 3) {
    s_fam_list_s[i] <- s_fam_list[i]
    s_fam_name[i] <- names(s_fam_list[i])
  } else {
    NULL
  }
}
s_fam_list <- s_fam_list_s[-which(sapply(s_fam_list_s, is.null))]
s_fam_name <- s_fam_name[!is.na(s_fam_name)]
names(s_fam_list) <- s_fam_name
saveRDS(s_fam_list, "sponge_family_physeq")


# All 
a_fam_list_s <- list()
a_fam_name <- vector()
for(i in seq_along(a_fam_list)) {
  if (nrow(tax_table(a_fam_list[[i]])) > 9) {
    a_fam_list_s[i] <- a_fam_list[i]
    a_fam_name[i] <- names(a_fam_list[i])
  } else {
    NULL
  }
}
a_fam_list <- a_fam_list_s[-which(sapply(a_fam_list_s, is.null))]
a_fam_name <- a_fam_name[!is.na(a_fam_name)]
names(a_fam_list) <- a_fam_name
saveRDS(a_fam_list, "all_family_physeq")

#Get metadata
c_meta <- metadata[[1]]
o_meta <- metadata[[2]]
s_meta <- metadata[[3]]
a_meta <- metadata[[4]]


## PACO: setup model ----

## Using sensitivity analysis as per Youngblut 2019 ##

## Load subsampled host trees
host_tree_l <- readRDS("PACo_results/host_tree_list_100")
summary(host_tree_l)


# Function to create input data (Youngblut et al. 2019)
make_paco_input <- function(host_tree, physeq){
  # subsampled phyloseq object
  physeq_f <- prune_samples(sample_names(physeq) %in% host_tree$tip.label, 
                            physeq) 
  
  # microbial distance matrix
  micro_D <- as.matrix(rescale_dist_mtx(cophenetic(phy_tree(physeq_f))))
  
  # host distance matrix
  host_D <- as.matrix(rescale_dist_mtx(cophenetic(host_tree)))
  
  # otu abundance matrix
  otu <- t(as.data.frame(as.matrix(otu_table(physeq_f)))) %>%
    apply(2, function(x) ifelse(x > 0, 1, 0)) %>% as.matrix
  
  # checking overlap
  x <- length(intersect(rownames(otu), rownames(host_D)))
  y <- length(union(rownames(otu), rownames(host_D)))
  stopifnot(x == y)
  
  # preparing paco data
  D <- prepare_paco_data(H=host_D, P=micro_D, HP=otu)
  D <- add_pcoord(D, correction='cailliez')
  return(D)
}

## Create paco objects
doParallel::registerDoParallel(2)

## Genus ##

#Load files
c_gen_list <- readRDS("coral_genera_physeq")
o_gen_list <- readRDS("octocoral_genera_physeq")
s_gen_list <- readRDS("sponge_genera_physeq")
a_gen_list <- readRDS("all_genera_physeq")

# Coral
paco_gen_c_l <- list()
for(i in seq_along(c_gen_list)) {
  paco_gen_c_l[[i]] <- plyr::llply(host_tree_l[[1]], make_paco_input, physeq=c_gen_list[[i]], .parallel=TRUE)
  print(paco_gen_c_l[[i]] %>% length)
}
saveRDS(paco_gen_c_l, "PACo_results/paco_list_dataFile_coral_gen")

# Octocoral
paco_gen_o_l <- list()
for(i in seq_along(o_gen_list)) {
  paco_gen_o_l[[i]] <- plyr::llply(host_tree_l[[2]], make_paco_input, physeq=o_gen_list[[i]], .parallel=TRUE)
  print(paco_gen_o_l[[i]] %>% length)
}
saveRDS(paco_gen_o_l, "PACo_results/paco_list_dataFile_octocoral_gen")

# Sponge
paco_gen_s_l <- list()
for(i in seq_along(s_gen_list)) {
  paco_gen_s_l[[i]] <- plyr::llply(host_tree_l[[3]], make_paco_input, physeq=s_gen_list[[i]], .parallel=TRUE)
  print(paco_gen_s_l[[i]] %>% length)
}
saveRDS(paco_gen_s_l, "PACo_results/paco_list_dataFile_sponge_gen")

# All 
paco_gen_a_l <- list()
for(i in seq_along(a_gen_list)) {
  paco_gen_a_l[[i]] <- plyr::llply(host_tree_l[[4]], make_paco_input, physeq=a_gen_list[[i]], .parallel=TRUE)
  print(paco_gen_a_l[[i]] %>% length)
}
saveRDS(paco_gen_a_l, "PACo_results/paco_list_dataFile_all_gen")


## Family ##

#Load files
c_fam_list <- readRDS("coral_family_physeq")
o_fam_list <- readRDS("octocoral_family_physeq")
s_fam_list <- readRDS("sponge_family_physeq")
a_fam_list <- readRDS("all_family_physeq")

# Coral
paco_fam_c_l <- list()
for(i in seq_along(c_fam_list)) {
  paco_fam_c_l[[i]] <- plyr::llply(host_tree_l[[1]], make_paco_input, physeq=c_fam_list[[i]], .parallel=TRUE)
  print(paco_fam_c_l[[i]] %>% length)
}
saveRDS(paco_fam_c_l, "PACo_results/paco_list_dataFile_coral_fam")

# Octocoral
paco_fam_o_l <- list()
for(i in seq_along(o_fam_list)) {
  paco_fam_o_l[[i]] <- plyr::llply(host_tree_l[[2]], make_paco_input, physeq=o_fam_list[[i]], .parallel=TRUE)
  print(paco_fam_o_l[[i]] %>% length)
}
saveRDS(paco_fam_o_l, "PACo_results/paco_list_dataFile_octocoral_fam")

# Sponge
paco_fam_s_l <- list()
for(i in seq_along(s_fam_list)) {
  paco_fam_s_l[[i]] <- plyr::llply(host_tree_l[[3]], make_paco_input, physeq=s_fam_list[[i]], .parallel=TRUE)
  print(paco_fam_s_l[[i]] %>% length)
}
saveRDS(paco_fam_s_l, "PACo_results/paco_list_dataFile_sponge_fam")

# All
paco_fam_a_l <- list()
for(i in seq_along(a_fam_list)) {
  paco_fam_a_l[[i]] <- plyr::llply(host_tree_l[[4]], make_paco_input, physeq=a_fam_list[[i]], .parallel=TRUE)
  print(paco_fam_a_l[[i]] %>% length)
}
saveRDS(paco_fam_a_l, "PACo_results/paco_list_dataFile_all_fam")



## Run model ----

# Function to run each paco model in parallel (Youngblut et al. 2019)
paco_each <- function(D, nperm=99, seed=8355){
  PACo(D, nperm=nperm, seed=seed, method='quasiswap', symmetric=TRUE)
}


## Genus ##

## Load data
paco_gen_c_l <- readRDS("PACo_results/paco_list_dataFile_coral_gen")
paco_gen_o_l <- readRDS("PACo_results/paco_list_dataFile_octocoral_gen")
paco_gen_s_l <- readRDS("PACo_results/paco_list_dataFile_sponge_gen")
paco_gen_a_l <- readRDS("PACo_results/paco_list_dataFile_all_gen")

# Coral
paco_res_l_c_gen <- list()
for(i in seq_along(paco_gen_c_l)) {
  paco_res_l_c_gen[[i]] <- future_lapply(paco_gen_c_l[[i]], paco_each, nperm=999, future.packages=c('paco'))
  print(paste(i, "complete"))
}
saveRDS(paco_res_l_c_gen, "PACo_results/GOF_test_results/paco_res_l_coral_gen")


# Octocoral
paco_res_l_o_gen <- list()
for(i in seq_along(paco_gen_o_l)) {
  paco_res_l_o_gen[[i]] <- future_lapply(paco_gen_o_l[[i]], paco_each, nperm=999, future.packages=c('paco'))
  print(paste(i, "complete"))
}
saveRDS(paco_res_l_o_gen, "PACo_results/GOF_test_results/paco_res_l_octocoral_gen")


# Sponge
paco_res_l_s_gen <- list()
for(i in seq_along(paco_gen_s_l)) {
  paco_res_l_s_gen[[i]] <- future_lapply(paco_gen_s_l[[i]], paco_each, nperm=999, future.packages=c('paco'))
  print(paste(i, "complete"))
}
saveRDS(paco_res_l_s_gen, "PACo_results/GOF_test_results/paco_res_l_sponge_gen")


# All
paco_res_l_a_gen <- list()
for(i in seq_along(paco_gen_a_l)) {
  paco_res_l_a_gen[[i]] <- future_lapply(paco_gen_a_l[[i]], paco_each, nperm=999, future.packages=c('paco'))
  print(paste(i, "complete"))
}
saveRDS(paco_res_l_a_gen, "PACo_results/GOF_test_results/paco_res_l_all_gen")


## Family ##

paco_fam_c_l <- readRDS("PACo_results/paco_list_dataFile_coral_fam")
paco_fam_o_l <- readRDS("PACo_results/paco_list_dataFile_octocoral_fam")
paco_fam_s_l <- readRDS("PACo_results/paco_list_dataFile_sponge_fam")
paco_fam_a_l <- readRDS("PACo_results/paco_list_dataFile_all_fam")

# Coral
paco_res_l_c_fam <- list()
for(i in seq_along(paco_fam_c_l)) {
  paco_res_l_c_fam[[i]] <- future_lapply(paco_fam_c_l[[i]], paco_each, nperm=999, future.packages=c('paco'))
  print(paste(i, "complete"))
}
saveRDS(paco_res_l_c_fam, "PACo_results/GOF_test_results/paco_res_l_coral_fam")


# Octocoral
paco_res_l_o_fam <- list()
for(i in seq_along(paco_fam_o_l)) {
  paco_res_l_o_fam[[i]] <- future_lapply(paco_fam_o_l[[i]], paco_each, nperm=999, future.packages=c('paco'))
  print(paste(i, "complete"))
}
saveRDS(paco_res_l_o_fam, "PACo_results/GOF_test_results/paco_res_l_octocoral_fam")


# Sponge
paco_res_l_s_fam <- list()
for(i in seq_along(paco_fam_s_l)) {
  paco_res_l_s_fam[[i]] <- future_lapply(paco_fam_s_l[[i]], paco_each, nperm=999, future.packages=c('paco'))
  print(paste(i, "complete"))
}
saveRDS(paco_res_l_s_fam, "PACo_results/GOF_test_results/paco_res_l_sponge_fam")


# All
paco_res_l_a_fam <- list()
for(i in seq_along(paco_fam_a_l)) {
  paco_res_l_a_fam[[i]] <- future_lapply(paco_fam_a_l[[i]], paco_each, nperm=999, future.packages=c('paco'))
  print(paste(i, "complete"))
}
saveRDS(paco_res_l_a_fam, "PACo_results/GOF_test_results/paco_res_l_all_fam")


## Format output ----

# Combine GOF list results

# Genus all
paco_res_l_a_gen <- readRDS("PACo_results/GOF_test_results/paco_res_l_all_gen")
gof_a_gen <- list()
for(i in seq_along(paco_res_l_a_gen)) {
  gof_a_gen[[i]] <- lapply(paco_res_l_a_gen[[i]], function(x) as.vector(x$gof)) %>%
    do.call(rbind, .) %>%
    as.data.frame
}
lapply(gof_a_gen, function(x) summary(as.numeric(x$p)))
gof_a_gen[[i]]$p %>% as.numeric %>% summary %>% print
gof_a_gen[[i]]$ss %>% as.numeric %>% summary %>% print

# adjust p
for(i in seq_along(gof_a_gen)) {
  gof_a_gen[[i]]$q <- p.adjust(gof_a_gen[[i]]$p, method = 'BH')
  gof_a_gen[[i]]$q %>% as.numeric %>% summary %>% print
}

# Family all
paco_res_l_a_fam <- readRDS("PACo_results/GOF_test_results/paco_res_l_all_fam")
gof_a_fam <- list()
for(i in seq_along(paco_res_l_a_fam)) {
  gof_a_fam[[i]] <- lapply(paco_res_l_a_fam[[i]], function(x) as.vector(x$gof)) %>%
    do.call(rbind, .) %>%
    as.data.frame
}
lapply(gof_a_fam, function(x) summary(as.numeric(x$p)))
lapply(gof_a_fam, function(x) summary(as.numeric(x$ss)))

# adjust p
for(i in seq_along(gof_a_fam)) {
  gof_a_fam[[i]]$q <- p.adjust(gof_a_fam[[i]]$p, method = 'BH')
  gof_a_fam[[i]]$q %>% as.numeric %>% summary %>% print
}


# Create GOF dataframe for each taxa


# Genus all
mean_pval <- vector()
mean_qval <- vector()
sd_pval <- vector()
mean_ss <- vector()
sd_ss <- vector()
for(i in seq_along(gof_a_gen)) {
  # get mean pval
  mean_pval[i] <- mean(as.numeric(gof_a_gen[[i]]$p))
  mean_qval[i] <- mean(as.numeric(gof_a_gen[[i]]$q))
  sd_pval[i] <- stats::sd(as.numeric(gof_a_gen[[i]]$p))
  
  # get mean sums of squares
  mean_ss[i] <- mean(as.numeric(gof_a_gen[[i]]$ss))
  sd_ss[i] <- stats::sd(as.numeric(gof_a_gen[[i]]$ss))
}

# Create dataframe
gof_a_gen_s <- data.frame(names(a_gen_list), mean_pval, mean_qval, sd_pval, mean_ss, sd_ss)
head(gof_a_gen_s)

# Remove uncultured as contain multiple phyla 
gof_a_gen_s <- gof_a_gen_s[-which(gof_a_gen_s$names.a_gen_list.=="uncultured"),]
gof_a_gen_s <- gof_a_gen_s[-which(gof_a_gen_s$names.a_gen_list.=="uncultured bacterium"),]

# order dataframe
gof_a_gen_s <- gof_a_gen_s[order(gof_a_gen_s$mean_qval, gof_a_gen_s$mean_ss),]
gof_a_gen_s

options(scipen = 999) # to get regular number instead of exponential 
write.table(gof_a_gen_s, "all_genus_sens_GOF.tsv", sep = "\t", row.names = FALSE)


# Family all
mean_pval <- vector()
mean_qval <- vector()
sd_pval <- vector()
mean_ss <- vector()
sd_ss <- vector()
for(i in seq_along(gof_a_fam)) {
  # get mean pval
  mean_pval[i] <- mean(as.numeric(gof_a_fam[[i]]$p))
  mean_qval[i] <- mean(as.numeric(gof_a_fam[[i]]$q))
  sd_pval[i] <- stats::sd(as.numeric(gof_a_fam[[i]]$p))
  
  # get mean sums of squares
  mean_ss[i] <- mean(as.numeric(gof_a_fam[[i]]$ss))
  sd_ss[i] <- stats::sd(as.numeric(gof_a_fam[[i]]$ss))
}

# Create dataframe
gof_a_fam_s <- data.frame(names(a_fam_list), mean_pval, mean_qval, sd_pval, mean_ss, sd_ss)
head(gof_a_fam_s)

# Order dataframe
gof_a_fam_s <- gof_a_fam_s[order(gof_a_fam_s$mean_qval, gof_a_fam_s$mean_ss),]
write.table(gof_a_fam_s, "all_family_sens_GOF.tsv", sep = "\t", row.names = FALSE)


