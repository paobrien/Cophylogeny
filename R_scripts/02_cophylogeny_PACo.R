
### Codivergence analysis adapted following the method of Youngblut et al. 2019. Nat Comms ###

# Paul A. O'Brien
# E. paul.obrien@my.jcu.edu.au


# 02 Run PACo model on core data for cophylogenetic analysis 


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

library(stringr)


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

# Core phyloseq data - all core filtered phyloseq data
physeq <- readRDS("physeq_core_all")

# Hot phylogeny 
tree_files <- list.files(pattern = "tree_final")[c(2,3,4,1)]
host_trees <- list()
for (f in seq_along(tree_files)) {
  host_trees[[f]] <- read.tree(tree_files[[f]])
}
summary(host_trees)

# Metadata and taxonomy
metadata <- readRDS("metadata_files")
taxonomy <- readRDS("asv_taxonomy_files")



## Subsample host tree for sensitivity analysis ----

# Function to randomly select one sample per host (From Youngblut et al. 2019)
# Recommended when unequal sample sizes and phylogenetic uncertainty - also less likely to get false positives

ntrees = 100

tree_subsample <- function(L, df, tree){                  
  # get subsample (note: subsampling within each species)
  to_keep = df %>% 
    group_by(Species) %>% 
    sample_n(1) %>%
    .$sample
  # subsampling tree
  to_rm = setdiff(tree$tip.label, to_keep)
  tree = drop.tip(tree, to_rm)
  return(tree)
}

# Execute function
registerDoParallel(3)

df <- list()
for(i in seq_along(metadata)){
  df[[i]] <- metadata[[i]] %>%
    mutate(sample = rownames(metadata[[i]])) %>% select(sample, Species) 
}

host_tree_l <- list()
for(i in seq_along(host_trees)) {
  host_tree_l[[i]] <- plyr::llply(as.list(1:ntrees), 
                                  function(x) tree_subsample(x, df[[i]], host_trees[[i]]), # x being a list of trees (1:ntrees)
                                  .parallel=TRUE)
}

# Check tree lengths
for (i in seq_along(host_tree_l)) {
  print(length(host_tree_l[[i]]))
  print(lapply(host_tree_l[[i]], function(x) x$tip.label %>% length) %>% unlist %>% summary)
}

# Check overlap of labels
for(i in seq_along(host_tree_l)) {
  print(setdiff(host_tree_l[[i]]$tip.label, rownames(metadata[[i]])) %>% length)
}

saveRDS(host_tree_l, "host_tree_list_100")



## PACo model ----

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


# Create list of paco objects
doParallel::registerDoParallel(3)
paco_l_all <- list()
for(i in seq_along(host_tree_l)) {
  paco_l_all[[i]] <- plyr::llply(host_tree_l[[i]], make_paco_input, physeq=physeq[[i]], .parallel=TRUE)
  print(paco_l_all[[i]] %>% length)
}

saveRDS(paco_l_all, "paco_list_dataFile_all")

## Run model ----

# Note: From Hutchinson et al. 2017, on symmetric = TRUE;
# Tests the dependence of both phylogenies on each other. May be important when condsidering diffuse interactions 
# (i.e., interactions differ depending on a third species), where it is not clear if one clade depends on the evolution of the other


############################
# Run the following on HPC #
############################


# Function to run each paco model in parallel (this time have 100 subsets of paco formatted data)
paco_each <- function(D, nperm=99, seed=8355){
  PACo(D, nperm=nperm, seed=seed, method='quasiswap', symmetric=TRUE)
}

# Use future env to then run each model in parallel using lapply
plan("multisession")

paco_res_l_coral <- future_lapply(paco_l_all[[1]], paco_each, nperm=999, future.packages=c('paco'))
length(paco_res_l_coral)

paco_res_l_octocoral <- future_lapply(paco_l_all[[2]], paco_each, nperm=999, future.packages=c('paco'))
length(paco_res_l_octocoral)

paco_res_l_sponge <- future_lapply(paco_l_all[[3]], paco_each, nperm=999, future.packages=c('paco'))
length(paco_res_l_sponge)

paco_res_l_all <- future_lapply(paco_l_all[[4]], paco_each, nperm=999, future.packages=c('paco'))
length(paco_res_l_all)

warnings()

# Save results
saveRDS(paco_res_l_coral, "paco_res_l_coral")
saveRDS(paco_res_l_octocoral, "paco_res_l_octocoral")
saveRDS(paco_res_l_sponge, "paco_res_l_sponge")
saveRDS(paco_res_l_all, "paco_res_l_all")

# Calculate goodness of fit
gof_c <- lapply(paco_res_l_coral, function(x) as.vector(x$gof)) %>%
  do.call(rbind, .) %>%
  as.data.frame

gof_o <- lapply(paco_res_l_octocoral, function(x) as.vector(x$gof)) %>%
  do.call(rbind, .) %>%
  as.data.frame

gof_s <- lapply(paco_res_l_sponge, function(x) as.vector(x$gof)) %>%
  do.call(rbind, .) %>%
  as.data.frame

gof_a <- lapply(paco_res_l_all, function(x) as.vector(x$gof)) %>%
  do.call(rbind, .) %>%
  as.data.frame

gof_results <- list(gof_c, gof_o, gof_s, gof_a)
names(gof_results) <- c("coral", "octocoral", "sponge", "all")

# Adjust p-value
for(i in seq_along(gof_results)) {
  gof_results[[i]]$q <- p.adjust(gof_results[[i]]$p, method = 'BH')
  gof_results[[i]]$q %>% as.numeric %>% summary %>% print
}

for(i in seq_along(gof_results)) {
  print(names(gof_results[i]))
  gof_results[[i]]$p %>% as.numeric %>% summary %>% print
  gof_results[[i]]$q %>% as.numeric %>% summary %>% print
  gof_results[[i]]$ss %>% as.numeric %>% summary %>% print
  
}

# Save results
saveRDS(gof_results, "gof_all")


##########################
# Transfer data to local #
##########################


############################
# Run the following on HPC #
############################


# Run model
plan("multisession")

paco_coral_links_l <- future_lapply(paco_res_l_coral, paco_links, future.packages=c('paco'))
length(paco_coral_links_l)

paco_octocoral_links_l <- future_lapply(paco_res_l_octocoral, paco_links, future.packages=c('paco'))
length(paco_octocoral_links_l)

paco_sponge_links_l <- future_lapply(paco_res_l_sponge, paco_links, future.packages=c('paco'))
length(paco_sponge_links_l)

paco_all_links_l <- future_lapply(paco_res_l_all, paco_links, future.packages=c('paco'))
length(paco_all_links_l)


warnings()


saveRDS(paco_coral_links_l, "paco_coral_links_l")
saveRDS(paco_octocoral_links_l, "paco_octocoral_links_l")
saveRDS(paco_sponge_links_l, "paco_sponge_links_l")
saveRDS(paco_all_links_l, "paco_all_links_l")


##########################
# Transfer data to local #
##########################


## Load data ##
links_files <- list.files("PACo_results/Individual_links", pattern = "links_l")[c(2,3,4,1)] # change order so coral first in list
paco_links_all <- list()
for(i in seq_along(links_files)) {
  paco_links_all[[i]] <- readRDS(paste0("PACo_results/Individual_links/", links_files[i]))
}

## Format output ----

# Functions from Youngblut et al. 2019
# added line to functions to correct Correct Br2-1 names ("-" causing problems)

# residuals
get_residuals <- function(rep, D_links_l){
  # residuals
  res = residuals_paco(D_links_l[[rep]]$proc) %>% as.data.frame
  
  # correct row names
  row.names(res) <- gsub("Br2-1", "Br2_1", row.names(res))
  row.names(res) <- gsub("Br2-2", "Br2_2", row.names(res))
  row.names(res) <- gsub("Br2-3", "Br2_3", row.names(res))
  row.names(res) <- gsub("Br2-4", "Br2_4", row.names(res))
  row.names(res) <- gsub("Br2-5", "Br2_5", row.names(res))
  
  colnames(res) = 'residuals'
  res = res %>%
    mutate(comparison = rownames(.),
           subsample_rep = rep) %>%
    separate(comparison, c('host', 'microbe'), sep='-') 
  
  # jackknife
  D_links_jk = do.call(rbind, D_links_l[[rep]]$jackknife) %>%
    t %>% as.data.frame
  
  # correct row names
  row.names(D_links_jk) <- gsub("Br2-1", "Br2_1", row.names(D_links_jk))
  row.names(D_links_jk) <- gsub("Br2-2", "Br2_2", row.names(D_links_jk))
  row.names(D_links_jk) <- gsub("Br2-3", "Br2_3", row.names(D_links_jk))
  row.names(D_links_jk) <- gsub("Br2-4", "Br2_4", row.names(D_links_jk))
  row.names(D_links_jk) <- gsub("Br2-5", "Br2_5", row.names(D_links_jk))
  
  D_links_jk = D_links_jk %>% mutate(comparison = rownames(.)) %>%
    separate(comparison, c('host', 'microbe'), sep='-') %>%
    inner_join(res, c('host'='host', 'microbe'='microbe'))
  
  # return
  return(D_links_jk)
}


links_all <- list()
for(i in seq_along(paco_links_all)) {
  links_all[[i]] <- lapply(as.list(1:length(paco_links_all[[i]])), get_residuals, D_links_l=paco_links_all[[i]])
  links_all[[i]] <- do.call(rbind, links_all[[i]])
}

lapply(links_all, head)


for(i in seq_along(links_all)) {
  links_all[[i]] <- links_all[[i]] %>%
    group_by(host, microbe) %>%
    summarize(mean_resid = mean(residuals),
              median_resid = median(residuals),
              sd_resid = sd(residuals),
              CV_resid = sd_resid / mean_resid * 100) %>%
    ungroup()
  print(head(links_all[[i]]))
}

lapply(links_all, dim)


# Add host taxonomy
# First correct Br- row names
for(i in seq_along(metadata)) {
  metadata[[i]]$SampleID <- gsub("Br2-1", "Br2_1", metadata[[i]]$SampleID) 
  metadata[[i]]$SampleID <- gsub("Br2-2", "Br2_2", metadata[[i]]$SampleID)
  metadata[[i]]$SampleID <- gsub("Br2-3", "Br2_3", metadata[[i]]$SampleID)
  metadata[[i]]$SampleID <- gsub("Br2-4", "Br2_4", metadata[[i]]$SampleID)
  metadata[[i]]$SampleID <- gsub("Br2-5", "Br2_5", metadata[[i]]$SampleID)
}


tmp <- list()
for(i in seq_along(links_all)) {
  # Get host taxonomy
  tmp[[i]] <- metadata[[i]] %>% 
    dplyr::select(Taxonomy, Species, Reef, Zone, SampleID)
  rownames(tmp[[i]]) <- 1:nrow(tmp[[i]])
  # Add to links dataframe
  links_all[[i]] <- inner_join(links_all[[i]], tmp[[i]], c('host'='SampleID')) 
  print(head(links_all[[i]]))
  print(dim(links_all[[i]]))
}

# Add microbe taxonomy
for(i in seq_along(taxonomy)) {
  taxonomy[[i]] <- taxonomy[[i]] %>% mutate(microbe = rownames(.))
  links_all[[i]] <- inner_join(links_all[[i]], taxonomy[[i]], c('microbe'))
  print(head(links_all[[i]]))
  print(dim(links_all[[i]]))
}




## Summarise results by host ----

# Plot results summary - look at relationship between mean and SD
p <- ggplot(links_all[[4]], aes(mean_resid, sd_resid, color=Taxonomy)) +
  geom_point(alpha=0.5) +
  theme_bw() 
p


# Coefficient of Variation
# Increase CV suggests that intra-species variation may influence residuals

for(i in seq_along(links_all)) {
  print(summary(links_all[[i]]$CV_resid))
}

p <- ggplot(links_all[[4]], aes(CV_resid)) +
  geom_histogram(binwidth=0.25) +
  labs(x='intra-species CV in PACo resimduals', y='# of host-microbe\ncomparisons') +
  facet_wrap(~ Species.x, scales='free_y') +
  theme_bw()
p


# Plot by host taxonomy
# First reorder factors
tmp <- list()
for(i in seq_along(links_all)) {
  tmp[[i]] <- links_all[[i]] %>% 
    group_by(Species.x) %>%
    mutate(median_resid = median(mean_resid)) %>%
    ungroup() %>%
    mutate(Species.x = Species.x %>% reorder(-median_resid)) 
  tmp[[i]]$Phylum <- tmp[[i]]$Phylum %>% as.character %>% replace_na("Unknown") 
}


# Plot by host group
p <- ggplot(tmp[[4]], aes(Taxonomy, mean_resid)) +
  geom_boxplot() +
  labs(x='', y='Residuals') +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle=45, hjust=1)
  )
p 
#+ geom_text(data=resids.groups.df, aes(x = Taxonomy, y = Ymax+0.005, label = letters),vjust=0) # for significance


# Anova on residuals (host group)
tmp[[4]]$fourth <-(tmp[[4]]$mean_resid ^ (1/4)) # Transform data
resids.lm <- lm(fourth ~ Taxonomy, tmp[[4]])
plot(resids.lm)

summary(resids.lm) # summarise linear model.
anova(resids.lm) # for anova

#post hoc test
glht.resids.sum <- summary(glht(resids.lm, linfct = mcp(Taxonomy = "Tukey")), test = adjusted("bonferroni"))
resids.groups <- cld(glht.resids.sum)
resids.groups.df <- fortify(resids.groups)
colnames(resids.groups.df) <- c("Taxonomy", "letters")

ymax <- tapply(tmp[[4]]$mean_resid, tmp[[4]]$Taxonomy, max)
resids.groups.df$Ymax <- ymax # add to plot above


# plot by host species
p <- ggplot(tmp[[1]]) +
  geom_boxplot(aes(Species.x, mean_resid), outlier.shape=NA) +
  geom_jitter(aes(Species.x, mean_resid, colour = Phylum), position=position_jitter(width=.2, height=0)) +
  labs(x='Species', y='Residuals') +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle=45, hjust=1)
  )
p

# Plot species grouped by host taxa
p <- ggplot(tmp[[4]]) +
  geom_boxplot(aes(Species.x, mean_resid)) +
  geom_jitter(aes(Species.x, mean_resid, colour = Phylum, alpha = 0.1), position=position_jitter(width=.2, height=0)) +
  labs(x='Species', y='Residuals') +
  theme_bw() +
  facet_grid(. ~ Taxonomy, scales='free_x', space='free_x')+
  guides(alpha = FALSE, colour=guide_legend(ncol = 1)) +
  theme(
    axis.text.x = element_text(angle=65, hjust=1)
  )
p 

ggsave("all_species_residuals_microbes.svg", device = "svg", width = 30, height = 17, units = "cm", dpi = 300)



# anova on residuals (host species)
# by host species
resids.lm <- lm(fourth ~ Species.x, tmp[[4]])
plot(resids.lm)

summary(resids.lm) # summarise linear model.
anova(resids.lm) # for anova


# Plot by host and zone
p <- ggplot(tmp[[4]]) +
  geom_boxplot(aes(Taxonomy, mean_resid, color=Zone)) +
  scale_color_discrete('Taxa') +
  labs(x='', y='Residuals') +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle=45, hjust=1)
  )
p


## Summarise by microbes ----
# By genus
tmp <- list()
for(i in seq_along(links_all)) {
  tmp[[i]] <- links_all[[i]] %>%
    filter(!is.na(mean_resid),
           Genus != '', 
           Genus != 'uncultured',
           Genus != 'unclassified',
           Genus != 'Candidatus Tectomicrobia bacterium RIFCSPLOWO2_12_FULL_69_37',
           Genus != 'alpha proteobacterium endosymbiont 2 of Inanidrilus makropetalos',
           Genus != 'uncultured Syntrophobacterales bacterium',
           Genus != 'uncultured delta proteobacterium',
           Genus != 'metagenome',
           Genus != 'uncultured archaeon',
           Genus != 'uncultured bacterium',
           Genus != 'uncultured Deferribacteres bacterium',
           Genus != 'uncultured marine bacterium',
           Genus != 'uncultured gamma proteobacterium') %>%
    group_by(Genus) %>%
    mutate(mean_res = mean(mean_resid, na.rm=TRUE),
           n = n() %>% as.numeric) %>%
    ungroup() %>%
    mutate(Genus = Genus %>% reorder(mean_res),
           Genus = Genus %>% reorder(Phylum %>% as.factor %>% as.numeric))
}


# By family
tmp <- list()
for(i in seq_along(links_all)) {
  tmp[[i]] <- links_all[[i]] %>%
    filter(!is.na(mean_resid),
           Family != '', 
           Family != 'uncultured',
           Family != 'unclassified',
           Family != 'Candidatus Tectomicrobia bacterium RIFCSPLOWO2_12_FULL_69_37',
           Family != 'alpha proteobacterium endosymbiont 2 of Inanidrilus makropetalos',
           Family != 'uncultured Syntrophobacterales bacterium',
           Family != 'uncultured delta proteobacterium',
           Family != 'metagenome',
           Family != 'uncultured archaeon',
           Family != 'uncultured bacterium',
           Family != 'uncultured Deferribacteres bacterium',
           Family != 'uncultured marine bacterium',
           Family != 'uncultured gamma proteobacterium') %>%
    group_by(Family) %>%
    mutate(mean_res = mean(mean_resid, na.rm=TRUE),
           n = n() %>% as.numeric) %>%
    ungroup() %>%
    mutate(Family = Family %>% reorder(mean_res),
           Family = Family %>% reorder(Phylum %>% as.factor %>% as.numeric))
}


p <- ggplot(tmp[[4]], aes(Family, mean_resid, color=Phylum)) +
  geom_boxplot() +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)
  )
p

p <- ggplot(tmp[[4]], aes(Family, n, fill=Phylum)) +
  geom_bar(stat='identity') +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)
  )
p



