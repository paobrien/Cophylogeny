### Codivergence analysis adapted following the method of Youngblut et al. 2019. Nat Comms ###

# Paul A. O'Brien
# E. paul.obrien@my.jcu.edu.au

# 04 Use Moran's I (Phylosignal) to find OTUs clustered by host species

setwd("~/Documents/R/Ramaciotti_16S_dataset/Codivergence_and_Neutral/Codivergence/Phylosignal")

library(dplyr)
library(tidyr)
library(ggplot2)
library(ape)
library(phyloseq)
library(phylosignal)
library(doParallel)

## Functions ----

# From https://github.com/leylabmpi/animal_gut_16S-uni (Youngblut et al. 2019)

## Phylosignal analysis on each tree subsample
phyloSignal_each <- function(i, tree_4d, methods='all', reps=99){
  # phylosignal
  physig_res = phyloSignal(tree_4d[[i]], methods=methods, reps=reps)
  
  # formatting output
  tmp1 = physig_res$stat 
  tmp1$OTU = rownames(tmp1)
  
  # calculating qvalue on results
  tmp2 = physig_res$pvalue
  tmp3 = apply(tmp2, 2, function(x) p.adjust(x, method='BH')) %>% as.data.frame
  tmp2$OTU = rownames(tmp2)
  tmp3$OTU = rownames(tmp2)
  
  tmp1 = tmp1 %>%
    gather(method, coef, -OTU)
  tmp2 = tmp2 %>%
    gather(method, pvalue, -OTU)
  tmp3 = tmp3 %>%
    gather(method, qvalue, -OTU)
  
  tmp1 %>%
    inner_join(tmp2, c('OTU', 'method')) %>%
    inner_join(tmp3, c('OTU', 'method')) %>%
    mutate(subsample_rep = i)
}

## Randomly select one sample per group
phylo4d_subsample <- function(L, df, otu, tree){
  # get subsample (one sample per species)
  df = df %>% 
    group_by(Species) %>% 
    sample_n(1)
  #df %>% head %>% print
  # getting OTU
  otu = otu[,df$sample] %>% t 
  # subsampling tree
  to_rm = setdiff(tree$tip, rownames(otu))
  tree = drop.tip(tree, to_rm)
  # creating phylo4d 
  tree_4d = phylobase::phylo4d(tree, tip.data=otu)
  return(tree_4d)
}


# Lipamoran analysis on each trait
lipamoran_each <- function(i, tree_4d, traits, reps=99){
  # phylosignal
  doParallel::registerDoParallel(threads)
  lipa_res = plyr::llply(as.list(traits), lipaMoran_per_OTU, tree_4d=tree_4d[[i]], reps=reps, .parallel=TRUE)
  # formatting
  lipa_res = do.call(rbind, lipa_res)
  lipa_res$subsample_rep = i
  return(lipa_res)
}


lipaMoran_per_OTU = function(trait, tree_4d, reps=9999){
  res = lipaMoran(tree_4d, trait=trait, reps=reps, prox.phylo = "nNodes")
  z = colnames(res$lipa)[1]
  x = res$lipa
  colnames(x) = c('coef')
  y = res$p.value 
  colnames(y) = c('pvalue')
  df = cbind(x,y) %>% as.data.frame
  df$OTU = z
  df$host = rownames(df)
  rownames(df) = 1:nrow(df)
  return(df)
}

## Load and format data ----

filepath <- "~/Documents/R/Ramaciotti_16S_dataset/Codivergence_and_Neutral/Codivergence/"

# Phyloseq file - core asvs
physeq <- readRDS(paste0(filepath, "physeq_core_all"))

# Host phylogeny - final trees from Niko
tree_files <- list.files(path = filepath, pattern = "tree_final")[c(2,3,4,1)]
host_trees <- list()
for (f in seq_along(tree_files)) {
  host_trees[[f]] <- read.tree(paste0(filepath, file = tree_files[[f]]))
}
summary(host_trees)

# ASV table
asv <- list()
for (i in seq_along(physeq)) {
  asv[[i]] <- as.data.frame(as.matrix(otu_table(physeq[[i]])))
  print(dim(asv[[i]]))
}

# Metadata
metadata <- list()
for (i in seq_along(physeq)) {
  metadata[[i]] <- as.data.frame(as.matrix(sample_data(physeq[[i]])))
}

# Ensure asv table and metadata are consistent
for(i in seq_along(asv)) {
  print(setdiff(colnames(asv[[i]]), rownames(metadata[[i]])))
  print(setdiff(rownames(metadata[[i]]), colnames(asv[[i]])))
}

# Taxonomy
tax <- list()
for(i in seq_along(physeq)) {
  tax[[i]] <- as.data.frame(as.matrix(tax_table(physeq[[i]])))
}

for (i in seq_along(tax)) {
  tax[[i]]$ASV <- gsub('^([0-9])', 'X\\1', rownames(tax[[i]]))  # add new ASV column. Add 'X' to any ASV that starts with 0-9 
  rownames(tax[[i]]) <- 1:nrow(tax[[i]])                                  # needed for downstream analyis
  rownames(tax[[i]]) <- 1:nrow(tax[[i]])
}

# Check ASV ID overlap
for(i in seq_along(asv)) {
  setdiff(rownames(asv[[i]]), rownames(otu_table(physeq[[i]]))) %>% 
    length %>% print
}

## Calculate residuals ----
# Scale residuals
asv_res <- list()
for(i in seq_along(asv)) {
  asv_res[[i]] <- as.data.frame(scale(asv[[i]]))
}


for(i in seq_along(asv_res)) {
  boxplot(asv_res[[i]][,1:20])
  print(summary(rowSums(asv_res[[i]])))
  print(summary(colSums(asv_res[[i]])))
}

# Save data
saveRDS(asv_res, "asv_residuals_all")


## Run phylosignal ----

# Sensitivity analysis - subsample 1 sample per species (Youngblut et al. 2019)

# Subsample tree
df <- list()
for(i in seq_along(metadata)) {
  df[[i]] <- metadata[[i]] %>%
    mutate(sample = rownames(.)) %>%
    dplyr::select(sample, Species) 
}

doParallel::registerDoParallel(2)
host_tree_4d <- list()
for (i in seq_along(df)) {
  host_tree_4d[[i]] <- plyr::llply(as.list(1:100), 
                                   function(x) phylo4d_subsample(x, df[[i]], asv_res[[i]], host_trees[[i]]),
                                   .parallel=TRUE)
  print(length(host_tree_4d[[i]]))
}

# Phylosignal calculation
methods <- c("I")   # just moran's I
doParallel::registerDoParallel(2)

physeq_sub_res <- list()
for(j in seq_along(host_tree_4d)) {
  physeq_sub_res[[j]] <- plyr::llply(as.list(1:length(host_tree_4d[[j]])),
                                     function(i) phyloSignal_each(i, host_tree_4d[[j]], methods=methods, reps=999),
                                     .parallel=TRUE)
  print(length(physeq_sub_res[[j]]))
}

for(i in seq_along(physeq_sub_res)) {
  physeq_sub_res[[i]] <- do.call(rbind, physeq_sub_res[[i]])
  print(head(physeq_sub_res[[i]]))
  print(dim(physeq_sub_res[[i]]))
}

# Add taxonomy
for(i in seq_along(physeq_sub_res)) {
  colnames(physeq_sub_res[[i]])[1] <- "ASV" # edit column name
  physeq_sub_res[[i]] <- inner_join(physeq_sub_res[[i]], tax[[i]], c('ASV')) 
  print(head(physeq_sub_res[[i]]))
}

## Find significant ASVs ----

# ASVs with sig q-value in >= 50% of reps
physeq_sub_res_s <-list()
for(i in seq_along(physeq_sub_res)) {
  physeq_sub_res_s[[i]] <- physeq_sub_res[[i]] %>%
    group_by(ASV, method, Kingdom, Phylum, Class, Order, Family, Genus) %>%
    summarize(median_coef = median(coef),
              mean_coef = mean(coef),
              sd_coef = sd(coef),
              pvalue = (length(pvalue) - sum(pvalue < 0.05)) / length(pvalue),      ## Gives proportion of p/qvalues above threshold within those replicates
              qvalue = (length(qvalue) - sum(qvalue < 0.50)) / length(qvalue)) %>%  # i.e., what proporation of ASVs had a significant value in 95% (0.05) of iterations
    ungroup()
  print(head(physeq_sub_res_s[[i]]))
}

# number of significant OTUs
n_sig_ASVs <- list()
for(i in seq_along(physeq_sub_res_s)) {
  n_sig_ASVs[[i]] <- physeq_sub_res_s[[i]] %>%
    filter(qvalue < 0.05) %>%          # Filters out the ASVs with qvalues below the threshold in the code above (^) 
    distinct(ASV) %>%                  # within a certain percentage of replicates (i.e., 0.5 = 50% of replicates[ie subsamples])
    nrow
  
  cat('Number of globally sig. ASVs:', n_sig_ASVs[[i]], '\n') 
}

# Plot
p <- ggplot(physeq_sub_res_s[[4]] %>% 
              filter(qvalue < 0.05) %>%
              mutate(Genus = Genus %>% reorder(Phylum %>% as.factor %>% as.numeric)), 
            aes(Family, fill=Phylum)) +
  geom_bar() +
  labs(y='# of OTUs') +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle=45, hjust=1)
  )
p

# Get significant ASVs
sig_asvs <- physeq_sub_res_s[[4]] %>% 
  filter(qvalue < 0.05)

# Save
saveRDS(sig_asvs, "Results/significant_phylosig_asv_sens_noBin_q50")
saveRDS(physeq_sub_res_s[[4]], "Results/phylosig_asv_all_sens_results")

sig_asvs <- readRDS("Results/significant_phylosig_asv_sens_noBin_q50")
physeq_sub_res_s <- readRDS("Results/phylosig_asv_all_sens_results")



## LipaMoran sensitivity analysism ----

# Next: the following will find which ASVs are significant locally - i.e., which host are they found in

# Subsample trees (1 sample per species)
df <- list()
for(i in seq_along(metadata)) {
  df[[i]] <- metadata[[i]] %>%
    mutate(sample = SampleID) %>%
    dplyr::select(sample, Species) 
}
doParallel::registerDoParallel(2)


# Subsample tree, 1 per species
host_tree_4d <- list()
for(i in seq_along(df)) {
  host_tree_4d[[i]] <-  plyr::llply(as.list(1:100), 
                                    function(x) phylo4d_subsample(x, df[[i]], asv[[i]], host_trees[[i]]),
                                    .parallel=TRUE)
  length(host_tree_4d[[i]]) %>% print
}

## Phylosignal calculation 

threads = 3
traits <- unique(sig_asvs$ASV) # significant ASVs from above
host_tree_4d <- host_tree_4d[[4]]

lipa_sub_res <- plyr::llply(as.list(1:length(host_tree_4d)), 
                            function(x) lipamoran_each(x, host_tree_4d, traits=traits, reps=999))

lipa_sub_res <- do.call(rbind, lipa_sub_res)
dim(lipa_sub_res) 
head(lipa_sub_res)

saveRDS(lipa_sub_res, "Results/lipa_sub_res")
lipa_sub_res <- readRDS("Results/lipa_sub_res")

# Add host data
lipa_sub_res <- lipa_sub_res %>%
  inner_join(metadata[[4]] %>% dplyr::select(SampleID, Species, Taxonomy, Zone),
             c('host'='SampleID')) 
head(lipa_sub_res)

## Find significant ASVs ----
# Sample size per group
lipa_sub_res %>%
  group_by(OTU, Species) %>%
  summarize(n=n()) %>%
  .$n %>% summary

# q-value adjustment
lipa_sub_res %>%
  group_by(subsample_rep) %>%
  mutate(pvalue = pvalue %>% as.numeric,
         qvalue = p.adjust(pvalue, method='BH')) %>%
  ungroup() %>%
  .$qvalue %>% summary

lipa_sub_res <- lipa_sub_res %>%
  group_by(subsample_rep) %>%
  mutate(pvalue = pvalue %>% as.numeric,
         qvalue = p.adjust(pvalue, method='BH')) %>%
  ungroup() 
head(lipa_sub_res)

## Find all OTUs with sig q-value in >= 50% of reps
lipa_sub_res_s <- lipa_sub_res %>%
  group_by(OTU, Species, Taxonomy) %>%
  summarize(median_coef = median(coef),
            mean_coef = mean(coef),
            sd_coef = sd(coef),
            pvalue = (length(pvalue) - sum(pvalue < 0.05)) / length(pvalue),
            qvalue = (length(qvalue) - sum(qvalue < 0.5)) / length(qvalue)) %>%
  ungroup() 
head(lipa_sub_res_s) 
dim(lipa_sub_res_s)

lipa_sub_res_s$pvalue %>% summary
lipa_sub_res_s$qvalue %>% summary
lipa_sub_res_s$OTU %>% unique %>% length

# filtering to just significant values
lipa_sub_res_s_f <- lipa_sub_res_s %>%
  filter(qvalue < 0.05) 

n_sig_OTUs <- lipa_sub_res_s_f$OTU %>% unique %>% length
cat('Number of sig. OTUs:', n_sig_OTUs, '\n')

saveRDS(lipa_sub_res_s_f, "Results/lipa_sig_asvs_all_sens_q50")
lipa_sub_res_s_f <- readRDS("Results/lipa_sig_asvs_all_sens_q50")

# Add taxonomy
names(lipa_sub_res_s_f)[names(lipa_sub_res_s_f) == "OTU"] <- "ASV"
lipa_sub_res_s_j <- lipa_sub_res_s_f %>%
  inner_join(tax[[4]], c('ASV')) 

length(unique(lipa_sub_res_s_j$Genus))

# Summary plot by host species/ bact Phylum
lipa_sub_res_s_j_s <- lipa_sub_res_s_j %>%
  distinct(Species.x, Phylum, Family, ASV, Taxonomy) %>%
  group_by(Species.x, Phylum, Family, Taxonomy) %>%
  summarize(n=n()) %>%
  ungroup()

p <- ggplot(lipa_sub_res_s_j_s, aes(Taxonomy, n, fill=Family)) +
  geom_bar(stat='identity') +
  labs(x='Host Taxonomy', y='# of OTUs') +
  facet_wrap(~ Phylum) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)
  )
p

ggsave(filename = "LIPA_asvs_by_host.eps", device = "eps", width = 25, height = 22, units = "cm", dpi = "print")








