### Codivergence analysis adapted following the method of Youngblut et al. 2019. Nat Comms ###

# Paul A. O'Brien
# E. paul.obrien@my.jcu.edu.au

# 05 Find microbial co-occurence patterns using a network analysis and identify moran ASVs from script 04

setwd("~/Documents/R/Ramaciotti_16S_dataset/Codivergence_and_Neutral/Codivergence")

library(dplyr) #
library(tidyr) #
library(ggplot2) #
library(tidygraph) #
library(igraph) #
library(phyloseq) #
library(networkD3)
library(cooccur) #
library(ggraph) #
library(cowplot)

## Load and format data ----
# phyloseq core data
physeq <- readRDS("physeq_core_all")

# Metadata
metadata <- list()
for (i in seq_along(physeq)) {
  metadata[[i]] <- as.data.frame(as.matrix(sample_data(physeq[[i]])))
}
# Get presence / abssence otu table
otu_bin <- list()
for(i in seq_along(physeq)) {
  otu_bin[[i]] <- as.data.frame(as.matrix(otu_table(physeq[[i]]))) %>% 
    apply(2, function(x) ifelse(x > 0, 1, 0)) %>% as.matrix
}

# Quick look at data
heatmap(otu_bin[[4]] %>% as.matrix %>% t, 
        cexCol = 0.001, cexRow = 0.5,
        margins=c(0.5, 8), scale='none')

## Run cooccur analysis ----

cooccur_res_l <- list()
for (i in seq_along(otu_bin)) {
  cooccur_res_l[[i]] <- cooccur(mat = otu_bin[[i]],
                                type="spp_site",
                                thresh=TRUE,
                                spp_names=TRUE)
}

saveRDS(cooccur_res_l, "cooccur_list")


## Summarise results ----

# Check data
lapply(cooccur_res_l, summary)

# Plot observed and expected cooccurences
obs.v.exp(cooccur_res_l[[1]])
obs.v.exp(cooccur_res_l[[2]])
obs.v.exp(cooccur_res_l[[3]])
obs.v.exp(cooccur_res_l[[4]])

# Extract pairwise effect sizes
ef_size <- list()
for(i in seq_along(cooccur_res_l)) {
  ef_size[[i]] <- effect.sizes(cooccur_res_l[[i]], standardized = TRUE)
  print(paste(i, "done"))
}

lapply(ef_size, dim)
lapply(ef_size, head)

saveRDS(ef_size, "effect_size_l")

# Get taxonomy
tax <- list()
for(i in seq_along(physeq)) {
  tax[[i]] <- physeq[[i]] %>%
    tax_table %>%
    as.matrix %>%
    as.data.frame %>%
    mutate(OTU = rownames(.))
}

# Combine data - only include signficant co-occurring interactions
qval_cutoff = 0.01
cooccur_res_s <- list()
for(i in seq_along(cooccur_res_l)) {
  cooccur_res_s[[i]] <- cooccur_res_l[[i]] %>% 
    prob.table %>%
    mutate(p_gt = p_gt %>% as.numeric,
           p_lt = p_lt %>% as.numeric,
           q_gt = p.adjust(p_gt, method='BH'),
           q_lt = p.adjust(p_lt, method='BH'),
           sign = ifelse(q_gt < qval_cutoff, 1, 0) + ifelse(q_lt < qval_cutoff, -1, 0),
           sign = sign %>% as.character) %>%
    inner_join(tax[[i]], c('sp1_name'='OTU')) %>%
    inner_join(tax[[i]], c('sp2_name'='OTU')) %>%
    inner_join(ef_size[[i]], c('sp1_name'='sp1',
                               'sp2_name'='sp2'))
}

# Format to view network
for(i in seq_along(cooccur_res_s)) {
  cooccur_res_s[[i]] <-  cooccur_res_s[[i]][-c(1:4)] %>%
    select(sp1_name, sp2_name, everything())
}

saveRDS(cooccur_res_s, "Co-occurence/cooccur_network_file")



## Get network metrics ----

# Load  files
cooccur_res_s <- readRDS("Co-occurence/cooccur_network_file") 
tax <- readRDS("asv_taxonomy_files")

# Format network data
ntwk <- lapply(cooccur_res_s, as_tbl_graph)

# Add taxonony to node data
for(i in seq_along(ntwk)) {
  tax[[i]] <- mutate(tax[[i]], OTU = rownames(tax[[i]]))
  ntwk[[i]] <- ntwk[[i]] %>%
    activate(nodes) %>%
    inner_join(tax[[i]], c('name'='OTU'))
} 

# Calculate weights based on effects
for(i in seq_along(ntwk)) {
  ntwk[[i]] <- ntwk[[i]] %>%
    activate(edges) %>%
    mutate(weight = ifelse(effects == 0, 1e-5, abs(effects)))
}

# Calculate centrality
for(i in seq_along(ntwk)) {
  ntwk[[i]] <- ntwk[[i]] %>%
    activate(nodes) %>%
    mutate(centrality_btw = centrality_betweenness(weights=weight),
           centrality_hub = centrality_hub(weights=weight),
           centrality_degree = centrality_degree(weights=weight))
}

# Summary
ntwk[[4]] %>% 
  activate(nodes) %>%
  as.data.frame %>%
  group_by(Kingdom, Phylum) %>%
  summarize(mean_centrality_btw = mean(centrality_btw),
            mean_centrality_hub = mean(centrality_hub),
            mean_centrality_degree = mean(centrality_degree)) %>%
  ungroup() %>%
  arrange(-mean_centrality_btw)



## Plot network ----

# Filter out random associations
ntwk_f <- list()
for(i in seq_along(ntwk)) {
  ntwk_f[[i]] <- ntwk[[i]] %>%
    activate(edges) %>%
    filter(sign != 0) %>%
    activate(nodes) %>%
    mutate(community = as.factor(group_walktrap(weights=weight))) %>%
    mutate(centrality_btw = centrality_betweenness(weights=weight),
           centrality_hub = centrality_hub(weights=weight),
           centrality_degree = centrality_degree(weights=weight),
           d = local_ave_degree()) %>%
    filter(!is.na(d))
}
ntwk_f[[4]]
saveRDS(ntwk_f, "network_filtered")

# Plot network using ggraph
p <- ntwk_f[[4]] %>%
  ggraph(layout = 'nicely') + 
  geom_edge_fan(aes(edge_colour=sign), show.legend=FALSE, alpha=0.25) + 
  geom_node_point(aes(size=centrality_btw, fill=Phylum, shape=community)) +
  scale_shape_manual('Sub-network', values=c(21, 22, 23, 24)) +
  scale_color_discrete('Phylum') +
  scale_size_continuous('Centrality betweenness', range=c(2.75,10)) +
  scale_edge_color_manual(values=c('black')) + 
  guides(fill=guide_legend(override.aes=list(shape=21, size = 3))) + # change shape in legend to something that can 'fill'
  theme_graph(base_family = 'Helvetica')
p
ggsave("co-occur_network_core_all_nicely.pdf", width = 22, height = 17, units = 'cm', dpi = 300, device = cairo_pdf)



## Format and plot metrics ----

# By microbe taxonomy
df <- list()
for(i in seq_along(ntwk_f)) {
  df[[i]] <- ntwk_f[[i]] %>%
    activate(nodes) %>%
    as.data.frame %>%
    dplyr::select(Phylum, Class, Order, Family, Genus, name, community, centrality_btw) %>%
    arrange(-centrality_btw)
}

# Plot node betweenness
p <- ggplot(df[[4]], aes(Family, centrality_btw, fill=community)) +
  geom_bar(stat='identity') +
  labs(y='Node betweenness') +
  coord_flip() +
  theme_bw() +
  theme(
    axis.text.y = element_text(size=8)
  )
p


# Taxonomy by community
df <- list()
for(i in seq_along(ntwk_f)) {
  df[[i]] <- ntwk_f[[i]] %>%
    activate(nodes) %>%
    as.data.frame %>%
    group_by(community) %>%
    mutate(total_OTUs = n()) %>%
    ungroup() %>%
    group_by(community, Phylum, Family) %>%
    summarize(n_OTUs = n() / first(total_OTUs) * 100) %>%
    group_by(Phylum) %>%
    mutate(sum_n_OTUs = sum(n_OTUs)) %>%
    ungroup() %>%
    mutate(Phylum = Phylum %>% reorder(sum_n_OTUs)) 
}

# Plot nodes by sub-network
p <- ggplot(df[[4]], aes(Family, n_OTUs, fill=Phylum)) +
  geom_bar(stat='identity') +
  labs(x='Node taxonomy', y='% of nodes in sub-network') +
  facet_grid(. ~ community) +
  coord_flip() +
  theme_bw() +
  theme(
    axis.text.y = element_text(size=8)
  )
p


# Sub-networks by host species
# Get node data as df
df_nodes <- list()
for(i in seq_along(ntwk_f)) {
  df_nodes[[i]] <- ntwk_f[[i]] %>%
    activate(nodes) %>%
    as.data.frame
  print(df_nodes[[i]] %>% head)
}

# Combine otu table with metadata & node data
otu <- list()
for(i in seq_along(physeq)) {
  otu[[i]] <- physeq[[i]] %>%
    otu_table %>%
    as.matrix %>%
    as.data.frame %>%
    mutate(OTU = rownames(.)) %>%
    gather(Sample, Count, -OTU) %>%
    group_by(Sample) %>%
    mutate(TOTAL_COUNT = sum(Count)) %>%
    ungroup() %>%
    mutate(Count = Count / TOTAL_COUNT * 100)  
}

otu_nodes <- list()
for (i in seq_along(physeq)) {
  otu_nodes[[i]] <- inner_join(otu[[i]], df_nodes[[i]], c('OTU'='name')) %>%
    inner_join(metadata[[i]], c('Sample'='SampleID'))
}

# Count distribution
# View sub-network by host group 
p <- ggplot(otu_nodes[[4]], aes(Taxonomy, Count, color=community)) +
  geom_boxplot() +
  scale_y_log10() +
  theme_bw()
p


# Grouped by host alternative 
otu_s <- list()
for(i in seq_along(otu_nodes)) {
  otu_s[[i]] <- otu_nodes[[i]] %>%
    mutate(PA = ifelse(Count > 0, 1, 0)) %>%
    group_by(community, Taxonomy, Species.y, OTU) %>%
    summarize(frac_pres = sum(PA) / n() * 100) %>%
    group_by(community, Taxonomy) %>%
    summarize(mean_frac_pres = mean(frac_pres),
              sd_frac_pres = sd(frac_pres) / sqrt(length(frac_pres))) %>%
    ungroup()
}


p <- ggplot(otu_s[[4]], aes(Taxonomy, mean_frac_pres, color=community)) +
  geom_linerange(aes(ymin=mean_frac_pres-sd_frac_pres,
                     ymax=mean_frac_pres+sd_frac_pres),
                 size=1, alpha=0.5) +
  geom_point() +
  labs(x='Host group', y='% presence among species') +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle=45, hjust=1)
  )
p

# Plot by host species
otu_s <- list()
for(i in seq_along(otu_nodes)) {
  otu_s[[i]] <- otu_nodes[[i]] %>%
    mutate(PA = ifelse(Count > 0, 1, 0)) %>%
    group_by(community, Taxonomy, Species.y, OTU) %>%
    summarize(frac_pres = sum(PA) / n() * 100) %>%
    group_by(community, Taxonomy, Species.y) %>%
    summarize(mean_frac_pres = mean(frac_pres),
              sd_frac_pres = sd(frac_pres) / sqrt(length(frac_pres))) %>%
    ungroup()
}

p2 <- ggplot(otu_s[[4]], aes(Species.y, ifelse(mean_frac_pres==0, NA, mean_frac_pres), color=community)) +
  geom_linerange(aes(ymin=mean_frac_pres-sd_frac_pres,
                     ymax=mean_frac_pres+sd_frac_pres),
                 size=0.75) +
  geom_point(aes(size = 1.25)) +
  coord_flip() +
  scale_colour_manual('sub-\nnetwork', values = c("#D55E00", "#56B4E9", "#CC79A7", "#009E73")) +
  #scale_color_discrete('sub-\nnetwork') +
  facet_grid(Taxonomy ~ ., scales='free_y', space='free_y') +
  labs(x='', y='Mean presence within species (%)') +
  theme_bw() 
p2
ggsave("network_community_by_host_species.eps",  width = 18, height = 20, units = 'cm', dpi = 300)

# Combine above plot with network plot

p3 <- ggdraw() +
  draw_plot(p, x = 0, y = .01, width = .60 , height = .99) + # change 'p' to netork plot above
  draw_plot(p2, x = .61, y = .05, width = .35, height = .90)
p3
ggsave("co-occurence_combined_plot.pdf", width = 35, height = 20, units = 'cm', dpi = 300, device = cairo_pdf)



## Subnetwork stats ----

# density = total edges / all possible edges
# Get sub-network IDs
comms <- list()
for(i in seq_along(ntwk_f)) {
  comms[[i]] = ntwk_f[[i]] %>%
    activate(nodes) %>%
    as.data.frame %>%
    .$community %>%
    unique %>% as.numeric 
}

# Calculate network properties
calc_subnetwork = function(commID, g){
  x = g %>%
    activate(nodes) %>%
    filter(community == commID) %>%
    mutate(
      adhesion = graph_adhesion(),
      phylum_assort = graph_assortativity(Phylum, directed = FALSE),
      class_assort = graph_assortativity(Class, directed = FALSE),
      genus_assort = graph_assortativity(Genus, directed = FALSE),
      component_count = graph_component_count(),
      min_cut = graph_min_cut(),
      mean_dist = graph_mean_dist(),
      graph_size = graph_size(),
      graph_order = graph_order(),
      graph_density = graph_size() / (graph_order() * graph_order()),
      edge_density = graph_size() / graph_order(),
      rel_neighbors = mean(centrality_degree() / graph_order()),
      max_centrality = max(centrality_degree()),
      rel_btw = mean(centrality_betweenness() / graph_order())) %>%
    as.data.frame() %>%
    distinct(adhesion, phylum_assort, class_assort, genus_assort, min_cut, 
             component_count, mean_dist, graph_size, graph_order, 
             graph_density, edge_density, rel_neighbors, max_centrality, rel_btw)
  x$community = commID
  return(x)
}

# summary by community
comm_stats <- list()
for(i in seq_along(comms)) {
  comm_stats[[i]] = lapply(as.list(comms[[i]]), calc_subnetwork, g=ntwk_f[[i]]) %>%
    do.call(rbind, .) %>%
    arrange(community)
}
comm_stats[[4]]

# formatting
comm_stats_s <- list()
for(i in seq_along(comm_stats)) {
  comm_stats_s[[i]] = comm_stats[[i]] %>%
    dplyr::select(community, graph_order, graph_density, max_centrality)
}
comm_stats_s[[4]]

write.table(comm_stats_s[[1]], file="coral_subntwk_stats.tsv", sep='\t', quote=FALSE, row.names=FALSE)
write.table(comm_stats_s[[2]], file="octocoral_subntwk_stats.tsv", sep='\t', quote=FALSE, row.names=FALSE)
write.table(comm_stats_s[[3]], file="sponge_subntwk_stats.tsv", sep='\t', quote=FALSE, row.names=FALSE)
write.table(comm_stats_s[[4]], file="all_subntwk_stats.tsv", sep='\t', quote=FALSE, row.names=FALSE)


## Moran OTUs ---- 

# load network file
ntwk_f <- readRDS("Co-occurence/network_filtered")

# Load lipa file with significant OTUs
moran_asvs <- readRDS("~/Documents/R/Ramaciotti_16S_dataset/Codivergence_and_Neutral/Codivergence/Phylosignal/Results/lipa_sig_asvs_all_sens_q50")
moran_asvs

# Remove X's from beginning of ASV (needed in earlier script)
lipa_asvs$ASV <- gsub('^X', '', lipa_asvs$ASV)

# Grab OTU names
sig_OTUs <- lipa_asvs$ASV %>% unique

# Combine with network file
ntwk_f_moran <- ntwk_f[[4]] %>%
  activate(nodes) %>%
  as.data.frame %>%
  group_by(community) %>%
  mutate(total_OTUs = name %>% unique %>% length) %>%
  ungroup() %>%
  mutate(Moran = name %in% sig_OTUs) # Change to lipa or physig 
ntwk_f_moran %>% head

# Summarise
ntwk_sum_tbl <- ntwk_f_moran %>%
  group_by(community) %>%
  summarize(comm_n_OTUs = name %>% unique %>% length,
            Moran_n_OTUs = sum(Moran)) %>%
  ungroup() %>%
  mutate(perc_Moran = Moran_n_OTUs / comm_n_OTUs * 100)
ntwk_sum_tbl
write.csv(ntwk_sum_tbl, "summary_table_network_analysis.csv")











