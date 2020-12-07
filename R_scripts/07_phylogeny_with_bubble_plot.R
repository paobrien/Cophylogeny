### Plot codivergence figure ###

# Paul A. O'Brien
# E. paul.obrien@my.jcu.edu.au

# 07 Plot microbe phylogeny with relative abundance bubble plot

setwd("~/Documents/R/Ramaciotti_16S_dataset/Codivergence_and_Neutral/Codivergence")

library(ggplot2)
library(ggtree)
library(reshape2)
library(dplyr)
library(cowplot)


## Note ----

# This script uses the ggtree function 'facet_widths' which can be found at:
# https://rdrr.io/github/GuangchuangYu/ggtree/src/R/facet-utilities.R


## Import and format data ----

# load data
taxonomy <- readRDS("asv_taxonomy_files")
physeq <- readRDS("asv_sw_rarefied_phyloseq_all") # with seawater
metadata <- readRDS("metadata_sw") # with seawater

# Subset trees of interest
endo <- subset_taxa(physeq, Family=="Endozoicomonadaceae") %>% phy_tree
nitro <- subset_taxa(physeq, Family=="Nitrosopumilaceae") %>% phy_tree
thermo <- subset_taxa(physeq, Family=="Thermoanaerobaculaceae") %>% phy_tree
sp_mi <- subset_taxa(physeq, Family=="Spirochaetaceae"| Family=="Microtrichaceae") %>% phy_tree ## combine spiro and micro for figure
trees <- list(endo, nitro, thermo, sp_mi)

# Get tip labels
tip_labs <- list()
for(i in seq_along(trees)) {
  tip_labs[[i]] <- trees[[i]]$tip.label
}

# Subset asv data to tip labels
asv <- otu_table(physeq) %>% as.matrix %>% as.data.frame %>% t
asv_sub <- list()
for(i in seq_along(tip_labs)) {
  asv_sub[[i]] <- asv[,which(colnames(asv) %in% tip_labs[[i]])]
  print(ncol(asv_sub[[i]]))
}


## Calculate and plot relative abundance ----

# Calculate relative abundance 
asv_rel <- list()
for (i in seq_along(asv_sub)){
  asv_rel[[i]] <- asv_sub[[i]][,1:ncol(asv_sub[[i]])] / rowSums(asv_sub[[i]][,1:ncol(asv_sub[[i]])]) * 100 
  asv_rel[[i]][is.nan(asv_rel[[i]])] <- 0 # Change NaNs to 0  (caused by dividing 0 by 0)
  print(asv_rel[[i]][1:5, 1:5])
}

# Add metadata
asv_rel <- lapply(asv_rel, cbind, metadata)

# reshape to long format
asv_long <- list()
for (i in seq_along(asv_rel)){
  asv_long[[i]] <- melt(asv_rel[[i]], variable.name = "ASV")
  colnames(asv_long[[i]])[colnames(asv_long[[i]]) == "variable"] <- "ASV"  
  print(head(asv_long[[i]]))
}

# Group by taxonomy and summarise
group_taxon <- list()
for (i in seq_along(asv_long)){
  group_taxon[[i]] <- asv_long[[i]] %>% group_by(Species, Taxonomy, ASV)
  print(head(group_taxon[[i]]))
}

taxon_summary <- list()
for (i in seq_along(group_taxon)){
  taxon_summary[[i]] <- summarise(group_taxon[[i]], Abundance = mean(value, na.rm = T))
  print(head(taxon_summary[[i]]))
}

tax_sum.df <- lapply(taxon_summary, as.data.frame)  
lapply(tax_sum.df, head)


# Relevel factors for plot (by host phylogeny)
unique(tax_sum.df[[1]]$Species) # for species names

for(i in seq_along(tax_sum.df)) {
  tax_sum.df[[i]]$Species <- factor(tax_sum.df[[i]]$Species, levels = c("S. hystrix", "S. hystrix_CS", "S. pistillata", "P. damicornis", "P. verrucosa", "D. heliopora", 
                                                                        "E. mammiformis", "P. cactus", "P. speciosa", "A. hyacinthus", "A. formosa", 
                                                                        "P. cylindrica", "P. cylindrica_CS", "P. massive", "Sarcophyton sp", "Sarcophyton sp_CS", 
                                                                        "Sinularia sp", "Sinularia sp2", "Sinularia sp_CS", "Cladiella sp", "Pinnigorgia sp", "I. hippuris", 
                                                                        "Briareum sp", "Briareum sp2", "Heteroxenia sp", "Clavularia sp", 
                                                                        "I. ramosa_CS", "I. ramosa", "Coscinoderma sp", "Ircinia sp", "C. foliascens", "C. singaporensis", "Seawater"))
}


for(i in seq_along(tax_sum.df)) {
  tax_sum.df[[i]]$Taxonomy <- factor(tax_sum.df[[i]]$Taxonomy, levels = c("Coral", "Octocoral", "Sponge", "Seawater"))
}

# Check bubble plot
bubble <- ggplot(tax_sum.df[[1]], aes(x = Species, y = ASV, color = Taxonomy)) +
  geom_point(aes(size = ifelse(Abundance==0, NA, Abundance))) + 
  scale_size(range = range(taxon_summary$Perentage)) +
  scale_size_area(max_size = 10) + 
  scale_colour_brewer(palette = "Dark2") +
  labs(x = "\nTaxa", y = "Relative Abundance (Phyla)\n") +
  guides(colour = guide_legend(title = "Taxa")) +
  theme_bw() +
  theme(axis.text.y=element_blank(),axis.title.y=element_blank(),  axis.ticks.y=element_blank(),
        legend.text=element_text(size=11), legend.title=element_text(size=13))
bubble



## Plot tree with bubble plot ----

# Note: ASV ids MUST be in first column
for(i in seq_along(tax_sum.df)) {
  tax_sum.df[[i]] <-tax_sum.df[[i]][,c(3,1,2,4)] 
}

# set colours to colour blind friendly
clr <- c("Coral" = "#0072B2", "Octocoral" = "#009E73", "Sponge" = "#E69F00", "Seawater" = "#CC79A7")

# Get Moran asvs to highlight in figure
moran_asvs <- readRDS("~/Documents/R/Ramaciotti_16S_dataset/Codivergence_and_Neutral/Codivergence/Phylosignal/Results/lipa_sig_asvs_all_sens_q50")
moran_asvs$ASV <- gsub('X([0-9])', '\\1', moran_asvs$ASV) # remove X's added in other script
moran_asvs <- moran_asvs$ASV %>% unique


# Plot spiro / micro
p <- ggtree(sp_mi, branch.length = "none", aes(color = label %in% moran_asvs_asvs)) + geom_hilight(node=1, fill="steelblue", alpha=.6) # plot microbe phylogeny
p2 <- facet_plot(p + xlim_tree(0), panel = "Relative Abundance", data = tax_sum.df[[4]], geom = geom_point,                      # add bubble plot
                 mapping = aes(x = Species, size = ifelse(Abundance==0, NA, Abundance), color = Taxonomy)) + 
  scale_colour_manual(values = c("#0072B2", "black", "#009E73", "#CC79A7", "#E69F00", "red")) +
  scale_x_discrete() + 
  scale_y_discrete() +
  theme_bw() +
  theme(axis.title.y=element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "none") + labs(size="Relative Abundance")
p2 <- facet_labeller(p2, c(Tree = "Cladogram")) %>% facet_widths(c(1, 1.5))                                                     # Change facet plot widths
p2

# Plot endo
p <- ggtree(endo, branch.length = "none", aes(color = label %in% moran_asvs))
p3 <- facet_plot(p + xlim_tree(0), panel = "Relative Abundance", data = tax_sum.df[[1]], geom = geom_point,
                 mapping = aes(x = Species, size = ifelse(Abundance==0, NA, Abundance), color = Taxonomy)) + 
  scale_colour_manual(values = c("#0072B2", "black", "#009E73", "#CC79A7", "#E69F00", "red")) +
  scale_x_discrete() + 
  scale_y_discrete() +
  guides(color=guide_legend(title="Host")) +
  theme_bw() +
  theme(axis.title.y=element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.text=element_text(size=11), legend.title=element_text(size=12)) + labs(size="Relative Abundance")
p3 <- facet_labeller(p3, c(Tree = "Cladogram")) %>% facet_widths(c(1, 1.5)) 
p3

# Combine to one plot
p4 <- ggdraw() +
  draw_plot(p2, x = 0, y = .01, width = .41 , height = .99) +
  draw_plot(p3, x = .41, y = .01, width = .59, height = .99)
p4

# Save and finalise in illustrator
ggsave(filename = "spiro_micro_endo_den_bubble.eps", device = "eps", width = 30, height = 30, units = "cm", dpi = "print")


# Plot nitro
p <- ggtree(nitro, branch.length = "none", aes(color = label %in% moran_asvs))
p2 <- facet_plot(p + xlim_tree(0), panel = "Relative Abundance", data = tax_sum.df[[2]], geom = geom_point,
                 mapping = aes(x = Species, size = ifelse(Abundance==0, NA, Abundance), color = Taxonomy)) + 
  scale_colour_manual(values = c("#0072B2", "black", "#009E73", "#CC79A7", "#E69F00", "red")) +
  scale_x_discrete() + 
  scale_y_discrete() +
  guides(color=guide_legend(title="Host")) +
  theme_bw() +
  theme(axis.title.y=element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "none") + labs(size="Relative Abundance")
p2 <- facet_labeller(p2, c(Tree = "Cladogram")) %>% facet_widths(c(1, 1.5)) 
p2

# Plot thermo
p <- ggtree(thermo, branch.length = "none", aes(color = label %in% moran_asvs))
p3 <- facet_plot(p + xlim_tree(0), panel = "Relative Abundance", data = tax_sum.df[[3]], geom = geom_point,
                 mapping = aes(x = Species, size = ifelse(Abundance==0, NA, Abundance), color = Taxonomy)) + 
  scale_colour_manual(values = c("#0072B2", "black", "#009E73", "#CC79A7", "#E69F00", "red")) +
  scale_x_discrete() + 
  scale_y_discrete() +
  guides(color=guide_legend(title="Host")) +
  theme_bw() +
  theme(axis.title.y=element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.text=element_text(size=11), legend.title=element_text(size=12)) + labs(size="Relative Abundance")
p3 <- facet_labeller(p3, c(Tree = "Cladogram")) %>% facet_widths(c(1, 1.5)) 
p3

# Combine to one plot
p4 <- ggdraw() +
  draw_plot(p2, x = 0, y = .01, width = .41 , height = .99) +
  draw_plot(p3, x = .41, y = .01, width = .59, height = .99)
p4

# save and finalise in illustrator 
ggsave(filename = "nitro_thermo_den_bubble.eps", device = "eps", width = 30, height = 30, units = "cm", dpi = "print")







