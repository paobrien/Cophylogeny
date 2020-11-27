### Plot codivergence figure ###

# Paul A. O'Brien
# E. paul.obrien@my.jcu.edu.au

# 06 Plot presence/absence heatmap of microbial OTUs in their hosts. Rows arranged by OTU phylogeny and columns by host phylogeny

setwd("~/Documents/R/Ramaciotti_16S_dataset/Codivergence_and_Neutral/Codivergence")

library(phyloseq)
library(dendextend)
library(ape)

## Load and format data ----

# Phyloseq data
physeq <- readRDS("physeq_core_all") # Core ASVs
physeq <- physeq[[4]] # all samples

# ASV taxonomy
tax <- readRDS("asv_taxonomy_files")
tax[[4]] <- tax[[4]][-c(250),] # Removed as not in tree / analysis


# Metadata
metadata <- readRDS("metadata_files")

# Host phylogeny 
host_tree <- read.tree("all_samples_tree_final")

# ASV phylogeny
bact_tree <- phy_tree(physeq)


## Plot heatmap of core ASVs vs Host Phylogeny ----

# Make host dendrogram
host_den <- cophenetic(host_tree) %>%
  as.dist %>% hclust %>% as.dendrogram
plot(host_den, horiz = TRUE)


# Make bacteria dendrogram

bact_den <- cophenetic(bact_tree) %>%
  as.dist %>% hclust %>% as.dendrogram
plot(bact_den, horiz = T)

# Get asv table and change to binary
asv <- as.data.frame(as.matrix(otu_table(physeq)))
asv_bin <- ifelse(asv > 0, 1, 0)

# Reorder matrix according to phylogeny
asv_ord <- asv_bin[bact_tree$tip.label, host_tree$tip.label] # reorder row by bact phylogeny & column by host phylogeny

# Highlight Moran ASVs
moran_asvs <- readRDS("~/Documents/R/Ramaciotti_16S_dataset/Codivergence_and_Neutral/Codivergence/Phylosignal/Results/lipa_sig_asvs_all_sens_q50") # import asvs
moran_asvs$ASV <- gsub('X([0-9])', '\\1', moran_asvs$ASV) # remove X's added in other script
moran_asvs <- moran_asvs$ASV %>% unique # get unique
bact_den <- bact_den %>% set("by_labels_branches_col",  value = moran_asvs, type = "any") # set branch colour
plot(bact_den)

# Get host dendrogram colours
host_den_c <- color_branches(host_den, k = 3) # For coloured host
plot(host_den_c)

# Check bacteria taxomony
tmp <- tax[[4]]
levels(tmp$Phylum) # check levels

# Add new factor level for unclassified
tmp$Phylum <- factor(tmp$Phylum, levels=c(levels(tmp$Phylum), "Unclassified"))

# Convert all NA's to Unclassified
tmp$Phylum[is.na(tmp$Phylum)] = "Unclassified"

# Check levels again
levels(tmp$Phylum)
tmp$Phylum %>% unique()

# Change label colours
labels_colors(bact_den) <- rainbow_hcl(n = 21, start = 100, end = 350)[sort_levels_values(
  as.numeric(tmp[,2])[order.dendrogram(bact_den)])]

# Change labels to phylum classification
labels(bact_den) <- paste(as.character(tmp[,2])[order.dendrogram(bact_den)],
                          "(",labels(bact_den),")", 
                          sep = "")

# Check dendrogram
phylum_leg <- tmp[,2][order.dendrogram(bact_den)] %>% unique
plot(bact_den, horiz = T)
legend("topleft", legend = phylum_leg, fill = rainbow_hcl(n = 21, start = 100, end = 350))

# Get the color codes to add colour strip
#host
col_labels_h <- get_leaves_branches_col(host_den_c)
col_labels_h <- col_labels_h[order(order.dendrogram(host_den_c))]

# bact
col_labels_b <- labels_colors(bact_den)
col_labels_b <- col_labels_b[order(order.dendrogram(bact_den))]

# Plot heatmap - binary
heatmap(as.matrix(asv_ord), scale = "row",  Rowv = bact_den, Colv =  host_den, labRow = "", labCol = "",
        col = c('black', 'orange'), ColSideColors = col_labels_h, RowSideColors = col_labels_b)

# Add heatmap legend
legend(x="topleft", legend=c("Present", "Absent"), 
       fill=c('orange', 'black'))

# add bacteria legend
legend("topright", legend = rev(phylum_leg), fill = rev(rainbow_hcl(n = 21, start = 100, end = 350)))

# Save plot
save("heatmap_with_dendrogram.eps", device = "eps", width = 20, height = 20, units = "cm", dpi = "print")

