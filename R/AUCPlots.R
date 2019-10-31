library(plotROC)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tmod)
source("./R/GeneSetBuilders.R")

adult_target_region <- "insula"
inclusionPattern <- adult_target_region
neoCortexOnly <- FALSE
geneSetName = "DisGeNet"
targetGroup <- "Cocaine-Related Disorders"

baseFolder <- paste0(inclusionPattern, ".neocortex.", neoCortexOnly)

matrixAdult <- read_tsv(paste0("./data/processed/", baseFolder, "/allen_HBA_brainarea_vs_genes_exp_qcNames.tsv"), guess_max=10000)
matrixFetal <- read_tsv(paste0("./data/processed/", baseFolder, "/allen_human_fetal_brain_brainarea_vs_genes_exp_qcNames.tsv"), guess_max = 10000) 

matrixFetal %<>% select(gene_symbol, contains(inclusionPattern)) %>% 
  gather(key="region", value="rank", contains(inclusionPattern)) %>% 
  mutate(source = "Fetal")
matrixAdult %<>% select(gene_symbol, contains(inclusionPattern)) %>% 
  gather(key="region", value="rank", contains(inclusionPattern)) %>% 
  mutate(source = "Adult")
combined <- bind_rows(matrixAdult, matrixFetal)

#redo rank to end at max gene count (no negative values)
combined %<>% group_by(source, region) %>% mutate(rank = rank(rank)) %>% arrange(-rank)

if (geneSetName == "PhenoCarta") {
  geneSetsToUse <- loadPhenocarta("human", allGenes)
  filterGenes <- T
} else if (geneSetName == "CYP") {
  geneSetsToUse <- loadCypGenes(allGenes)
  filterGenes <- F
} else if (geneSetName == "Custom") {
  geneSetsToUse <- loadFileSets(prefix="Custom")
  filterGenes <- F
} else if (geneSetName == "Darmanis") {
  geneSetsToUse <- loadFileSets(prefix="Darmanis")
  filterGenes <- F
} else if (geneSetName == "Extra") {
  geneSetsToUse <- loadFileSets(prefix="Extra")
  filterGenes <- F
} else if (geneSetName == "DisGeNet") {
  geneSetsToUse <- loadDisGeNetSets(allGenes)
  filterGenes <- T
} else if (geneSetName == "GO") {
  if (exists("geneSetsGO") && length(geneSetsGO$MODULES2GENES) > 1000 ) { #assume it's already loaded - needs a fix to see if the variable is declared
  } else {
    geneSetsGO <- loadGOSets(allGenes)
  }
  geneSetsToUse <- geneSetsGO
  filterGenes <- T
}

targetGroupID <- dplyr::filter(tbl_df(geneSetsToUse$MODULES), Title == targetGroup)$ID
genesOfInterest <- unlist(geneSetsToUse$MODULES2GENES[targetGroupID])

combined %<>% mutate(present = gene_symbol %in% genesOfInterest * 1)
combined %<>% mutate(combinedLabel = paste(source, region))

cbbPalette <- c("#E69F00", "#000000", "#000000", "#E69F00")

(AUCPlot <- ggplot(combined, aes(d = present, m = rank, color=combinedLabel)) + ylab("") + 
    #geom_roc(n.cuts=0, linetype = "dashed") + 
    geom_roc(data = combined %>% filter(source == "Fetal"), n.cuts=0, linetype = "dotted") + 
    geom_roc(data = combined %>% filter(source == "Adult"), n.cuts=0, linetype = "solid") + 
    style_roc() + coord_cartesian(expand=F) +
    theme(legend.position = c(1,0), legend.justification = c(1, 0), legend.background= element_rect(fill = "transparent", colour = "transparent"), plot.margin=unit(c(.5,.5,.5,.5),"cm")) + 
    labs(color='Gene Group')  + 
    #facet_grid(dummy ~ ., switch="y")  + ylab("") +
    theme(strip.background = element_blank(), strip.placement = "inside", strip.text = element_blank()) +
    geom_abline(slope = 1, intercept = 0, colour = "grey90", size = 0.2) +
    scale_fill_manual(values=cbbPalette) + scale_colour_manual(values=cbbPalette) +
    guides(color = guide_legend(override.aes = list(linetype = c("solid","solid", "dotted","dotted")))) +
    theme(legend.key.width = unit(2, "line"))
    )
#save as 6x6 pdf


throwawayAUC <- combined
throwawayAUC$rank <- -1*throwawayAUC$rank
throwawayAUC$combinedLabel <- factor(throwawayAUC$combinedLabel, levels= c("Adult long insular gyri", "Adult short insular gyri",  "Fetal dysgranular insular cortex", "Fetal granular insular cortex"))

throwawayAUC$linetype <-"dotted"

(rasterPlot <- ggplot(throwawayAUC, aes(x = rank, y = present, color= combinedLabel)) + 
    geom_blank() + 
    geom_vline(data = filter(throwawayAUC, present == 1), aes(xintercept=rank, color=combinedLabel, linetype=combinedLabel)) + #,color="black") + #, size=0.07) + 
    theme_bw()+coord_cartesian(expand=F) +
    ylab("Transcriptomic cell type") + 
    facet_wrap(~combinedLabel, strip.position="top",ncol=1) + #, switch = "both"
    theme(strip.background = element_blank(), strip.placement = "inside") + #, strip.text.y = element_text(angle = 180)) +
    theme(axis.title.y = element_blank(),  axis.text.y=element_blank(), axis.ticks.y=element_blank(),axis.ticks.x=element_blank()) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_x_continuous(name = paste0("Insula Specific Expression (",length(unique(throwawayAUC$gene_symbol))," genes)"), breaks= c(min(throwawayAUC$rank)+700, max(throwawayAUC$rank)-700), labels = c("Enriched", "Depleted")) +
    scale_colour_manual(values=cbbPalette) +
    scale_linetype_manual(values = c("solid","solid", "dashed","dashed")) + 
    guides(color=FALSE, linetype=FALSE)
)

(bothPlots <- plot_grid(AUCPlot, rasterPlot, nrow = 2, align = "v", rel_heights=c(1,0.6),scale = 0.95, labels = c("A", "B"))) #add labels = c("A", "B"), for manuscript
#save as 6x8 PDF  -remove cords fixed?
#poster 5x9 PDF
(bothPlots <- plot_grid(AUCPlot, rasterPlot, nrow = 2, align = "v", rel_heights=c(1,0.45),scale = 1)) #add labels = c("A", "B"), for manuscript
