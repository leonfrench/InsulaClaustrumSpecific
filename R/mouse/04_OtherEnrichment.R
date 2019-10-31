source("./R/GeneSetBuilders.R")
library(readr)
library(ggplot2)
library(magrittr)
library(readr)
library(tmod)
library(homologene)

detach("package:dplyr", unload=TRUE)
library(dplyr)


#geneSetName <- "PhenoCarta" 
#geneSetName <- "DisGeNet"
#geneSetName <- "Custom" 
geneSetName <- "Extra" 
#geneSetName <- "Darmanis" 

agranular_specific <- read_csv("./results/mouse/Agranular insular area.zScores.csv")
#agranular_specific <- read_csv("./results/mouse/Agranular insular area.zScores.two_experiments.csv") %>% rename(expression_energy_z = expression_energy_z_mean)

#try with sagittal only to start
print("Restricting to sagittal")
agranular_specific %<>% filter(plane == "sagittal")
agranular_specific %<>% arrange(-expression_energy_z)

human_map <- as_tibble(mouse2human(agranular_specific$entrezid))

human_agranular_specific <- left_join(agranular_specific, human_map, by=c("entrezid" = "mouseID"))

human_agranular_specific %<>% filter(!is.na(humanGene))

human_agranular_specific %>% select(humanGene, expression_energy_z) %>% group_by(humanGene) %>% summarize(n=dplyr::n()) %>% arrange(-n)
human_agranular_specific %>% select(humanGene, expression_energy_z) %>% group_by(humanGene) %>% summarize(n=dplyr::n()) %>% group_by(n) %>% summarize(count=dplyr::n())

human_agranular_specific %>% filter(humanGene == "PVALB")
#for multiple hits, pick one with highest ABS z score
human_agranular_specific %<>% select(humanGene, expression_energy_z) %>% group_by(humanGene) %>% arrange(-abs(expression_energy_z)) %>% summarize_all(first)
sortedGenes <- human_agranular_specific %>% arrange(-expression_energy_z) %>% .$humanGene


if (geneSetName == "PhenoCarta") {
  geneSetsToUse <- loadPhenocarta("human", sortedGenes)
  filterGenes <- T
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
  geneSetsToUse <- loadDisGeNetSets(sortedGenes, min=5, max=400)
  filterGenes <- T
}

result <- tbl_df(tmodUtest(c(sortedGenes), mset=geneSetsToUse, qval = 1.1, filter = T))
result %<>% rowwise() %>% mutate(P.Value = if_else(AUC < 0.5, 1 - P.Value, P.Value)) %>% ungroup() %>% mutate(adj.P.Val=p.adjust(P.Value, method="fdr")) #tmod runs one-sided tests

result %<>% arrange(P.Value)
#print top 20 at top and bottom
print(head(filter(result, AUC > 0.5) %>% dplyr::select(-ID), n=20))

write_csv(result, paste0("./results/mouse/Agranular insular area.zScores-", geneSetName, ".csv"))

