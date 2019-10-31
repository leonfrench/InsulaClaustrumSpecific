library(umap)
library(cowplot)
library(reshape2)
library(dplyr)
library(magrittr)
library(cowplot)
library(ggplot2)
library(readr)

inclusionPattern <- "insula"
neoCortexOnly <- FALSE

baseFolder <- paste0(inclusionPattern, ".neocortex.", neoCortexOnly)

set.seed(1)
#load fetal and human - files written by limma.allRegions.R
expressionMatrixAdult <- read_tsv(paste0("./data/processed/",baseFolder,"/allen_HBA_sampleMatrix_qcNames.tsv"))
#cut down to 2000 probes?

expressionMatrixFetal <- read_tsv(paste0("./data/processed/",baseFolder,"/allen_human_fetal_brain_sampleMatrix_qcNames.tsv"))

sampleAnnotAdult <- read_tsv(paste0("./data/processed/",baseFolder,"/allen_HBA_sampleMatrix_colAnnotations_qcNames.tsv"))

sampleAnnotFetal <- read_tsv(paste0("./data/processed/",baseFolder,"/allen_human_fetal_brain_sampleMatrix_colAnnotations_qcNames.tsv"))

sampleAnnotAdult %<>% filter(grepl("claustrum|insula", structure_name_left_right_stripped))
sampleAnnotFetal %<>% filter(grepl("claustrum|insula", structure_name_left_right_stripped))

expressionMatrixFetal <- expressionMatrixFetal[, c("probe_name", sampleAnnotFetal$uniqueID)]
dim(expressionMatrixFetal)

expressionMatrixAdult <- expressionMatrixAdult[, c("probe_name", sampleAnnotAdult$uniqueID)]
dim(expressionMatrixAdult)

#combine fetal and adult?
#expressionMatrixForPlot <- inner_join(expressionMatrixAdult, expressionMatrixFetal)
#sampleAnnotForPlot <- bind_rows(sampleAnnotFetal, sampleAnnotAdult)

#just adult
expressionMatrixForPlot <- expressionMatrixAdult
sampleAnnotForPlot <- sampleAnnotAdult

#just fetal
#expressionMatrixForPlot <- expressionMatrixFetal
#sampleAnnotForPlot <- sampleAnnotFetal

#generate UMAP plot
expressionMatrixForPlot %<>% select(-probe_name)
umap_embeddings <- umap(t(as.matrix(expressionMatrixForPlot)))

forPlot <- data.frame(uniqueID = colnames(expressionMatrixForPlot), stringsAsFactors = F) %>% mutate(UMAP1 = umap_embeddings$layout[,1], UMAP2 = umap_embeddings$layout[,2]) %>% as_tibble()
forPlot <- inner_join(forPlot, sampleAnnotForPlot %>% select(uniqueID, structure_name_left_right_stripped, donorID))
forPlot %<>% rename(Region = structure_name_left_right_stripped)

ggplot(forPlot, aes(x=UMAP1, y=UMAP2, color = Region)) + 
  geom_point() +theme_bw() + theme(legend.position="bottom")

forPlot %>% filter(UMAP1 < 0) %>% select(-donorID)
forPlot %>% filter(UMAP1 < 0) %>% select(donorID)
forPlot %>% filter(donorID == "normalized_microarray_donor10021")
