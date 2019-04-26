library(GO.db)
library(org.Hs.eg.db)
library(annotate)
library(readr)
library(tmod)
library(magrittr)
library(foreach)
library(doMC)
registerDoMC(cores=7)
detach("package:dplyr", unload=TRUE)
library(dplyr)

source("./R/GeneSetBuilders.R")

#sourceExpression <- "adult"
sourceExpression <- "fetal"

#geneSetName <- "CYP" 
#geneSetName <- "GO" 
#geneSetName <- "PhenoCarta" 
#geneSetName <- "DisGeNet"
#geneSetName <- "Custom" 
geneSetName <- "Extra" 
#geneSetName <- "Darmanis" 

#inclusionPattern <- "insula"
inclusionPattern <- "claustrum"
neoCortexOnly <- FALSE

baseFolder <- paste0(inclusionPattern, ".neocortex.", neoCortexOnly)

#automatically get the target regions
if (sourceExpression == "fetal") {
  targetRegions <- gsub(".allen_human_fetal_brain.csv","",list.files(paste0("./results/", baseFolder, "/limma/"), pattern = ".allen_human_fetal_brain.csv"))
} else {
  targetRegions <- gsub(".allen_HBA.csv","",list.files(paste0("./results/", baseFolder, "/limma/"), pattern = ".allen_HBA.csv"))
}
#targetRegions <- c("medial habenular nucleus", "lateral habenular nucleus")

if (sourceExpression == "adult") {
  matrix <- read_tsv(paste0("./data/processed/", baseFolder, "/allen_HBA_brainarea_vs_genes_exp_qcNames.tsv"), guess_max=10000)
} else if (sourceExpression == "fetal") {
  matrix <- read_tsv(paste0("./data/processed/", baseFolder, "/allen_human_fetal_brain_brainarea_vs_genes_exp_qcNames.tsv"), guess_max = 10000) 
}

allGenes <- unique(matrix$gene_symbol)

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

print(paste("Number of annotated genes:", length(unique(geneSetsToUse$GENES$ID))))
#complex to handle collapsing
print(paste("Number of gene sets:", 
            nrow(as_tibble(data.frame(ID = names(geneSetsToUse$MODULES2GENES))) %>% rowwise() %>% 
                   mutate(genes = paste(sort(unlist(geneSetsToUse$MODULES2GENES[ID])), collapse = " ")) %>% 
                   select(genes) %>% distinct())))

#get AUC values for every region and gene set
regionNames <- setdiff(colnames(matrix), "gene_symbol")

#multithreaded
allResults<-foreach(region=regionNames,.combine=rbind) %dopar% {
  matrix %<>% arrange_(paste0("desc(`",region,"`)"))
  sortedGenes <- matrix$gene_symbol
  
  result <- tbl_df(tmodUtest(c(sortedGenes), mset=geneSetsToUse, qval = 1, filter = filterGenes))
  result %<>% select(ID, AUC) %>% mutate(region = region)
  result
}

#take the AUC's for the region, count how many regions have higher AUC values 
#targetRegion <- "lateral habenular nucleus"
#targetRegion <- "medial habenular nucleus"

for(targetRegion in targetRegions) {
  #add column for higher AUC than the target regions
  targetResults <- allResults %>% filter(region == targetRegion) %>% select(ID, targetAUC = AUC)
  targetResults <- inner_join(allResults,targetResults)
  targetResults %<>% mutate(betterAUC = AUC >= targetAUC)
  summary <- targetResults %>% 
    filter(region != targetRegion) %>% 
    group_by(ID) %>% 
    summarize(betterAUCCount = sum(betterAUC), n=dplyr::n())
  summary %<>% mutate(specificity = betterAUCCount/n) %>% arrange(specificity)
  #summary %<>% filter(betterAUCCount < 1)
  
  #join with original
  matrix %<>% arrange_(paste0("desc(`", targetRegion,"`)"))
  sortedGenes <- matrix$gene_symbol
  result <- tbl_df(tmodUtest(c(sortedGenes), mset=geneSetsToUse, qval = 1.01, filter = filterGenes))

  
  #running as a one sided test
  result %<>% mutate(P.Value = if_else(AUC < 0.5, 1 - P.Value, P.Value)) %>% mutate(adj.P.Val = p.adjust(P.Value, method="fdr")) %>% arrange(P.Value)
  result %<>% rowwise() %>% mutate(aspect = Ontology(ID))
  result <- inner_join(result, summary %>% select(-n)) 
  
  
  #merge same GO groups using the genes in the set
  result %<>% rowwise() %>% mutate(genes = paste(sort(unlist(geneSetsToUse$MODULES2GENES[ID])), collapse = " "))
  
  result %<>% ungroup() %>% group_by(genes, N1) %>% arrange(Title) %>% 
    summarize(MainTitle = first(Title),  ID=paste(ID, collapse=","), AUC = first(AUC), P.Value= first(P.Value), betterAUCCount = first(betterAUCCount), specificity = first(specificity), aspect= first(aspect), otherNames = if_else(dplyr::n() > 1, paste(Title[2:length(Title)], collapse=", "), ""))
  result %<>% ungroup() %>% select(-genes)
  
  #adjust again
  result %<>% ungroup() %>% mutate(adj.P.Value=p.adjust(P.Value, method = "fdr")) %>% arrange(P.Value)
  result %<>% mutate(rank = rank(P.Value))
  result %<>% dplyr::select(MainTitle, geneCount = N1, AUC, P.Value, adj.P.Value, betterAUCCount, everything()) 
  
  write_tsv(result, paste0("./results/", baseFolder, "/", sourceExpression, "-", targetRegion, "-", geneSetName, ".tsv"))
}


# AUCDiff <- allResults %>% filter(ID %in% c("CTRA.up", "CTRA.down"))
# AUCDiff <- tidyr::spread(AUCDiff, ID, AUC) 
# AUCDiff %<>% mutate(diff = CTRA.down - CTRA.up)
# AUCDiff %>% arrange(diff) %>% head()
# AUCDiff %>% arrange(diff) %>% tail()
# #as.data.frame(AUCDiff %>% arrange(-diff) )
# AUCDiff %<>% mutate(rank = rank(diff), rankDown = rank(CTRA.down))
# AUCDiff %>% filter(grepl("habenula", region))


