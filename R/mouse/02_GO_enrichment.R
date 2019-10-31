library(readr)
library(ggplot2)
library(magrittr)
#library(AnnotationDbi)
library(readr)
library(tmod)
library(annotate)
library(GO.db)
library(org.Mm.eg.db)
detach("package:dplyr", unload=TRUE)
library(dplyr)

agranular_specific <- read_csv("./results/mouse/Agranular insular area.zScores.csv")
#use confident hits 
#agranular_specific <- read_csv("./results/mouse/Agranular insular area.zScores.two_experiments.csv") %>% rename(expression_energy_z = expression_energy_z_mean)

#try with sagittal only to start
print("Restricting to sagittal")
agranular_specific %<>% filter(plane == "sagittal")

#deal with duplicates by taking the most extreme zscore
agranular_specific %<>% group_by(geneSymbol) %>% arrange(-abs(expression_energy_z)) %>% summarize_all(first)

agranular_specific %<>% arrange(-expression_energy_z)

#parameters to set:
species <- "mouse"
maxGOgroupSize <- 400
minGOgroupSize <- 5

print(paste("GO Database date:", tbl_df(GO_dbInfo()) %>% filter(name == "GOSOURCEDATE") %>% .$value))

goSource <- 'org.Mm.eg'

sortedGenes <- agranular_specific$geneSymbol 

#select GO groups
if (exists("geneSetsGO") && length(geneSetsGO$MODULES2GENES) > 1000 ) { #assume it's already loaded - needs a fix to see if the variable is declared
} else {
  go_object <- as.list(org.Mm.egGO2ALLEGS)

  symbolsInGO <- getSYMBOL(unique(unlist(go_object)), data=goSource)
  
  #build GO sets for tmod -slow
  tmodNames <- data.frame()
  modules2genes <- list()
  goGroupName <- names(go_object)[1]
  showMethods(Term)
  
  goCount <- length(go_object)
  count <- 1
  for(goGroupName in names(go_object)) {
    if (count %% 1000 == 0) print(paste(count, "of", goCount))
    count <- count + 1
    
    goGroup <- go_object[goGroupName]
    geneIDs <- unique(unlist(goGroup, use.names=F))  #discard evidence codes
    genesymbols <- unique(getSYMBOL(geneIDs, data=goSource))
    
    genesymbols <- intersect(genesymbols, sortedGenes) #get size after intersecting with our full gene set
    if (!(length(genesymbols) >= minGOgroupSize & length(genesymbols) <= maxGOgroupSize)) next();
    
    modules2genes[goGroupName] <- list(genesymbols)
    
    tmodNames <- rbind(tmodNames, data.frame(ID=goGroupName, Title = Term(goGroupName)))
  }
  geneSetsGO <- makeTmod(modules = tmodNames, modules2genes = modules2genes)
}

result <- tbl_df(tmodUtest(c(sortedGenes), mset=geneSetsGO, qval = 1, filter = T))
result %<>% rowwise() %>% mutate(P.Value = if_else(AUC < 0.5, 1 - P.Value, P.Value)) %>% ungroup() %>% mutate(adj.P.Val=p.adjust(P.Value, method="fdr")) #tmod runs one-sided tests
result %<>% rowwise() %>% mutate(aspect = Ontology(ID)) #add the source ontology (could be filterd for just biological process)

#collapse genesets that have the exact same set of genes
#don't collapse as we need to merge it 
#result %<>% rowwise() %>% mutate(genes = paste(sort(unlist(geneSetsGO$MODULES2GENES[ID])), collapse = " "))
#result %<>% ungroup() %>% group_by(genes, N1) %>% arrange(Title) %>% 
#  summarize(MainTitle = dplyr::first(Title),  ID=paste(ID, collapse=","), AUC = dplyr::first(AUC), P.Value= dplyr::first(P.Value), aspect= dplyr::first(aspect), otherNames = if_else(dplyr::n() > 1, paste(Title[2:length(Title)], collapse=", "), ""))
#result %<>% ungroup() %>% mutate(adj.P.Val=p.adjust(P.Value, method="fdr"))
#result %<>% dplyr::select(-genes)
result %<>% arrange(P.Value)
#print top 20 at top and bottom
print(head(filter(result, AUC > 0.5) %>% dplyr::select(-ID), n=20))

write_csv(result, "./results/mouse/Agranular insular area.zScores-GO.csv")
