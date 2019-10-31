library(cowplot)
library(reshape2)
detach("package:dplyr", unload=TRUE)
library(dplyr)
library(magrittr)
library(readr)

inclusionPattern <- "insula"
#inclusionPattern <- "claustrum"
#inclusionPattern <- "superior colliculus"

neoCortexOnly <- FALSE

baseFolder <- paste0(inclusionPattern, ".neocortex.", neoCortexOnly)

allGOResults <- NULL
for (file in list.files(path = paste0("./results/",baseFolder), pattern = "-GO.tsv")) {
  result <- read_tsv(paste0(paste0("./results/",baseFolder),"/", file)) %>% mutate(source = file)
  allGOResults <- bind_rows(result, allGOResults)
}
allGOResults %<>% tidyr::separate(source, sep="[-]", into=c("age", "region")) #this warning is fine

#check for non matches
allGOResults %>% group_by(ID) %>% summarize(n=dplyr::n()) %>% filter(n!=4)

#unique groups:
length(unique(allGOResults$ID))

#how many significant?
allGOResults %>% filter(age == "adult") %>% group_by(region) %>% filter(adj.P.Value < 0.05) %>% summarize(n=dplyr::n())
allGOResults %>% filter(age == "fetal") %>% group_by(region) %>% filter(adj.P.Value < 0.05) %>% summarize(n=dplyr::n())

#how many significant and best AUC
allGOResults %>% filter(age == "adult") %>% group_by(region) %>% filter(betterAUCCount == 0) %>% summarize(n=dplyr::n())

nicotineGroups <- allGOResults %>% filter(grepl( "nicotine", MainTitle)) %>% arrange(P.Value) %>% group_by(age, region) %>% mutate(adj.P.Value = p.adjust(P.Value, method="fdr"))
nicotineGroups %<>% filter(age=="adult")
nicotineGroups %>% ungroup() %>% select(MainTitle) %>% distinct()
nicotineGroups %>% filter(adj.P.Value < 0.05)

alcoholGroups <- allGOResults %>% filter(grepl( "alcohol", MainTitle)) %>% arrange(P.Value) %>% group_by(age, region) %>% mutate(adj.P.Value = p.adjust(P.Value, method="fdr"))
alcoholGroups %<>% filter(age=="adult")
alcoholGroups %>% ungroup() %>% select(MainTitle) %>% distinct()
alcoholGroups %>% filter(adj.P.Value < 0.05) 
alcoholGroups %>% select(MainTitle, ID)
alcoholGroups %>% filter(MainTitle == "response to alcohol")
alcoholGroups %>% filter(adj.P.Value < 0.05) %>% select(MainTitle, ID) #look at which genes
alcoholGroups %>% filter(MainTitle == "primary alcohol catabolic process")
alcoholGroups %>% filter(MainTitle == "secondary alcohol biosynthetic process")


allGOResults %>% filter(grepl( "cocaine", MainTitle)) %>% arrange(P.Value) %>% group_by(age, region) %>% mutate(adj.P.Value = p.adjust(P.Value, method="fdr"))
allGOResults %>% filter(grepl( "cocaine", MainTitle)) %>% arrange(P.Value) %>% group_by(age, region) %>% mutate(adj.P.Value = p.adjust(P.Value, method="fdr")) %>% select(MainTitle, P.Value, adj.P.Value)

allGOResults %>% filter(grepl( "morphine", otherNames) | grepl( "morphine", MainTitle)) %>% arrange(P.Value) %>% group_by(age, region) %>% mutate(adj.P.Value = p.adjust(P.Value, method="fdr"))
allGOResults %>% filter(grepl( "morphine", otherNames) | grepl( "morphine", MainTitle)) %>% arrange(P.Value) %>% group_by(age, region) %>% mutate(adj.P.Value = p.adjust(P.Value, method="fdr")) %>% select(otherNames)

allGOResults %>% filter(grepl( "amphet", MainTitle)) %>% arrange(P.Value) %>% group_by(age, region) %>% mutate(adj.P.Value = p.adjust(P.Value, method="fdr"))

allGOResults %>% filter(MainTitle == "response to xenobiotic stimulus") %>% arrange(P.Value) %>% filter(age=="adult")
allGOResults %>% filter(MainTitle == "xenobiotic metabolic process") %>% arrange(P.Value) %>% filter(age=="adult")

unique(allGOResults$region)
for (target_region in unique(allGOResults %>% filter(age=="adult") %>% .$region)) {
  top15 <- allGOResults %>% filter(age=="adult", region == target_region) %>% arrange(P.Value) %>% head(15) %>% select(-age, -region)
  top15 %>% mutate(AUC = signif(AUC, digits=3), adj.P.Value = signif(adj.P.Value, digits=3)) %>%
    select(Name = MainTitle, Genes = geneCount, AUROC=AUC, `Specificity Rank` = betterAUCCount, `p-valueFDR` = adj.P.Value) %>%
    write_tsv(paste0("./results/", baseFolder,"/adult-", target_region,"-GO.top15.tsv"))
}
  
#filter for specific and significant in fetal
sigAndSpec <- allGOResults %>% group_by(ID, region) %>% filter(age=="adult", adj.P.Value < 0.05, betterAUCCount <= 1) %>% arrange(rank)
sigAndSpec %>% group_by(region) %>% summarize(n=dplyr::n()) 
sigAndSpec %<>% select(MainTitle, region, everything(), -age) %>% arrange(P.Value)
########### 
write_tsv(sigAndSpec, paste0("./results/",baseFolder, "/adult-", target_region, "-GO.SigAndSpec.tsv"))

fetalEnriched <- allGOResults %>% filter(age=="fetal", P.Value < 0.05) %>% select(ID, region, fetalrank = rank, fetalAUC= AUC, fetalP = P.Value)


#map the two regions between adult/fetal
if (inclusionPattern == "insula") {
  fetalEnriched %<>% mutate(mapped_adult_region = if_else(region == "granular insular cortex", "long insular gyri", "short insular gyri"))
} else {
  fetalEnriched %<>% mutate(mapped_adult_region = region)
}

#print out spec and sig (twice)
inner_join(sigAndSpec, fetalEnriched, by = c("ID", "region" = "mapped_adult_region")) %>% group_by(region) %>% summarize(n=dplyr::n()) 
inner_join(sigAndSpec, fetalEnriched, by = c("ID", "region" = "mapped_adult_region")) %>% group_by(region) %>% select(MainTitle, geneCount, region, AUC,  adj.P.Value, fetalAUC, fetalP)

unique(sigAndSpec$region)
for (target_region in unique(sigAndSpec$region)) {
  bothSigAndSpec <- inner_join(sigAndSpec, fetalEnriched, by = c("ID", "region" = "mapped_adult_region")) %>% filter(region==target_region) 
  bothSigAndSpec %>% ungroup() %>% mutate(AUC = signif(AUC, digits=3), adj.P.Value = signif(adj.P.Value, digits=3), fetalAUC = signif(fetalAUC, digits=3), fetalP = signif(fetalP, digits=3)) %>%
    select(Name = MainTitle, ID, Genes = geneCount, AUROC=AUC, `Specificity Rank` = betterAUCCount, `p-valueFDR` = adj.P.Value, `Fetal AUROC` = fetalAUC,  `Fetal p-value` = fetalP) %>%
    write_tsv(paste0("./results/", baseFolder,"/adult-", target_region,"-GO.fetalvalidated.SigAndSpec.tsv"))
  print(target_region)
  print(binom.test(nrow(bothSigAndSpec), nrow(sigAndSpec %>% filter(region==target_region)), p = 0.05, alternative = c("greater"), conf.level = 0.95))
}


#filter on better hits only - only works for same region names in adult/fetal
allGOResults %>% group_by(MainTitle, region) %>% filter(betterAUCCount < 3) %>% summarize(n = dplyr::n()) %>% filter(n>1)

allGOResults %>% group_by(MainTitle, region) %>% filter(betterAUCCount < 4) %>% summarize(n = dplyr::n()) %>% filter(n>1)



