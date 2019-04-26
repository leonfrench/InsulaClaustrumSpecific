library(cowplot)
library(reshape2)
detach("package:dplyr", unload=TRUE)
library(dplyr)
library(magrittr)
library(readr)

inclusionPattern <- "insula"
#inclusionPattern <- "claustrum"
neoCortexOnly <- FALSE

baseFolder <- paste0(inclusionPattern, ".neocortex.", neoCortexOnly)


allDisResults <- NULL
for (file in list.files(path = paste0("./results/", baseFolder), pattern = "-DisGeNet.tsv")) {
  result <- read_tsv(paste0("./results/",baseFolder,"/", file)) %>% mutate(source = file)
  allDisResults <- bind_rows(result, allDisResults)
}
allDisResults %<>% tidyr::separate(source, sep="[-]", into=c("age", "region")) #this warning is fine

#check for non matches
allDisResults %>% group_by(ID) %>% summarize(n=dplyr::n()) %>% filter(n!=4)

#unique groups:
length(unique(allDisResults$ID))

#how many significant?
allDisResults %>% filter(age == "adult") %>% group_by(region) %>% filter(adj.P.Value < 0.05) %>% summarize(n=dplyr::n())

allDisResults %>% filter(age == "adult") %>% group_by(region) %>% filter(betterAUCCount == 0) %>% summarize(n=dplyr::n())
allDisResults %>% filter(age == "adult") %>% group_by(region) %>% filter(adj.P.Value < 0.05) %>% select(MainTitle)

allDisResults %>% filter(age == "fetal") %>% group_by(region) %>% filter(adj.P.Value < 0.05) 
allDisResults %>% filter(age == "fetal") %>% group_by(region) %>% filter(adj.P.Value < 0.05) %>% select(MainTitle)
allDisResults %>% filter(age == "fetal") %>% group_by(region) %>% filter(betterAUCCount == 0)

alcoholGroups <- allDisResults %>% filter(grepl( "[Aa]lcohol", MainTitle)) %>% arrange(P.Value) %>% group_by(age, region) %>% mutate(adj.P.Value = p.adjust(P.Value, method="fdr"))
alcoholGroups %<>% filter(age=="adult") %>% group_by(age, region) %>% mutate(adj.P.Value = p.adjust(P.Value, method="fdr"))
alcoholGroups %>% ungroup() %>% select(MainTitle) %>% distinct()
alcoholGroups %>% arrange(P.Value) %>% select(MainTitle, region, everything(), -otherNames, -otherNames, -aspect)

allDisResults %>% filter(grepl( "cocaine", MainTitle, ignore.case = T)) %>% filter(age=="adult") %>% arrange(P.Value) %>% group_by(age, region) %>% mutate(adj.P.Value = p.adjust(P.Value, method="fdr"))
allDisResults %>% filter(grepl( "amphet", MainTitle, ignore.case = T)) %>%  arrange(P.Value) %>% group_by(age, region) %>% mutate(adj.P.Value = p.adjust(P.Value, method="fdr"))
allDisResults %>% filter(grepl( "Marijuana", MainTitle, ignore.case = T)) %>% filter(age=="adult") %>% arrange(P.Value) %>% group_by(age, region) %>% mutate(adj.P.Value = p.adjust(P.Value, method="fdr"))

depressionGroups <- allDisResults %>% filter(MainTitle %in% c("Depressive Symptoms", "Depressed mood", "Drug-induced depressive state", "Depression, Bipolar")) %>% arrange(P.Value) %>% group_by(age, region) %>% mutate(adj.P.Val = p.adjust(P.Value, method="fdr"))
depressionGroups %<>% filter(age=="adult") %>% group_by(age, region) %>% mutate(adj.P.Value = p.adjust(P.Value, method="fdr"))
depressionGroups %>% ungroup() %>% select(MainTitle) %>% distinct()
depressionGroups %>% arrange(P.Value) %>% select(MainTitle, region, everything(), -aspect, -otherNames)

#add fetal p-values to top ten lists
fetalDisResults <- allDisResults %>% filter(age=="fetal") %>% select(ID, region, fetalrank = rank, fetalAUC= AUC, fetalP = P.Value)

#map the two regions between adult/fetal
if (inclusionPattern == "insula") {
  fetalDisResults %<>% mutate(mapped_adult_region = if_else(region == "granular insular cortex", "long insular gyri", "short insular gyri"))
} else {
  fetalDisResults %<>% mutate(mapped_adult_region = region)
}

adultDisResults <- inner_join(allDisResults %>% filter(age=="adult"), fetalDisResults, by = c("ID", "region" = "mapped_adult_region"))

#write out the adult region top ten lists
unique(adultDisResults$region)
for (target_region in unique(adultDisResults$region)) {
  top10 <- adultDisResults %>% filter(region == target_region) %>% arrange(P.Value) %>% head(10) %>% select(-age, -region)
  top10 %<>% mutate(AUC = signif(AUC, digits=3), adj.P.Value = signif(adj.P.Value, digits=3)) %>%
    select(Name = MainTitle, Genes = geneCount, AUROC=AUC, `Specificity Rank` = betterAUCCount, `p-valueFDR` = adj.P.Value, fetalAUC, fetalP)
  write_tsv(top10, paste0("./results/", baseFolder, "/adult-top10-", target_region, "-DisGeNet.top10.tsv"))
}

#hits with signifcant adult and fetal enrichment
adultDisResults %>% filter(adj.P.Value < 0.05, fetalP < 0.05) %>% select(MainTitle)
