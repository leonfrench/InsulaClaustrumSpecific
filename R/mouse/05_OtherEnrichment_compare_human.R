detach("package:dplyr", unload=TRUE)
library(readr)
library(magrittr)
library(dplyr)

#geneSetName <- "DisGeNet"
#geneSetName <- "Custom" 
geneSetName <- "Extra" 

mouse_enrichments <- read_csv(paste0("./results/mouse/Agranular insular area.zScores-", geneSetName, ".csv"), guess_max = 10000)
human_long_enrichments <- read_tsv(paste0("./results/insula.neocortex.FALSE/adult-long insular gyri-", geneSetName, ".tsv"))
human_short_enrichments <- read_tsv(paste0("./results/insula.neocortex.FALSE/adult-short insular gyri-", geneSetName, ".tsv"))

mouse_slim <- mouse_enrichments %>% select(ID, mouse.P.Value = P.Value, mouse.geneCount = N1, mouse.AUC = AUC)

human_long_enrichments %<>% left_join(mouse_slim) %>% select(MainTitle, P.Value, adj.P.Value, mouse.P.Value, betterAUCCount,geneCount, mouse.geneCount, AUC, mouse.AUC, everything())
human_short_enrichments %<>% left_join(mouse_slim) %>% select(MainTitle, P.Value, adj.P.Value, mouse.P.Value, betterAUCCount,geneCount, mouse.geneCount,AUC, mouse.AUC, everything())

(human_long_enrichments_filtered <- human_long_enrichments %>% filter(mouse.P.Value < 0.1, adj.P.Value < 0.05) )
(human_short_enrichments_filtered <- human_short_enrichments %>% filter(mouse.P.Value < 0.1, adj.P.Value < 0.05) )

write_tsv(human_long_enrichments, paste0("./results/insula.neocortex.FALSE/adult-long insular gyri-", geneSetName, "-plusMouse.tsv"))
write_tsv(human_short_enrichments, paste0("./results/insula.neocortex.FALSE/adult-short insular gyri-", geneSetName, "-plusMouse.tsv"))
