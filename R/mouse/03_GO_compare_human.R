detach("package:dplyr", unload=TRUE)
library(readr)
library(magrittr)
library(dplyr)

mouse_go <- read_csv("./results/mouse/Agranular insular area.zScores-GO.csv", guess_max = 10000)
human_long_go <- read_tsv("./results/insula.neocortex.FALSE/adult-long insular gyri-GO.SigAndSpec.tsv")
human_short_go <- read_tsv("./results/insula.neocortex.FALSE/adult-short insular gyri-GO.SigAndSpec.tsv")

human_short_go %<>% rowwise() %>% mutate(singleID = strsplit(x=ID, split=",")[[1]][1])
human_long_go %<>% rowwise() %>% mutate(singleID = strsplit(x=ID, split=",")[[1]][1])

mouse_slim <- mouse_go %>% mutate(P.Value = signif(P.Value, digits=3)) %>% select(singleID = ID, mouse.P.Value = P.Value, mouse.geneCount = N1)
mouse_slim
#human_long_go %<>% left_join(mouse_slim) %>% select(MainTitle, P.Value, adj.P.Value, mouse.P.Value, betterAUCCount,geneCount, mouse.geneCount, everything())
#human_short_go %<>% left_join(mouse_slim) %>% select(MainTitle, P.Value, adj.P.Value, mouse.P.Value, betterAUCCount,geneCount, mouse.geneCount, everything())
human_long_go %<>% left_join(mouse_slim)
human_short_go %<>% left_join(mouse_slim)

human_long_go %>% write_tsv("./results/insula.neocortex.FALSE/adult-long insular gyri-GO.SigAndSpec.plusMouse.tsv")
human_short_go %>% write_tsv("./results/insula.neocortex.FALSE/adult-short insular gyri-GO.SigAndSpec.plusMouse.tsv")

#########################
nrow(human_long_go %>% filter(adj.P.Value < 0.05))
nrow(human_long_go_filtered <- human_long_go %>% filter(mouse.P.Value < 0.05, adj.P.Value < 0.05) )
nrow(human_short_go %>% filter(adj.P.Value < 0.05))
nrow(human_short_go_filtered <- human_short_go %>% filter(mouse.P.Value < 0.05, adj.P.Value < 0.05) )

write_tsv(human_long_go, "./results/insula.neocortex.FALSE/adult-long insular gyri-GO-plusMouse.tsv")
write_tsv(human_short_go, "./results/insula.neocortex.FALSE/adult-short insular gyri-GO-plusMouse.tsv")

write_tsv(human_long_go_filtered, "./results/insula.neocortex.FALSE/adult-long insular gyri-GO-plusMouse.filtered.tsv")
write_tsv(human_short_go_filtered, "./results/insula.neocortex.FALSE/adult-short insular gyri-GO-plusMouse.filtered.tsv")
