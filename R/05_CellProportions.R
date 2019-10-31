library(lazyeval)
library(tidyr)
library(readr)
library(dplyr)
library(magrittr)

inclusionPattern <- "insula"
fullinclusionPattern <- "insula|claustrum"

neoCortexOnly <- FALSE

baseFolder <- paste0(inclusionPattern, ".neocortex.", neoCortexOnly)

cellProp <- read_csv(paste0("./results/", baseFolder, "/Adult.Regional.cellEstimates.csv"))
cellProp %<>% mutate_if(is.numeric, funs(. * -1)) %>% mutate_if(is.numeric, rank)

cellProp %<>% filter(grepl(fullinclusionPattern, structure_name_left_right_stripped))
cellProp <- melt(cellProp)
cellProp <- as_tibble(cellProp %>% arrange(value))
cellProp %<>% spread(structure_name_left_right_stripped, value)
cellProp %<>% mutate(min = apply(cellProp, 1, min)) %>% arrange(min) %>% select(-min)

nm1 <- colnames(cellProp)[2:ncol(cellProp)]
cellProp %>% mutate_(min = interp(~pmin(v1), v1= as.name(nm1))) %>% arrange(min) %>% select(-min)


#arrange by the lower of the two regions - needs to be generalized
cellProp %<>% mutate(min = pmin(`lateral habenular nucleus`, `medial habenular nucleus`)) %>% arrange(min) %>% select(-min)

write_csv(cellProp, paste0("./results/", baseFolder, "/Adult.relative.cellEstimates.csv"))

cellProp %>% filter(grepl("NeuroExpresso", variable)) %>% arrange(`long insular gyri`)
cellProp %>% filter(grepl("Darm", variable)) %>% mutate(variable = gsub("Darmanis.", "", variable)) %>% write_csv(paste0("./results/", baseFolder, "/Adult.Darmanis.relative.cellEstimates.csv"))
