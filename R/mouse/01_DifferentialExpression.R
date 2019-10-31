library(magrittr)
library(data.tree)
library(readr)
detach("package:dplyr", unload=TRUE)
library(dplyr) 

#TODO - generalize to other regions

df <- read_tsv("./data/raw/mouse_ABA/allMerged.tsv.zip")
exp <- read_tsv("./data/raw/mouse_ABA/experiments.tsv")
structures <-read_tsv("./data/raw/mouse_ABA/structures.tsv")

exp %>% group_by(geneSymbol, plane) %>% summarize(dplyr::n()) %>% group_by(plane) %>% summarize(dplyr::n())

#uncomment this to remove low expressing genes
#nonExpressors <- readLines("./data/other gene lists/mouse/ZieselAndAllenNonExpressors.txt")
#length(intersect(exp$geneSymbol, nonExpressors))
#exp %<>% filter(!geneSymbol %in% nonExpressors)

df %>% select(regionName) %>% distinct() %>% filter(grepl('sula', regionName)) 

structures %>% select(name) %>% distinct() %>% filter(grepl('sula', name)) 

#create tree of structures
structures <- inner_join(structures %>% select(parent_structure_id = ID, parent_name=name) %>% distinct(), structures)
structures %<>% select(from = parent_name, to = name)
FromDataFrameNetwork(structures)
structures %>% filter(from == "root")
structures %<>% mutate(from = if_else(from == "root", "all terms", from))
structure_tree <- FromDataFrameNetwork(structures)

#get height - Agranular insular area is height = 3
structure_tree$`Basic cell groups and regions`$Cerebrum$`Cerebral cortex`$`Cortical plate`$Isocortex$`Agranular insular area`$height

#get regions at the same tree height for comparison
regionsAtLevel3 <- structure_tree$Get("name", filterFun = function(node) { node$height == 3 })
regionsAtLevel3 <- names(regionsAtLevel3)

length(regionsAtLevel3)
regionsAtLevel3 <- intersect(regionsAtLevel3,  df$regionName)
length(regionsAtLevel3)


#reduce to level3 regions
df %<>% filter(regionName %in% regionsAtLevel3) 
scale_this <- function(x){
  (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
}
df %<>% group_by(label, plane) %>% mutate(expression_energy_z = scale_this(expression_energy))

#genes sorted by insula z-score
agranular_specific <- df %>% filter(regionName == "Agranular insular area") %>% arrange(-expression_energy_z)

agranular_specific <- inner_join(agranular_specific, exp) %>% select(geneSymbol, plane, expression_energy_z, everything())
agranular_specific
write_csv(agranular_specific, "./results/mouse/Agranular insular area.zScores.csv")

#look at average of two image sets
two_experiments <- agranular_specific %>% group_by(geneSymbol) %>% summarize(expression_energy_z_mean = mean(expression_energy_z), n=n(), entrezid= first(entrezid), plane="sagittal") %>% filter(n>1) %>% arrange(-expression_energy_z_mean)

write_csv(two_experiments, "./results/mouse/Agranular insular area.zScores.two_experiments.csv")
