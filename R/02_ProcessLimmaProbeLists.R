library(org.Hs.eg.db)
library(annotate)
library(readr)
library(magrittr)
detach("package:dplyr", unload=TRUE)
library(dplyr)

#these don't line up so well
#adult_target <- "long"
#fetal_target <- "granular insular cortex"

#adult_target <- "short"
#fetal_target <- "dysgranular insular cortex"

adult_target <- "claustrum"
fetal_target <- "claustrum"

inclusionPattern <- "claustrum"
neoCortexOnly <- FALSE

baseFolder <- paste0(inclusionPattern, ".neocortex.", neoCortexOnly)

(adultFilename <- list.files(paste0("./results/",baseFolder,"/limma/"), pattern = paste0("^",adult_target,".*allen_HBA.geneSummary.csv"), full.names = T)[1])
(fetalFilename <- list.files(paste0("./results/",baseFolder,"/limma/"), pattern = paste0("^",fetal_target,".*allen_human_fetal_brain.geneSummary.csv"), full.names = T)[1] )

#load fetal and adult gene-wise data
adult <- read_csv(adultFilename, guess_max = 20000)
fetal <- read_csv(fetalFilename, guess_max = 20000)

adultProbes <- read_csv(gsub("geneSummary.","", adultFilename))
fetalProbes <- read_csv(gsub("geneSummary.","", fetalFilename))

#check counts
length(unique(adultProbes$gene_symbol))
length(unique(fetalProbes$gene_symbol))
length(intersect(adultProbes$gene_symbol, fetalProbes$gene_symbol))

length(setdiff(adultProbes$gene_symbol, fetalProbes$gene_symbol))
length(setdiff(fetalProbes$gene_symbol, adultProbes$gene_symbol))
sort(setdiff(fetalProbes$gene_symbol, adultProbes$gene_symbol))
sort(setdiff(adultProbes$gene_symbol, fetalProbes$gene_symbol))

length(unique(adultProbes$probe_name))
length(unique(fetalProbes$probe_name))
length(intersect(adultProbes$probe_name, fetalProbes$probe_name))

adultProbes %<>% mutate(adj.p.value = p.adjust(p.value, method="fdr")) %>% arrange(p.value) %>% 
  mutate(sigUpRegulated = adj.p.value < 0.05 & t > 0) %>% print()
fetalProbes %<>% mutate(adj.p.value = p.adjust(p.value, method="fdr")) %>% arrange(p.value) %>% 
  mutate(sigUpRegulated = adj.p.value < 0.05 & t > 0) %>% print()

supplementTable <- inner_join(
  adultProbes %>% select(probe_name, gene_symbol, p.value, t, adj.p.value), 
  fetalProbes %>% select(probe_name, p.value, t, adj.p.value), suffix = c(".adult", ".fetal")
  , by="probe_name")

paste("Significantly enriched probes for adult", adult_target, ":", nrow(adultProbes %>% filter(sigUpRegulated)))
paste("Significantly enriched probes for fetal", fetal_target, ":", nrow(fetalProbes %>% filter(sigUpRegulated)))

#get the probe level corrected max p-value
threshold <- max((adultProbes %>% filter(sigUpRegulated))$p.value)
adult %<>% arrange(-pValueWithDirection) %>% mutate(is_probe_corrected_significant = p.value <= threshold & direction > 0) %>% print()
threshold <- max((fetalProbes %>% filter(sigUpRegulated))$p.value)
fetal %<>% arrange(-pValueWithDirection) %>% mutate(is_probe_corrected_significant = p.value <= threshold & direction > 0) %>% print()

sigAdultGenes <- unique(adult %>% filter(is_probe_corrected_significant) %>% .$gene_symbol)
sigFetalGenes <- unique(fetal %>% filter(is_probe_corrected_significant) %>% .$gene_symbol)

cat("Significantly enriched genes for adult", adult_target, ":",  length(sigAdultGenes),
    "\nSignificantly enriched genes for fetal", fetal_target, ":",  length(sigFetalGenes),
    "\nIntersecting genes", ":",  length(intersect(sigAdultGenes, sigFetalGenes)))

adult %<>% mutate(rank = rank(-pValueWithDirection))
fetal %<>% mutate(rank = rank(-pValueWithDirection))

geneToName <- adultProbes %>% select(gene_symbol, entrez_id) %>% distinct() %>% filter(!is.na(entrez_id))

#add gene names
geneToName %<>% mutate(name = unlist(lookUp(as.character(entrez_id), "org.Hs.eg", "GENENAME"))) %>% distinct() %>% print()
geneToName %<>% filter(!is.na(name))

adult <- left_join(adult, geneToName) %>% select(gene_symbol, name, entrez_id, everything())
fetal <- left_join(fetal, geneToName) %>% select(gene_symbol, name, entrez_id, everything())
supplementTable <- left_join(supplementTable, geneToName) %>% select(gene_symbol, name, entrez_id, everything())


#top 20
get_top_20 <- function(probe_table, gene_table) {
  #symbol, name, number of significant probes, p-value, fetal rank
  sigProbeCount <- probe_table %>% filter(sigUpRegulated) %>% group_by(gene_symbol) %>% summarize(sigProbes = dplyr::n())
  top20 <- gene_table %>% arrange(rank) %>% head(20)
  top20 %<>% inner_join(sigProbeCount)
  top20 %<>% mutate(p.value = signif(p.value, digits=3)) %>% 
    select(`Gene Symbol` = gene_symbol, `Name` = name, `Significant probes` = sigProbes, `p-value` = p.value) 
  top20
}

adult20 <- get_top_20(adultProbes, adult)
fetal20 <- get_top_20(fetalProbes, fetal)
#show overlapping
#compare adult and fetal top 20
inner_join(adult20, fetal20, by = "Gene Symbol")

write_csv(adult, paste0(adultFilename, ".addedStats.csv") )
write_csv(fetal, paste0(fetalFilename, ".addedStats.csv") )
write_csv(supplementTable, paste0(adultFilename, ".supplementTable.csv") )
write_csv(adult20, paste0(adultFilename, ".top20.csv") )
write_csv(fetal20, paste0(fetalFilename, ".top20.csv") )


