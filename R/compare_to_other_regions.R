library(doMC)
registerDoMC(cores=6)
library(here)
library(readr)
library(magrittr)
library(dplyr)
source(here("R", "AUCFunction.R"))

adult_target_region <- "claustrum"
inclusionPattern <- adult_target_region
neoCortexOnly <- FALSE

baseFolder <- paste0(inclusionPattern, ".neocortex.", neoCortexOnly)

unique(fetal %>% filter(is_probe_corrected_significant) %>% .$gene_symbol)
#claustrum
for_top_genes <- read_csv(here("results", baseFolder, "limma", paste0(inclusionPattern, ".allen_HBA.geneSummary.csv.addedStats.csv")), guess_max=10000)
sigAdultGenes <- unique(for_top_genes %>% filter(is_probe_corrected_significant) %>% .$gene_symbol)

adult_matrix <- read_tsv(paste0("./data/processed/", baseFolder, "/allen_HBA_brainarea_vs_genes_exp_qcNames.tsv"), guess_max=10000)
fetal_matrix <- read_tsv(paste0("./data/processed/", baseFolder, "/allen_human_fetal_brain_brainarea_vs_genes_exp_qcNames.tsv"), guess_max = 10000) 

#adult_matrix %<>% head(5000)
#fetal_matrix %<>% head(5000)

#add adult_target_region to fetal matrix
fetal_matrix <- inner_join(adult_matrix %>% select(gene_symbol, adult_target_region), fetal_matrix, by="gene_symbol", suffix=c('.adult',''))

#adult_cor_matrix <- cor(adult_matrix %>% select(-gene_symbol),m='s')
#adult_cor_matrix %<>% as_tibble(rownames = "region")
fetal_cor_matrix <- cor(fetal_matrix %>% select(-gene_symbol),m='p')
fetal_cor_matrix %<>% as_tibble(rownames = "region")

target_col <- paste0(adult_target_region, ".adult")
fetal_cor_matrix %>% select(region, target_col) %>% arrange_(target_col) %>% tail(20)
fetal_cor_matrix %>% select(region, target_col) %>% arrange_(target_col) %>% filter(grepl("insula", region))

#run AUC on the sigAdultGenes versus the fetal data

#adult
forIndices <- adult_matrix %>% select(gene_symbol)
forIndices %<>% mutate(isTargetGene = gene_symbol %in% head(sigAdultGenes, n=20))
#forIndices %<>% mutate(isTargetGene = gene_symbol %in% sigAdultGenes)
targetIndices <- forIndices$isTargetGene
wilcoxTests <- foreach(oneCol=iter(as.data.frame(adult_matrix %>% select(-gene_symbol)), by='col'), .combine=rbind) %dopar% {
  data.frame(auc = auroc_analytic(oneCol, as.numeric(targetIndices)), 
             pValue=wilcox.test(oneCol[targetIndices], oneCol[!targetIndices], conf.int = F)$p.value)
}
wilcoxTests$Region <- colnames(adult_matrix %>% select(-gene_symbol))
wilcoxTests %<>% as_tibble()
#wilcoxTests %<>% filter(Region != adult_target_region)
wilcoxTests %<>% mutate(Rank = rank(-auc))
wilcoxTests %<>% arrange(-auc) %>% print()
wilcoxTests_adult <- wilcoxTests

#fetal
forIndices <- fetal_matrix %>% select(gene_symbol)
forIndices %<>% mutate(isTargetGene = gene_symbol %in% head(sigAdultGenes, n=20))
#forIndices %<>% mutate(isTargetGene = gene_symbol %in% c("MIR133A1"))
#forIndices %<>% mutate(isTargetGene = gene_symbol %in% c("NR4A2"))
targetIndices <- forIndices$isTargetGene

wilcoxTests <- foreach(oneCol=iter(as.data.frame(fetal_matrix %>% select(-gene_symbol)), by='col'), .combine=rbind) %dopar% {
  data.frame(auc = auroc_analytic(oneCol, as.numeric(targetIndices)), 
             pValue=wilcox.test(oneCol[targetIndices], oneCol[!targetIndices], conf.int = F)$p.value)
}
wilcoxTests$Region <- colnames(fetal_matrix %>% select(-gene_symbol))
wilcoxTests %<>% as_tibble()
#wilcoxTests %<>% filter(Region != adult_target_region)
wilcoxTests %<>% filter(Region != paste0(adult_target_region, ".adult"))
wilcoxTests %<>% mutate(Rank = rank(-auc))
#handle acronyms
wilcoxTests %<>% mutate(Region = gsub("IZ", "Intermediate zone", Region))
wilcoxTests %<>% mutate(Region = gsub("MZ", "Marginal zone", Region))
wilcoxTests %<>% arrange(-auc) %>% print()
wilcoxTests_fetal <- wilcoxTests

#wilcoxTests_fetal %>% filter(grepl("insula", Region)) %>% select(Region) %>% distinct() %>% arrange(Region)
#wilcoxTests_fetal %>% filter(grepl("palli", Region))
#wilcoxTests_fetal %>% filter(grepl("utamen", Region))

combined_table <- inner_join(wilcoxTests_adult %>% select(Rank, Region, AUROC=auc), wilcoxTests_fetal %>% select(Rank, Region, AUROC=auc), by="Rank", suffix = c(" Adult"," Fetal"))
combined_table %>% write_csv(here("results", baseFolder, "Nearest_region_table.csv"))

combined_table %>% head(n=10) %>% mutate_if(is.numeric, round, 3) %>% write_csv(here("results", baseFolder, "Nearest_region_table.top10.csv"))
