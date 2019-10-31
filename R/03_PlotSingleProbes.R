library(cowplot)
library(reshape2)
library(dplyr)
library(magrittr)
library(cowplot)
library(ggplot2)
library(readr)

#inclusionPattern <- "habenula"
#inclusionPattern <- "claustrum"
inclusionPattern <- "insula"
#inclusionPattern <- "superior colliculus"
neoCortexOnly <- FALSE

baseFolder <- paste0(inclusionPattern, ".neocortex.", neoCortexOnly)

#load fetal and human - files written by limma.allRegions.R
expressionMatrixAdult <- read_tsv(paste0("./data/processed/",baseFolder,"/allen_HBA_sampleMatrix_qcNames.tsv"))
expressionMatrixFetal <- read_tsv(paste0("./data/processed/",baseFolder,"/allen_human_fetal_brain_sampleMatrix_qcNames.tsv"))

sampleAnnotAdult <- read_tsv(paste0("./data/processed/",baseFolder,"/allen_HBA_sampleMatrix_colAnnotations_qcNames.tsv"))
sampleAnnotFetal <- read_tsv(paste0("./data/processed/",baseFolder,"/allen_human_fetal_brain_sampleMatrix_colAnnotations_qcNames.tsv"))

#grab it from the limma results file
#just grab one probe file
adultProbeInfoFilename <- list.files(paste0("./results/", baseFolder, "/limma/"), pattern = ".allen_HBA.csv", full.names = T)[1]
probesInfo <- read_csv(adultProbeInfoFilename) %>% select(-t, -p.value, -is_qc_pass)

#probesInfo <- read_csv(paste0("./results/limma/lateral habenular nucleus.allen_HBA.csv")) %>% select(-t, -p.value, -is_qc_pass)
#genes_of_interest <- c("GPR151")
#genes_of_interest <- c("UGT2B15")
#genes_of_interest <- c("TRPM8")


#genes_of_interest <- c("CYP3A7")

#genes_of_interest <- c("CYP3A4", "CYP3A43", "CYP3A5","CYP3A7")
genes_of_interest <- c("CYP3A4", "CYP3A5","CYP3A7")
genes_of_interest <- c("ANXA1")
#genes_of_interest <- c("HTR2C")
genes_of_interest <- c("MIR133A1")
genes_of_interest <- c("NR4A2")
genes_of_interest <- c("UCMA")





#filter for probes for gene of interest
probesInfo %<>% filter(gene_symbol %in% genes_of_interest)

fetalFocus <- expressionMatrixFetal %>% filter(probe_name %in% probesInfo$probe_name)
adultFocus <- expressionMatrixAdult %>% filter(probe_name %in% probesInfo$probe_name)

fetalFocus <- as_tibble(melt(fetalFocus, value.name = "Expression", variable.name = "uniqueID"))
adultFocus <- as_tibble(melt(adultFocus, value.name = "Expression", variable.name = "uniqueID"))

if (length(genes_of_interest) > 1) {
  probesInfo %<>% mutate(joinedName = paste0(probe_name, " (", gene_symbol, ")"))
  adultFocus <- inner_join(adultFocus, probesInfo) %>% mutate(probe_name = joinedName)
  fetalFocus <- inner_join(fetalFocus, probesInfo) %>% mutate(probe_name = joinedName)
  fetalFocus$probe_name <- factor(fetalFocus$probe_name, levels = unique((fetalFocus %>% arrange(gene_symbol))$probe_name))
  adultFocus$probe_name <- factor(adultFocus$probe_name, levels = unique((adultFocus %>% arrange(gene_symbol))$probe_name))
}


fetalFocus %<>% inner_join(sampleAnnotFetal %>% select(uniqueID, Region = structure_name_left_right_stripped, Donor=donorID))
adultFocus %<>% inner_join(sampleAnnotAdult %>% select(uniqueID, Region = structure_name_left_right_stripped, Donor=donorID))

adultFocus %<>% mutate(Donor = gsub("normalized_microarray_donor", "", Donor))
fetalFocus %<>% mutate(Donor = gsub("lmd_matrix_", "", Donor))

adultFocus %<>% mutate(Region = if_else(grepl(inclusionPattern, Region), Region, "remaining structures"))
adultFocus %>% group_by(Region) %>% summarize(n=dplyr::n())
fetalFocus %<>% mutate(Region = if_else(grepl(inclusionPattern, Region), Region, "remaining structures"))
fetalFocus %>% group_by(Region) %>% summarize(n=dplyr::n())

#are they seperable?
adultFocus %>% group_by(probe_name, Region) %>% summarize(maxE = max(Expression))

adultFocus %<>% mutate(Donor = paste0("ID", Donor))
fetalFocus %<>% mutate(Donor = paste0("ID", Donor))

ggplot(adultFocus, aes(y=Expression, x=probe_name, color=Region)) + geom_point() + facet_wrap( ~Donor)
ggplot(fetalFocus, aes(y=Expression, x=probe_name, color=Region)) + geom_point() + facet_wrap( ~Donor)

#color bind freindly palette
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
if (inclusionPattern == "claustrum") { #colors for two regions
  cbbPalette <- c("#000000", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
}
unique(adultFocus$Region)

if (inclusionPattern == "insula") {
  adultFocus %<>% mutate(Region = factor(Region, levels = c("short insular gyri", "long insular gyri", "remaining structures")))
  fetalFocus %<>% mutate(Region = factor(Region, levels = c("dysgranular insular cortex", "granular insular cortex", "remaining structures")))  
}

(adultPlot <- ggplot(adultFocus, aes(y=Expression, x=Donor, fill = Region, color = Region)) + 
  geom_violin(data = filter(adultFocus, Region == "remaining structures")) + theme_bw() +
  geom_dotplot(data = filter(adultFocus, Region != "remaining structures"), binaxis = "y", stackdir = "center", dotsize=.62) +
  #geom_point(data = filter(adultFocus, Region != "remaining structures")) + 
  facet_wrap( ~probe_name, nrow=1) + #, scales="free_y", nrow=1) +
  scale_fill_manual(values=cbbPalette, limits= levels(adultFocus$Region)) + scale_colour_manual(values=cbbPalette, limits= levels(adultFocus$Region)) + xlab("")
)
if (inclusionPattern == "insula") {
  (adultPlot <- adultPlot + theme(legend.position="bottom"))
} else {
  (adultPlot <- adultPlot + guides(fill=FALSE, color=FALSE))
}


(fetalPlot <- ggplot(fetalFocus, aes(y=Expression, x=Donor, fill = Region, color=Region)) + 
  geom_violin(data = filter(fetalFocus, Region == "remaining structures")) + theme_bw() +
  #geom_point(data = filter(fetalFocus, Region != "remaining structures")) + 
  geom_dotplot(data = filter(fetalFocus, Region != "remaining structures"), binaxis = "y", stackdir = "center", dotsize=.5) +
  facet_wrap( ~probe_name, nrow=1) + #, scales="free_y", nrow=1) +
  scale_fill_manual(values=cbbPalette) + scale_colour_manual(values=cbbPalette) + theme(legend.position="bottom")
)

fetalPlot <- fetalPlot + theme(plot.margin = margin(.85, 0.1, 0, 0.1, "cm"))
adultPlot <- adultPlot + theme(plot.margin = margin(.85, 0.1, 0, 0.1, "cm"))

#plot_grid(adultPlot, fetalPlot, rel_heights = c(8.5,10), labels = c("A (adult brain)", "B (fetal brain)"), nrow=2, label_x=0.5, hjust=0.5)
plot_grid(adultPlot, fetalPlot, rel_heights = c(10,10), labels = c("A) Adult brain", "B) Fetal brain"), nrow=2, hjust=0)
#save as for GPR151 7x6.6 PDF
#for CYP 11x6.6
#for UCMA use 9x9
#for NR4A2 use 11x9
