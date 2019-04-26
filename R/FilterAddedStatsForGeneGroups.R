sourceExpression <- "adult"
#sourceExpression <- "fetal"
#inclusionPattern <- "insula"
inclusionPattern <- "claustrum"
source("./R/GeneSetBuilders.R")

neoCortexOnly <- FALSE

baseFolder <- paste0(inclusionPattern, ".neocortex.", neoCortexOnly)

#geneSetName <- "CYP" 
#geneSetName <- "GO" 
#geneSetName <- "PhenoCarta" 
geneSetName <- "DisGeNet"
#geneSetName <- "Custom" 
#geneSetName <- "Extra" 

targetGroup <- "positive regulation of hemopoiesis"
targetGroup <- "Drug-induced depressive state"
targetGroup <- "glycosaminoglycan binding"
targetGroup <- "regulation of odontogenesis of dentin-containing tooth"
targetGroup <- "icosanoid secretion"
targetGroup <- "leukocyte chemotaxis"
targetGroup <- "acetylcholine binding"
targetGroup <- "immunoglobulin binding"
targetGroup <- "neurotransmitter transporter activity"
targetGroup <- "positive regulation of hemopoiesis"
targetGroup <- "Uremia"
targetGroup <-"positive regulation of alcohol biosynthetic process"
targetGroup <-"FTO Obesity cluster"
targetGroup <- "associative learning"
targetGroup <- "regulation of muscle contraction"
targetGroup <- "oxygen transporter activity"
targetGroup <- "dendritic spine membrane"
targetGroup <- "Substance Withdrawal Syndrome"
targetGroup <- "cardiac muscle contraction"
targetGroup <- "oxygen transport"
targetGroup <- "regulation of macroautophagy"
targetGroup <- "Cocaine-Related Disorders"
targetGroup <- "Li et al 23andme NonTRD vrs TRD"
targetGroup <- "Claustrum markers"
targetGroup <- "oxygen transport"
targetGroup <- "dopamine receptor signaling pathway"
targetGroup <- "nucleotide-sugar metabolic process"
targetGroup <-"Severe mental retardation (I.Q. 20-34)"





if (sourceExpression == "adult") {
  pattern <- ".allen_HBA.geneSummary.csv.addedStats.csv"
  #  medial <- read_csv("./results/limma/medial habenular nucleus.allen_HBA.geneSummary.csv.addedStats.csv", guess_max = 130000)
#  lateral <- read_csv("./results/limma/lateral habenular nucleus.allen_HBA.geneSummary.csv.addedStats.csv", guess_max = 130000)
} else if (sourceExpression == "fetal") {
  pattern <- ".allen_human_fetal_brain.geneSummary.csv.addedStats.csv"
#  medial <- read_csv("./results/limma/medial habenular nucleus.allen_human_fetal_brain.geneSummary.csv.addedStats.csv", guess_max = 130000)
#  lateral <- read_csv("./results/limma/lateral habenular nucleus.allen_human_fetal_brain.geneSummary.csv.addedStats.csv", guess_max = 130000)
}
targetRegions <- gsub(pattern,"",list.files(paste0("./results/", baseFolder, "/limma/"), pattern = pattern))

allGenes <- read_tsv(paste0("./data/processed/", baseFolder, "/allen_HBA_brainarea_vs_genes_exp_qcNames.tsv"), guess_max=10000) %>% .$gene_symbol

if (geneSetName == "PhenoCarta") {
  geneSetsToUse <- loadPhenocarta("human", allGenes)
} else if (geneSetName == "CYP") {
  geneSetsToUse <- loadCypGenes(allGenes)
} else if (geneSetName == "Custom") {
  geneSetsToUse <- loadFileSets(geneSetName)
} else if (geneSetName == "Extra") {
  geneSetsToUse <- loadFileSets(geneSetName)
} else if (geneSetName == "DisGeNet") {
  geneSetsToUse <- loadDisGeNetSets(allGenes)
} else if (geneSetName == "GO") {
  if (exists("geneSetsGO") && length(geneSetsGO$MODULES2GENES) > 1000 ) { #assume it's already loaded - needs a fix to see if the variable is declared
  } else {
    geneSetsGO <- loadGOSets(allGenes)
  }
  geneSetsToUse <- geneSetsGO
}

dir.create(paste0("./results/", baseFolder, "/limma/filtered/"))
dir.create(paste0("./results/", baseFolder, "/limma/filtered/", geneSetName))
baseFilename <- paste0("./results/", baseFolder, "/limma/filtered/", geneSetName, "/")

targetGroupID <- dplyr::filter(tbl_df(geneSetsToUse$MODULES), Title == targetGroup)$ID
for (divisionFileName in targetRegions) {
  toFilter <- read_csv(paste0("./results/", baseFolder, "/limma/", divisionFileName, pattern), guess_max = 130000)
  filter(toFilter, gene_symbol %in% unlist(geneSetsToUse$MODULES2GENES[targetGroupID])) %>% 
    dplyr::arrange(desc(pValueWithDirection)) %>% 
    write_tsv(paste0(baseFilename, divisionFileName,".", targetGroup, ".",sourceExpression,".addedStats.tsv"))
}

