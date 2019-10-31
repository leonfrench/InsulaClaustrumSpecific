########################
#Venn Diagram
#install.packages('venneuler')
library(here)
library(venneuler)
library(ggrepel)

#input is a tmod result table and the number of top groups - don't send the whole result
getEulerDiagram <- function(goResult, geneSets) {
  combined <- NULL
  for(group in goResult$ID) {
    groupName <- as.character(filter(geneSets$MODULES, ID == group)$Title[1])
    combined <- rbind(combined, data.frame(elements= unlist(geneSets$MODULES2GENES[group]), sets=groupName))
  }
  
  v <- venneuler(combined)
  plot(v)
  
  (vForGGplot <- tbl_df(data.frame(diameter=v$diameters, v$centers, color=v$colors, goName=v$labels, stringsAsFactors = F)))
  
  #credited to user5061 and baptiste via stackoverflow.com
  circularise <- function(d, n=360){
    angle <- seq(-pi, pi, length = n)
    make_circle <- function(x,y,r,goName){
      data.frame(x=x+r*cos(angle), y=y+r*sin(angle), goName)
    }
    lmat <- mapply(make_circle, goName = d[,"goName"], 
                   x = d[,"x"], y=d[,"y"], r=d[,"diameter"]/2, SIMPLIFY = FALSE)
    do.call(rbind, lmat)
  }
  
  circles <- circularise(vForGGplot)
  vForGGplot <- as.data.frame(vForGGplot)
  diagram <- ggplot(vForGGplot) + geom_blank(aes(x, y)) + 
    geom_polygon(aes(x,y, group=goName, fill=goName), data=circles, alpha=0.4) +
    coord_fixed() + theme_void(base_size = 16) + theme(legend.position="none") +
    geom_label_repel(aes(x, y, fill=goName, label = goName), color="black", segment.alpha=0, box.padding = unit(0.2, "lines")) 
  diagram
}

if (exists("geneSetsGO") && length(geneSetsGO$MODULES2GENES) > 1000 ) { #assume it's already loaded - needs a fix to see if the variable is declared
} else {
  geneSetsGO <- loadGOSets(allGenes)
}
geneSetsToUse <- geneSetsGO


results <- read_tsv(here("results", "insula.neocortex.FALSE", "adult-short insular gyri-GO.SigAndSpec.plusMouse.tsv"))

getEulerDiagram(results, geneSetsToUse)

