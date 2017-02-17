######################################
#
#Recursive Hypergeometric Enrichment
#Crosstalk Network 
#
######################################

#RHEX-NW

#Inputs: GeneBunch, list of pathways / genesets , gene universe, enrichment p-value
#
#1)For every pathway, hypergeometric (GeneBunch, Pathway)
##List of EnrichedPathways
#
#2) find intersections between EnrichedPathways
##List of IntersectionSets
#
#3)For every IntersectionSet, hypergeometric (GeneBunch, IntersectionSet)
##List of EnrichedIntersectionSets
#
#4) Network construction: Populate with EnrichedPathways
#
#5) If EnrichedIntersectionSets[i,j] != cutoff, link EnrichedPathways[i,j]
##Return Network

#Output: a network, composed of genesets enriched in a GeneBunch,
#linked if they intersect in a set of genes that is also enriched 
#enrichment evaluated through hypergeometric test (Genebunch in Geneset)

########################################

########################################
#Libraries
########################################
library(HTSanalyzeR)
########################################
#Inputs
########################################
file_list<-list.files(pattern = "KEGG")
Pathways<-list()
for(i in file_list){
  x<-as.matrix(read.table(i, header=TRUE))
  assign(paste(colnames(x)), i) ##nombre de columna = nombre variable
  y<-(read.table(i, header=TRUE, stringsAsFactors = FALSE))
  Pathways<-c(Pathways, y)
  rm(x)
  rm(i)
}
rm(list = grep(pattern = "KEGG", x = ls(), value = TRUE))

#list of pathways 
short_pathways = Pathways[1:50]

# universe 
FilePath = "placeholder"
nombres_genes=read.table(file = FilePath, stringsAsFactors = FALSE)
nombres_genes = nombres_genes$V1
nombres_genes = nombres_genes[-grep(pattern = "_at", x = nombres_genes)]

#GeneBunch
GeneBunch = sample(nombres_genes, size = 100, replace = FALSE)
GeneBunch = unique(c(short_pathways$KEGG_ABC_TRANSPORTERS[1:20], short_pathways$KEGG_ADHERENS_JUNCTION[1:20], short_pathways$KEGG_ALLOGRAFT_REJECTION[1:20]))

#1)For every pathway, hypergeometric (GeneBunch, Pathway)

adjpvalue = 0.01

enrichment = multiHyperGeoTest(collectionOfGeneSets = short_pathways, 
                  universe = nombres_genes, 
                  hits = GeneBunch, 
                  minGeneSetSize = 1, 
                  pAdjustMethod = "BH")
enrichment = as.data.frame(enrichment)

##List of EnrichedPathways
EnrichedPathways = rownames(enrichment)[which(enrichment$Adjusted.Pvalue<adjpvalue)]

#2) find intersections between EnrichedPathways

EnrichedPathways_List = short_pathways[EnrichedPathways]
ListIntersectionSets = list()
for(i in seq_along(EnrichedPathways_List)){
  ListIntersectionSets[[i]] <-lapply(X = tail(EnrichedPathways_List, 
                                            n = length(EnrichedPathways_List) - i), 
                                   FUN = function(x) intersect(x, EnrichedPathways_List[[i]]))
  names(ListIntersectionSets)[i]<-names(EnrichedPathways_List)[i]
  }

##List of IntersectionSets

##List of NONEMPTY IntersectionSets
ListIntersectionSets_NONEMPTY = Filter(f = length, unlist(ListIntersectionSets, recursive = FALSE))

#3)For every IntersectionSet, hypergeometric (GeneBunch, IntersectionSet)
enrichment_intersections = multiHyperGeoTest(collectionOfGeneSets = ListIntersectionSets_NONEMPTY, 
                               universe = nombres_genes, 
                               hits = GeneBunch, 
                               minGeneSetSize = 1, 
                               pAdjustMethod = "BH")
enrichment_intersections = as.data.frame(enrichment_intersections)

EnrichedIntersectionSets = rownames(enrichment_intersections)[which(enrichment_intersections$Adjusted.Pvalue<adjpvalue)]
EnrichedIntersectionSets=strsplit(x = EnrichedIntersectionSets, split = "\\.")

##List of EnrichedIntersectionSets

#4) Network construction: Populate with EnrichedPathways
#
#5) If EnrichedIntersectionSets[i,j] != cutoff, link EnrichedPathways[i,j]
##Return Network
DF <- data.frame(do.call(rbind, EnrichedIntersectionSets)) #contains only significant edges
library(igraph)
g = graph_from_edgelist(as.matrix(DF))
#make sure to add unconnected pathways that may be enriched
#

if(length(setdiff(EnrichedPathways, V(g)$name)) != 0) {
  add.vertices(graph = g, 
               nv = length(setdiff(EnrichedPathways, V(g)$name)), 
               name = setdiff(EnrichedPathways, V(g)$name))
}

#Test
# gg= g
# eptt = c(EnrichedPathways, "hugo", "paco", "luis")
# if(length(setdiff(eptt, V(gg)$name)) != 0) {
#   gg = add.vertices(graph = gg, 
#                nv = length(setdiff(eptt, V(gg)$name)), 
#                name = setdiff(eptt, V(gg)$name))
# }
# plot(gg)
