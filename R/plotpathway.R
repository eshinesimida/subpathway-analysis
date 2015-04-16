#plot the signature network
plotpathway_Node<-function(mapkG3,de=DE_Colorectal,all=ALL_Colorectal){
  library(Rgraphviz)
  library(RColorBrewer)
  library(org.Hs.eg.db)
  library(RBGL)
  library(grid)
  mapkGnodedata <- getKEGGnodeData(mapkG3)
  allID<-names(mapkGnodedata)
  diffID<-getdiffNode(mapkG3,de)
  logcol <- c() 
  logcol[allID]<-"gray"
  logcol[match(diffID,allID)] <- "red"
  names(logcol) <- allID
  nA <- makeNodeAttrs(mapkG3,fillcolor=logcol, width=2, height=2)
  plot(mapkG3, "dot", nodeAttrs=nA)
}

#plot the gene network
plotpathway_Exgene<-function(mapkG2,de=DE_Colorectal,all=ALL_Colorectal){
  library(Rgraphviz)
  library(RColorBrewer)
  library(org.Hs.eg.db)
  library(RBGL)
  library(grid)
  deKID <- names(de) #differential genes
  allKID <- all#all genes
  geneID<-sapply(strsplit(nodes(mapkG2),":"),function(x) {x[2]})
  isDiffExp <- geneID %in% deKID
  logfcs <- DE_Colorectal[match(geneID, deKID)]
  names(logfcs) <- geneID
  logfcs[is.na(logfcs)] <- 0
  undetected <- !geneID %in% allKID
  logcol <- c() 
  logcol[logfcs==0] <- "darkgrey"
  logcol[undetected] <- "yellow"
  logcol[logfcs>0] <- "red"
  logcol[logfcs<0] <- "blue"
  names(logcol) <- names(logfcs)
  nA <- makeNodeAttrs(mapkG2,fillcolor=logcol, width=2, height=2)
  plot(mapkG2, "dot", nodeAttrs=nA)
}
