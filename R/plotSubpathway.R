#plot the signature network of subpathway
plotsubpathway_Node<-function(sub,mapkG3,de=DE_Colorectal,all=ALL_Colorectal){
  library(Rgraphviz)
  library(RColorBrewer)
  library(org.Hs.eg.db)
  library(RBGL)
  library(grid)
  
  allID<-nodes(sub)
  diffID<-getdiffNode(mapkG3,de)
  diffs<-allID[allID %in% diffID]
  logcol <- c() 
  logcol[allID]<-"gray"
  logcol[match(diffs,allID)] <- "red"
  names(logcol) <- allID
  nA <- makeNodeAttrs(sub,fillcolor=logcol, width=2, height=2)
  plot(sub, "dot", nodeAttrs=nA)
}


#change the signature network of subpathway inte the gene newtork of subpathway
KEGGsubpathway2subGraph<-function (pathway,diffgenes, genesOnly = TRUE, expandGenes = TRUE) 
{
  stopifnot(is(pathway, "KEGGPathway"))
  #pathway <- splitKEGGgroup(pathway)
  pathway <- splitKEGGgroup_compound(pathway)
  if (expandGenes) {
     pathway <- expandKEGGsubPathway_gene(pathway,diffgenes)
  }
  knodes <- nodes(pathway)
  kedges <- unique(attributes(pathway)$edges)
  node.entryIDs <- getEntryID(knodes)
  edge.entryIDs <- getEntryID(kedges)
  V <- node.entryIDs
  edL <- vector("list", length = length(V))
  names(edL) <- V
  if (is.null(nrow(edge.entryIDs))) {
    for (i in seq(along = edL)) {
      edL[[i]] <- list()
    }
  }
  else {
    for (i in 1:length(V)) {
      id <- node.entryIDs[i]
      hasRelation <- id == edge.entryIDs[, "Entry1ID"]
      if (!any(hasRelation)) {
        edL[[i]] <- list(edges = NULL)
      }
      else {
        entry2 <- unname(edge.entryIDs[hasRelation, "Entry2ID"])
        edL[[i]] <- list(edges = entry2)
      }
    }
  }
  gR <- new("graphNEL", nodes = V, edgeL = edL, edgemode = "directed")
  names(kedges) <- sapply(kedges, function(x) paste(getEntryID(x), 
                                                    collapse = "~"))
  env.node <- new.env()
  env.edge <- new.env()
  assign("nodes", knodes, envir = env.node)
  assign("edges", kedges, envir = env.edge)
  nodeDataDefaults(gR, "KEGGNode") <- env.node
  edgeDataDefaults(gR, "KEGGEdge") <- env.edge
  if (genesOnly) {
    gR <- subGraphByNodeType(gR, "gene")
  }
  return(gR)
}
