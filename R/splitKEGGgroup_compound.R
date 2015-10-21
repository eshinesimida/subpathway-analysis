splitKEGGgroup_compound<-function (pathway){
  pathway <- splitKEGGgroup(pathway)
  pnodes <- nodes(pathway)
  
  pedges <- attributes(pathway)$edges
  if (length(pedges) == 0) 
    return(pathway)
  types <- sapply(pnodes, getType)
  if (any(types == "compound")) {
    isCompound <- names(pnodes)[types == "compound"]
    edgeEntry <- sapply(pedges, getEntryID)
    compound<-c()
    for(i in 1:length(isCompound)){
      compoundAsID<- isCompound[i] %in% edgeEntry[1L,] && isCompound[i] %in% edgeEntry[2L,]
      
      compound<-c(compound,compoundAsID)
    }
    isCompound[compound]
    if(length(isCompound[compound])==0){
      return(pathway)
    }
    dim(edgeEntry)
    edgeRelation<-t(edgeEntry)
    compoundAsID <- edgeEntry[1L, ] %in% isCompound[compound] | edgeEntry[2L, 
                                                                          ] %in% isCompound[compound]
    newly <- list()
    for (i in 1:length(isCompound[compound])) {
      e<-pedges[which(edgeRelation[,2]==isCompound[compound][i])]
      node1<-edgeRelation[which(edgeRelation[,2]==isCompound[compound][i]),1]
      node2<-edgeRelation[which(edgeRelation[,1]==isCompound[compound][i]),2]
      expandmodel <- expand.grid(node1, node2)
      
      enews <- list()
      for (j in 1:nrow(expandmodel)) {
        enews[[j]] <-e[[1]]
        enews[[j]]@entry1ID<-as.character(expandmodel[j,1])
        enews[[j]]@entry2ID<-as.character(expandmodel[j,2])
        enews[[j]]@type<-"PPrel"
        names(enews)<-"relation"
      }
      newly<-c(newly,enews)
    }
    newEdges <- pedges[!compoundAsID]
    newEdges <- c(newEdges, newly)
    edges(pathway) <- newEdges
  }
  return(pathway)
}
