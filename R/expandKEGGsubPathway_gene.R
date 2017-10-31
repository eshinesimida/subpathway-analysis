expandKEGGsubPathway_gene<-function (pathway,subNode=diffgenes) 
{
  nodes.old <- nodes(pathway)
  nodes.old <- subNode
  edges.old <- attributes(pathway)$edges
  nodes.new <- list()
  entryMap <- list()
  for (i in seq(along = nodes.old)) {
    expanded <- expandKEGGNode(nodes.old[[i]])
    newEntryIDs <- sapply(expanded, getEntryID)
    names(expanded) <- newEntryIDs
    nodes.new[[i]] <- expanded
    oldEntryID <- getEntryID(nodes.old[[i]])
    entryMap[[i]] <- data.frame(oldEntryID = I(oldEntryID), 
                                newEntryID = I(newEntryIDs))
  }
  nodes.new <- unlist(nodes.new)
  entryMap <- do.call(rbind, entryMap)
  isDuplicatedNode <- duplicated(sapply(nodes.new, getEntryID))
  nodes.new <- nodes.new[!isDuplicatedNode]
  edges.new <- list()
  for (i in seq(along = edges.old)) {
    edge.old <- edges.old[[i]]
    entryIDs.old <- getEntryID(edge.old)
    entry1ID.new <- with(entryMap, newEntryID[oldEntryID == 
                                                entryIDs.old[1L]])
    entry2ID.new <- with(entryMap, newEntryID[oldEntryID == 
                                                entryIDs.old[2L]])
    if (!(length(entry1ID.new) >= 1 & length(entry2ID.new) >= 
          1)) {
      warning("Missing entries detected in the KGML file. If it is not a KO file, please check its integrity\n")
      next
    }
    expand <- expand.grid(entry1ID.new, entry2ID.new)
    edge.new <- list()
    tmp <- edge.old
    for (j in 1:nrow(expand)) {
      entryID(tmp) <- c(as.character(expand[j, 1]), as.character(expand[j, 
                                                                        2]))
      edge.new[[j]] <- tmp
    }
    edges.new[[i]] <- edge.new
  }
  edges.new <- unlist(edges.new)
  pathway.new <- pathway
  nodes(pathway.new) <- as.list(nodes.new)
  edges(pathway.new) <- as.list(edges.new)
  return(pathway.new)
}