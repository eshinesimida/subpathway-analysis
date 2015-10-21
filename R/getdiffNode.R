getdiffNode<-function(mapkG3,de=DE_Colorectal){
  mapkGnodedata <- getKEGGnodeData(mapkG3)
  node_pathway<-list()
  for(i in 1:length(mapkGnodedata)){
    node_pathway[[i]]<-sapply(strsplit(mapkGnodedata[[i]]@name,":"),function(x) {x[2]})
  }
  names(node_pathway)<-names(mapkGnodedata)
  node_ID<-c()
  for(i in 1:length(node_pathway)){
    if(any(node_pathway[[i]] %in% names(de))){
      node_ID<-c(node_ID,names(node_pathway[i]))
    }
  }
  if(length(node_ID)==0){
    return(0)
  }
  return(node_ID)
}
