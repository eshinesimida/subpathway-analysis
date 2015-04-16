#Generate subpathway based on MST alogrithm.

subpathway_set<-function(g_u,mapkpathway,de=DE_Colorectal,n=4){
  library(RBGL)
  diffNode<-getdiffNode(g_u,de)
  distance<-c()
  g_igraph<-igraph.from.graphNEL(g_u)
  a<-clusters(g_igraph)
  b<-which(a$csize>=5)
  if(length(b)==0)
    return(NULL)
  G<-list()
  for(i in 1:length(b)){
    index<-which(a$membership==b[i])
    G[[i]]<-V(g_igraph)$name[index]
  }
  short_path<-list()
  for(i in 1:length(b)){
    g<-subGraph(G[[i]],g_u)
    diffs<-diffNode[diffNode %in% nodes(g)]
    if(length(diffs)<2) next
    distance<-c()
    for(j in 1:length(diffs)){
      for(k in 1:length(diffs)){
        distance<-c(distance,sp.between(g,diffs[j],diffs[k])[[1]]$length)
      }
    }
    distance_m<-matrix(distance,nrow=length(diffs),ncol=length(diffs),byrow=T)
    rownames(distance_m)<-diffs
    colnames(distance_m)<-diffs
    diag(distance_m)<-0
    g1<-graph.adjacency(distance_m,mode="undirected",weighted=T)
    g2<-igraph.to.graphNEL(g1)
    m<-mstree.kruskal(g2)
    shortest_path_set<-c()
    if(max(m$weights)-1<=n+1){
      for(j in 1:length(m$weights)){
        c<-sp.between(g_u,m$edgeList[1,j],m$edgeList[2,j])[[1]]$path_detail
        shortest_path_set<-c(shortest_path_set,c)
      }
    }
    shortest_path_set<-unique(shortest_path_set)
    if(length(shortest_path_set)<5) next
    short_path[[i]]<-shortest_path_set
    
  }
  if(length(short_path)==0)
    return(NULL)
  short_path<-short_path[sapply(short_path,function(x) {!is.null(x)})]
  subpathID<-c()
  for(i in 1:length(short_path)){
    subpathID<-c(subpathID,paste(attributes(mapkpathway)$pathwayInfo@number,"_",i,sep=""))
  }
  names(short_path)<-subpathID
  return(short_path)
}