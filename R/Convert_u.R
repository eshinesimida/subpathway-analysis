#Change directed graph into undirected graph

Convert_u<-function(g){
  library(igraph)
  g_igraph<-igraph.from.graphNEL(g)
  g_u<-as.undirected(g_igraph)
  g_d<-igraph.to.graphNEL(g_u)
}