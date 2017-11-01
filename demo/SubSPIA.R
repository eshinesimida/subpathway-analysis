library(graphite)
library(Rgraphviz)
library(SPIA)
library(XML)

data(colorectalcancer)
library(hgu133plus2.db)
x<-hgu133plus2ENTREZID
top$ENTREZ<-unlist(as.list(x[top$ID]))
top<-top[!is.na(top$ENTREZ),]
top<-top[!duplicated(top$ENTREZ),]
tg1<-top[top$adj.P.Val<0.05,]
DE_Colorectal=tg1$logFC
names(DE_Colorectal)<-as.vector(tg1$ENTREZ)
ALL_Colorectal<-top$ENTREZ

#kgml.path=system.file("extdata/keggxml/hsa",package="SPIA")
#mdir=kgml.path
#paths<-dir(mdir,pattern=".xml")
#mapkpathway<-try(parseKGML(paste(mdir,paths[2],sep="/")),TRUE)
L<-list()

rel<-c("activation",
       "compound",
       "binding/association",
       "expression",
       "inhibition",
       "activation_phosphorylation",
       "phosphorylation",
       "inhibition_phosphorylation",
       "inhibition_dephosphorylation",
       "dissociation",
       "dephosphorylation",
       "activation_dephosphorylation",
       "state change",
       "activation_indirect effect",
       "inhibition_ubiquination",
       "ubiquination",
       "expression_indirect effect",
       "inhibition_indirect effect",
       "repression",
       "dissociation_phosphorylation",
       "indirect effect_phosphorylation",
       "activation_binding/association",
       "indirect effect",
       "activation_compound",
       "activation_ubiquination"
)
betas=c(1,0,0,1,-1,1,0,-1,-1,0,0,1,0,1,-1,0,1,-1,-1,0,0,1,0,1,1)
names(betas)=rel

L[[1]]<-NA
info<-NULL
info<-rbind(info,c(NA,NA))
out.path="."
organism="hsa"


kgml.path=system.file("extdata/keggxml/hsa",package="SPIA")
mdir=kgml.path
paths<-dir(mdir,pattern=".xml")
info<-NULL
alr<-NULL
L_all<-list()
for(p in 1:length(paths)){
  L<-list()
  info<-NULL
  #  info<-rbind(info,c(NA,NA))
  mapkpathway<-try(parseKGML(paste(mdir,paths[p],sep="/")),TRUE)
  #mapkG3<-KEGGpathway2Graph(mapkpathway,expandGenes=F)
  mapkG3<-KEGGpathway2Graph_c(mapkpathway,expandGenes=F)
  
  #mapkG3 transform to undirected graph
  g_u<-Convert_u(mapkG3)
  #find subpathway
  
  subpathwaylist<-subpathway_set(g_u,mapkpathway,de=DE_Colorectal,n=4)
  if(length(subpathwaylist)==0)
    next
  for(i in 1:length(subpathwaylist)){
    info1<-c(names(subpathwaylist)[i],mapkpathway@pathwayInfo@title)
    info<-rbind(info,info1)
  }
  
  for(q in 1:length(subpathwaylist)){
    L[[q]]<-NA
    sub<-subGraph(subpathwaylist[[q]],mapkG3)
    diffNodes<-nodes(sub)
    allNodes <- getKEGGnodeData(mapkG3)
    diffgenes<-allNodes[match(diffNodes,names(allNodes))]
    mapk<-splitKEGGgroup(mapkpathway)
    #edg<-mapkpathway@edges
    edg<-attributes(mapk)$edges
    edg_id<-getEntryID(edg)
    #pathway <- splitKEGGgroup_compound(mapkpathway)
    
    refty<-unlist(lapply(diffgenes,function(x){x@type}))
    
    index<-names(sub@edgeData@data)
    
    unlist(strsplit(index[1],"[|]"))
    outdegree<-sapply(strsplit(index,split="[|]"),function(x){x[1]})
    indegree<-sapply(strsplit(index,split="[|]"),function(x){x[2]})
    sub_edge<-data.frame(start=outdegree,end=indegree)
    index_sub<-c()
    for(j in 1:length(rownames(sub_edge))){
      for(i in 1:length(edg)){
        #if(edg[i]$relation@entry1ID==sub_edge[j,1] && edg[i]$relation@entry2ID==sub_edge[j,2])
        if(edg_id[i,1]==sub_edge[j,1] && edg_id[i,2]==sub_edge[j,2]) 
          index_sub<-c(index_sub,i)
      }
    }
    if (length(index_sub)==0)
      next
    edg_sub<-edg[index_sub]# edge set of subpathway
    
    
    G<-list()
    for(kk in 1:length(diffgenes)){
      if((diffgenes[[kk]])@type%in%c("gene","group")){
        if((diffgenes[[kk]])@type=="gene"){
          tmpp<-(diffgenes[[kk]])@name 
          G[[kk]]<-unlist(lapply(strsplit(tmpp,split=":"),function(x){x[2]}))
        } 
        if((diffgenes[[kk]])@type=="group"){
          gc=na.omit((diffgenes[[kk]])@component[refty[(diffgenes[[kk]])@component]=="gene"])
          if(length(gc)>=1){
            tmpp<- unlist(lapply((diffgenes[gc]),function(x){x@name}))
            G[[kk]]<-unlist(lapply(strsplit(tmpp,split=":"),function(x){x[2]}))
          }else{G[[kk]]<-"other"}
        } 
        
      }else{G[[kk]]<-"other"}
    } 
    
    nods<-names(diffgenes)
    names(G)<-nods
    
    
    N<-length(nods)
    
    LT<-list()
    for ( j in 1:length(rel)){
      LT[[j]]<-matrix(0,N,N)
      rownames(LT[[j]])<-colnames(LT[[j]])<-nods
    }
    names(LT)<-rel
    
    for(i in 1:length(edg_sub)){
      #tmp<-edg_sub[i]$rel
      tmp<-edg_sub[[i]]
      #tmp2<-edg_sub[i]$rel@subtype
      tmp2<-tmp@subtype
      if(length(tmp2)>0){
        zz<-NULL
        
        for(k in 1:length(tmp2)){
          zz<-c(zz,tmp2[k]$subtype@name)
        }
        zz<-sort(zz)
        zz<-zz[!zz%in%c(""," ")]
        zz=paste(zz,collapse="_")
        
        if(zz=="dephosphorylation_inhibition"){zz="inhibition_dephosphorylation"}
        if(zz=="indirect effect_inhibition"){zz="inhibition_indirect effect"}
        if(zz=="activation_indirect effect_phosphorylation"){zz=c("activation_indirect effect","activation_phosphorylation")}
        
        for(jj in 1:length(rel)){
          LT[[jj]][tmp@entry2ID,tmp@entry1ID]<-ifelse(rel[jj]%in%zz,1,0)
        }
        alr <-c(alr,paste(zz,collapse="_"))
      }#end for edges
      
      
      LT[["nodes"]]<-G
      LT[["title"]]<-mapkpathway@pathwayInfo@title
      LT[["NumberOfReactions"]]<-length(mapkpathway@reactions)
      
    }
    
    L[[q]]<-LT
  }
  nok<-is.na(L)
  L[nok]<-NULL
  info<-info[!nok,]
  info=rbind(info)
  
  names(L)<-info[,1]
  L_all<-c(L_all,L)
}
L<-list()
L<-L_all
nL<-L

for (ll in 1:length(L)){
  
  cornodes<-L[[ll]]$nodes[!is.na(L[[ll]]$nodes) & L[[ll]]$nodes!="other"]
  allgns<- unique(unlist(cornodes))
  
  for(re in rel){
    L[[ll]][[re]]<-L[[ll]][[re]][names(cornodes),names(cornodes)]
    nL[[ll]][[re]]<-matrix(0,length(allgns),length(allgns))
    rownames(nL[[ll]][[re]])<-colnames(nL[[ll]][[re]])<-allgns
    for(n in 1:length(cornodes)){
      if(cornodes[n]!="other"){
        nd<- names(cornodes)[n]
        s=names(L[[ll]][[re]][,nd][L[[ll]][[re]][,nd]!=0])
        length(s)
        
        if (length(s)>=1){
          for (jb in 1:length(s)){
            frm=unlist(cornodes[n]) 
            to=unlist(cornodes[s[jb]])
            if(!is.null(frm)& !is.null(to)){
              IND<-cbind(rep(frm,each=length(to)),rep(to,length(frm)))
              for(kb in 1:dim(IND)[1]){
                nL[[ll]][[re]][IND[kb,2],IND[kb,1]]<-1
              }
            }#end if
          }
        }#if downstream genes
      }#if valid node
    }#for nodes
  }#for reltion
  nL[[ll]]$nodes<-allgns
  
}

sumrel<-unlist(lapply(nL,function(x){sum(unlist(lapply(x[names(betas[betas!=0])],sum)))>0}))

path.info<-nL[unlist(lapply(nL,function(x){x$NumberOfReactions}))<1 & sumrel]

save(path.info,file=paste(out.path,"/",organism,"SPIA.RData",sep=""))

res_c2<-spia(de=DE_Colorectal, all=ALL_Colorectal, organism="hsa",data.dir="./")
res_c[,-12]

