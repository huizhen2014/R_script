##clusterprofiler results
#test_mf <- clusterProfiler_result_mf
##获取test_mf 中go term的children/parents go term，且要存在test_mf中
library(GO.db)
library(igraph)
links_mf <- data.frame()
i<-0
for(id in test_mf$ID){
  if(!is.na(GOMFOFFSPRING[[id]])){
    for(child in GOMFOFFSPRING[[id]]){
      print (child)
      if(child %in% test_mf$ID){
        i <- i+1
        links_mf[i,"from"] <- id
        links_mf[i,"to"] <- child
      }
    }
  }else if(!is.na(GOMFANCESTOR[[id]])){
    for(parent in GOMFANCESTOR[[id]]){
      if(parent %in% test_mf$ID){
        i <- i+1
        links_mf[i,"from"] <- parent
        links_mf[i,"to"] <- id
      }
    }
  }
}
nodes_mf <-test_mf[,c("ID","qvalue","Description","Count")]

net_mf <- graph_from_data_frame(d=links_mf,vertices = nodes_mf,
                                directed=T)
col_num <- length(unique(nodes_mf$qvalue))
col <- colorRampPalette(c("green","red"))(col_num)
col_value <- match(nodes_mf$qvalue,sort(unique(nodes_mf$qvalue)))
net_mf <- simplify(net_mf,remove.multiple=T)
V(net_mf)$color <- col[col_value]
V(net_mf)$size <- V(net_mf)$Count*0.7
V(net_mf)$label <- nodes_mf$Description
V(net_mf)$label.size <- V(net_mf)$Count * 0.0005
V(net_mf)$vertex.label.dist <- 2
V(net_mf)$label.cex <- 0.5
E(net_mf)$arrow.size <- 0.01
l <- layout_with_dn(net_mf)
plot(net_mf,layout=l)

##hubs and authorities
hs <- hub_score(net_mf,weight=NA)$vector
as <- authority_score(net_mf,weight=NA)$vector
plot(net_mf,vertex.size=hs*50,main="hubs")
plot(net_mf,vertex.size=as*30,main="Authorities")

##layout
V(net_mf)$size <- V(net_mf)$Count
l <- layout_with_graphopt(net_mf)
plot(net_mf,layout=l)
l <- layout_with_kk(net_mf)
plot(net_mf,layout=l)
plot(net_mf,layout=layout_with_lgl)



