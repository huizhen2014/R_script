##Enrichment map organizes enriched terms into a network with edges 
##connecting overlapping gene sets. In this way, mutually overlapping
##gene sets are tend to cluster together, making it easy to identify
##functional module. The results obtained from Hytpergeometric test.
##挑选p.adjust值(q值)小于等于0.05的GO terms或挑选前30个GO terms，

##MF BP CC only print the qvalue <= 0.05
##Relationships between GO terms
library(igraph)
library(reshape2)
library(ggplot2)
up_igraph_results_table <- list()
down_igraph_results_table <- list()
##DE up
for(i in 1:3){
  num <- 0
  vertices <- up_go_results_table[[i]][,c(1,2,4,7)]
  colnames(vertices) <- c("GO.ID","Term","Significant","qvalue")
  for(j in 1:(nrow(up_go_results_table[[i]])-1)){
    from <- up_go_results_table[[i]][j,]$GO.ID
    from_genes <- unlist( strsplit( up_go_results_table[[i]][j,]$Sig_Genes,","))
    for(k in (j+1):nrow(up_go_results_table[[i]])){
      end <- up_go_results_table[[i]][k,]$GO.ID
      end_genes <- unlist( strsplit( up_go_results_table[[i]][k,]$Sig_Genes,","))
      from_end_num <- sum(from_genes %in% end_genes)
      num <- num+1
      links[num,"from"] <- from
      links[num,"end"] <- end
      links[num,"count"] <- from_end_num
    }
  }
  links <- links[links$count>0,]
  up_igraph_results_table[[i]] <- list(vertices,links)
}
##draw up the netwroks
for(n in 1:3){
  tmp <- up_igraph_results_table[[n]]
  vertices <- tmp[[1]]
  d <- tmp[[2]]
  net <- graph_from_data_frame(
    d=d,vertices=vertices,directed = F)
  rbPal <- colorRampPalette(c("yellow","red"))
  V(net)$color <- rbPal(10)[as.numeric(cut(-log10(vertices$qvalue),breaks = 10))]
  V(net)$label <- vertices$Term
  V(net)$label.family <- "Times"
  V(net)$label.cex <- 0.6
  V(net)$label.dist <- 0.5
  E(net)$width <- d$count*0.15
  
  pdf(paste0("28vs21_up_",go_type[n],"_star_network.pdf"))
  plot(net,layout=layout_as_star,main=paste0(
    "28vs21_up_",go_type[n],"_star_network"))
  dev.off()
  
  pdf(paste0("28vs21_up_",go_type[n],"_network.pdf"))
  plot(net,main=paste0(
    "28vs21_up_",go_type[n],"_network"))
  dev.off()
}

##DE down    
for(i in 1:3){
  num <- 0
  vertices <- data.frame()
  links <- data.frame()
  vertices <- down_go_results_table[[i]][,c(1,2,4,7)]
  colnames(vertices) <- c("GO.ID","Term","Significant","qvalue")
  for(j in 1:(nrow(down_go_results_table[[i]])-1)){
    from <- down_go_results_table[[i]][j,]$GO.ID
    from_genes <- unlist( strsplit( down_go_results_table[[i]][j,]$Sig_Genes,","))
    for(k in (j+1):nrow(down_go_results_table[[i]])){
      end <- down_go_results_table[[i]][k,]$GO.ID
      end_genes <- unlist( strsplit( down_go_results_table[[i]][k,]$Sig_Genes,","))
      from_end_num <- sum(from_genes %in% end_genes)
      num <- num+1
      links[num,"from"] <- from
      links[num,"end"] <- end
      links[num,"count"] <- from_end_num
    }
  }
  links <- links[links$count>0,]
  down_igraph_results_table[[i]] <- list(vertices,links)
}

##draw down the netwroks
for(n in 1:3){
  tmp <- down_igraph_results_table[[n]]
  vertices <- tmp[[1]]
  d <- tmp[[2]]
  net <- graph_from_data_frame(
    d=d,vertices=vertices,directed = F)
  rbPal <- colorRampPalette(c("yellow","red"))
  V(net)$color <- rbPal(10)[as.numeric(cut(-log10(vertices$qvalue),breaks = 10))]
  V(net)$label <- vertices$Term
  V(net)$label.family <- "Times"
  V(net)$label.cex <- 0.6
  V(net)$label.dist <- 0.5
  E(net)$width <- d$count*0.15
  
  pdf(paste0("28vs21_down_",go_type[n],"_star_network.pdf"))
  plot(net,layout=layout_as_star,main=paste0(
    "28vs21_down_",go_type[n],"_star_network"))
  dev.off()
  
  pdf(paste0("28vs21_down_",go_type[n],"_network.pdf"))
  plot(net,main=paste0(
    "28vs21_down_",go_type[n],"_network"))
  dev.off()
}

##up heatmap of GO terms with genes
for(m in 1:3){
  genes_up <- vector()
  tmp=up_go_results_table[[m]]
  for(i in 1:nrow(tmp)){
    genes_up <- append(genes_up, unlist(strsplit(tmp[i,]$Sig_Genes,",")))
  }
  genes_up <- sort(unique(genes_up))
  
  Data <- data.frame(matrix(1:length(genes_up),nrow=1))
  for(j in 1:nrow(tmp)){
    Data[j,] <- as.integer(genes_up %in% unlist(strsplit(tmp[j,]$Sig_Genes,",")))
  }
  colnames(Data) <- factor(genes_up,levels=genes_up)
  rownames(Data) <- factor(tmp$GO.ID,levels=rev(tmp$GO.ID))
  x1 <- vector()
  x2 <- vector()
  y1 <- vector()
  y2 <- vector()
  q <- vector()
  for(k in 1:nrow(Data)){
    for(n in 1:ncol(Data)){
      x1 <- append(x1,Data[k,n]*(n-0.45))
      x2 <- append(x2,Data[k,n]*(n+0.45))
      y1 <- append(y1, Data[k,n]*(k-0.45))
      y2 <- append(y2, Data[k,n]*(k+0.45))
    }
    q <- append(q,rep(tmp[k,]$qvalue,length(genes_up)))
  }
  d <- data.frame(x1=x1,x2=x2,y1=y1,y2=y2,q=q)
  
p <- ggplot() + theme_bw()+ geom_rect(
  data=d,mapping=aes(xmin=x1,xmax=x2,ymin=y1,ymax=y2,fill=q))+
  scale_fill_gradient(low="red",high="blue")+
  scale_y_continuous(breaks=seq(1,length(tmp$GO.ID)),labels=tmp$GO.ID,expand = c(0,0))+
  scale_x_continuous(breaks=seq(1,length(genes_up)),labels=genes_up,expand = c(0,0))+
    theme(axis.text.x=element_text(angle=60,vjust=1,hjust=1,size=6.5),
          plot.title=element_text(hjust = 0.5))+
  labs(title=paste0(go_type[m],"_Heatmap"),y="GO Terms",x="DE Genes",fill="Qvalue")
ggsave(paste0("./GO_enrichment_results/","Kp_28vs21_Up_",go_type[m],"_heapmap.pdf"),
       plot=p,width = 28,height=18,units = "cm")
}

##down heatmap of GO terms with genes
for(m in 1:3){
  genes_down <- vector()
  tmp=down_go_results_table[[m]]
  for(i in 1:nrow(tmp)){
    genes_down <- append(genes_down, unlist(strsplit(tmp[i,]$Sig_Genes,",")))
  }
  genes_down <- sort(unique(genes_down))
  
  Data <- data.frame(matrix(1:length(genes_down),nrow=1))
  for(j in 1:nrow(tmp)){
    Data[j,] <- as.integer(genes_down %in% unlist(strsplit(tmp[j,]$Sig_Genes,",")))
  }
  colnames(Data) <- factor(genes_down,levels=genes_down)
  rownames(Data) <- factor(tmp$GO.ID,levels=rev(tmp$GO.ID))
  x1 <- vector()
  x2 <- vector()
  y1 <- vector()
  y2 <- vector()
  q <- vector()
  for(k in 1:nrow(Data)){
    for(n in 1:ncol(Data)){
      x1 <- append(x1,Data[k,n]*(n-0.45))
      x2 <- append(x2,Data[k,n]*(n+0.45))
      y1 <- append(y1, Data[k,n]*(k-0.45))
      y2 <- append(y2, Data[k,n]*(k+0.45))
    }
    q <- append(q,rep(tmp[k,]$qvalue,length(genes_down)))
  }
  d <- data.frame(x1=x1,x2=x2,y1=y1,y2=y2,q=q)
  
  p <- ggplot() + theme_bw()+ geom_rect(
    data=d,mapping=aes(xmin=x1,xmax=x2,ymin=y1,ymax=y2,fill=q))+
    scale_fill_gradient(low="red",high="blue")+
    scale_y_continuous(breaks=seq(1,length(tmp$GO.ID)),labels=tmp$GO.ID,expand = c(0,0))+
    scale_x_continuous(breaks=seq(1,length(genes_down)),labels=genes_down,expand = c(0,0))+
    theme(axis.text.x=element_text(angle=60,vjust=1,hjust=1,size=6.5),
          plot.title=element_text(hjust = 0.5))+
    labs(title=paste0(go_type[m],"_Heatmap"),y="GO Terms",x="DE Genes",fill="Qvalue")
  ggsave(paste0("./GO_enrichment_results/","Kp_28vs21_Down_",go_type[m],"_heapmap.pdf"),
         plot=p,width = 28,height=18,units = "cm")
}


