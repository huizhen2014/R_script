
group_name <- unique(group$site)
dir.create("anosim_two",recursive=T)
anosim_result_two <- NULL

for( i in 1:(length(group_name)-1)){
  for(j in (i+1):length(group_name)){
    group_ij <- subset(group,site %in% c(group_name[i],group_name[j]))
    otu_ij <- otu[group_ij$names,]
    anosim_result_otu_ij <- anosim(otu_ij,group_ij$site,permutations = 999,distance = "bray")
    
    if(anosim_result_otu_ij$signif <= 0.001)Sig <- "***"
    else if(anosim_result_otu_ij$signif <= 0.01)Sig <- "**"
    else if(anosim_result_otu_ij$signif <= 0.05)Sig <- "*"
    else Sig <- NA
    
    anosim_result_two <- rbind(anosim_result_two,
                               c(paste(group_name[i],group_name[j],sep="/"),
                                 "Bray-Curtis",anosim_result_otu_ij$statistic,
                                 anosim_result_otu_ij$signif,Sig))
    pdf(paste("anosim_two/anosim.",group_name[i],"_",group_name[j],".pdf",sep=""),
        width=7,height=5)
    plot(anosim_result_otu_ij,col=c("gray","red","blue"))
    dev.off()
  }
}

anosim_result_two <- data.frame(anosim_result_two,stringsAsFactors = F)
names(anosim_result_two) <- c("group","distance","R","P_value","Sig")






