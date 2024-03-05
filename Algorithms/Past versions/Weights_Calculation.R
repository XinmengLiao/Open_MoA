library(igraph)

### Weighted Network Construction
jnk_scorenetwork <- read.csv('HepG2_specific.csv',header = T)
jnk_scorenetwork <- unique(jnk_scorenetwork)
which(jnk_scorenetwork$common_score != 1)
jnk_fdr<- read.table('TGFB1.csv',sep = ",",header = T)
jnk_fdr <- jnk_fdr[which(jnk_fdr$FDR<1),]
jnk_fdr <- jnk_fdr[which(jnk_fdr$FDR != "NA"),]
jnk_fdr <- unique(jnk_fdr)
jnk_fdr <- jnk_fdr[which(jnk_fdr$Gene %in% jnk_scorenetwork$Vertex.A | jnk_fdr$Gene %in% jnk_scorenetwork$Vertex.B),]
length(which(duplicated(jnk_fdr$Gene) == T))

### Delete the TGFB1 perturbed genes which do not have regulatory relationship
fdr <- NULL
for(i in 1:nrow(jnk_fdr)){
  print(i)
  index1 <- jnk_scorenetwork[which(jnk_scorenetwork$Vertex.B == jnk_fdr$Gene[i]),]
  index2 <- index1[which(index1$Relationship == 'Regulatory'),5]
  if(length(index2 > 0)){
  fdr <- rbind(fdr,jnk_fdr[i,])
  } else{
  fdr <- fdr
  }
}
jnk_fdr <- fdr

edges = data.frame(from = jnk_scorenetwork$Vertex.A,
                   to = jnk_scorenetwork$Vertex.B)
combine <- graph_from_data_frame(edges, directed = T)


### Penalty Scores calculation

for (j in 1:nrow(jnk_fdr)){
  print(j) 
  jGene <- jnk_fdr[j,1] 
  jRMinds <- which(jnk_scorenetwork$Vertex.B == jGene & jnk_scorenetwork$Relationship != "Regulatory")
  if (length(jRMinds)>0) {
    jEdges <- data.frame(from = jnk_scorenetwork$Vertex.A[-jRMinds],
                         to = jnk_scorenetwork$Vertex.B[-jRMinds])
  } else {
    jEdges <- data.frame(from = jnk_scorenetwork$Vertex.A,
                         to = jnk_scorenetwork$Vertex.B)
  }
  jCombine <- graph_from_data_frame(jEdges, directed = T)
  j_result <- get.shortest.paths(jCombine,from = 'TGFB1',to = jGene,output = 'both',
                                 mode='out')
  j_tf_subnetwork <- edges[j_result$epath[[1]],]
  j_index <- nrow(j_tf_subnetwork)
  if (j_index>0) {
    j_all <- all_simple_paths(graph = jCombine, from = 'TGFB1',to = jGene,
                              cutoff = j_index,mode = 'out' )
    for (k in 1:length(j_all)) {
      each_a <- data.frame(unlist(j_all[k],use.name=T))
      index <- as.numeric(nrow(each_a))
      each_a$name <- rownames(each_a)
      rownames(each_a) <- NULL
      each_b <- data.frame(each_a[-1,-1])
      each_a <- data.frame(each_a[-index,-1])
      merge_1 <- cbind(each_a,each_b)
      colnames(merge_1) <- c('Vertex.A','Vertex.B')
      merge_1$common_score <- jnk_fdr[j,2]^(1/length(j_all)) 
      for (i in 1:nrow(merge_1)) {
        i_index1 <- which(jnk_scorenetwork$Vertex.A == merge_1[i,1] & jnk_scorenetwork$Vertex.B == merge_1[i,2])
        jnk_scorenetwork[i_index1,6] <- jnk_scorenetwork[i_index1,6] * merge_1[i,3]
      }
    }
  }
}


write.csv(jnk_scorenetwork,'TGFb1_specific_weighted_subnetwork.csv',row.names = F)
