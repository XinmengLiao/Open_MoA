library(xlsx)
library(igraph)

### Weighted Network Construction
jnk_scorenetwork <- read.csv('~/Desktop/Sweden/Scorenetwork 2.0/Datapreparation/Hepg2_Network_top10_HpaCcle_TPMmax>1.csv',header = T)
which(jnk_scorenetwork$common_score != 1)
jnk_fdr <- read.csv('TGFB1_trt_sh_1.csv')
metformin <- read.table('Metformin/HEPG2_gene_padj_final.txt',header = T)
jnk_fdr <- jnk_fdr[which(jnk_fdr$FDR<1),]
jnk_fdr <- jnk_fdr[which(jnk_fdr$FDR != "NA"),]
jnk <- metformin[,c(1,4)]
jnk <- jnk[which(jnk[,2] < 1 & jnk[,2] != 'NA'),]
jnk_fdr <- jnk
colnames(jnk_fdr)[1] <- 'Gene'
jnk_fdr <- jnk_fdr[which(jnk_fdr$Gene %in% jnk_scorenetwork$Vertex.A | jnk_fdr$Gene %in% jnk_scorenetwork$Vertex.B),]

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
  j_result <- get.shortest.paths(jCombine,from = 'MAPK10',to = jGene,output = 'both',
                                 mode='out')
  j_tf_subnetwork <- edges[j_result$epath[[1]],]
  j_index <- nrow(j_tf_subnetwork)
  if (j_index>0) {
    j_all <- all_simple_paths(graph = jCombine, from = 'MAPK10',to = jGene,
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


jnk_index4 <- which(jnk_scorenetwork$common_score != 1)
jnk_verify <- jnk_scorenetwork[jnk_index4,]
verify99 <- jnk_verify[which(jnk_verify$common_score < quantile(jnk_verify$common_score,0.01)),]
write.xlsx(verify99,'MAPK10_verify99.xlsx',row.names = F)
write.csv(jnk_scorenetwork,'JNK_scorenetwork_TPMmax>1.csv',row.names = F)
write.csv(jnk_verify,'Metformin_verify_top10.csv',row.names = F)


### Find shortest path ###
jnkRMinds <- which(jnk_scorenetwork$Vertex.B == 'PKLR' & jnk_scorenetwork$Relationship != "Regulatory")
if (length(jnkRMinds)>0) {
  jnkEdges <- data.frame(from = jnk_scorenetwork$Vertex.A[-jnkRMinds],
                         to = jnk_scorenetwork$Vertex.B[-jnkRMinds])
  jnkCommonScore = jnk_scorenetwork$common_score[-jnkRMinds]
} else {
  jnkEdges <- data.frame(from = jnk_scorenetwork$Vertex.A,
                         to = jnk_scorenetwork$Vertex.B)
  jnkCommonScore = jnk_scorenetwork$common_score
}

jnk_combine <- graph_from_data_frame(jnkEdges, directed = T)
jnk_result <- get.shortest.paths(jnk_combine,from = 'MAPK10', to = 'PKLR', mode = 'out',
                                 output = 'both', weights = jnkCommonScore)
jnk_subnetwork <- jnkEdges[jnk_result$epath[[1]],]


## labeling shortest path results
ver.A1 <- jnk_subnetwork[,1]
ver.B1 <- jnk_subnetwork[,2]

result1 <- NULL
for (i in 1:length(ver.A1)) {
  index <- match(jnk_scorenetwork$Vertex.A,ver.A1[i])
  each <- jnk_scorenetwork[which(!is.na(index)),]
  result1 <- rbind(result1, each)
}

result2 <- NULL
for (j in 1:length(ver.B1)) {
  index2 <- match(result1$Vertex.B, ver.B1[j])
  each2 <- result1[which(!is.na(index2)),]
  result2 <- rbind(result2,each2)
}
index3 <- which(result2$Vertex.A == result2$Vertex.B)
jnk_shortest <- result2[-index3,]
jnk_shortest <- jnk_shortest[c(-3),]
write.xlsx(jnk_shortest[1:2,],'TGFB1_SMAD3_CREBBP.xlsx',row.names = F)

merge_file <- NULL
for (k in 1:length(result)) {
  each_a <- data.frame(unlist(result[k],use.name=T))
  index <- as.numeric(nrow(each_a))
  each_a$name <- rownames(each_a)
  rownames(each_a) <- NULL
  each_b <- data.frame(each_a[-1,-1])
  each_a <- data.frame(each_a[-index,-1])
  merge_1 <- cbind(each_a,each_b)
  merge_file <- rbind(merge_file,merge_1)
}


jnk_shortest <- NULL
for (i in 1:nrow(merge_file)) {
  print(i)
  index1 <- match(scorenetwork$Vertex.A,merge_file[i,1])
  index2 <- match(scorenetwork$Vertex.B,merge_file[i,2])
  index <- which(!is.na(index1) & !is.na(index2))
  jnk_shortest <- rbind(jnk_shortest,scorenetwork[index,])
}
jnk_shortest <- jnk_shortest[!duplicated(jnk_shortest),]


write.xlsx(jnk_shortest,file = '~/Desktop/Sweden/MAPK10/all_shortest(unique).xlsx',row.names = F)

