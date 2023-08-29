library(xlsx)
library(igraph)

### Weighted Network Construction
jnk_scorenetwork <- read.csv('~/Desktop/Sweden/Scorenetwork_3.0/HepG2_specific.csv',header = T)
jnk_scorenetwork <- unique(jnk_scorenetwork)
which(jnk_scorenetwork$common_score != 1)
#jnk_fdr <- read.csv('~/Desktop/Sweden/Scorenetwork 2.0/Metformin/HEPG2_gene_padj_final.txt',sep = "\t")
jnk_fdr<- read.table('WNT1/Wnt1_sh4.csv',sep = ",",header = T)
jnk_fdr <- jnk_fdr[,c(1,2)]
colnames(jnk_fdr)[2] <- 'FDR'
colnames(jnk_fdr)[1] <- 'Gene'
jnk_fdr <- jnk_fdr[which(jnk_fdr$FDR<1),]
jnk_fdr <- jnk_fdr[which(jnk_fdr$FDR != "NA"),]
jnk_fdr <- unique(jnk_fdr)



#jnk_fdr <- read.csv('Wnt1_sh4.csv',sep = ",")
#jnk_fdr <- jnk_fdr[which(jnk_fdr$FDR<1),]
#jnk_fdr <- jnk_fdr[which(jnk_fdr$FDR != "NA"),]

jnk_fdr <- jnk_fdr[which(jnk_fdr$Gene %in% jnk_scorenetwork$Vertex.A | jnk_fdr$Gene %in% jnk_scorenetwork$Vertex.B),]
length(which(duplicated(jnk_fdr$Gene) == T))
#jnk_fdr <- jnk_fdr[which(jnk_fdr$Gene %in% jnk_scorenetwork$Vertex.A | jnk_fdr$Gene %in% jnk_scorenetwork$Vertex.B),]


###删掉jnk_fdr中无regulatory关系的点
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

####################
### for loop new ###
####################
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
  j_result <- get.shortest.paths(jCombine,from = 'WNT1',to = jGene,output = 'both',
                                 mode='out')
  j_tf_subnetwork <- edges[j_result$epath[[1]],]
  j_index <- nrow(j_tf_subnetwork)
  if (j_index>0) {
    j_all <- all_simple_paths(graph = jCombine, from = 'WNT1',to = jGene,
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


original_scorenetwork <- read.csv('~/Desktop/Sweden/Scorenetwork 2.0/Wnt1/WNT_scorenetwork_TPMmax>1.csv',header = T)
original_verify <- read.csv('~/Desktop/Sweden/Scorenetwork 2.0/MAPK10/scorenetwork4.0/JNK_verify.csv')

jnk_verify <- jnk_scorenetwork[which(jnk_scorenetwork$common_score != 1),]
write.csv(jnk_scorenetwork,'TGFB1/TGFb1_scorenetwork_top10_TPMmax>1.csv',row.names = F)
write.csv(jnk_verify,'TGFB1/TGFb_verify_top10.csv',row.names = F)


## import data 

setwd("~/Desktop/Sweden/Scorenetwork_3.0")
jnk_scorenetwork <- read.csv('WNT1/wnt1_scorenetwork_high_TPMmax>1.csv')

#将所有连接PKLR点的PPI关系删掉,重新指定最后的一个关系是regulation
jnkRMinds <- which(jnk_scorenetwork$Vertex.B == 'CCND3' & jnk_scorenetwork$Relationship != "Regulatory")
#这里删除了Vertex B是PKLR且关系是PPI的所有点
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
jnk_result <- get.shortest.paths(jnk_combine,from = 'WNT1', to = 'CCND3', mode = 'out',
                                 output = 'both', weights = jnkCommonScore)
jnk_subnetwork <- jnkEdges[jnk_result$epath[[1]],]



## 给shortest path 加 label 
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
result2 <- result2[c(1,4),]
result2$PoA <- 1-result2$common_score


write.xlsx(result2,'WNT1/WNT1_shortest.xlsx',row.names = F,sheetName = 'WNT1_CCND3',append = T)


##去掉label中重复的行和内容
a <- data.frame(jnk_verify[,1])
colnames(a) <- 'a'
b <- data.frame(jnk_verify[,3])
colnames(b) <- 'a'
c <- rbind(a,b)
c <- c[!duplicated(c),]


### get all shortest path
## all shortst path
## MAPK10
scorenetwork <- read.csv('WNT1/wnt1_scorenetwork_high_TPMmax>1.csv',header = T)

jnkRMinds <- which(scorenetwork$Vertex.B == 'CCND1' & scorenetwork$Relationship != "Regulation")
#这里删除了Vertex B是PKLR且关系是PPI的所有点
if (length(jnkRMinds)>0) {
  jnkEdges <- data.frame(from = scorenetwork$Vertex.A[-jnkRMinds],
                         to = scorenetwork$Vertex.B[-jnkRMinds])
  jnkCommonScore = scorenetwork$common_score[-jnkRMinds]
} else {
  jnkEdges <- data.frame(from = scorenetwork$Vertex.A,
                         to = scorenetwork$Vertex.B)
  jnkCommonScore = scorenetwork$common_score
}

jnk_combine <- graph_from_data_frame(jnkEdges, directed = T)

result <- all_simple_paths(jnk_combine, from = 'CTNNB1', to = 'CCND1',mode = 'out', cutoff = 2)

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


#### functional analysis
verify99 <- jnk_verify[which(jnk_verify$common_score < quantile(jnk_verify$common_score,0.01)),]
verify995 <- jnk_verify[which(jnk_verify$common_score < quantile(jnk_verify$common_score,0.005)),]
verify99$Rank <- 297-rank(verify99$common_score)
write.xlsx(verify995,'Metformin-verify995.xlsx',row.names = F)
write.xlsx(verify99,'WNT_verify99.xlsx',row.names = F)


length(which(verify99$Relationship == 'PPI'))
length(which(verify99$Relationship == 'Regulatory'))
length(which(verify99$Vertex.A == 'SMAD3'))

gene99 <- verify99
gene1 <- data.frame(gene99[,1])
colnames(gene1) <- 'a'
gene2 <- data.frame(gene99[,3])
colnames(gene2) <- 'a'
gene3 <- rbind(gene1, gene2)
gene3 <- data.frame(unique(gene3))

write.table(gene3,file = 'WNT_gene99.txt',row.names = F,quote = F)

a <- read.csv
a <- a[which(a$Relationship != 'Drug_Gene'),]
drugabn <- jnk_scorenetwork[which(jnk_scorenetwork$Relationship  == 'Drug_Gene'),]

jnk_scorenetwork <- read.csv('wnt1_hepg2_scorenetwork_top10.csv')

a <- jnk_verify$Vertex.A
b <- jnk_verify$Vertex.B
c <- unique(append(a,b))
length(c)


a <- verify99$Vertex.A
b <- verify99$Vertex.B
c <- unique(append(a,b))
length(c)


### calculate nodes and edges numbers in the subnetwork 
jnk_scorenetwork <- read.csv('WNT1/wnt1_scorenetwork_high_TPMmax>1.csv')
verify <- jnk_scorenetwork[which(jnk_scorenetwork$common_score != 1),]
verify99 <- verify[which(verify$common_score < quantile(verify$common_score, 0.01)),]

node1 <- as.matrix(verify99$Vertex.A)
node2 <- as.matrix(verify99$Vertex.B)
node <- as.matrix(rbind(node1, node2))
length(unique(node))

write.table(node, 'nodes99.txt',quote = F,row.names = F, col.names = F)
