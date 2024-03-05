
if (!require("igraph", quietly = TRUE)) {
  install.packages("igraph")
  library(igraph)
} else {
  library(igraph)
}

###-------------------------------###
### Weighted Network Construction
###-------------------------------###

weighted.network = function(unweighted.network, starting.point, DEGs){
  
  # Create igraph object
  edges = data.frame(from = unweighted.network$Vertex.A,
                     to = unweighted.network$Vertex.B)
  combine <- graph_from_data_frame(edges, directed = T)
  
  # penalty score calculation 
  for (j in 1:length(DEGs)){
    print(j)
    jGene <- unlist(DEGs)[j]
    jRMinds <- which(unweighted.network$Vertex.B == jGene & unweighted.network$Relationship != "Regulatory")
    if (length(jRMinds)>0) {
      jEdges <- data.frame(from = unweighted.network$Vertex.A[-jRMinds],
                           to = unweighted.network$Vertex.B[-jRMinds])
    } else {
      jEdges <- data.frame(from = unweighted.network$Vertex.A,
                           to = unweighted.network$Vertex.B)
    }
    jCombine <- graph_from_data_frame(jEdges, directed = T)
    j_result <- get.shortest.paths(jCombine,from = starting.point,to = jGene,output = 'both',
                                   mode='out')
    j_tf_subnetwork <- edges[j_result$epath[[1]],]
    j_index <- nrow(j_tf_subnetwork)
    if (j_index>0) {
      j_all <- all_simple_paths(graph = jCombine, from = starting.point ,to = jGene,
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
        merge_1$Penalty_Score <- jnk_fdr$FDR[j]^(1/length(j_all)) 
        for (i in 1:nrow(merge_1)) {
          i_index1 <- which(unweighted.network$Vertex.A == merge_1$Vertex.A[i] & unweighted.network$Vertex.B == merge_1$Vertex.B[i])
          unweighted.network$Penalty_Score[i_index1] <- unweighted.network$Penalty_Score[i_index1] * merge_1$Penalty_Score[i]
        }
      }
    }
  }
  
  # Covert penalty score to confidence score
  unweighted.network$Confidence_Score <-  1 - unweighted.network$Penalty_Score
  weighted.network <- unweighted.network
  
  return(weighted.network)
}



###-------------------------------###
### Potential pathway prediction
### (Shortest path identification)
###-------------------------------###


pathway.prediction <- function(weighted.network,starting.point, endpoint){
  jnkRMinds <- which(weighted.network$Vertex.B == starting.point & weighted.network$Relationship != "Regulatory")
  if (length(jnkRMinds)>0) {
    jnkEdges <- data.frame(from = weighted.network$Vertex.A[-jnkRMinds],
                           to = weighted.network$Vertex.B[-jnkRMinds])
    jnkCommonScore = weighted.network$Penalty_Score[-jnkRMinds]
  } else {
    jnkEdges <- data.frame(from = weighted.network$Vertex.A,
                           to = weighted.network$Vertex.B)
    jnkCommonScore = weighted.network$Penalty_Score
  }
  
  jnk_combine <- graph_from_data_frame(jnkEdges, directed = T)
  jnk_result <- get.shortest.paths(jnk_combine,from = starting.point, to = endpoint, mode = 'out',
                                   output = 'both', weights = jnkCommonScore)
  jnk_shortest <- jnkEdges[jnk_result$epath[[1]],]
  
  
  ver.A <- jnk_shortest[,1]
  ver.B <- jnk_shortest[,2]
  
  result1 <- NULL
  for (i in 1:length(ver.A)) {
    index <- match(jnk_scorenetwork$Vertex.A,ver.A[i])
    each <- jnk_scorenetwork[which(!is.na(index)),]
    result1 <- rbind(result1, each)
  }
  
  result2 <- NULL
  for (j in 1:length(ver.B)) {
    index2 <- match(result1$Vertex.B, ver.B[j])
    each2 <- result1[which(!is.na(index2)),]
    result2 <- rbind(result2,each2)
  }
  
  dup <- which(result2$Vertex.A == result2$Vertex.B)
  result2 <- result2[-dup,]
  prediction <- result2
  
  return(prediction)
}





