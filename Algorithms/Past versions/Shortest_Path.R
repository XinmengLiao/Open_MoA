library(igraph)

jnk_scorenetwork <- read.csv('TGFb1_specific_weighted_subnetwork.csv', header = T)

### the potential core pathway prediction (the shortest path)

jnkRMinds <- which(jnk_scorenetwork$Vertex.B == 'EP300' & jnk_scorenetwork$Relationship != "Regulatory")
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
jnk_result <- get.shortest.paths(jnk_combine,from = 'TGFB1', to = 'EP300', mode = 'out',
                                 output = 'both', weights = jnkCommonScore)
jnk_shortest <- jnkEdges[jnk_result$epath[[1]],]


ver.A <- jnk_subnetwork[,1]
ver.B <- jnk_subnetwork[,2]

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
result2$PoA <- 1-result2$common_score

dup <- which(result2$Vertex.A == result2$Vertex.B)
result2 <- result2[-dup,]
shortest_path <- result2


write.csv(shortest_path,'TGFB1_EP300.csv',row.names = F)
