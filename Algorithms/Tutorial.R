
### Notices:
### In the tutorial, we are going to predict the Mechanism of Action (MoA) of JNK-IN-5A. 
### Therefore, we will use the HepG2 cell line specific Integrated Network (Datasets/HepG2_specific.csv) and
### the JNK-IN-5A data (Test/JNK-IN-5A.csv) for constructing the weighted network. 
### Notably, MAPK10, the downstream target of JNK-IN-5A, will be used as the starting point 
### since JNK-IN-5A is not included in our drug databases.
### Here, we are going to figure out how JNK-IN-5A affects gene PKLR expression.
### Hence, MAPK10 and PKLR will be set as the starting point and endpoint, respectively. 


### Locations of files and scripts:
### 1) HepG2 specific Integrated Network: Datasets/HepG2_specific.csv
### 2) JNK-IN-5A transcriptomic data: Test/JNK-IN-5A.csv
### 3) Open MoA functions: Algorithms/OpenMoA Functions.R


if (!require("igraph", quietly = TRUE)) {
  install.packages("igraph")
  library(igraph)
} else {
  library(igraph)
}


## Read HepG2 cell line specific Integrated Network as jnk_scorenetwork
jnk_scorenetwork <- read.csv('HepG2_specific.csv',header = T)
jnk_scorenetwork <- unique(jnk_scorenetwork)

# Make sure all the penalty scores are 1.00 at the beginning 
which(jnk_scorenetwork$Penalty_Score != 1) 

## Read JNK-IN-5A.csv file as jnk_fdr
a = read.xlsx("btad666_supplementary_data/Supplementary Table 1.xlsx",sheetIndex = 2) %>% select(-2) 
names(a) = c("Gene","FDR")
a = jnk_fdr
jnk_fdr<- read.table('JNK-IN-5A.csv',sep = ",",header = T)
jnk_fdr <- jnk_fdr[which(jnk_fdr$FDR<1),]
jnk_fdr <- jnk_fdr[which(jnk_fdr$FDR != "NA"),]
jnk_fdr <- unique(jnk_fdr)
jnk_fdr <- jnk_fdr[which(jnk_fdr$Gene %in% jnk_scorenetwork$Vertex.A | jnk_fdr$Gene %in% jnk_scorenetwork$Vertex.B),]
length(which(duplicated(jnk_fdr$Gene) == T))


## Filtered out the JNK-IN-5A perturbed genes which do not have a regulatory relationship
fdr <- NULL
for(i in 1:nrow(jnk_fdr)){
  #print(i)
  index1 <- jnk_scorenetwork[which(jnk_scorenetwork$Vertex.B == jnk_fdr$Gene[i]),]
  index2 <- index1[which(index1$Relationship == 'Regulatory'),5]
  if(length(index2 > 0)){
  fdr <- rbind(fdr,jnk_fdr[i,])
  } else{
  fdr <- fdr
  }
}
jnk_fdr <- fdr
rm(fdr, index1, index2)


## Weighted network construction 
jnk_weightnetwork <-  weighted.network(unweighted.network = jnk_scorenetwork,
                                     starting.point = "MAPK10",
                                     DEGs = jnk_fdr$Gene)

## Pathway prediction
pathway <- pathway.prediction(weighted.network = jnk_weightnetwork,
                              starting.point = "MAPK10",
                              endpoint = "PKLR")

## Show results
vertex <- paste(unique(pathway$Vertex.B),collapse = " -> ")
print(paste0("The predicted pathway is: ","MAPK10"," -> " ,vertex))
